#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <openacc.h>

#include "./include/Common.h"
#include "./include/Stencil.h"
#include "./include/Check.h"
#include "./include/Force.h"
#include "./include/Viscosity.h"
#include "./include/Stokes.h"
#include "./include/Config.h"
#include "./include/IOBin.h"
#include "./include/State.h"
#include "./include/Velocity.h"

// ============================================================================
// 生产路径：配置驱动 + .fpd 断点重启
//
// 输入输出【只有】.fpd 一种格式：init 文件、检查点、重启文件都是它，
// 文件里的 step 决定从哪继续。可视化由 tools/fpd2vtk.py 离线转换。
//
// 自检与验证路径（--check / --lambda / --noise / --equipart）已迁移到
// 独立可执行 fpd_check（src/CheckMain.cpp + CheckStencil.cpp + CheckNoise.cpp）。
// ============================================================================
static int run_production(const FpdConfig& c)
{
    NS_Config cfg = make_ns_config(c);
    PhiParams pp  = make_phi_params(c);

    std::string err;
    if (!ensure_dir(c.out_dir, err)) { std::cerr << "错误: " << err << "\n"; return 2; }

    // --- 决定初始状态：来自 .fpd 还是内置默认 ---
    CkptHeader h;
    h.Nx = cfg.Nx; h.Ny = cfg.Ny; h.Nz = cfg.Nz;
    h.N = 1; h.step = 0;
    h.dt = c.dt; h.kT = c.kT; h.radius = c.radius; h.xi = c.xi; h.ratio_eta = c.ratio_eta;
    h.noise_on = c.noise_on; h.seed = c.seed; h.rng_draws = 0;

    if (!c.init_file.empty())
    {
        CkptHeader fh;
        if (!read_ckpt_header(c.init_file.c_str(), fh, err))
        { std::cerr << "错误: " << err << "\n"; return 2; }

        // 维度冲突必须报错，不静默采信任何一方
        if (fh.Nx != cfg.Nx || fh.Ny != cfg.Ny || fh.Nz != cfg.Nz)
        {
            std::cerr << "错误: " << c.init_file << " 的网格是 "
                      << fh.Nx << "x" << fh.Ny << "x" << fh.Nz
                      << "，配置文件要求 " << cfg.Nx << "x" << cfg.Ny << "x" << cfg.Nz << "\n";
            return 2;
        }
        h = fh;

        // 配置里的 seed 优先 —— 否则 --set seed=... 会被 init 文件里的值静默吞掉。
        // 逐位重启要求原运行与续跑用同一配置，因此「配置优先」正是所需语义。
        if (h.seed != c.seed)
        {
            std::cout << "[注意] " << c.init_file << " 记录的 seed 是 " << h.seed
                      << "，改用配置里的 " << c.seed << "\n";
            h.seed = c.seed;
        }
    }

    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N    = h.N;
    const long start_step = (long)h.step;

    if (start_step >= c.n_steps)
    {
        std::cout << "起始步 " << start_step << " 已达到目标 " << c.n_steps << "，无事可做\n";
        return 0;
    }

    // 状态：分配 + 装配指针 + 设备映射 + plan/gen。ST_FULL = 生产全套。
    // 映射用运行时 API（见 State.h），RNG 统一 Philox（见 State.cpp）。
    FpdState st;
    st.init(cfg, N, ST_FULL, h.seed);

    CkptArrays arr = st.ckpt_arrays();

    if (!c.init_file.empty())
    {
        CkptHeader tmp;
        if (!load_checkpoint(c.init_file.c_str(), tmp, arr, err))
        { std::cerr << "错误: " << err << "\n"; st.finish(); return 2; }
        std::cout << "从 " << c.init_file << " 载入（step=" << h.step
                  << "，N=" << N << "，rng_draws=" << h.rng_draws << "）\n";
    }
    else
    {
        // 内置默认：单粒子放盒心的一般亚格点位置
        st.Rx[0] = cfg.Nx / 2.0 + 0.3; st.Ry[0] = cfg.Ny / 2.0 + 0.7; st.Rz[0] = cfg.Nz / 2.0 + 0.5;
        st.Rux[0] = st.Rx[0]; st.Ruy[0] = st.Ry[0]; st.Ruz[0] = st.Rz[0];
        std::cout << "未指定 init_file，使用内置默认（单粒子放盒心）\n";
    }

    // st.init 在 load 之前 copyin 了初始 0；现在主机端有 load/默认 的值，推上去
    st.upload(ST_PARTICLE | ST_VELOCITY);

    std::cout << "开始模拟：" << cfg.Nx << "x" << cfg.Ny << "x" << cfg.Nz
              << "  N=" << N << "  dt=" << cfg.dt << "  kT=" << c.kT
              << "  W=" << cfg.W << "  a=" << pp.radius
              << "  eta_c/eta_l=" << pp.ratio_eta
              << "  步数 " << start_step << " -> " << c.n_steps << std::endl;

    const auto t_start = std::chrono::steady_clock::now();
    int rc = 0;

    // 只做文件 IO。arr 每次 write 现场从 st 取，不持有陈旧指针。
    struct CkptWriter
    {
        const FpdConfig& c;
        const CkptHeader& h;
        FpdState& st;

        bool write(long step) const
        {
            CkptHeader oh = h;
            oh.step = step;
            // 仅供人读：实际使用的 offset 由 step 直接算出，读回时不采信此值
            oh.rng_draws = 2ULL * (unsigned long long)step
                         * 3ULL * (unsigned long long)h.size();
            oh.flags = c.save_pressure ? FPD_FLAG_PRESSURE : 0u;

            CkptArrays oa = st.ckpt_arrays();
            if (!c.save_pressure) { oa.p = 0; }

            const std::string path = ckpt_path(c.out_dir, c.run_name, step);
            std::string e2;
            if (!save_checkpoint(path.c_str(), oh, oa, e2))
            { std::cerr << "错误: " << e2 << "\n"; return false; }
            return true;
        }
    };
    const CkptWriter writer{c, h, st};

    for (long step = start_step; step < c.n_steps; step++)
    {
        if (step % c.interval_ckpt == 0)
        {
            // 拉全检查点需要的数组（含 Fx/Fy/Fz —— 原来漏过，力恒 0 时无害，
            // 加粒子间力后会写出陈旧主机端零值；现在一步拉全，不再有遗漏）
            st.download(ST_VELOCITY | ST_PARTICLE);
            if (!writer.write(step)) { rc = 2; break; }
        }

        if (step % c.interval_log == 0)
        {
            st.download(ST_VELOCITY | ST_PARTICLE);
            bool bad = false;
            for (int i = 0; i < size; i++)
            { if (!std::isfinite(st.vx[i])) { bad = true; break; } }
            if (!std::isfinite(st.Rx[0])) { bad = true; }

            const double el = std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - t_start).count();
            const long done = step - start_step;
            std::cout << "step " << std::setw(9) << step
                      << "   R = (" << std::fixed << std::setprecision(4)
                      << st.Rx[0] << ", " << st.Ry[0] << ", " << st.Rz[0] << ")"
                      << "   " << std::setprecision(1) << (done > 0 ? done / el : 0.0) << " 步/秒"
                      << (bad ? "   *** 检测到 NaN/Inf，中止 ***" : "")
                      << std::endl;

            // 发散时先落盘再走正常清理 —— 原来直接 return 2 会跳过所有 exit data
            if (bad)
            {
                writer.write(step);
                rc = 3;
                break;
            }
        }

        update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz, st.fx, st.fy, st.fz);
        step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p, st.fx, st.fy, st.fz,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz, st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan, st.gen, st.randD, st.randN,
                           st.tmp_fx, st.tmp_fy, st.tmp_fz, step);
        update_particle_velocity(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                 st.sum_phix, st.sum_phiy, st.sum_phiz,
                                 st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);
        update_particle_position(cfg, N, st.Rx, st.Ry, st.Rz,
                                 st.Rux, st.Ruy, st.Ruz, st.Vx, st.Vy, st.Vz);
    }

    // 终态检查点（目标步），失败路径已在上面落过盘
    if (rc == 0)
    {
        st.download(ST_VELOCITY | ST_PARTICLE);
        if (!writer.write(c.n_steps)) { rc = 2; }
    }

    st.finish();

    std::cout << (rc == 0 ? "完成" : (rc == 3 ? "因发散中止（已落盘）" : "因错误中止")) << std::endl;
    return rc;
}

int main(int argc, char** argv)
{
    // 自检与验证路径已迁移到 fpd_check。保留识别旧命令，但只打印迁移提示
    // —— 绝不静默改变行为（静默是 CLAUDE.md 反复点名的事故类别）。
    if (argc > 1 && (std::strcmp(argv[1], "--check") == 0 ||
                     std::strcmp(argv[1], "--lambda") == 0 ||
                     std::strcmp(argv[1], "--noise") == 0 ||
                     std::strcmp(argv[1], "--equipart") == 0))
    {
        std::cerr << "自检/验证已迁移：请用 ./build/fpd_check " << argv[1] << " ...\n";
        return 2;
    }

    // --- 生产路径：<配置文件> [--set key=value ...] ---
    if (argc < 2)
    {
        std::cerr << "用法: fpd <配置文件> [--set key=value ...]\n\n";
        print_config_help();
        return 2;
    }

    FpdConfig fc;
    std::string err;
    if (!load_config(argv[1], fc, err)) { std::cerr << "错误: " << err << "\n"; return 2; }

    for (int i = 2; i < argc; i++)
    {
        if (std::strcmp(argv[i], "--set") == 0 && i + 1 < argc)
        {
            if (!apply_override(fc, argv[++i], err))
            { std::cerr << "错误: " << err << "\n"; return 2; }
        }
        else
        {
            std::cerr << "错误: 无法识别的参数 '" << argv[i] << "'\n";
            return 2;
        }
    }

    if (!validate_config(fc, err)) { std::cerr << "错误: " << err << "\n"; return 2; }

    if (!ensure_dir(fc.out_dir, err)) { std::cerr << "错误: " << err << "\n"; return 2; }
    // 用 run_name 前缀，避免同一 out_dir 跑两个 run_name 时互相覆盖
    const std::string used = fc.out_dir + "/" + fc.run_name + ".config.used";
    if (!dump_config(fc, used.c_str(), err)) { std::cerr << "错误: " << err << "\n"; return 2; }

    return run_production(fc);
}
