#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <openacc.h>

#include "Common.h"
#include "Stencil.h"
#include "Check.h"
#include "Force.h"
#include "Viscosity.h"
#include "Stokes.h"
#include "Config.h"
#include "IOBin.h"
#include "State.h"
#include "Potential.h"
#include "Velocity.h"

// ============================================================================
// 生产路径：配置驱动 + .fpd 断点重启
//
// 输入输出【只有】.fpd 一种格式：init 文件、检查点、重启文件都是它，
// 文件里的 step 决定从哪继续。可视化由 tools/fpd2vtk.py 离线转换。
//
// 自检与验证路径（--check / --lambda 等）已迁移到独立可执行 fpd_check
// （tests/Check*.cpp）。
// ============================================================================
static int run_production(const FpdConfig& c)
{
    NS_Config cfg = make_ns_config(c);
    PhiParams pp  = make_phi_params(c);

    std::string err;

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
        // 旧的 z 周期存档由版本号（FPD_VERSION 1 → 2）在 read_ckpt_header 里
        // 直接拒绝，不需要额外的标志位检查 —— 见 IOBin.h 的版本说明。
        // 壁面不变量的【直接】检查在载入数组之后做（见下）。
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

    // 壁面不变量【直接】检查：上壁那一层 vz 必须逐位为 0。
    //
    // 这取代了原来的 FPD_FLAG_WALL_Z 标志位 —— 标志位从来只是这个不变量的【代理】，
    // 而直接查不变量还能抓住手工构造、被别的工具改写、或将来损坏的 v2 文件，
    // 那些情况标志位是拦不住的。破坏它的后果是泊松相容性 Σ_k b̂_k = 0 失效，
    // 表现为压力的线性漂移，且看不出来。
    {
        // IDX 是宏，内部引用裸的 Nx/Ny，必须靠同名局部变量代入（照 Wall.h 的约定）
        const int Nx = cfg.Nx, Ny = cfg.Ny;
        int nbad = 0;
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                if (st.vz[IDX(i, j, cfg.Nz - 1)] != 0.0) { nbad++; }
            }
        }
        if (nbad > 0)
        {
            std::cerr << "错误: " << (c.init_file.empty() ? std::string("内置默认初态") : c.init_file)
                      << " 的 vz[:,:," << cfg.Nz - 1 << "] 有 " << nbad
                      << " 个非零值。上壁法向速度必须恒 0，否则破坏泊松的相容性条件。\n";
            st.finish();
            return 2;
        }
    }

    // st.init 在 load 之前 copyin 了初始 0；现在主机端有 load/默认 的值，推上去
    st.upload(ST_PARTICLE | ST_VELOCITY);

    // 粒子间势、壁面势与外场
    const PotentialParams pot  = make_potential_params(c);
    const WallParams     wall = make_wall_params(c);
    const ExternalField  ext  = make_external_field(c);

    // 背景力密度补偿恒为 0：无滑移壁面本身就是真实的动量汇，均匀外场的反冲由
    // 壁面吸收，不需要 bg = −ΣF/size 去抵消「整盒漂移」。
    // ⚠️ 管道保留（update_force_field 的 bg 形参不动，见 Force.cpp 的力路径约束），
    //    只是传 0。删形参会动到力路径，代价远大于收益。
    const double bgx = 0.0, bgy = 0.0, bgz = 0.0;

    if (c.potential != "none")
    {
        std::cout << "势: " << c.potential << "   shift=" << c.pot_shift
                  << "   rcut=" << pot.rcut << std::endl;
    }
    if (c.wallpotential != "none")
    {
        std::cout << "壁面势: " << c.wallpotential << "   shift=" << c.wall_shift
                  << "   rcut=" << wall.pot.rcut << "   a=" << wall.radius
                  << "（自变量是表面间隙 h = z_壁 - a）" << std::endl;
    }
    if (ext.gx != 0.0 || ext.gy != 0.0 || ext.gz != 0.0)
    {
        std::cout << "外场: (" << ext.gx << ", " << ext.gy << ", " << ext.gz
                  << ")   背景补偿 = 0（壁面即动量汇）" << std::endl;
    }

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
            oh.flags = (c.save_pressure ? FPD_FLAG_PRESSURE : 0u);

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

    // 粒子越界检查。只在已经 download 过 ST_PARTICLE 的地方调用 ——
    // 不做每步的 device→host 同步（那会把 GPU 流水线打断）。
    // ⚠️ 代价是检测延迟到下一个检查点/日志间隔。壁面势会把粒子在碰到壁之前弹回，
    //    所以这个检查正常情况下是纯保险 —— 但保险必须留着：壁面势配软了（残余力比
    //    报警那条）就压不住，重力会把粒子按穿，那时这里必须中止而不是继续跑。
    // ⚠️ 壁面势开启时，合法区间【收紧】到「粒子表面不越过壁面」：z ∈ [-1/2+a, Nz-1/2-a]。
    //    这不是洁癖。壁面势的自变量是间隙 h，而 pair_force_over_r 对负自变量
    //    会给出【符号反向】的力（h<0 时 h² 仍为正，力把粒子往壁里推）——
    //    一条无诊断的错物理路径。所以前置条件 h > 0 必须在这里兜住。
    //    关闭壁面势时保持原判据（中心在盒内），不改变既有运行的语义。
    const bool wall_active = (wall.pot.type != POT_NONE);
    const double zlo = wall_active ? (-0.5 + wall.radius) : -0.5;
    const double zhi = wall_active ? ((double)cfg.Nz - 0.5 - wall.radius)
                                   : ((double)cfg.Nz - 0.5);

    auto check_wall_bounds = [&](long step) -> bool
    {
        for (int n = 0; n < N; n++)
        {
            if (st.Rz[n] < zlo || st.Rz[n] > zhi)
            {
                std::cerr << "错误: step " << step << " 粒子 " << n
                          << " 的 Rz = " << st.Rz[n] << " 越出允许区间 [" << zlo
                          << ", " << zhi << "]"
                          << (wall_active ? "（壁面势要求粒子表面不穿墙，h > 0）"
                                          : "（粒子中心穿墙）")
                          << "。中止（不 clamp）\n";
                return true;
            }
        }
        return false;
    };

    for (long step = start_step; step < c.n_steps; step++)
    {
        // 力：用当前 device 上的 R 算 F（GPU 全程 device，无 host 往返）。
        // potential=none 且外场为零时 F 恒 0，与改动前数值路径逐位一致。
        compute_particle_forces(cfg, pot, ext, wall, N, st.Rx, st.Ry, st.Rz,
                                st.Fx, st.Fy, st.Fz);

        if (step % c.interval_ckpt == 0)
        {
            // 拉全检查点需要的数组（含 Fx/Fy/Fz —— 原来漏过，力恒 0 时无害，
            // 加粒子间力后会写出陈旧主机端零值；现在一步拉全，不再有遗漏）
            st.download(ST_VELOCITY | ST_PARTICLE);
            if (check_wall_bounds(step)) { writer.write(step); rc = 3; break; }
            if (!writer.write(step)) { rc = 2; break; }
        }

        if (step % c.interval_log == 0)
        {
            st.download(ST_VELOCITY | ST_PARTICLE);
            bool bad = false;
            for (int i = 0; i < size; i++)
            { if (!std::isfinite(st.vx[i])) { bad = true; break; } }
            if (!std::isfinite(st.Rx[0])) { bad = true; }
            if (check_wall_bounds(step)) { bad = true; }

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
                           st.sum_phix, st.sum_phiy, st.sum_phiz, bgx, bgy, bgz,
                           st.fx, st.fy, st.fz);
        step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p, st.fx, st.fy, st.fz,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz, st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan_xy, st.tri_w, st.diag,
                           st.gen, st.randD, st.randN,
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
        // 循环里最后一次 update_particle_position 更新了 R，而 F 只在循环开头算。
        // 补算一次，保证「文件里的 F 对应文件里的 R」。
        compute_particle_forces(cfg, pot, ext, wall, N, st.Rx, st.Ry, st.Rz,
                                st.Fx, st.Fy, st.Fz);
        st.download(ST_VELOCITY | ST_PARTICLE);
        if (!writer.write(c.n_steps)) { rc = 2; }
    }

    st.finish();

    // rc=3 是「物理上不该发生的事」：NaN/Inf 发散，或壁面模式下粒子穿墙。
    // 具体原因已经在上面的错误行里打印过，这里不再猜。
    std::cout << (rc == 0 ? "完成"
                          : (rc == 3 ? "中止（NaN/发散或粒子穿墙，已落盘）"
                                     : "因错误中止")) << std::endl;
    return rc;
}

int main(int argc, char** argv)
{
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
