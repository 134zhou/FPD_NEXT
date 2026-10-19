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
#include "Config.h"
#include "IOBin.h"
#include "State.h"
#include "Potential.h"
#include "Poisson.h"
#include "Step.h"

// ============================================================================
// 生产路径：配置参数 + .fpd 初态/断点重启
//
// 输入输出【只有】.fpd 一种格式：init 文件、检查点、重启文件都是它，
// 文件里的 step 决定从哪继续。可视化由 tools/fpd2vtk.py 离线转换。
//
// 自检与验证路径（--check / --lambda 等）已迁移到独立可执行 fpd_check
// （tests/Check*.cpp）。
// ============================================================================
static int run_production(const FpdConfig& c)
{
    std::string err;
    CkptHeader h;
    CkptReader reader;
    if (!open_checkpoint(c.init_file.c_str(), h, reader, err))
    { std::cerr << "错误: " << err << "\n"; return 2; }

    FpdConfig effective = c;
    effective.Nx = h.Nx; effective.Ny = h.Ny; effective.Nz = h.Nz;
    const StepParams sp = make_step_params(effective);
    const NS_Config cfg = sp.cfg;
    const PhiParams pp  = sp.pp;
    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N = h.N;
    const long start_step = (long)h.step;
    if (start_step >= c.n_steps)
    {
        close_checkpoint(reader);
        std::cout << "起始步 " << start_step << " 已达到目标 " << c.n_steps << "，无事可做\n";
        return 0;
    }

    FpdState st;
    st.init(cfg, N, ST_FULL, c.seed);
    CkptArrays arr = st.ckpt_arrays();
    if (!read_checkpoint(reader, h, arr, err))
    { std::cerr << "错误: " << err << "\n"; st.finish(); return 2; }
    std::cout << "从 " << c.init_file << " 载入（step=" << h.step
              << "，N=" << N << "）\n";

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
            std::cerr << "错误: " << c.init_file
                      << " 的 vz[:,:," << cfg.Nz - 1 << "] 有 " << nbad
                      << " 个非零值。上壁法向速度必须恒 0，否则破坏泊松的相容性条件。\n";
            st.finish();
            return 2;
        }
    }

    // st.init 在载入前映射了初始 0；现在把文件数组上传到设备
    st.upload(ST_PARTICLE | ST_VELOCITY);

    // 粒子间势、壁面势与重力加速度。重力是体力 rho*g，净体密度 (rho-rho_l)*g，
    // 由 step_fields 里的 apply_gravity 进入面心力场，不进粒子力（见 Step.h、Mass.h）。
    const PotentialParams& pot  = sp.pot;
    const WallParams&      wall = sp.wall;
    const ExternalField&   ext  = sp.ext;

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
        std::cout << "重力加速度: (" << ext.gx << ", " << ext.gy << ", " << ext.gz
                  << ")   rho_p/rho_l = " << c.ratio_rho
                  << "   净体密度 = (rho-rho_l)*g" << std::endl;
    }
    if (cfg.density_on)
    {
        std::cout << "密度对比: rho_p/rho_l = " << c.ratio_rho
                  << "   压力方程 div((1/rho)grad p)=b，FFT 预条件 PCG（容差 " << cfg.pcg_tol
                  << "，上限 " << PCG_MAX_ITER << " 次）" << std::endl;
    }

    std::cout << "开始模拟：" << cfg.Nx << "x" << cfg.Ny << "x" << cfg.Nz
              << "  N=" << N << "  dt=" << cfg.dt << "  kT=" << c.kT
              << "  W=" << cfg.W << "  a=" << pp.radius
              << "  eta_c/eta_l=" << pp.ratio_eta
              << "  rho_p/rho_l=" << c.ratio_rho
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
            CkptArrays oa = st.ckpt_arrays();

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
    // 检测延迟到下一个检查点/日志间隔；穿墙时中止。
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

    // 变系数压力求解的迭代数统计（只在 density_on 时有意义）
    long cg_sum = 0, cg_max = 0, cg_cnt = 0;

    for (long step = start_step; step < c.n_steps; step++)
    {
        // 力：用当前 device 上的 R 算 F（GPU 全程 device，无 host 往返）。
        step_forces(sp, st);

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
            st.download(ST_SOLVER);             // 只含 diag；变系数路径看 PCG 残差
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
                      << "   " << std::setprecision(1) << (done > 0 ? done / el : 0.0) << " 步/秒";
            if (cfg.density_on)
            {
                std::cout << "   cg=" << (cg_cnt > 0 ? (double)cg_sum / (double)cg_cnt : 0.0)
                          << "/" << cg_max << "   |Σb|=" << std::scientific
                          << std::setprecision(1) << st.diag[0] << "   残差=" << st.diag[1]
                          << "   真残差=" << st.diag[2] << std::fixed;
            }
            std::cout << (bad ? "   *** 检测到 NaN/Inf 或穿墙，中止 ***" : "")
                      << std::endl;

            cg_sum = 0; cg_max = 0; cg_cnt = 0;

            // 发散时先落盘再走正常清理 —— 原来直接 return 2 会跳过所有 exit data
            if (bad)
            {
                writer.write(step);
                rc = 3;
                break;
            }
        }

        step_fields(sp, st);
        const int cg_iter = step_fluid(cfg, st, step);
        // PCG 未收敛会让 div v 悄悄不为 0：这是算法级失败。每步都查 ——
        // solve_pressure_var 在主机侧写好了 diag[1]（见 Poisson.h），不需要额外同步。
        // 收敛 ⟺ diag[1] <= pcg_tol；NaN 也会落进「未收敛」。
        bool pcg_failed = false;
        if (cfg.density_on)
        {
            cg_sum += cg_iter; cg_cnt++;
            if (cg_iter > cg_max) { cg_max = cg_iter; }
            pcg_failed = !(st.diag[1] <= cfg.pcg_tol);
        }
        step_particles(sp, st);

        // 未收敛时把这一步走完再落盘，文件里的 v 与 R 属于同一时刻 step+1。
        if (pcg_failed)
        {
            std::cerr << "错误: step " << step << " PCG 未收敛：迭代 " << cg_iter
                      << " 次，相对残差 " << std::scientific << st.diag[1]
                      << "（容差 " << cfg.pcg_tol << "）。落盘后中止\n" << std::fixed;
            step_forces(sp, st);
            st.download(ST_VELOCITY | ST_PARTICLE);
            writer.write(step + 1);
            rc = 3;
            break;
        }
    }

    // 终态检查点（目标步），失败路径已在上面落过盘
    if (rc == 0)
    {
        // 循环里最后一次 update_particle_position 更新了 R，而 F 只在循环开头算。
        // 补算一次，保证「文件里的 F 对应文件里的 R」。
        step_forces(sp, st);
        st.download(ST_VELOCITY | ST_PARTICLE);
        if (!writer.write(c.n_steps)) { rc = 2; }
    }

    st.finish();

    // rc=3 是「物理上不该发生的事」：NaN/Inf 发散、粒子穿墙，或 PCG 未收敛。
    // 具体原因已经在上面的错误行里打印过，这里不再猜。
    std::cout << (rc == 0 ? "完成"
                          : (rc == 3 ? "中止（NaN/发散、粒子穿墙或 PCG 未收敛，已落盘）"
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
