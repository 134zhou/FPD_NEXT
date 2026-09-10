#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>

#include "./include/Tests.h"
#include "./include/Stencil.h"
#include "./include/Check.h"
#include "./include/Analysis.h"
#include "./include/Force.h"
#include "./include/Viscosity.h"
#include "./include/Velocity.h"
#include "./include/Stokes.h"
#include "./include/State.h"

// ============================================================================
// 测试 3A：纯流体噪声谱（无粒子）
//
// 不可压流体每个 k 模式只有 2 个横向自由度（纵向被投影算子消掉），
// 能量均分给出   <|v|^2> = 2 kT / (rho * dV) = 2 kT   （rho = dV = 1）
//
// 这个测试【完全绕开】相场、力投影、速度平均的全部机制，直接检验：
//   (a) W = sqrt(2kT/dt) 的标定是否正确
//   (b) 随机应力缺少 -2/3 delta delta 迹项的影响（不可压投影应当消掉迹部分）
//   (c) FFT 投影本身是否正确
//
// 【时间尺度】最慢模式弛豫 tau = (L/2pi)^2 / nu。必须用小立方盒，
// 否则 L=128 时要 40 万步才平衡。无粒子时 eta=1，稳定性限制 dt < dx^2/(2*3*eta) = 1/6，
// 可以用比生产算例大得多的步长。
// ============================================================================
int run_noise_check(int L, double dt, double kT, long n_steps)
{
    NS_Config cfg = make_ns_config(L, L, L, dt, kT, true);
    PhiParams pp  = make_phi_params(3.2, 1.0, 50.0);   // 无粒子，仅占位

    const int size = L * L * L;
    const int N    = 0;                                 // 关键：无粒子

    const double nu       = 1.0;                        // eta_l / rho
    const double tau_slow = (L/(2.0*M_PI))*(L/(2.0*M_PI))/nu;   // 最慢模式弛豫时间
    const long   n_equil  = (long)(5.0 * tau_slow / dt);        // 弃掉 5 个弛豫时间

    std::cout << "=== 测试 3A：纯流体噪声谱 ===\n"
              << "  盒子 " << L << "^3   dt = " << dt << "   kT = " << kT
              << "   W = " << cfg.W << "\n"
              << "  最慢模式 tau = " << std::fixed << std::setprecision(1) << tau_slow
              << " (= " << (long)(tau_slow/dt) << " 步)"
              << "   弃掉前 " << n_equil << " 步\n";

    if (n_steps <= n_equil)
    {
        std::cout << "  ** 总步数不足以平衡，至少需要 " << (n_equil*2) << " 步 **\n";
        return 1;
    }

    // N=0：nalloc = max(0,1) = 1，粒子数组按 1 个哑元分配
    FpdState st;
    st.init(cfg, N, ST_FULL, 20260905ULL);

    double acc_v2 = 0.0;
    long   n_samp = 0;

    for (long step = 0; step < n_steps; step++)
    {
        // N=0：eta 被重置为 1（背景值），f 被清零。正是纯流体所需。
        update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                           st.fx, st.fy, st.fz);
        step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p, st.fx, st.fy, st.fz,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz, st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan, st.plan_xy, st.tri_w, st.diag,
                           st.gen, st.randD, st.randN,
                           st.tmp_fx, st.tmp_fy, st.tmp_fz, step);

        // 平衡后每 100 步采一次样（去关联）
        if (step >= n_equil && step % 100 == 0)
        {
            double s = 0.0;
            #pragma acc parallel loop reduction(+:s) present(st.vx, st.vy, st.vz)
            for (int i = 0; i < size; i++)
            {
                s += st.vx[i]*st.vx[i] + st.vy[i]*st.vy[i] + st.vz[i]*st.vz[i];
            }
            acc_v2 += s / (double)size;
            n_samp++;
        }

        if (step % 10000 == 0)
        {
            std::cout << "  step " << std::setw(7) << step
                      << (n_samp > 0
                          ? ("   <|v|^2>/2kT = " + std::to_string(acc_v2/n_samp/(2.0*kT)))
                          : std::string("   (平衡中)"))
                      << std::endl;
        }
    }

    const double v2       = acc_v2 / (double)n_samp;
    const double v2_theory = 2.0 * kT;
    const double ratio    = v2 / v2_theory;

    std::cout << std::scientific << std::setprecision(6)
              << "  <|v|^2> 实测 = " << v2 << "   (" << n_samp << " 个样本)\n"
              << "  <|v|^2> 理论 = " << v2_theory << "  = 2kT\n"
              << std::fixed << std::setprecision(4)
              << "  比值 = " << ratio << "   偏差 " << std::showpos
              << (ratio-1.0)*100.0 << "%" << std::noshowpos
              << (std::fabs(ratio-1.0) < 0.02 ? "   PASS" : "   FAIL") << "\n";

    st.finish();

    return (std::fabs(ratio - 1.0) < 0.02) ? 0 : 1;
}

// ============================================================================
// 测试 3B：粒子能量均分   <|V|^2> = 2kT/M_i - 2kT/(rho L^3)
//
// 三个配置逐级加一个变量，出问题能直接定位：
//   ghost  eta_c=1，位置冻结  —— 纯运动学 + 离散化
//   frozen eta_c=50，位置冻结 —— 加入变粘度噪声
//   moving eta_c=50，位置更新 —— 完整物理
//
// 【不需要任何代码分支】Viscosity.cpp 的 d_eta = ratio_eta - 1（C3 修复的产物）
// 让 ghost 靠 ratio_eta=1.0 零成本实现，且走【完全相同】的代码路径
// （sum_phi 照算、stencil_point 照跑），只是 eta 场是平的。
//
// 核心问题：3A 测得的 O(dt) 偏差系数 335 是在 eta=1 下的。粒子内 eta_c=50 让
// 扩散数 nu*dt/dx^2 大 50 倍。若偏差按此放大，dt=0.002 下比值约 1.34；
// 若是流体式，约 1.007。10% 精度即可干净判别。
// ============================================================================
int run_equipartition_check(EquipartMode mode, int L, double dt, double kT,
                            long n_steps, unsigned long long seed)
{
    const double ratio_eta = (mode == EQ_GHOST) ? 1.0 : 50.0;
    const bool   moving    = (mode == EQ_MOVING);
    const char*  mname     = (mode == EQ_GHOST) ? "ghost"
                           : (mode == EQ_FROZEN) ? "frozen" : "moving";

    NS_Config cfg = make_ns_config(L, L, L, dt, kT, true);
    PhiParams pp  = make_phi_params(3.2, 1.0, ratio_eta);

    const int  N    = 1;
    const double nu = 1.0;                                   // 溶剂才是最慢模式
    const double tau_slow = (L/(2.0*M_PI))*(L/(2.0*M_PI))/nu;
    const long n_equil = (long)(5.0 * tau_slow / dt);

    // 采样间隔在【物理时间】上固定为 0.5 —— 否则 dt 越小样本间越相关，
    // 同样的物理时间会得到更多但更相关的样本，分块平台反而更难出现。
    const int samp_every = (int)(0.5 / dt + 0.5) > 0 ? (int)(0.5 / dt + 0.5) : 1;

    // 常量表与目标值
    FpdConstants fc = compute_fpd_constants(pp, cfg,
                                            L/2.0 + 0.3, L/2.0 + 0.7, L/2.0 + 0.5, 4);
    const double target = equipartition_target(fc, kT, cfg);

    std::printf("=== 测试 3B：粒子能量均分  [%s] ===\n", mname);
    std::printf("  盒子 %d^3   dt = %g   kT = %g   W = %.4f   eta_c/eta_l = %.1f   seed = %llu\n",
                L, dt, kT, cfg.W, ratio_eta, seed);
    std::printf("  M_i = %.4f (亚格点散布 %.4f%%)   目标 <|V|^2> = %.6e\n",
                fc.mass_mean, fc.spread_subgrid*100.0, target);
    std::printf("  最慢模式 tau = %.1f (%ld 步)   弃掉前 %ld 步   每 %d 步采样 (dt_samp = %.2f)\n",
                tau_slow, (long)(tau_slow/dt), n_equil, samp_every, samp_every*dt);
    std::printf("  采样物理时长 T = %.0f   预期样本数 %ld\n",
                (n_steps - n_equil)*dt, (n_steps - n_equil)/samp_every);

    if (n_steps <= n_equil + 600L*samp_every)
    {
        std::printf("  ** 总步数不足（分块平台至少需 ~600 个样本），建议 >= %ld 步 **\n",
                    n_equil + 2000L*samp_every);
        return 1;
    }

    FpdState st;
    st.init(cfg, N, ST_FULL, seed);

    st.Rx[0] = L/2.0 + 0.3; st.Ry[0] = L/2.0 + 0.7; st.Rz[0] = L/2.0 + 0.5;
    st.Rux[0] = st.Rx[0]; st.Ruy[0] = st.Ry[0]; st.Ruz[0] = st.Rz[0];
    st.upload(ST_PARTICLE);

    std::vector<double> s_v2, s_vx2, s_vy2, s_vz2;
    s_v2.reserve((n_steps - n_equil)/samp_every + 8);

    const auto t0 = std::chrono::steady_clock::now();

    for (long step = 0; step < n_steps; step++)
    {
        // 每步都调，让测试走生产路径（frozen 下 R 不变，开销可忽略）
        update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                           st.fx, st.fy, st.fz);
        step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p, st.fx, st.fy, st.fz,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz, st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan, st.plan_xy, st.tri_w, st.diag,
                           st.gen, st.randD, st.randN,
                           st.tmp_fx, st.tmp_fy, st.tmp_fz, step);
        update_particle_velocity(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                 st.sum_phix, st.sum_phiy, st.sum_phiz,
                                 st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);
        if (moving)
        {
            update_particle_position(cfg, N, st.Rx, st.Ry, st.Rz,
                                     st.Rux, st.Ruy, st.Ruz, st.Vx, st.Vy, st.Vz);
        }

        if (step >= n_equil && step % samp_every == 0)
        {
            st.download(ST_PARTICLE);
            s_vx2.push_back(st.Vx[0]*st.Vx[0]);
            s_vy2.push_back(st.Vy[0]*st.Vy[0]);
            s_vz2.push_back(st.Vz[0]*st.Vz[0]);
            s_v2.push_back(st.Vx[0]*st.Vx[0] + st.Vy[0]*st.Vy[0] + st.Vz[0]*st.Vz[0]);
        }

        if (step % 50000 == 0)
        {
            const double el = std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - t0).count();
            std::printf("  step %8ld   %6.1f 步/秒   样本 %zu\n",
                        step, step > 0 ? step/el : 0.0, s_v2.size());
            std::fflush(stdout);
        }
    }

    const double el = std::chrono::duration<double>(
                          std::chrono::steady_clock::now() - t0).count();

    // --- 分块平均 ---
    const int n = (int)s_v2.size();
    BlockStat bs[32];
    const int nbs = blocking_analysis(s_v2.data(), n, bs, 32);
    const BlockStat pk = pick_plateau(bs, nbs, 20);
    std::printf("\n");
    print_blocking(bs, nbs, pk);

    const double mean = pk.mean;
    const double se   = pk.stderr_mean;
    const double ratio = mean / target;
    const double rerr  = se / target;

    // 各向同性
    double mx=0, my=0, mz=0;
    for (int i = 0; i < n; i++) { mx += s_vx2[i]; my += s_vy2[i]; mz += s_vz2[i]; }
    mx/=n; my/=n; mz/=n;

    std::printf("\n  <|V|^2> 实测 = %.6e  +/- %.2e   (%d 样本, 有效独立 %.0f)\n",
                mean, se, n, pk.block_len > 0 ? (double)pk.n_block : 0.0);
    std::printf("  <|V|^2> 目标 = %.6e   = 2kT/M_i - 2kT/(rho L^3)\n", target);
    std::printf("  比值 = %.4f +/- %.4f   偏差 %+.2f%% +/- %.2f%%   [%.1f sigma]\n",
                ratio, rerr, (ratio-1.0)*100.0, rerr*100.0,
                se > 0 ? std::fabs(mean-target)/se : 0.0);
    std::printf("  各分量 <V_a^2>/(target/3) = %.4f / %.4f / %.4f  (各向同性检查)\n",
                mx/(target/3), my/(target/3), mz/(target/3));
    std::printf("  用时 %.1f 秒   %.1f 步/秒\n", el, n_steps/el);

    const bool ok = (pk.block_len > 0) && (std::fabs(mean-target) < 3.0*se);
    std::printf("  %s\n", ok ? "PASS" : (pk.block_len < 0 ? "FAIL（无平台）" : "FAIL（超 3 sigma）"));

    st.finish();

    return ok ? 0 : 1;
}
