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
#include "./include/IO.h"
#include "./include/Velocity.h"

// ============================================================================
// Phase 1 自检：力守恒 + 亚格点不变性
//
// 判据：给定粒子受力 F，摊布到网格后必须满足 Σ_grid f_α == F_α，
// 且该等式与粒子在格胞内的分数位置无关。
//
// ⚠️ 这个判据检验的是 stencil_point 的【自洽性】，对它的内部公式错误是盲的
//    （投影和归一化会一起错）。独立的数值对照见 spike/。
// ============================================================================
static int run_force_conservation_check(NS_Config cfg, PhiParams pp)
{
    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N    = 1;

    std::vector<double> Rx(N), Ry(N), Rz(N);
    std::vector<double> Fx(N), Fy(N), Fz(N);
    std::vector<double> sum_phix(N), sum_phiy(N), sum_phiz(N);
    std::vector<double> fx(size), fy(size), fz(size);
    std::vector<double> eta(size), etaXY(size), etaYZ(size), etaZX(size);

    double *d_Rx = Rx.data(), *d_Ry = Ry.data(), *d_Rz = Rz.data();
    double *d_Fx = Fx.data(), *d_Fy = Fy.data(), *d_Fz = Fz.data();
    double *d_spx = sum_phix.data(), *d_spy = sum_phiy.data(), *d_spz = sum_phiz.data();
    double *d_fx = fx.data(), *d_fy = fy.data(), *d_fz = fz.data();
    double *d_eta = eta.data(), *d_eXY = etaXY.data(), *d_eYZ = etaYZ.data(), *d_eZX = etaZX.data();

    #pragma acc enter data create(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])
    #pragma acc enter data create(d_spx[0:N], d_spy[0:N], d_spz[0:N])
    #pragma acc enter data create(d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc enter data create(d_eta[0:size], d_eXY[0:size], d_eYZ[0:size], d_eZX[0:size])

    // 三个亚格点位置：整数、任意分数、半格
    const double px[3] = {64.0, 64.3, 64.5};
    const double py[3] = {32.0, 32.7, 32.5};
    const double pz[3] = {16.0, 16.5, 16.5};

    const double TOL = 1e-12;
    int failures = 0;

    std::cout << "=== Phase 1 自检：力守恒 + 亚格点不变性 ===\n";
    std::cout << std::scientific << std::setprecision(3);

    for (int t = 0; t < 3; t++)
    {
        Rx[0] = px[t]; Ry[0] = py[t]; Rz[0] = pz[t];
        Fx[0] = 1.0;   Fy[0] = 2.0;   Fz[0] = 3.0;

        #pragma acc update device(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])

        update_viscosity_fields(cfg, pp, N, d_Rx, d_Ry, d_Rz,
                                d_spx, d_spy, d_spz, d_eta, d_eXY, d_eYZ, d_eZX);
        update_force_field(cfg, pp, N, d_Rx, d_Ry, d_Rz, d_Fx, d_Fy, d_Fz,
                           d_spx, d_spy, d_spz, d_fx, d_fy, d_fz);

        #pragma acc update host(d_fx[0:size], d_fy[0:size], d_fz[0:size])
        #pragma acc update host(d_spx[0:N], d_spy[0:N], d_spz[0:N])

        double sx = 0.0, sy = 0.0, sz = 0.0;
        for (int i = 0; i < size; i++) { sx += fx[i]; sy += fy[i]; sz += fz[i]; }

        const double ex = std::fabs(sx - Fx[0]);
        const double ey = std::fabs(sy - Fy[0]);
        const double ez = std::fabs(sz - Fz[0]);
        const bool ok = (ex < TOL && ey < TOL && ez < TOL);
        if (!ok) { failures++; }

        std::cout << "  R = (" << std::fixed << std::setprecision(1)
                  << Rx[0] << ", " << Ry[0] << ", " << Rz[0] << ")"
                  << std::scientific << std::setprecision(3)
                  << "   |err| = " << ex << " / " << ey << " / " << ez
                  << "   " << (ok ? "PASS" : "FAIL") << "\n";

        // 各向同性 sanity check：三个方向的归一化因子应几乎相等
        const double smax = std::fmax(sum_phix[0], std::fmax(sum_phiy[0], sum_phiz[0]));
        const double smin = std::fmin(sum_phix[0], std::fmin(sum_phiy[0], sum_phiz[0]));
        const double spread = (smax - smin) / smin;
        const bool iso_ok = (spread < 0.01);
        if (!iso_ok) { failures++; }

        std::cout << "      sum_phi = " << std::fixed << std::setprecision(6)
                  << sum_phix[0] << " / " << sum_phiy[0] << " / " << sum_phiz[0]
                  << "   离散度 " << std::setprecision(3) << spread * 100.0 << "%   "
                  << (iso_ok ? "PASS" : "FAIL") << "\n";
    }

    #pragma acc exit data delete(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])
    #pragma acc exit data delete(d_spx[0:N], d_spy[0:N], d_spz[0:N])
    #pragma acc exit data delete(d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc exit data delete(d_eta[0:size], d_eXY[0:size], d_eYZ[0:size], d_eZX[0:size])

    std::cout << (failures == 0 ? "全部通过\n" : "有失败项\n");
    return failures;
}

// ============================================================================
// Phase 1 自检：多粒子重叠
//
// C2（速度平均误用全场 Σφ_n 而非粒子自己的 φ_i）在 N=1 时数学上恰好抵消，
// 曾骗过 20 万步的模拟。必须用重叠粒子才能暴露。
//
// 两个判据：
//  (a) 力守恒推广：Σ_grid f_α == Σ_n F_α[n]
//  (b) 均匀流场检验：v ≡ const 时，V_i = ∫vφ_i/∫φ_i = v 对每个粒子恒成立，
//      与重叠与否无关。C2 会让重叠处的 V 显著大于 v。
// ============================================================================
static int run_overlap_check(NS_Config cfg, PhiParams pp)
{
    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N    = 2;

    std::vector<double> Rx(N), Ry(N), Rz(N);
    std::vector<double> Fx(N), Fy(N), Fz(N);
    std::vector<double> Vx(N), Vy(N), Vz(N);
    std::vector<double> sum_phix(N), sum_phiy(N), sum_phiz(N);
    std::vector<double> fx(size), fy(size), fz(size);
    std::vector<double> eta(size), etaXY(size), etaYZ(size), etaZX(size);
    std::vector<double> vx(size), vy(size), vz(size);

    // 均匀流场
    const double V_UNIFORM = 1.0;
    for (int i = 0; i < size; i++) { vx[i] = V_UNIFORM; vy[i] = V_UNIFORM; vz[i] = V_UNIFORM; }

    double *d_Rx = Rx.data(), *d_Ry = Ry.data(), *d_Rz = Rz.data();
    double *d_Fx = Fx.data(), *d_Fy = Fy.data(), *d_Fz = Fz.data();
    double *d_Vx = Vx.data(), *d_Vy = Vy.data(), *d_Vz = Vz.data();
    double *d_spx = sum_phix.data(), *d_spy = sum_phiy.data(), *d_spz = sum_phiz.data();
    double *d_fx = fx.data(), *d_fy = fy.data(), *d_fz = fz.data();
    double *d_vx = vx.data(), *d_vy = vy.data(), *d_vz = vz.data();
    double *d_eta = eta.data(), *d_eXY = etaXY.data(), *d_eYZ = etaYZ.data(), *d_eZX = etaZX.data();

    #pragma acc enter data create(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])
    #pragma acc enter data create(d_Vx[0:N], d_Vy[0:N], d_Vz[0:N], d_spx[0:N], d_spy[0:N], d_spz[0:N])
    #pragma acc enter data create(d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc enter data create(d_eta[0:size], d_eXY[0:size], d_eYZ[0:size], d_eZX[0:size])
    #pragma acc enter data copyin(d_vx[0:size], d_vy[0:size], d_vz[0:size])

    // 间距 4.0 < 2a = 6.4，两球显著重叠
    const double SEP[3] = {4.0, 6.4, 20.0};   // 重叠 / 刚接触 / 分离
    const char*  DESC[3] = {"重叠", "刚接触", "分离"};

    int failures = 0;
    std::cout << "=== Phase 1 自检：多粒子重叠（N=2）===\n";

    for (int t = 0; t < 3; t++)
    {
        Rx[0] = 64.0 - 0.5*SEP[t]; Ry[0] = 32.3; Rz[0] = 16.0;
        Rx[1] = 64.0 + 0.5*SEP[t]; Ry[1] = 32.3; Rz[1] = 16.0;
        Fx[0] = 1.0; Fy[0] = 2.0; Fz[0] = 3.0;
        Fx[1] = -0.5; Fy[1] = 0.25; Fz[1] = 1.5;

        #pragma acc update device(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])

        update_viscosity_fields(cfg, pp, N, d_Rx, d_Ry, d_Rz,
                                d_spx, d_spy, d_spz, d_eta, d_eXY, d_eYZ, d_eZX);
        update_force_field(cfg, pp, N, d_Rx, d_Ry, d_Rz, d_Fx, d_Fy, d_Fz,
                           d_spx, d_spy, d_spz, d_fx, d_fy, d_fz);
        update_particle_velocity(cfg, pp, N, d_Rx, d_Ry, d_Rz,
                                 d_spx, d_spy, d_spz, d_vx, d_vy, d_vz,
                                 d_Vx, d_Vy, d_Vz);

        #pragma acc update host(d_fx[0:size], d_fy[0:size], d_fz[0:size])
        #pragma acc update host(d_Vx[0:N], d_Vy[0:N], d_Vz[0:N], d_eta[0:size])

        // (a) 力守恒推广到多粒子
        double sx = 0.0, sy = 0.0, sz = 0.0;
        for (int i = 0; i < size; i++) { sx += fx[i]; sy += fy[i]; sz += fz[i]; }
        const double ex = std::fabs(sx - (Fx[0] + Fx[1]));
        const double ey = std::fabs(sy - (Fy[0] + Fy[1]));
        const double ez = std::fabs(sz - (Fz[0] + Fz[1]));
        const bool f_ok = (ex < 1e-12 && ey < 1e-12 && ez < 1e-12);

        // (b) 均匀流场下 V_i 必须等于 v，与重叠无关
        double vmax_err = 0.0;
        for (int n = 0; n < N; n++)
        {
            vmax_err = std::fmax(vmax_err, std::fabs(Vx[n] - V_UNIFORM));
            vmax_err = std::fmax(vmax_err, std::fabs(Vy[n] - V_UNIFORM));
            vmax_err = std::fmax(vmax_err, std::fabs(Vz[n] - V_UNIFORM));
        }
        const bool v_ok = (vmax_err < 1e-12);

        // 粘度峰值：重叠区 Σφ 可能 > 1，η 会超过 η_c，这是 FPD 的已知行为
        double eta_max = 0.0;
        for (int i = 0; i < size; i++) { eta_max = std::fmax(eta_max, eta[i]); }

        if (!f_ok || !v_ok) { failures++; }

        std::cout << "  间距 " << std::fixed << std::setprecision(1) << SEP[t]
                  << " (" << DESC[t] << ")"
                  << std::scientific << std::setprecision(3)
                  << "   力守恒 |err| = " << std::fmax(ex, std::fmax(ey, ez))
                  << "   V 偏差 = " << vmax_err
                  << std::fixed << std::setprecision(2)
                  << "   eta_max = " << eta_max
                  << "   " << ((f_ok && v_ok) ? "PASS" : "FAIL") << "\n";
    }

    #pragma acc exit data delete(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])
    #pragma acc exit data delete(d_Vx[0:N], d_Vy[0:N], d_Vz[0:N], d_spx[0:N], d_spy[0:N], d_spz[0:N])
    #pragma acc exit data delete(d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc exit data delete(d_eta[0:size], d_eXY[0:size], d_eYZ[0:size], d_eZX[0:size])
    #pragma acc exit data delete(d_vx[0:size], d_vy[0:size], d_vz[0:size])

    std::cout << (failures == 0 ? "全部通过\n" : "有失败项\n");
    return failures;
}

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
static int run_noise_check(int L, double dt, double kT, long n_steps)
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

    std::vector<double> vx(size,0.0), vy(size,0.0), vz(size,0.0), p(size,0.0);
    std::vector<double> fx(size,0.0), fy(size,0.0), fz(size,0.0);
    std::vector<double> eta(size), etaXY(size), etaYZ(size), etaZX(size);
    std::vector<double> fft_data(size*2), randD(size*3), randN(size*3);
    std::vector<double> tfx(size), tfy(size), tfz(size);
    std::vector<double> pdx(size), pdy(size), pdz(size), pnx(size), pny(size), pnz(size);
    std::vector<double> Rx(1,0.0), Ry(1,0.0), Rz(1,0.0), Fx(1,0.0), Fy(1,0.0), Fz(1,0.0);
    std::vector<double> spx(1,1.0), spy(1,1.0), spz(1,1.0);

    double *d_vx=vx.data(), *d_vy=vy.data(), *d_vz=vz.data(), *d_p=p.data();
    double *d_fx=fx.data(), *d_fy=fy.data(), *d_fz=fz.data();
    double *d_eta=eta.data(), *d_eXY=etaXY.data(), *d_eYZ=etaYZ.data(), *d_eZX=etaZX.data();
    double *d_fft=fft_data.data(), *d_rD=randD.data(), *d_rN=randN.data();
    double *d_tfx=tfx.data(), *d_tfy=tfy.data(), *d_tfz=tfz.data();
    double *d_pdx=pdx.data(), *d_pdy=pdy.data(), *d_pdz=pdz.data();
    double *d_pnx=pnx.data(), *d_pny=pny.data(), *d_pnz=pnz.data();
    double *d_Rx=Rx.data(), *d_Ry=Ry.data(), *d_Rz=Rz.data();
    double *d_Fx=Fx.data(), *d_Fy=Fy.data(), *d_Fz=Fz.data();
    double *d_spx=spx.data(), *d_spy=spy.data(), *d_spz=spz.data();

    cufftHandle plan;
    CUFFT_CHECK(cufftPlan3d(&plan, cfg.Nz, cfg.Ny, cfg.Nx, CUFFT_Z2Z));
    curandGenerator_t gen;
    CURAND_CHECK(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
    CURAND_CHECK(curandSetPseudoRandomGeneratorSeed(gen, 20260905ULL));

    #pragma acc enter data copyin(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size])
    #pragma acc enter data copyin(d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc enter data copyin(d_Rx[0:1], d_Ry[0:1], d_Rz[0:1], d_Fx[0:1], d_Fy[0:1], d_Fz[0:1])
    #pragma acc enter data copyin(d_spx[0:1], d_spy[0:1], d_spz[0:1])
    #pragma acc enter data create(d_eta[0:size], d_eXY[0:size], d_eYZ[0:size], d_eZX[0:size])
    #pragma acc enter data create(d_fft[0:size*2], d_rD[0:size*3], d_rN[0:size*3])
    #pragma acc enter data create(d_tfx[0:size], d_tfy[0:size], d_tfz[0:size])
    #pragma acc enter data create(d_pdx[0:size], d_pdy[0:size], d_pdz[0:size])
    #pragma acc enter data create(d_pnx[0:size], d_pny[0:size], d_pnz[0:size])

    double acc_v2 = 0.0;
    long   n_samp = 0;

    for (long step = 0; step < n_steps; step++)
    {
        // N=0：eta 被重置为 1（背景值），f 被清零。正是纯流体所需。
        update_viscosity_fields(cfg, pp, N, d_Rx, d_Ry, d_Rz,
                                d_spx, d_spy, d_spz, d_eta, d_eXY, d_eYZ, d_eZX);
        update_force_field(cfg, pp, N, d_Rx, d_Ry, d_Rz, d_Fx, d_Fy, d_Fz,
                           d_spx, d_spy, d_spz, d_fx, d_fy, d_fz);
        step_navier_stokes(cfg, d_vx, d_vy, d_vz, d_p, d_fx, d_fy, d_fz,
                           d_eta, d_eXY, d_eYZ, d_eZX,
                           d_pdx, d_pdy, d_pdz, d_pnx, d_pny, d_pnz,
                           d_fft, plan, gen, d_rD, d_rN, d_tfx, d_tfy, d_tfz);

        // 平衡后每 100 步采一次样（去关联）
        if (step >= n_equil && step % 100 == 0)
        {
            double s = 0.0;
            #pragma acc parallel loop reduction(+:s) present(d_vx, d_vy, d_vz)
            for (int i = 0; i < size; i++)
            {
                s += d_vx[i]*d_vx[i] + d_vy[i]*d_vy[i] + d_vz[i]*d_vz[i];
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

    #pragma acc exit data delete(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size])
    #pragma acc exit data delete(d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc exit data delete(d_Rx[0:1], d_Ry[0:1], d_Rz[0:1], d_Fx[0:1], d_Fy[0:1], d_Fz[0:1])
    #pragma acc exit data delete(d_spx[0:1], d_spy[0:1], d_spz[0:1])
    #pragma acc exit data delete(d_eta[0:size], d_eXY[0:size], d_eYZ[0:size], d_eZX[0:size])
    #pragma acc exit data delete(d_fft[0:size*2], d_rD[0:size*3], d_rN[0:size*3])
    #pragma acc exit data delete(d_tfx[0:size], d_tfy[0:size], d_tfz[0:size])
    #pragma acc exit data delete(d_pdx[0:size], d_pdy[0:size], d_pdz[0:size])
    #pragma acc exit data delete(d_pnx[0:size], d_pny[0:size], d_pnz[0:size])
    CUFFT_CHECK(cufftDestroy(plan));
    CURAND_CHECK(curandDestroyGenerator(gen));

    return (std::fabs(ratio - 1.0) < 0.02) ? 0 : 1;
}

int main(int argc, char** argv)
{
    // Phase 2 会把这些搬进配置文件
    const int    Nx = 128, Ny = 64, Nz = 32;
    const double dt        = 0.001;
    const double kT        = 0.05;   // 对应旧代码的 W=10（见 README C6，尚未标定）
    const bool   noise_on  = true;
    const double radius    = 3.2;
    const double xi        = 1.0;
    const double ratio_eta = 50.0;
    const long   n_steps   = 200000;

    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, kT, noise_on);
    PhiParams pp  = make_phi_params(radius, xi, ratio_eta);

    if (argc > 1 && std::strcmp(argv[1], "--check") == 0)
    {
        int fails = run_force_conservation_check(cfg, pp);
        std::cout << "\n";
        fails += run_overlap_check(cfg, pp);
        return fails;
    }

    // 测试 3A：--noise [L] [dt] [kT] [steps]
    if (argc > 1 && std::strcmp(argv[1], "--noise") == 0)
    {
        const int    L_    = (argc > 2) ? std::atoi(argv[2]) : 32;
        const double dt_   = (argc > 3) ? std::atof(argv[3]) : 0.01;
        const double kT_   = (argc > 4) ? std::atof(argv[4]) : 1.0;
        const long   nst_  = (argc > 5) ? std::atol(argv[5]) : 200000;
        return run_noise_check(L_, dt_, kT_, nst_);
    }

    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N    = 1;

    std::vector<double> vx(size, 0.0), vy(size, 0.0), vz(size, 0.0), p(size, 0.0);
    std::vector<double> fx(size, 0.0), fy(size, 0.0), fz(size, 0.0);

    std::vector<double> Rx(N, 0.0), Ry(N, 0.0), Rz(N, 0.0);      // 粒子位置
    std::vector<double> Rux(N, 0.0), Ruy(N, 0.0), Ruz(N, 0.0);   // 不折叠位置（MSD 用）
    std::vector<double> Vx(N, 0.0), Vy(N, 0.0), Vz(N, 0.0);
    std::vector<double> Fx(N, 0.0), Fy(N, 0.0), Fz(N, 0.0);

    std::vector<double> sum_phix(N, 0.0), sum_phiy(N, 0.0), sum_phiz(N, 0.0);
    std::vector<double> eta(size), etaXY(size), etaYZ(size), etaZX(size);
    std::vector<double> fft_data(size * 2);
    std::vector<double> randD(size * 3), randN(size * 3);

    std::vector<double> h_tmp_fx(size), h_tmp_fy(size), h_tmp_fz(size);
    std::vector<double> pi_dx(size), pi_dy(size), pi_dz(size);
    std::vector<double> pi_nx(size), pi_ny(size), pi_nz(size);

    Rx[0] = 64.0; Ry[0] = 32.0; Rz[0] = 16.0;
    Rux[0] = Rx[0]; Ruy[0] = Ry[0]; Ruz[0] = Rz[0];

    double* d_vx = vx.data(); double* d_vy = vy.data(); double* d_vz = vz.data();
    double* d_p  = p.data();
    double* d_fx = fx.data(); double* d_fy = fy.data(); double* d_fz = fz.data();

    double* d_Rx = Rx.data(); double* d_Ry = Ry.data(); double* d_Rz = Rz.data();
    double* d_Rux = Rux.data(); double* d_Ruy = Ruy.data(); double* d_Ruz = Ruz.data();
    double* d_Vx = Vx.data(); double* d_Vy = Vy.data(); double* d_Vz = Vz.data();
    double* d_Fx = Fx.data(); double* d_Fy = Fy.data(); double* d_Fz = Fz.data();

    double* d_spx = sum_phix.data();
    double* d_spy = sum_phiy.data();
    double* d_spz = sum_phiz.data();
    double* d_eta = eta.data();
    double* d_etaXY = etaXY.data();
    double* d_etaYZ = etaYZ.data();
    double* d_etaZX = etaZX.data();
    double* d_fft = fft_data.data();

    double* tmp_fx = h_tmp_fx.data();
    double* tmp_fy = h_tmp_fy.data();
    double* tmp_fz = h_tmp_fz.data();

    double* d_pi_dx = pi_dx.data(); double* d_pi_dy = pi_dy.data(); double* d_pi_dz = pi_dz.data();
    double* d_pi_nx = pi_nx.data(); double* d_pi_ny = pi_ny.data(); double* d_pi_nz = pi_nz.data();

    double* d_randD = randD.data(); double* d_randN = randN.data();

    cufftHandle plan;
    CUFFT_CHECK(cufftPlan3d(&plan, cfg.Nz, cfg.Ny, cfg.Nx, CUFFT_Z2Z));
    curandGenerator_t gen;
    CURAND_CHECK(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
    CURAND_CHECK(curandSetPseudoRandomGeneratorSeed(gen, 1234ULL));

    #pragma acc enter data copyin(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size], d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc enter data copyin(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Rux[0:N], d_Ruy[0:N], d_Ruz[0:N])
    #pragma acc enter data copyin(d_Vx[0:N], d_Vy[0:N], d_Vz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])

    #pragma acc enter data create(d_spx[0:N], d_spy[0:N], d_spz[0:N])
    #pragma acc enter data create(d_eta[0:size], d_etaXY[0:size], d_etaYZ[0:size], d_etaZX[0:size])
    #pragma acc enter data create(d_fft[0:size*2])

    #pragma acc enter data create(tmp_fx[0:size], tmp_fy[0:size], tmp_fz[0:size])
    #pragma acc enter data create(d_pi_dx[0:size], d_pi_dy[0:size], d_pi_dz[0:size], d_pi_nx[0:size], d_pi_ny[0:size], d_pi_nz[0:size], d_randD[0:size*3], d_randN[0:size*3])

    std::cout << "开始模拟：" << cfg.Nx << "x" << cfg.Ny << "x" << cfg.Nz
              << "  N=" << N << "  dt=" << cfg.dt << "  kT=" << kT
              << "  W=" << cfg.W << "  a=" << pp.radius
              << "  eta_c/eta_l=" << pp.ratio_eta << std::endl;

    const auto t_start = std::chrono::steady_clock::now();

    for (long step = 0; step < n_steps; step++)
    {
        update_viscosity_fields(cfg, pp, N, d_Rx, d_Ry, d_Rz,
                                d_spx, d_spy, d_spz,
                                d_eta, d_etaXY, d_etaYZ, d_etaZX);

        update_force_field(cfg, pp, N, d_Rx, d_Ry, d_Rz,
                           d_Fx, d_Fy, d_Fz,
                           d_spx, d_spy, d_spz,
                           d_fx, d_fy, d_fz);

        step_navier_stokes(cfg, d_vx, d_vy, d_vz, d_p, d_fx, d_fy, d_fz,
                           d_eta, d_etaXY, d_etaYZ, d_etaZX,
                           d_pi_dx, d_pi_dy, d_pi_dz, d_pi_nx, d_pi_ny, d_pi_nz,
                           d_fft, plan, gen, d_randD, d_randN,
                           tmp_fx, tmp_fy, tmp_fz);

        update_particle_velocity(cfg, pp, N, d_Rx, d_Ry, d_Rz,
                                 d_spx, d_spy, d_spz,
                                 d_vx, d_vy, d_vz,
                                 d_Vx, d_Vy, d_Vz);

        update_particle_position(cfg, N, d_Rx, d_Ry, d_Rz,
                                 d_Rux, d_Ruy, d_Ruz,
                                 d_Vx, d_Vy, d_Vz);

        if (step % 10000 == 0)
        {
            #pragma acc update host(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size])
            #pragma acc update host(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N])
            #pragma acc update host(d_Vx[0:N], d_Vy[0:N], d_Vz[0:N])
            save_fluid_vtk("result", (int)step, cfg, vx, vy, vz, p);
            save_particles_vtk("result", (int)step, N, Rx, Ry, Rz, Vx, Vy, Vz, Fx, Fy, Fz);

            // 发散早期预警：NaN/Inf 一出现立即停，不要跑满几小时才发现
            bool bad = false;
            for (int i = 0; i < size; i++)
            {
                if (!std::isfinite(vx[i]) || !std::isfinite(vy[i]) ||
                    !std::isfinite(vz[i]) || !std::isfinite(p[i])) { bad = true; break; }
            }
            if (!std::isfinite(Rx[0]) || !std::isfinite(Vx[0])) { bad = true; }

            const auto now = std::chrono::steady_clock::now();
            const double el = std::chrono::duration<double>(now - t_start).count();
            std::cout << "step " << std::setw(8) << step
                      << "   R = (" << std::fixed << std::setprecision(4)
                      << Rx[0] << ", " << Ry[0] << ", " << Rz[0] << ")"
                      << "   " << std::setprecision(1) << (step > 0 ? step / el : 0.0) << " 步/秒"
                      << (bad ? "   *** 检测到 NaN/Inf，中止 ***" : "")
                      << std::endl;
            if (bad) { return 2; }
        }
    }

    #pragma acc exit data delete(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size], d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc exit data delete(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Rux[0:N], d_Ruy[0:N], d_Ruz[0:N])
    #pragma acc exit data delete(d_Vx[0:N], d_Vy[0:N], d_Vz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])

    #pragma acc exit data delete(d_spx[0:N], d_spy[0:N], d_spz[0:N])
    #pragma acc exit data delete(d_eta[0:size], d_etaXY[0:size], d_etaYZ[0:size], d_etaZX[0:size])
    #pragma acc exit data delete(d_fft[0:size*2])

    #pragma acc exit data delete(tmp_fx[0:size], tmp_fy[0:size], tmp_fz[0:size])
    #pragma acc exit data delete(d_pi_dx[0:size], d_pi_dy[0:size], d_pi_dz[0:size], d_pi_nx[0:size], d_pi_ny[0:size], d_pi_nz[0:size], d_randD[0:size*3], d_randN[0:size*3])

    CUFFT_CHECK(cufftDestroy(plan));
    CURAND_CHECK(curandDestroyGenerator(gen));

    std::cout << "\nDone. Final vx[0]: " << vx[0] << std::endl;
    return 0;
}
