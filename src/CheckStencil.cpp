#include <iostream>
#include <iomanip>
#include <cmath>

#include "./include/Tests.h"
#include "./include/Stencil.h"
#include "./include/Analysis.h"
#include "./include/Force.h"
#include "./include/Viscosity.h"
#include "./include/Velocity.h"
#include "./include/State.h"
#include "./include/Potential.h"

// ============================================================================
// Phase 1 自检：力守恒 + 亚格点不变性
//
// 判据：给定粒子受力 F，摊布到网格后必须满足 Σ_grid f_α == F_α，
// 且该等式与粒子在格胞内的分数位置无关。
//
// ⚠️ 这个判据检验的是 stencil_point 的【自洽性】，对它的内部公式错误是盲的
//    （投影和归一化会一起错）。独立的数值对照见 spike/。
// ============================================================================
int run_force_conservation_check(NS_Config cfg, PhiParams pp)
{
    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N    = 1;

    FpdState st;
    st.init(cfg, N, ST_STENCIL, 0);

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
        st.Rx[0] = px[t]; st.Ry[0] = py[t]; st.Rz[0] = pz[t];
        st.Fx[0] = 1.0;   st.Fy[0] = 2.0;   st.Fz[0] = 3.0;

        st.upload(ST_PARTICLE);

        update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                           st.fx, st.fy, st.fz);

        st.download(ST_PHI | ST_PARTICLE);

        double sx = 0.0, sy = 0.0, sz = 0.0;
        for (int i = 0; i < size; i++) { sx += st.fx[i]; sy += st.fy[i]; sz += st.fz[i]; }

        const double ex = std::fabs(sx - st.Fx[0]);
        const double ey = std::fabs(sy - st.Fy[0]);
        const double ez = std::fabs(sz - st.Fz[0]);
        const bool ok = (ex < TOL && ey < TOL && ez < TOL);
        if (!ok) { failures++; }

        std::cout << "  R = (" << std::fixed << std::setprecision(1)
                  << st.Rx[0] << ", " << st.Ry[0] << ", " << st.Rz[0] << ")"
                  << std::scientific << std::setprecision(3)
                  << "   |err| = " << ex << " / " << ey << " / " << ez
                  << "   " << (ok ? "PASS" : "FAIL") << "\n";

        // 各向同性 sanity check：三个方向的归一化因子应几乎相等
        const double smax = std::fmax(st.sum_phix[0], std::fmax(st.sum_phiy[0], st.sum_phiz[0]));
        const double smin = std::fmin(st.sum_phix[0], std::fmin(st.sum_phiy[0], st.sum_phiz[0]));
        const double spread = (smax - smin) / smin;
        const bool iso_ok = (spread < 0.01);
        if (!iso_ok) { failures++; }

        std::cout << "      sum_phi = " << std::fixed << std::setprecision(6)
                  << st.sum_phix[0] << " / " << st.sum_phiy[0] << " / " << st.sum_phiz[0]
                  << "   离散度 " << std::setprecision(3) << spread * 100.0 << "%   "
                  << (iso_ok ? "PASS" : "FAIL") << "\n";
    }

    st.finish();

    // --- 交错混叠体检 ---
    // 上面只检了 Σφ_α（那只是 k=0 分量）。由 Parseval，Σ_k|φ̂_α(k)|² = N_g·Σ_r φ_α²，
    // 所以 Σφ_α² 的三方向一致性才是混叠的直接判据 —— 这是曾经漏掉的那一半。
    {
        FpdConstants fc = compute_fpd_constants(pp, cfg, 64.3, 32.7, 16.5, 4);
        const bool a_ok = (fc.spread_int_phi2 < 1e-4);
        const bool s_ok = (fc.spread_subgrid  < 0.01);
        if (!a_ok) { failures++; }
        if (!s_ok) { failures++; }
        std::cout << std::scientific << std::setprecision(3)
                  << "  交错混叠  Σφ² 三方向散布 = " << fc.spread_int_phi2
                  << "   " << (a_ok ? "PASS" : "FAIL") << "\n"
                  << std::fixed << std::setprecision(4)
                  << "  系统误差下限  M_i 亚格点散布 = " << fc.spread_subgrid*100.0
                  << "%   " << (s_ok ? "PASS" : "FAIL") << "\n";
    }

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
int run_overlap_check(NS_Config cfg, PhiParams pp)
{
    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N    = 2;

    FpdState st;
    st.init(cfg, N, ST_KINEMATIC, 0);

    // 均匀流场
    const double V_UNIFORM = 1.0;
    for (int i = 0; i < size; i++) { st.vx[i] = V_UNIFORM; st.vy[i] = V_UNIFORM; st.vz[i] = V_UNIFORM; }
    st.upload(ST_VELOCITY);

    // 间距 4.0 < 2a = 6.4，两球显著重叠
    const double SEP[3] = {4.0, 6.4, 20.0};   // 重叠 / 刚接触 / 分离
    const char*  DESC[3] = {"重叠", "刚接触", "分离"};

    int failures = 0;
    std::cout << "=== Phase 1 自检：多粒子重叠（N=2）===\n";

    for (int t = 0; t < 3; t++)
    {
        st.Rx[0] = 64.0 - 0.5*SEP[t]; st.Ry[0] = 32.3; st.Rz[0] = 16.0;
        st.Rx[1] = 64.0 + 0.5*SEP[t]; st.Ry[1] = 32.3; st.Rz[1] = 16.0;
        st.Fx[0] = 1.0; st.Fy[0] = 2.0; st.Fz[0] = 3.0;
        st.Fx[1] = -0.5; st.Fy[1] = 0.25; st.Fz[1] = 1.5;

        st.upload(ST_PARTICLE);

        update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                           st.fx, st.fy, st.fz);
        update_particle_velocity(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                 st.sum_phix, st.sum_phiy, st.sum_phiz,
                                 st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);

        st.download(ST_PHI | ST_PARTICLE);

        // (a) 力守恒推广到多粒子
        double sx = 0.0, sy = 0.0, sz = 0.0;
        for (int i = 0; i < size; i++) { sx += st.fx[i]; sy += st.fy[i]; sz += st.fz[i]; }
        const double ex = std::fabs(sx - (st.Fx[0] + st.Fx[1]));
        const double ey = std::fabs(sy - (st.Fy[0] + st.Fy[1]));
        const double ez = std::fabs(sz - (st.Fz[0] + st.Fz[1]));
        const bool f_ok = (ex < 1e-12 && ey < 1e-12 && ez < 1e-12);

        // (b) 均匀流场下 V_i 必须等于 v，与重叠无关
        double vmax_err = 0.0;
        for (int n = 0; n < N; n++)
        {
            vmax_err = std::fmax(vmax_err, std::fabs(st.Vx[n] - V_UNIFORM));
            vmax_err = std::fmax(vmax_err, std::fabs(st.Vy[n] - V_UNIFORM));
            vmax_err = std::fmax(vmax_err, std::fabs(st.Vz[n] - V_UNIFORM));
        }
        const bool v_ok = (vmax_err < 1e-12);

        // 粘度峰值：重叠区 Σφ 可能 > 1，η 会超过 η_c，这是 FPD 的已知行为
        double eta_max = 0.0;
        for (int i = 0; i < size; i++) { eta_max = std::fmax(eta_max, st.eta[i]); }

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

    st.finish();

    std::cout << (failures == 0 ? "全部通过\n" : "有失败项\n");
    return failures;
}

// ============================================================================
// J5 端到端力链路 + GPU 力确定性
//
// 测的是「compute_particle_forces 的 GPU F 正确且确定 + 被 update_force_field
// 正确投影」。三条判据：
//   (a) GPU 全矩阵 vs CPU 半矩阵（两份独立实现）对照 —— 力的正确性
//   (b) GPU 版同一构型算两次逐位 —— reduction 确定性（J8 的核心目标）
//   (c) 力投影守恒 Σf == ΣF（bg=0），以及 bg 补偿后 Σf == 0（发现 1 的补偿逻辑）
// ============================================================================
int run_force_pipeline_check(NS_Config cfg, PhiParams pp)
{
    const int size = cfg.Nx * cfg.Ny * cfg.Nz;
    const int N = 4;

    PotentialParams pot = make_potential_params(POT_LJ126, SHIFT_ENERGY,
                                                57.1428571429, 7.4, 0, 0, 0, 15.0);
    ExternalField ext{1.0, -2.0, 3.0};

    FpdState st;
    st.init(cfg, N, ST_STENCIL, 0);

    // 构型：4 个粒子，间距 < rcut=15（有力）
    st.Rx[0]=64; st.Ry[0]=32; st.Rz[0]=16;
    st.Rx[1]=74; st.Ry[1]=32; st.Rz[1]=16;    // 与 0 间距 10
    st.Rx[2]=64; st.Ry[2]=42; st.Rz[2]=16;    // 与 0 间距 10
    st.Rx[3]=69; st.Ry[3]=37; st.Rz[3]=16;    // 与 0 间距 √50
    st.upload(ST_PARTICLE);

    int failures = 0;
    std::cout << "=== J5 端到端力链路 + GPU 力确定性 ===\n";
    std::cout << std::scientific << std::setprecision(3);

    // (a)+(b) GPU 算两次逐位 + 与 CPU 参考对照
    compute_particle_forces(cfg, pot, ext, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz);
    st.download(ST_PARTICLE);
    double g0x[N], g0y[N], g0z[N];
    for (int n = 0; n < N; n++) { g0x[n] = st.Fx[n]; g0y[n] = st.Fy[n]; g0z[n] = st.Fz[n]; }

    compute_particle_forces(cfg, pot, ext, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz);
    st.download(ST_PARTICLE);
    bool bitwise = true;
    for (int n = 0; n < N && bitwise; n++)
    {
        if (st.Fx[n] != g0x[n] || st.Fy[n] != g0y[n] || st.Fz[n] != g0z[n]) { bitwise = false; }
    }

    double cFx[N], cFy[N], cFz[N];
    compute_particle_forces_cpu(cfg, pot, ext, N, st.Rx, st.Ry, st.Rz, cFx, cFy, cFz);
    double fmax = 0.0, relmax = 0.0;
    for (int n = 0; n < N; n++)
    {
        fmax = std::fmax(fmax, std::fabs(st.Fx[n]));
        fmax = std::fmax(fmax, std::fabs(st.Fy[n]));
        fmax = std::fmax(fmax, std::fabs(st.Fz[n]));
        relmax = std::fmax(relmax, std::fabs(st.Fx[n] - cFx[n]));
        relmax = std::fmax(relmax, std::fabs(st.Fy[n] - cFy[n]));
        relmax = std::fmax(relmax, std::fabs(st.Fz[n] - cFz[n]));
    }
    relmax = fmax > 0.0 ? relmax / fmax : 0.0;

    const bool a_ok = relmax < 1e-12;
    if (!a_ok) { failures++; }
    if (!bitwise) { failures++; }
    std::cout << "  GPU vs CPU 对照   max 相对差 = " << relmax
              << "   " << (a_ok ? "PASS" : "FAIL") << "\n";
    std::cout << "  GPU 算两次逐位   " << (bitwise ? "PASS" : "FAIL") << "\n";

    // (c) 力投影守恒：Σf == ΣF（bg=0）
    update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                            st.sum_phix, st.sum_phiy, st.sum_phiz,
                            st.eta, st.etaXY, st.etaYZ, st.etaZX);
    update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                       st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                       st.fx, st.fy, st.fz);
    st.download(ST_PHI | ST_PARTICLE);
    double sx = 0, sy = 0, sz = 0, sumFx = 0, sumFy = 0, sumFz = 0;
    for (int i = 0; i < size; i++) { sx += st.fx[i]; sy += st.fy[i]; sz += st.fz[i]; }
    for (int n = 0; n < N; n++) { sumFx += st.Fx[n]; sumFy += st.Fy[n]; sumFz += st.Fz[n]; }
    const double e1 = std::fmax(std::fabs(sx - sumFx), std::fmax(std::fabs(sy - sumFy), std::fabs(sz - sumFz)));
    const bool c0_ok = e1 < 1e-12;
    if (!c0_ok) { failures++; }
    std::cout << "  力守恒 Σf == ΣF  (bg=0)   |err| = " << e1
              << "   " << (c0_ok ? "PASS" : "FAIL") << "\n";

    // (c') bg 补偿后 Σf == 0
    const double bgx = -sumFx / (double)size;
    const double bgy = -sumFy / (double)size;
    const double bgz = -sumFz / (double)size;
    update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                       st.sum_phix, st.sum_phiy, st.sum_phiz, bgx, bgy, bgz,
                       st.fx, st.fy, st.fz);
    st.download(ST_PHI);
    sx = 0; sy = 0; sz = 0;
    for (int i = 0; i < size; i++) { sx += st.fx[i]; sy += st.fy[i]; sz += st.fz[i]; }
    const double e2 = std::fmax(std::fabs(sx), std::fmax(std::fabs(sy), std::fabs(sz)));
    const bool cc_ok = e2 < 1e-12;
    if (!cc_ok) { failures++; }
    std::cout << "  力守恒 Σf == 0  (bg 补偿)  |err| = " << e2
              << "   " << (cc_ok ? "PASS" : "FAIL") << "\n";

    st.finish();
    std::cout << (failures == 0 ? "全部通过\n" : "有失败项\n");
    return failures;
}
