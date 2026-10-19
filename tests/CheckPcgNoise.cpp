#include <cstdio>
#include <cmath>
#include <vector>

#include "State.h"
#include "Step.h"
#include "PcgCalib.h"
#include "Tests.h"

// ============================================================================
// --pcg-noise（需 GPU）：在内置合成构型上跑 fpd_calib 的同一核心（PcgCalib.cpp），
// 并给出硬判据 PN1：默认 pcg_tol = 1e-12 时 ε_V、ε_E <= 1e-3。
//
// 构型：32³、N=4（前两个支撑域与 WCA 截断都重叠）、sed.cfg 的物理参数
// （ratio_rho=2、ratio_eta=50、kT=0.25、dt=0.002、WCA 粒子间势与壁面势、g_z=-0.058718）。
// 生产配置的标定用 fpd_calib，这里只是常驻回归。
// ============================================================================
int run_pcg_noise(int Nx, int Ny, int Nz)
{
    if (Nx <= 0 || Ny <= 0 || Nz <= 0) { Nx = Ny = Nz = 32; }

    StepParams sp;
    sp.cfg = make_ns_config(Nx, Ny, Nz, 0.002, 0.25, /*noise_on=*/true);
    sp.cfg.density_on = 1;
    sp.pp  = make_phi_params(3.2, 1.0, 50.0);
    sp.pot = make_potential_params(POT_WCA, SHIFT_ENERGY, 10.0, 5.7018, 0, 0, 0, 0);
    sp.wall.pot = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0);
    sp.wall.radius = sp.pp.radius;
    sp.ext = {0.0, 0.0, -0.058718};
    sp.d_rho = 1.0;

    const int N = 4;
    const double px[N] = {0.38 * Nx, 0.38 * Nx + 6.2, 0.70 * Nx, 0.25 * Nx};
    const double py[N] = {0.50 * Ny, 0.50 * Ny + 0.6, 0.26 * Ny, 0.77 * Ny};
    const double pz[N] = {0.44 * Nz, 0.44 * Nz + 0.3, 0.64 * Nz, 0.56 * Nz};

    FpdState st;
    st.init(sp.cfg, N, ST_FULL, 20261019ULL);
    for (int n = 0; n < N; n++)
    {
        st.h_Rx[n] = px[n]; st.h_Ry[n] = py[n]; st.h_Rz[n] = pz[n];
        st.h_Rux[n] = px[n]; st.h_Ruy[n] = py[n]; st.h_Ruz[n] = pz[n];
    }
    st.upload(ST_PARTICLE);

    CalibOpts o;
    o.warm = 1000;
    o.gap  = 200;
    o.meas = 4;

    std::printf("=== PCG 误差 vs 热噪声（需 GPU）===\n");
    std::printf("  网格 %d x %d x %d   N=%d   ratio_rho=2   kT=0.25   dt=0.002   W=%.3f\n",
                Nx, Ny, Nz, N, sp.cfg.W);

    std::vector<CalibPoint> pts;
    int fails = calib_run(sp, st, 0, o, pts);
    st.finish();
    if (fails != 0) { return fails; }

    double tau_star = 0.0;
    calib_report(pts, o, 1e-12, tau_star);

    const bool ok = tau_star > 0.0;      // τ* 存在 ⟺ 1e-12 及其之前各档都满足门槛
    std::printf("  %-60s %s\n", "PN1 默认 pcg_tol=1e-12：ε_V、ε_E <= 1e-3", ok ? "PASS" : "FAIL");
    std::printf("\n%s\n", ok ? "完成" : "有失败项");
    return ok ? 0 : 1;
}
