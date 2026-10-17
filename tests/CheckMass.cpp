#include <cstdio>
#include <cmath>
#include <vector>

#include "Common.h"
#include "Mass.h"
#include "Potential.h"
#include "Stencil.h"
#include "Viscosity.h"
#include "State.h"
#include "Tests.h"
#include "Check.h"

// ============================================================================
// 质量密度场判据（--check-mass，需 GPU）
//
// 判据独立性的做法与 --check 的 W4/J3 一致：参考实现【不复用】stencil_point 的
// 模板盒，而是对全网格逐点、逐粒子直接算距离与权重。它和被测实现共享的只有
// 相场函数 order() 与截断半径 range2 —— 那两条是物理定义，不是实现假设。
//
//   MA1 三个面心的 a=1/rho 与暴力枚举逐点对照
//   MA2 sum_faces(1/a-1) = d_rho * sum_n S1（三方向，望远镜恒等式）
//   MA3 纯流体区 a=1.0 逐位；全网格最小 a 与 1/ratio 同量级
//   MA4 重力体力总和 = d_rho*g*sum_n S1；流体区逐位为 0
//   MA5 N=0：a 恒 1、重力体力恒 0
//   MA6 d_rho=0：a 恒 1.0（逐位）
//   MA7 壁面平衡：F_wall(z_rest) + d_rho*S1z(z_rest)*g_z = 0，
//       且 GPU 实际施加的重力总和就是 d_rho*S1z*g_z
// ============================================================================

static int g_fail = 0;
static void check(bool ok, const char* what)
{
    std::printf("  %-60s %s\n", what, ok ? "PASS" : "FAIL");
    if (!ok) { g_fail++; }
}

static inline double wrap_xy(double d, double L)
{
    return d - L * std::round(d / L);
}

// 独立参考：暴力枚举全网格。loc: 0=FACE_X, 1=FACE_Y, 2=FACE_Z。
// 与 stencil_point 的 z 截断一致：FACE_Z 的 k=Nz-1 不参与（保持 1/rho=1）。
static void brute_field(int loc, NS_Config cfg, PhiParams pp, double d_rho, int N,
                        const double* Rx, const double* Ry, const double* Rz,
                        std::vector<double>& a_out)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const double ox = (loc == 0) ? 0.5 : 0.0;
    const double oy = (loc == 1) ? 0.5 : 0.0;
    const double oz = (loc == 2) ? 0.5 : 0.0;
    a_out.assign((size_t)Nx * Ny * Nz, 1.0);

    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                if (loc == 2 && k == Nz - 1) { continue; }   // 上壁法向面
                const double x = (double)i + ox;
                const double y = (double)j + oy;
                const double z = (double)k + oz;
                double s = 0.0;
                for (int n = 0; n < N; n++)
                {
                    const double dx = wrap_xy(x - Rx[n], (double)Nx);
                    const double dy = wrap_xy(y - Ry[n], (double)Ny);
                    const double dz = z - Rz[n];
                    const double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 <= pp.range2)
                    {
                        s += order(dx, dy, dz, pp.radius, pp.inv_xi);
                    }
                }
                a_out[(size_t)IDX(i, j, k)] = 1.0 / (1.0 + d_rho * s);
            }
}

// 独立参考：一个粒子在高度 z 处、FACE_Z 面心的 sum(phi)（壁面 z 截断一致）。
static double brute_s1z(NS_Config cfg, PhiParams pp, double Rx, double Ry, double z)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    double s = 0.0;
    for (int k = 0; k < Nz - 1; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                const double dx = wrap_xy((double)i - Rx, (double)Nx);
                const double dy = wrap_xy((double)j - Ry, (double)Ny);
                const double dz = (double)k + 0.5 - z;
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 <= pp.range2)
                {
                    s += order(dx, dy, dz, pp.radius, pp.inv_xi);
                }
            }
    return s;
}

// 三个面心场对网格求和（GPU 结果已在主机）
static double sum_field(const std::vector<double>& v)
{
    double s = 0.0;
    for (size_t t = 0; t < v.size(); t++) { s += v[t]; }
    return s;
}

int run_check_mass(int Nx, int Ny, int Nz)
{
    g_fail = 0;
    std::printf("=== 质量密度场判据（需 GPU）===\n");
    if (Nx <= 0) { Nx = 32; Ny = 32; Nz = 32; }

    const double dt = 0.002, kT = 0.0;
    const bool noise_on = false;
    const double radius = 3.2, xi = 1.0, ratio_eta = 50.0, ratio_rho = 2.0;
    const double d_rho = ratio_rho - 1.0;

    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, kT, noise_on);
    cfg.density_on = 1;
    PhiParams pp = make_phi_params(radius, xi, ratio_eta);

    const int N = 3;
    const double RX[3] = {10.3, 20.2, 15.6};
    const double RY[3] = { 9.7, 21.4,  6.1};
    const double RZ[3] = { 9.5, 15.5, 21.0};

    FpdState st;
    st.init(cfg, N, ST_FULL, 20261017ULL);
    for (int n = 0; n < N; n++) { st.Rx[n] = RX[n]; st.Ry[n] = RY[n]; st.Rz[n] = RZ[n]; }
    st.upload(ST_PARTICLE);

    update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                            st.sum_phix, st.sum_phiy, st.sum_phiz,
                            st.eta, st.etaXY, st.etaYZ, st.etaZX);
    update_mass_fields(cfg, pp, d_rho, N, st.Rx, st.Ry, st.Rz,
                       st.inv_rho_x, st.inv_rho_y, st.inv_rho_z);
    st.download(ST_PHI | ST_PARTICLE);

    const size_t size = st.size;

    // --- MA1：逐点暴力枚举对照 ---
    {
        const double* gpu[3] = { st.inv_rho_x, st.inv_rho_y, st.inv_rho_z };
        double worst[3] = {0.0, 0.0, 0.0};
        for (int loc = 0; loc < 3; loc++)
        {
            std::vector<double> ref;
            brute_field(loc, cfg, pp, d_rho, N, st.Rx, st.Ry, st.Rz, ref);
            double mx = 0.0;
            for (size_t t = 0; t < size; t++)
            {
                const double d = std::fabs(gpu[loc][t] - ref[t]);
                if (d > mx) { mx = d; }
            }
            worst[loc] = mx;
        }
        check(worst[0] < 1e-12 && worst[1] < 1e-12 && worst[2] < 1e-12,
              "MA1 面心 1/rho 与独立暴力枚举逐点一致（三方向）");
        std::printf("      max|a_gpu - a_ref| = %.3e / %.3e / %.3e\n",
                    worst[0], worst[1], worst[2]);
    }

    // --- MA2：望远镜恒等式 ---
    {
        const double* gpu[3] = { st.inv_rho_x, st.inv_rho_y, st.inv_rho_z };
        const double* s1[3]  = { st.sum_phix, st.sum_phiy, st.sum_phiz };
        bool ok = true;
        for (int loc = 0; loc < 3; loc++)
        {
            double lhs = 0.0;
            for (size_t t = 0; t < size; t++) { lhs += (1.0 / gpu[loc][t] - 1.0); }
            double rhs = 0.0;
            for (int n = 0; n < N; n++) { rhs += s1[loc][n]; }
            rhs *= d_rho;
            const double rel = (std::fabs(rhs) > 0.0) ? std::fabs(lhs - rhs) / std::fabs(rhs) : std::fabs(lhs);
            std::printf("      %s: sum(1/a-1) = %.6e   d_rho*sum(S1) = %.6e   相对差 %.3e\n",
                        loc == 0 ? "x" : (loc == 1 ? "y" : "z"), lhs, rhs, rel);
            if (rel > 1e-13) { ok = false; }
        }
        check(ok, "MA2 sum_faces(1/a-1) = d_rho * sum_n S1（三方向，相对差 < 1e-13）");
    }

    // --- MA3：纯流体区 a=1.0 逐位；最小 a 与 1/ratio 同量级 ---
    {
        // 角落离三个粒子都 > range：三方向都应为 1.0
        const int c = IDX(0, 0, 0);
        bool exact = (st.inv_rho_x[c] == 1.0 && st.inv_rho_y[c] == 1.0 && st.inv_rho_z[c] == 1.0);
        double amin = st.inv_rho_x[0];
        for (size_t t = 1; t < size; t++)
        {
            if (st.inv_rho_x[t] < amin) { amin = st.inv_rho_x[t]; }
            if (st.inv_rho_y[t] < amin) { amin = st.inv_rho_y[t]; }
            if (st.inv_rho_z[t] < amin) { amin = st.inv_rho_z[t]; }
        }
        std::printf("      纯流体角点 a = (%.17g, %.17g, %.17g)   min a = %.6f   1/ratio = %.6f\n",
                    st.inv_rho_x[c], st.inv_rho_y[c], st.inv_rho_z[c], amin, 1.0 / ratio_rho);
        check(exact, "MA3 纯流体区 a 逐位为 1.0");
        check(std::fabs(amin - 1.0 / ratio_rho) < 1e-2,
              "MA3 全网格最小 a 与 1/ratio 同量级（< 1e-2）");
    }

    // --- MA4：重力体力 ---
    {
        const ExternalField g{0.3, -0.7, 0.2};
        // fx/fy/fz 在 init 时是 0；update_viscosity_fields 不碰它们。
        apply_gravity(cfg, g, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                      st.fx, st.fy, st.fz);
        st.download(ST_PHI);

        const double* f[3] = { st.fx, st.fy, st.fz };
        const double gc[3] = { g.gx, g.gy, g.gz };
        const double* s1[3] = { st.sum_phix, st.sum_phiy, st.sum_phiz };
        bool ok = true;
        for (int loc = 0; loc < 3; loc++)
        {
            double lhs = sum_field(std::vector<double>(f[loc], f[loc] + size));
            double rhs = 0.0;
            for (int n = 0; n < N; n++) { rhs += s1[loc][n]; }
            rhs *= d_rho * gc[loc];
            const double rel = (std::fabs(rhs) > 0.0) ? std::fabs(lhs - rhs) / std::fabs(rhs) : std::fabs(lhs);
            std::printf("      %s: sum(f_g) = %.6e   期望 %.6e   相对差 %.3e\n",
                        loc == 0 ? "x" : (loc == 1 ? "y" : "z"), lhs, rhs, rel);
            if (rel > 1e-13) { ok = false; }
        }
        check(ok, "MA4 sum(f_g) = d_rho * g * sum_n S1（三方向）");
        // 远离粒子的角落：三个面心都应为 0（逐位）
        const int c = IDX(0, 0, 0);
        check(st.fx[c] == 0.0 && st.fy[c] == 0.0 && st.fz[c] == 0.0,
              "MA4 流体区重力体力逐位为 0");
    }

    // --- MA6：d_rho=0 时 a 恒 1.0（逐位）---
    {
        update_mass_fields(cfg, pp, 0.0, N, st.Rx, st.Ry, st.Rz,
                           st.inv_rho_x, st.inv_rho_y, st.inv_rho_z);
        st.download(ST_PHI);
        bool exact = true;
        for (size_t t = 0; t < size && exact; t++)
        {
            if (st.inv_rho_x[t] != 1.0 || st.inv_rho_y[t] != 1.0 || st.inv_rho_z[t] != 1.0)
            { exact = false; }
        }
        check(exact, "MA6 d_rho=0（等密度）时三个面心 a 逐位为 1.0");
    }

    st.finish();

    // --- MA5：N=0 ---
    {
        FpdState s0;
        s0.init(cfg, 0, ST_FULL, 7ULL);
        update_mass_fields(cfg, pp, d_rho, 0, s0.Rx, s0.Ry, s0.Rz,
                           s0.inv_rho_x, s0.inv_rho_y, s0.inv_rho_z);
        const ExternalField g{0.3, -0.7, 0.2};
        apply_gravity(cfg, g, s0.inv_rho_x, s0.inv_rho_y, s0.inv_rho_z,
                      s0.fx, s0.fy, s0.fz);
        s0.download(ST_PHI);
        bool ok = true;
        for (size_t t = 0; t < s0.size; t++)
        {
            if (s0.inv_rho_x[t] != 1.0 || s0.inv_rho_y[t] != 1.0 || s0.inv_rho_z[t] != 1.0) { ok = false; }
            if (s0.fx[t] != 0.0 || s0.fy[t] != 0.0 || s0.fz[t] != 0.0) { ok = false; }
        }
        check(ok, "MA5 N=0：a 恒 1、重力体力恒 0（逐位）");
        s0.finish();
    }

    // --- MA7：壁面平衡 ---
    {
        WallParams w;
        w.pot    = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0);
        w.radius = radius;

        // 等价净力约 -10 的重力加速度：g = F / (d_rho * S1)，S1 取体相 170.306
        const double gz = -10.0 / (d_rho * 170.3061);
        const double Rx = 16.5, Ry = 16.5;

        // 独立二分：F_total(z) = F_wall(z) + d_rho*S1z(z)*gz
        auto F_total = [&](double z) {
            return wall_force_z(w, cfg, z) + d_rho * brute_s1z(cfg, pp, Rx, Ry, z) * gz;
        };
        double lo = -0.5 + radius + 1e-3;
        double hi = lo + 8.0;
        for (int it = 0; it < 80; it++)
        {
            const double mid = 0.5 * (lo + hi);
            if (F_total(mid) > 0.0) { lo = mid; } else { hi = mid; }
        }
        const double z_rest = 0.5 * (lo + hi);
        const double F_wall = wall_force_z(w, cfg, z_rest);
        const double F_g = d_rho * brute_s1z(cfg, pp, Rx, Ry, z_rest) * gz;
        std::printf("      z_rest = %.6f   F_wall = %.6e   F_g = %.6e   F_total = %.3e\n",
                    z_rest, F_wall, F_g, F_total(z_rest));
        check(std::fabs(F_wall + F_g) < 1e-10 * std::fabs(F_g),
              "MA7 独立二分出的 z_rest 上 壁面力 + 净重力 = 0");

        FpdState s1;
        s1.init(cfg, 1, ST_FULL, 11ULL);
        s1.Rx[0] = Rx; s1.Ry[0] = Ry; s1.Rz[0] = z_rest;
        s1.upload(ST_PARTICLE);
        update_mass_fields(cfg, pp, d_rho, 1, s1.Rx, s1.Ry, s1.Rz,
                           s1.inv_rho_x, s1.inv_rho_y, s1.inv_rho_z);
        const ExternalField g{0.0, 0.0, gz};
        apply_gravity(cfg, g, s1.inv_rho_x, s1.inv_rho_y, s1.inv_rho_z,
                      s1.fx, s1.fy, s1.fz);
        s1.download(ST_PHI);
        double fz = 0.0;
        for (size_t t = 0; t < s1.size; t++) { fz += s1.fz[t]; }
        std::printf("      GPU 实际施加的 sum(fz) = %.6e   d_rho*S1z*gz = %.6e\n", fz, F_g);
        const double rel = std::fabs(fz - F_g) / std::fabs(F_g);
        check(rel < 1e-12, "MA7 GPU 实际重力总和 = d_rho*S1z(z_rest)*gz");
        check(std::fabs(fz + F_wall) < 1e-10 * std::fabs(F_g),
              "MA7 GPU 重力 + 壁面力 = 0（沉降的力平衡）");
        s1.finish();
    }

    std::printf("\n%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}
