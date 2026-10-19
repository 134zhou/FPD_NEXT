#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <openacc.h>

#include "Poisson.h"
#include "State.h"
#include "Stokes.h"
#include "Viscosity.h"
#include "Force.h"
#include "Velocity.h"
#include "Mass.h"
#include "Potential.h"
#include "Tests.h"

// ============================================================================
// --pcg-noise（需 GPU）：变系数压力求解的误差与【同一步注入的热噪声】比大小。
//
// 投影后的速度误差是 δv = -dt·a·∇e（e = p_k - p_exact）。它的 ρ 加权能量
//   Σ_f ρ_f |δv_f|² = dt² eᵀ A e
// 拿来和本步噪声增量投影后的 ρ 加权能量 Σ_f ρ_f |P Δv*_f|² 比，得
//   ε_E(k) = sqrt( Σρ|δv|² / Σρ|PΔv*|² )
// 粒子层面：δV_i、ΔV_i 分别是 δv、PΔv* 的相场加权平均（生产 update_particle_velocity），
//   ε_V(k) = max_i |δV_i| / rms_i |ΔV_i|
//
// 噪声增量的取法：同一状态、同一 RNG slot 跑两次 step_navier_stokes，W 正常与 W=0。
// v* 对 W 是线性的，两次的 tmp_f*（= v*）之差就是本步噪声增量 Δv*；投影用参考解。
//
// 独立实现约定：右端 b、修正步、能量都在主机侧手写，只调用生产的 solve_pressure_var
// 作为被测对象（第 k 次迭代 = PcgTrace{tol=0, max_iter=k}，p0=0 固定所以可重放）。
// 只打印标量。
//
// 同一测量点另评估压力分裂法（cfg.split=1）：压力历史 pⁿ、pⁿ⁻¹ 取自参考 PCG 轨迹的前两步
// （最有利的历史），跑一步生产 step_navier_stokes(split_now=1)，误差 δv = v_split - v_exact。
// ============================================================================

static int g_fail = 0;
static void check(bool ok, const char* what)
{
    std::printf("  %-60s %s\n", what, ok ? "PASS" : "FAIL");
    if (!ok) { g_fail++; }
}

struct HostFace
{
    std::vector<double> x, y, z;
    void assign(size_t n) { x.assign(n, 0.0); y.assign(n, 0.0); z.assign(n, 0.0); }
};

// b = (1/dt)·div v*（壁面：k=0 的下邻居是恒 0 的下壁法向面）
static void host_rhs(NS_Config cfg, const HostFace& v, std::vector<double>& b)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    b.assign((size_t)Nx * Ny * Nz, 0.0);
    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                const int ijk = IDX(i, j, k);
                const int im = (i - 1 + Nx) % Nx;
                const int jm = (j - 1 + Ny) % Ny;
                const double vzb = (k > 0) ? v.z[IDX(i, j, k - 1)] : 0.0;
                b[ijk] = ((v.x[ijk] - v.x[IDX(im, j, k)])
                        + (v.y[ijk] - v.y[IDX(i, jm, k)])
                        + (v.z[ijk] - vzb)) / cfg.dt;
            }
}

// out = -dt·a·∇q（上壁法向面不修正，恒 0）
static void host_grad_corr(NS_Config cfg, const HostFace& a, const std::vector<double>& q,
                           HostFace& out)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    out.assign((size_t)Nx * Ny * Nz);
    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                const int ijk = IDX(i, j, k);
                out.x[ijk] = -cfg.dt * a.x[ijk] * (q[IDX((i + 1) % Nx, j, k)] - q[ijk]);
                out.y[ijk] = -cfg.dt * a.y[ijk] * (q[IDX(i, (j + 1) % Ny, k)] - q[ijk]);
                out.z[ijk] = (k < Nz - 1) ? -cfg.dt * a.z[ijk] * (q[IDX(i, j, k + 1)] - q[ijk]) : 0.0;
            }
}

// Σ_f ρ_f v_f²，ρ_f = 1/a_f；上壁法向面不是自由度
static double rho_energy(NS_Config cfg, const HostFace& a, const HostFace& v)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    double s = 0.0;
    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                const int ijk = IDX(i, j, k);
                s += v.x[ijk] * v.x[ijk] / a.x[ijk] + v.y[ijk] * v.y[ijk] / a.y[ijk];
                if (k < Nz - 1) { s += v.z[ijk] * v.z[ijk] / a.z[ijk]; }
            }
    return s;
}

// 用生产 solve_pressure_var 解 L p = b，固定迭代 max_iter 次（tol=0）。res 返回末次相对残差。
static void solve_k(FpdState& st, const std::vector<double>& b, int max_iter,
                    std::vector<double>& p, double& res)
{
    const size_t size = st.size;
    for (size_t t = 0; t < size; t++) { st.h_cg_r[t] = b[t]; }
    #pragma acc update device(st.cg_r[0:size])

    std::vector<double> hist(max_iter + 1, 0.0);
    PcgTrace tr;
    tr.res = hist.data();
    tr.cap = max_iter + 1;
    tr.tol = 0.0;
    tr.max_iter = max_iter;
    solve_pressure_var(st.cfg, st.cg_r, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                       st.fft, st.tri_w, st.plan_xy, st.p, st.diag,
                       st.cg_s, st.cg_zq, &tr);
    #pragma acc update self(st.p[0:size])
    p.assign(st.p, st.p + size);
    res = (tr.n > 0) ? hist[tr.n - 1] : 0.0;
}

// 速度场的粒子相场加权平均（生产 update_particle_velocity；sum_phi 已由调用方算好）
static void particle_avg(FpdState& st, PhiParams pp, const HostFace& v, std::vector<double>& V)
{
    const size_t size = st.size;
    for (size_t t = 0; t < size; t++) { st.h_vx[t] = v.x[t]; st.h_vy[t] = v.y[t]; st.h_vz[t] = v.z[t]; }
    #pragma acc update device(st.vx[0:size], st.vy[0:size], st.vz[0:size])
    update_particle_velocity(st.cfg, pp, st.N, st.Rx, st.Ry, st.Rz,
                             st.sum_phix, st.sum_phiy, st.sum_phiz,
                             st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);
    const int N = st.N;
    #pragma acc update self(st.Vx[0:N], st.Vy[0:N], st.Vz[0:N])
    V.assign(3 * (size_t)N, 0.0);
    for (int n = 0; n < N; n++) { V[3 * n] = st.Vx[n]; V[3 * n + 1] = st.Vy[n]; V[3 * n + 2] = st.Vz[n]; }
}

// 生产主循环一步之前的场准备：力 → 粘度 → 密度 → 力密度 → 重力
static void prep_fields(FpdState& st, PhiParams pp, PotentialParams pot, WallParams wall,
                        ExternalField g, double d_rho)
{
    const NS_Config cfg = st.cfg;
    const int N = st.N;
    compute_particle_forces(cfg, pot, wall, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz);
    update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                            st.sum_phix, st.sum_phiy, st.sum_phiz,
                            st.eta, st.etaXY, st.etaYZ, st.etaZX);
    update_mass_fields(cfg, pp, d_rho, N, st.Rx, st.Ry, st.Rz,
                       st.inv_rho_x, st.inv_rho_y, st.inv_rho_z);
    update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                       st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                       st.fx, st.fy, st.fz);
    apply_gravity(cfg, g, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z, st.fx, st.fy, st.fz);
}

static void ns_step(FpdState& st, NS_Config cfg, long step, int split_now = 0)
{
    step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p, st.p_prev, split_now,
                       st.fx, st.fy, st.fz,
                       st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                       st.cg_r, st.cg_s, st.cg_zq,
                       st.eta, st.etaXY, st.etaYZ, st.etaZX,
                       st.pi_dx, st.pi_dy, st.pi_dz,
                       st.pi_nx, st.pi_ny, st.pi_nz,
                       st.fft, st.plan_xy, st.tri_w, st.diag,
                       st.gen, st.randD, st.randN,
                       st.tmp_fx, st.tmp_fy, st.tmp_fz, step);
}

static void get_vstar(FpdState& st, HostFace& v)
{
    const size_t size = st.size;
    #pragma acc update self(st.tmp_fx[0:size], st.tmp_fy[0:size], st.tmp_fz[0:size])
    v.x.assign(st.tmp_fx, st.tmp_fx + size);
    v.y.assign(st.tmp_fy, st.tmp_fy + size);
    v.z.assign(st.tmp_fz, st.tmp_fz + size);
}

static void put_velocity(FpdState& st, const HostFace& v)
{
    const size_t size = st.size;
    for (size_t t = 0; t < size; t++) { st.h_vx[t] = v.x[t]; st.h_vy[t] = v.y[t]; st.h_vz[t] = v.z[t]; }
    #pragma acc update device(st.vx[0:size], st.vy[0:size], st.vz[0:size])
}

static double vec_norm3(const std::vector<double>& V, int n)
{
    return std::sqrt(V[3 * n] * V[3 * n] + V[3 * n + 1] * V[3 * n + 1] + V[3 * n + 2] * V[3 * n + 2]);
}

int run_pcg_noise(int Nx, int Ny, int Nz)
{
    g_fail = 0;
    if (Nx <= 0 || Ny <= 0 || Nz <= 0) { Nx = Ny = Nz = 32; }

    // sed.cfg 的物理参数：ratio_rho=2，WCA 粒子间势与壁面势，重力按体力
    const double dt = 0.002, kT = 0.25, ratio_rho = 2.0;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, kT, /*noise_on=*/true);
    cfg.density_on = 1;
    NS_Config cfg0 = cfg;
    cfg0.W = 0.0;
    NS_Config cfgS = cfg;           // 分裂法只在测量点用；预热轨迹始终是 PCG
    cfgS.split = 1;
    cfgS.split_rho0 = 1.0;          // min(1, ratio_rho)
    const PhiParams pp = make_phi_params(3.2, 1.0, 50.0);
    const PotentialParams pot = make_potential_params(POT_WCA, SHIFT_ENERGY, 10.0, 5.7018, 0, 0, 0, 0);
    WallParams wall;
    wall.pot = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0);
    wall.radius = pp.radius;
    const ExternalField g = {0.0, 0.0, -0.058718};
    const double d_rho = ratio_rho - 1.0;

    // 4 个粒子，离两壁 >= range；前两个支撑域与 WCA 截断都重叠
    const int N = 4;
    const double px[N] = {0.38 * Nx, 0.38 * Nx + 6.2, 0.70 * Nx, 0.25 * Nx};
    const double py[N] = {0.50 * Ny, 0.50 * Ny + 0.6, 0.26 * Ny, 0.77 * Ny};
    const double pz[N] = {0.44 * Nz, 0.44 * Nz + 0.3, 0.64 * Nz, 0.56 * Nz};

    FpdState st;
    st.init(cfgS, N, ST_FULL, 20261019ULL);     // cfgS：顺带分配 p_prev
    for (int n = 0; n < N; n++)
    {
        st.h_Rx[n] = px[n]; st.h_Ry[n] = py[n]; st.h_Rz[n] = pz[n];
        st.h_Rux[n] = px[n]; st.h_Ruy[n] = py[n]; st.h_Ruz[n] = pz[n];
    }
    st.upload(ST_PARTICLE);

    const size_t size = st.size;
    const double dim = (double)Nx * Ny * (2.0 * Nz - 1.0) + 1.0;
    const int KMAX = 30, KREF = 300;
    const long n_warm = 1000, n_gap = 200;
    const int n_meas = 4;

    std::printf("=== PCG 误差 vs 热噪声（需 GPU）===\n");
    std::printf("  网格 %d x %d x %d   N=%d   ratio_rho=%.1f   kT=%.2f   dt=%.3f   W=%.3f\n",
                Nx, Ny, Nz, N, ratio_rho, kT, dt, cfg.W);
    std::printf("  预热 %ld 步后每隔 %ld 步测一次，共 %d 次；参考解 %d 次迭代\n",
                n_warm, n_gap, n_meas, KREF);

    // 每个 k：各次测量的最大值
    std::vector<double> res_k(KMAX + 1, 0.0), eE_k(KMAX + 1, 0.0), eV_k(KMAX + 1, 0.0),
                        eT_k(KMAX + 1, 0.0);
    double bfrac_min = 1e300, bfrac_max = 0.0;
    double sE_max = 0.0, sV_max = 0.0, sdiv_max = 0.0;
    std::vector<double> p_nm1(size), p_n(size);

    long step = 0;
    for (int m = 0; m < n_meas; m++)
    {
        const long target = n_warm + (long)m * n_gap;
        for (; step < target; step++)
        {
            prep_fields(st, pp, pot, wall, g, d_rho);
            ns_step(st, cfg, step);
            update_particle_velocity(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                     st.sum_phix, st.sum_phiy, st.sum_phiz,
                                     st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);
            update_particle_position(cfg, N, st.Rx, st.Ry, st.Rz,
                                     st.Rux, st.Ruy, st.Ruz, st.Vx, st.Vy, st.Vz);
            // 测量点前两步的压力 = 分裂法所需的历史 pⁿ⁻¹、pⁿ
            if (step == target - 2 || step == target - 1)
            {
                #pragma acc update self(st.p[0:size])
                std::vector<double>& dst = (step == target - 2) ? p_nm1 : p_n;
                dst.assign(st.p, st.p + size);
            }
        }

        // --- 本步：场准备、两次显式步（W 正常 / W=0）---
        HostFace v0, vsA, vsB, a;
        #pragma acc update self(st.vx[0:size], st.vy[0:size], st.vz[0:size])
        v0.x.assign(st.vx, st.vx + size); v0.y.assign(st.vy, st.vy + size); v0.z.assign(st.vz, st.vz + size);

        prep_fields(st, pp, pot, wall, g, d_rho);
        #pragma acc update self(st.inv_rho_x[0:size], st.inv_rho_y[0:size], st.inv_rho_z[0:size])
        a.x.assign(st.inv_rho_x, st.inv_rho_x + size);
        a.y.assign(st.inv_rho_y, st.inv_rho_y + size);
        a.z.assign(st.inv_rho_z, st.inv_rho_z + size);

        ns_step(st, cfg, step);  get_vstar(st, vsA);
        put_velocity(st, v0);
        ns_step(st, cfg0, step); get_vstar(st, vsB);

        std::vector<double> bA, bB, bN(size);
        host_rhs(cfg, vsA, bA);
        host_rhs(cfg, vsB, bB);
        HostFace dvs;
        dvs.assign(size);
        for (size_t t = 0; t < size; t++)
        {
            bN[t] = bA[t] - bB[t];
            dvs.x[t] = vsA.x[t] - vsB.x[t];
            dvs.y[t] = vsA.y[t] - vsB.y[t];
            dvs.z[t] = vsA.z[t] - vsB.z[t];
        }
        double nbA = 0.0, nbN = 0.0;
        for (size_t t = 0; t < size; t++) { nbA += bA[t] * bA[t]; nbN += bN[t] * bN[t]; }
        const double bfrac = std::sqrt(nbN / nbA);
        bfrac_min = std::min(bfrac_min, bfrac);
        bfrac_max = std::max(bfrac_max, bfrac);

        // 噪声增量投影：PΔv* = Δv* - dt a ∇p_N
        std::vector<double> pref, pN;
        double rref = 0.0, rN = 0.0;
        solve_k(st, bA, KREF, pref, rref);
        solve_k(st, bN, KREF, pN, rN);
        HostFace corrN, Pdv;
        host_grad_corr(cfg, a, pN, corrN);
        Pdv.assign(size);
        for (size_t t = 0; t < size; t++)
        {
            Pdv.x[t] = dvs.x[t] + corrN.x[t];
            Pdv.y[t] = dvs.y[t] + corrN.y[t];
            Pdv.z[t] = dvs.z[t] + corrN.z[t];
        }
        const double E_noise = rho_energy(cfg, a, Pdv);
        std::vector<double> Vn;
        particle_avg(st, pp, Pdv, Vn);
        double vn2 = 0.0;
        for (int n = 0; n < N; n++) { vn2 += vec_norm3(Vn, n) * vec_norm3(Vn, n); }
        const double Vn_rms = std::sqrt(vn2 / N);

        std::printf("  测量 %d（step %ld）：|b_noise|/|b| = %.3f   Σρ|PΔv*|² = %.4e"
                    "（kT·dim = %.4e）   rms|ΔV_noise| = %.3e\n",
                    m + 1, step, bfrac, E_noise, kT * dim, Vn_rms);

        for (int k = 1; k <= KMAX; k++)
        {
            std::vector<double> pk, e(size);
            double rk = 0.0;
            solve_k(st, bA, k, pk, rk);
            for (size_t t = 0; t < size; t++) { e[t] = pk[t] - pref[t]; }
            HostFace dv;
            host_grad_corr(cfg, a, e, dv);
            const double Ee = rho_energy(cfg, a, dv);
            std::vector<double> dV;
            particle_avg(st, pp, dv, dV);
            double dVmax = 0.0;
            for (int n = 0; n < N; n++) { dVmax = std::max(dVmax, vec_norm3(dV, n)); }

            res_k[k] = std::max(res_k[k], rk);
            eE_k[k]  = std::max(eE_k[k], std::sqrt(Ee / E_noise));
            eV_k[k]  = std::max(eV_k[k], dVmax / Vn_rms);
            eT_k[k]  = std::max(eT_k[k], std::sqrt(Ee / (kT * dim)));
        }

        // --- 压力分裂法一步（W 正常，与 vsA 同一个 v*）---
        {
            for (size_t t = 0; t < size; t++) { st.h_p[t] = p_n[t]; st.h_p_prev[t] = p_nm1[t]; }
            #pragma acc update device(st.p[0:size], st.p_prev[0:size])
            put_velocity(st, v0);
            ns_step(st, cfgS, step, 1);
            #pragma acc update self(st.vx[0:size], st.vy[0:size], st.vz[0:size], st.diag[0:3])
            HostFace corr, dv;
            host_grad_corr(cfg, a, pref, corr);
            dv.assign(size);
            for (size_t t = 0; t < size; t++)
            {
                dv.x[t] = st.vx[t] - (vsA.x[t] + corr.x[t]);
                dv.y[t] = st.vy[t] - (vsA.y[t] + corr.y[t]);
                dv.z[t] = st.vz[t] - (vsA.z[t] + corr.z[t]);
            }
            const double sdiv = st.diag[2];
            const double Ee = rho_energy(cfg, a, dv);
            std::vector<double> dV;
            particle_avg(st, pp, dv, dV);
            double dVmax = 0.0;
            for (int n = 0; n < N; n++) { dVmax = std::max(dVmax, vec_norm3(dV, n)); }
            sE_max = std::max(sE_max, std::sqrt(Ee / E_noise));
            sV_max = std::max(sV_max, dVmax / Vn_rms);
            sdiv_max = std::max(sdiv_max, sdiv);
        }

        // 恢复本步之前的速度，交回生产轨迹（下一轮 prep+NS 从 v0 重新走这一步）
        put_velocity(st, v0);
    }

    std::printf("\n  逐迭代（%d 次测量的最大值）：\n", n_meas);
    std::printf("  %4s  %12s  %12s  %12s  %12s\n", "k", "||r||/||b||", "ε_E", "ε_V", "ε_热");
    for (int k = 1; k <= KMAX; k++)
    {
        std::printf("  %4d  %12.3e  %12.3e  %12.3e  %12.3e\n", k, res_k[k], eE_k[k], eV_k[k], eT_k[k]);
    }

    // 按容差档汇总：首个满足 res <= tol 的 k（各次测量取最大残差，故偏保守）
    std::printf("\n  按容差档（ε 门槛 1e-3）：\n");
    std::printf("  %8s  %4s  %12s  %12s  %s\n", "pcg_tol", "k", "ε_E", "ε_V", "");
    double loosest = -1.0;
    static const double tols[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12};
    for (double tol : tols)
    {
        int kk = 0;
        for (int k = 1; k <= KMAX; k++) { if (res_k[k] <= tol) { kk = k; break; } }
        if (kk == 0) { std::printf("  %8.0e  %4s\n", tol, "--"); continue; }
        const bool ok = eE_k[kk] <= 1e-3 && eV_k[kk] <= 1e-3;
        if (ok && loosest < 0.0) { loosest = tol; }
        std::printf("  %8.0e  %4d  %12.3e  %12.3e  %s\n", tol, kk, eE_k[kk], eV_k[kk],
                    ok ? "满足" : "不满足");
    }
    std::printf("  |b_noise|/|b| ∈ [%.3f, %.3f]\n", bfrac_min, bfrac_max);
    if (loosest > 0.0) { std::printf("  满足 ε<=1e-3 的最松容差：%.0e\n", loosest); }

    std::printf("\n  压力分裂法（历史取自参考轨迹，%d 次测量的最大值）：\n", n_meas);
    std::printf("    ε_E = %.3e   ε_V = %.3e   真残差 ||div v||/||div v*|| = %.3e   %s\n",
                sE_max, sV_max, sdiv_max,
                (sE_max <= 1e-3 && sV_max <= 1e-3) ? "满足 ε<=1e-3" : "不满足 ε<=1e-3");

    // 硬判据：生产默认容差下误差远小于噪声
    int k12 = 0;
    for (int k = 1; k <= KMAX; k++) { if (res_k[k] <= 1e-12) { k12 = k; break; } }
    check(k12 > 0 && eE_k[k12] <= 1e-3 && eV_k[k12] <= 1e-3,
          "PN1 默认 pcg_tol=1e-12：ε_E、ε_V <= 1e-3");

    st.finish();
    std::printf("\n%s\n", g_fail == 0 ? "完成" : "有失败项");
    return g_fail;
}
