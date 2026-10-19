#include "PcgCalib.h"

#include <cstdio>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <openacc.h>

#include "Poisson.h"
#include "Velocity.h"

// ---------------------------------------------------------------------------
// 主机侧的独立实现（不调用 Stokes.cpp 内部）
// ---------------------------------------------------------------------------
namespace
{

struct HostFace
{
    std::vector<double> x, y, z;
    void assign(size_t n) { x.assign(n, 0.0); y.assign(n, 0.0); z.assign(n, 0.0); }
};

// out = div v / dt（壁面：k=0 的下邻居是恒 0 的下壁法向面）
void host_div(NS_Config cfg, const HostFace& v, std::vector<double>& out)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    out.assign((size_t)Nx * Ny * Nz, 0.0);
    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                const int ijk = IDX(i, j, k);
                const double vzb = (k > 0) ? v.z[IDX(i, j, k - 1)] : 0.0;
                out[ijk] = ((v.x[ijk] - v.x[IDX((i - 1 + Nx) % Nx, j, k)])
                          + (v.y[ijk] - v.y[IDX(i, (j - 1 + Ny) % Ny, k)])
                          + (v.z[ijk] - vzb)) / cfg.dt;
            }
}

// out = -dt·a·∇q（上壁法向面不修正，恒 0）
void host_grad_corr(NS_Config cfg, const HostFace& a, const std::vector<double>& q, HostFace& out)
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
double rho_energy(NS_Config cfg, const HostFace& a, const HostFace& v)
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

double norm2(const std::vector<double>& b)
{
    double s = 0.0;
    for (size_t t = 0; t < b.size(); t++) { s += b[t] * b[t]; }
    return std::sqrt(s);
}

void get_face(const double* x, const double* y, const double* z, size_t size, HostFace& v)
{
    #pragma acc update self(x[0:size], y[0:size], z[0:size])
    v.x.assign(x, x + size); v.y.assign(y, y + size); v.z.assign(z, z + size);
}

void put_velocity(FpdState& st, const HostFace& v)
{
    const size_t size = st.size;
    std::copy(v.x.begin(), v.x.end(), st.h_vx.begin());
    std::copy(v.y.begin(), v.y.end(), st.h_vy.begin());
    std::copy(v.z.begin(), v.z.end(), st.h_vz.begin());
    #pragma acc update device(st.vx[0:size], st.vy[0:size], st.vz[0:size])
}

// 生产 PCG（被测对象）解 L p = b；tol/max_iter 由 PcgTrace 覆盖，初值 p0=0 与生产相同。
int solve_to(FpdState& st, const std::vector<double>& b, double tol, int max_iter,
             std::vector<double>& p, double& res)
{
    const size_t size = st.size;
    std::copy(b.begin(), b.end(), st.h_cg_r.begin());
    #pragma acc update device(st.cg_r[0:size])
    PcgTrace tr;
    tr.tol = tol;
    tr.max_iter = max_iter;
    const int it = solve_pressure_var(st.cfg, st.cg_r, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                                      st.fft, st.tri_w, st.plan_xy, st.p, st.diag,
                                      st.cg_s, st.cg_zq, &tr);
    #pragma acc update self(st.p[0:size])
    p.assign(st.p, st.p + size);
    res = st.diag[1];      // solve_pressure_var 在主机侧写好
    return it;
}

// 速度场的粒子相场加权平均（生产 update_particle_velocity；sum_phi 已由 step_fields 算好）。
// 会覆盖设备上的 v 与 V，调用方负责恢复。返回每个粒子 |V_n|。
std::vector<double> particle_speed(const StepParams& sp, FpdState& st, const HostFace& v)
{
    put_velocity(st, v);
    const int N = st.N;
    update_particle_velocity(sp.cfg, sp.pp, N, st.Rx, st.Ry, st.Rz,
                             st.sum_phix, st.sum_phiy, st.sum_phiz,
                             st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);
    #pragma acc update self(st.Vx[0:N], st.Vy[0:N], st.Vz[0:N])
    std::vector<double> s(N);
    for (int n = 0; n < N; n++)
    {
        s[n] = std::sqrt(st.Vx[n] * st.Vx[n] + st.Vy[n] * st.Vy[n] + st.Vz[n] * st.Vz[n]);
    }
    return s;
}

const double REF_TOL = 1e-15;
const int    REF_MAX = 300;

CalibPoint calib_measure(const StepParams& sp, FpdState& st, long step,
                         const std::vector<double>& levels)
{
    const size_t size = st.size;
    const int N = st.N;
    const NS_Config cfg = sp.cfg;
    NS_Config cfg0 = cfg;
    cfg0.W = 0.0;

    CalibPoint pt;
    pt.step = step;

    // --- 1. 测量前状态 ---
    HostFace v0, a, P, tmp;
    get_face(st.vx, st.vy, st.vz, size, v0);
    #pragma acc update self(st.Vx[0:N], st.Vy[0:N], st.Vz[0:N])
    const std::vector<double> V0x(st.Vx, st.Vx + N), V0y(st.Vy, st.Vy + N), V0z(st.Vz, st.Vz + N);

    step_forces(sp, st);
    step_fields(sp, st);
    get_face(st.inv_rho_x, st.inv_rho_y, st.inv_rho_z, size, a);

    // --- 2. 两次显式步：P ← v*_A，tmp ← v*_B，P ← Δv* ---
    std::vector<double> bA, bN, pref, pk;
    step_fluid(cfg, st, step);
    get_face(st.tmp_fx, st.tmp_fy, st.tmp_fz, size, P);
    host_div(cfg, P, bA);
    put_velocity(st, v0);
    step_fluid(cfg0, st, step);
    get_face(st.tmp_fx, st.tmp_fy, st.tmp_fz, size, tmp);
    for (size_t t = 0; t < size; t++)
    {
        P.x[t] -= tmp.x[t]; P.y[t] -= tmp.y[t]; P.z[t] -= tmp.z[t];
    }
    host_div(cfg, P, bN);
    pt.noise_frac = norm2(bN) / norm2(bA);

    // --- 3. 噪声增量投影：P ← Δv* - dt·a·∇p_N ---
    double r = 0.0;
    solve_to(st, bN, REF_TOL, REF_MAX, pk, r);
    host_grad_corr(cfg, a, pk, tmp);
    for (size_t t = 0; t < size; t++)
    {
        P.x[t] += tmp.x[t]; P.y[t] += tmp.y[t]; P.z[t] += tmp.z[t];
    }
    const double E_noise = rho_energy(cfg, a, P);
    {
        const std::vector<double> s = particle_speed(sp, st, P);
        double s2 = 0.0;
        for (int n = 0; n < N; n++) { s2 += s[n] * s[n]; }
        pt.Vn_rms = std::sqrt(s2 / N);
    }
    std::vector<double>().swap(bN);
    P.assign(0); P.x.shrink_to_fit(); P.y.shrink_to_fit(); P.z.shrink_to_fit();

    // --- 4. 参考压力与各档误差 ---
    solve_to(st, bA, REF_TOL, REF_MAX, pref, pt.ref_res);
    for (double tau : levels)
    {
        CalibLevel L;
        L.tol = tau;
        L.iters = solve_to(st, bA, tau, PCG_MAX_ITER, pk, L.res);
        for (size_t t = 0; t < size; t++) { pk[t] -= pref[t]; }
        host_grad_corr(cfg, a, pk, tmp);
        L.eE = std::sqrt(rho_energy(cfg, a, tmp) / E_noise);
        const std::vector<double> s = particle_speed(sp, st, tmp);
        L.eV = *std::max_element(s.begin(), s.end()) / pt.Vn_rms;
        pt.lv.push_back(L);
    }

    // --- 5. 恢复测量前的 v、V（F、场、p 由调用方的正常一步重算）---
    put_velocity(st, v0);
    std::copy(V0x.begin(), V0x.end(), st.h_Vx.begin());
    std::copy(V0y.begin(), V0y.end(), st.h_Vy.begin());
    std::copy(V0z.begin(), V0z.end(), st.h_Vz.begin());
    #pragma acc update device(st.Vx[0:N], st.Vy[0:N], st.Vz[0:N])
    return pt;
}

} // namespace

std::vector<double> calib_levels()
{
    std::vector<double> lv;
    for (int e = 2; e <= 12; e++)
    {
        lv.push_back(std::pow(10.0, -e));
        if (e < 12) { lv.push_back(std::pow(10.0, -e) * 0.3); }
    }
    return lv;
}

int calib_run(const StepParams& sp, FpdState& st, long start_step,
              const CalibOpts& o, std::vector<CalibPoint>& out)
{
    const std::vector<double> levels = calib_levels();
    const long total = o.warm + o.gap * (long)(o.meas - 1) + 1;
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    out.clear();
    double last_print = -1e300;
    for (long i = 0; i < total; i++)
    {
        const long step = start_step + i;
        if (i >= o.warm && (i - o.warm) % o.gap == 0)
        {
            const double tm = elapsed();
            out.push_back(calib_measure(sp, st, step, levels));
            const CalibPoint& pt = out.back();
            if (!std::isfinite(pt.Vn_rms) || pt.Vn_rms <= 0.0)
            {
                std::fprintf(stderr, "错误: step %ld 的噪声增量无效（Vn_rms = %g）\n", step, pt.Vn_rms);
                return 1;
            }
            std::printf("  测量 %zu/%d  step %ld   噪声占比 %.4f   rms|ΔV_noise| = %.3e   "
                        "参考残差 %.1e   用时 %.1f s\n",
                        out.size(), o.meas, step, pt.noise_frac, pt.Vn_rms, pt.ref_res, elapsed() - tm);
            std::fflush(stdout);
        }
        // 正常的生产一步（轨迹用配置里的 pcg_tol）
        step_forces(sp, st);
        step_fields(sp, st);
        step_fluid(sp.cfg, st, step);
        step_particles(sp, st);
        // 进度：第一步之后与此后每 30 秒一行（大网格时可据此估计总耗时）
        if (i == 0 || elapsed() - last_print >= 30.0)
        {
            last_print = elapsed();
            std::printf("  step %ld / %ld   已用 %.1f s（%.2f 步/秒）\n",
                        step + 1, start_step + total, last_print, (i + 1) / last_print);
            std::fflush(stdout);
        }
    }
    return 0;
}

double calib_report(const std::vector<CalibPoint>& pts, const CalibOpts& o,
                    double cur_tol, double& tau_star)
{
    const size_t nl = pts.empty() ? 0 : pts[0].lv.size();
    std::vector<CalibLevel> agg(nl);
    double nf_min = 1e300, nf_max = 0.0;
    for (size_t l = 0; l < nl; l++)
    {
        agg[l] = pts[0].lv[l];
        for (const CalibPoint& p : pts)
        {
            agg[l].iters = std::max(agg[l].iters, p.lv[l].iters);
            agg[l].res = std::max(agg[l].res, p.lv[l].res);
            agg[l].eV = std::max(agg[l].eV, p.lv[l].eV);
            agg[l].eE = std::max(agg[l].eE, p.lv[l].eE);
        }
    }
    for (const CalibPoint& p : pts)
    {
        nf_min = std::min(nf_min, p.noise_frac);
        nf_max = std::max(nf_max, p.noise_frac);
    }

    // 从最严一档往松走，连续满足门槛的最后一档就是 τ*
    tau_star = -1.0;
    for (size_t l = nl; l-- > 0;)
    {
        if (agg[l].eV <= o.target && agg[l].eE <= o.target) { tau_star = agg[l].tol; }
        else { break; }
    }

    std::printf("\n  各档（%zu 个测量点的最大值；ε 门槛 %.0e）：\n", pts.size(), o.target);
    std::printf("  %9s %5s %11s %11s %11s %9s\n", "pcg_tol", "迭代", "实际残差", "ε_V", "ε_E", "ε_V/残差");
    double cv_min = 1e300, cv_max = 0.0;
    for (size_t l = 0; l < nl; l++)
    {
        const CalibLevel& L = agg[l];
        const double cv = L.res > 1e-14 ? L.eV / L.res : 0.0;
        if (L.res > 1e-14) { cv_min = std::min(cv_min, cv); cv_max = std::max(cv_max, cv); }
        const bool ok = L.eV <= o.target && L.eE <= o.target;
        std::printf("  %9.0e %5d %11.3e %11.3e %11.3e %9.2f  %s\n",
                    L.tol, L.iters, L.res, L.eV, L.eE, cv, ok ? "满足" : "不满足");
    }
    std::printf("\n  噪声在右端的占比 ||b_noise||/||b|| ∈ [%.4f, %.4f]\n", nf_min, nf_max);
    if (cv_max > 0.0) { std::printf("  ε_V / 实际残差 ∈ [%.2f, %.2f]\n", cv_min, cv_max); }

    if (tau_star <= 0.0)
    {
        std::printf("  连 1e-12 也不满足门槛：噪声在右端占比可能极低，请检查配置（noise_on、kT、重力）。\n");
        return -1.0;
    }

    // 推荐值：τ*/margin，向下取到候选档
    double rec = agg[nl - 1].tol;
    int rec_it = agg[nl - 1].iters;
    for (size_t l = 0; l < nl; l++)
    {
        if (agg[l].tol <= tau_star / o.margin * (1.0 + 1e-9)) { rec = agg[l].tol; rec_it = agg[l].iters; break; }
    }
    int it12 = agg[nl - 1].iters;
    std::printf("  满足门槛的最松档 τ* = %.0e；余量 %.0f 倍 → 推荐 %.0e（%d 次/步，1e-12 为 %d 次/步）\n",
                tau_star, o.margin, rec, rec_it, it12);
    std::printf("  当前配置 pcg_tol = %.0e：%s\n", cur_tol,
                cur_tol <= tau_star * (1.0 + 1e-9) ? "满足门槛" : "不满足门槛");
    return rec;
}
