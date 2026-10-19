#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <openacc.h>

#include "Poisson.h"
#include "State.h"
#include "Stokes.h"
#include "Viscosity.h"
#include "Mass.h"
#include "Tests.h"
#include "Check.h"

// ============================================================================
// 压力泊松求解器自检。两个入口：
//   --check-tridiag  纯 CPU，无卡可跑：三对角 Thomas + build_tridiag_coeffs
//                    对照稠密 Gauss，含 (0,0) 奇异列的定规。
//   --check-poisson  需 GPU：算子往返判据（手写 div∘grad 作用出右端，求解器
//                    还原 p，两者只差一个常数）。
//
// 独立实现约定（与 Potential.h 的 GPU/CPU 双实现同一条规矩）：这里的稠密 Gauss、
// 手写 div∘grad 参考算子、主机 Thomas 都是【独立手写】，不调用生产求解器内部，
// 否则「判据」退化成自洽的、盲的。
// ============================================================================

static int g_fail = 0;
static void check(bool ok, const char* what)
{
    std::printf("  %-56s %s\n", what, ok ? "PASS" : "FAIL");
    if (!ok) { g_fail++; }
}

// 确定性伪随机（LCG，避免 <random> 跨平台差异）
static void fill_rand(std::vector<double>& v, unsigned seed)
{
    unsigned s = seed;
    for (size_t t = 0; t < v.size(); t++)
    {
        s = s * 1664525u + 1013904223u;
        v[t] = ((double)(s & 0xFFFFFF) / 16777216.0) * 2.0 - 1.0;
    }
}

// ============================================================================
// 独立参考 1：稠密 Gauss 消元（带部分选主元），解三对角 A(λ_xy) x = b。
// gauged：是否为 (0,0) 奇异列（d_0 = λ_xy - 2 而非 λ_xy - 1）。
// ============================================================================
static void dense_gauss(int Nz, double lam, bool gauged, const std::vector<double>& b,
                        std::vector<double>& x)
{
    std::vector<double> A((size_t)Nz * Nz, 0.0);
    for (int k = 0; k < Nz; k++)
    {
        double dk = (k == 0 || k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
        if (gauged && k == 0) { dk -= 1.0; }
        A[(size_t)k * Nz + k] = dk;
        if (k > 0)       { A[(size_t)k * Nz + (k-1)] = 1.0; }
        if (k < Nz - 1) { A[(size_t)k * Nz + (k+1)] = 1.0; }
    }
    std::vector<double> M = A;
    x = b;

    for (int col = 0; col < Nz; col++)
    {
        int piv = col;
        double best = std::fabs(M[(size_t)col * Nz + col]);
        for (int r = col + 1; r < Nz; r++)
        {
            const double v = std::fabs(M[(size_t)r * Nz + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (piv != col)
        {
            for (int c = 0; c < Nz; c++) { std::swap(M[(size_t)col*Nz+c], M[(size_t)piv*Nz+c]); }
            std::swap(x[col], x[piv]);
        }
        const double pv = M[(size_t)col * Nz + col];
        for (int r = col + 1; r < Nz; r++)
        {
            const double f = M[(size_t)r * Nz + col] / pv;
            for (int c = col; c < Nz; c++) { M[(size_t)r * Nz + c] -= f * M[(size_t)col * Nz + c]; }
            x[r] -= f * x[col];
        }
    }
    for (int r = Nz - 1; r >= 0; r--)
    {
        for (int c = r + 1; c < Nz; c++) { x[r] -= M[(size_t)r * Nz + c] * x[c]; }
        x[r] /= M[(size_t)r * Nz + r];
    }
}

// ============================================================================
// 独立参考 2：主机 Thomas，用 build_tridiag_coeffs 产出的 tri_w。
// w 指向第 (i,j) 列首（tri_w + i + j*Nx），stride = Nx*Ny。
// 与设备 thomas_sweep 同算法、独立代码路径（真正的独立参考是上面的 dense_gauss）。
// ============================================================================
static void host_thomas(const double* w, int Nz, int stride, std::vector<double>& b)
{
    b[0] *= w[0];
    for (int k = 1; k < Nz; k++)
    {
        b[k] = (b[k] - b[k-1]) * w[(size_t)k * stride];
    }
    for (int k = Nz - 2; k >= 0; k--)
    {
        b[k] = b[k] - w[(size_t)k * stride] * b[k+1];
    }
}

// 未定规 Neumann 算子 A0（λ_xy=0）的残差 |A0 x - b|：验证定规解满足原方程。
static double residual_neumann(int Nz, const std::vector<double>& x, const std::vector<double>& b)
{
    double mx = 0.0;
    for (int k = 0; k < Nz; k++)
    {
        const double dk = (k == 0 || k == Nz - 1) ? -1.0 : -2.0;
        double lhs = dk * x[k];
        if (k > 0)       { lhs += x[k-1]; }
        if (k < Nz - 1) { lhs += x[k+1]; }
        const double r = std::fabs(lhs - b[k]);
        if (r > mx) { mx = r; }
    }
    return mx;
}

// ============================================================================
// --check-tridiag（纯 CPU）：build_tridiag_coeffs + 主机 Thomas vs 稠密 Gauss。
// 覆盖 Nz ∈ {2,6,32,64}；对每个 (i,j) 列逐一对照；对 (0,0) 奇异列额外验
// 定规（x_0 ≈ 0 且满足未定规 Neumann 方程）。
// ============================================================================
int run_check_tridiag()
{
    g_fail = 0;
    std::printf("=== 泊松三对角自检（纯 CPU，无卡可跑）===\n");

    const int grids[][3] = { {4, 3, 2}, {4, 3, 6}, {4, 4, 32}, {2, 2, 64} };

    for (int g = 0; g < 4; g++)
    {
        const int Nx = grids[g][0], Ny = grids[g][1], Nz = grids[g][2];
        const int stride = Nx * Ny;
        NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false);

        std::vector<double> tri_w((size_t)Nx * Ny * Nz);
        build_tridiag_coeffs(cfg, tri_w.data());

        double worst_rel = 0.0;

        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                const double lam = 2.0 * (cos(2.0 * M_PI * i / Nx) - 1.0)
                                 + 2.0 * (cos(2.0 * M_PI * j / Ny) - 1.0);
                const bool singular = (i == 0 && j == 0);

                std::vector<double> b(Nz), b2;
                fill_rand(b, 9000u + (unsigned)(i * 97 + j * 13 + Nz));
                if (singular)
                {
                    // 相容投影：Σb = 0（精确算术下自动成立，浮点下显式减均值）
                    double mean = 0.0;
                    for (int k = 0; k < Nz; k++) { mean += b[k]; }
                    mean /= (double)Nz;
                    for (int k = 0; k < Nz; k++) { b[k] -= mean; }
                }
                b2 = b;

                // 稠密 Gauss（独立参考）
                std::vector<double> xg;
                dense_gauss(Nz, lam, singular, b, xg);

                // 主机 Thomas（用生产 build_tridiag_coeffs 的 tri_w）
                host_thomas(tri_w.data() + i + j * Nx, Nz, stride, b2);

                double maxmag = 0.0, maxd = 0.0;
                for (int k = 0; k < Nz; k++)
                {
                    const double m = std::fabs(xg[k]);
                    if (m > maxmag) { maxmag = m; }
                    const double d = std::fabs(b2[k] - xg[k]);
                    if (d > maxd) { maxd = d; }
                }
                const double rel = maxmag > 0.0 ? maxd / maxmag : maxd;
                if (rel > worst_rel) { worst_rel = rel; }
            }
        }

        char lbl[64];
        std::snprintf(lbl, sizeof(lbl), "三对角 Nz=%d：Thomas vs 稠密 Gauss（全 %d 列）", Nz, Nx * Ny);
        check(worst_rel < 1e-12, lbl);
        std::printf("      Nz=%2d  max 相对差 = %.3e\n", Nz, worst_rel);
    }

    // 奇异列定规的汇总断言（取所有 grid 中最坏）
    // 在循环里已累计 worst_x0 / worst_res，这里用最后一轮的值打印一次
    {
        // 重新跑一次 (0,0) 列的独立检查，得到可报告的数值
        const int Nx = 4, Ny = 4, Nz = 32, stride = Nx * Ny;
        NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false);
        std::vector<double> tri_w((size_t)Nx * Ny * Nz);
        build_tridiag_coeffs(cfg, tri_w.data());
        std::vector<double> b(Nz);
        fill_rand(b, 4242u);
        double mean = 0.0;
        for (int k = 0; k < Nz; k++) { mean += b[k]; }
        mean /= Nz;
        for (int k = 0; k < Nz; k++) { b[k] -= mean; }
        std::vector<double> b0 = b;
        std::vector<double> x = b;
        host_thomas(tri_w.data(), Nz, stride, x);
        const double res = residual_neumann(Nz, x, b0);
        const double x0 = std::fabs(x[0]);
        check(x0 < 1e-13 && res < 1e-13, "奇异列 (0,0) 定规：x_0 ≈ 0 且满足原 Neumann 方程");
        std::printf("      |x_0| = %.3e   |A0 x - b| = %.3e\n", x0, res);
    }

    std::printf("\n%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}

// ============================================================================
// 独立参考 3：手写离散算子 div∘grad（壁面：z 向 Neumann ghost）。
//   km=(k==0?0:k-1), kp=(k==Nz-1?Nz-1:k+1)     （Neumann ghost p[-1]=p[0] 等）
//   p[ip]+p[im]+p[jp]+p[jm]+p[kp]+p[km] - 6 p[ijk]
//
// ⚠️ 这是【手写的独立参考】，不许改成调用 Poisson.cpp 的代码 —— 那会把往返判据
//    变成自洽的、盲的（AGENTS.md：独立参考不得复用待测实现）。
// ============================================================================
static void apply_operator(int Nx, int Ny, int Nz, const double* p, double* b)
{
    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                const int ijk = IDX(i, j, k);
                const int ip = (i + 1) % Nx, im = (i - 1 + Nx) % Nx;
                const int jp = (j + 1) % Ny, jm = (j - 1 + Ny) % Ny;
                const int kp = (k == Nz - 1) ? Nz - 1 : k + 1;
                const int km = (k == 0)      ? 0      : k - 1;
                b[ijk] = p[IDX(ip, j, k)] + p[IDX(im, j, k)]
                       + p[IDX(i, jp, k)] + p[IDX(i, jm, k)]
                       + p[IDX(i, j, kp)] + p[IDX(i, j, km)]
                       - 6.0 * p[ijk];
            }
        }
    }
}

// 单个盒子的算子往返：随机 p_ref → 手写算子 → solve_pressure → 还原（只差常数）。
// 返回失败数。
static int roundtrip_grid(int Nx, int Ny, int Nz, const char* tag)
{
    const int size = Nx * Ny * Nz;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 12345ULL);

    // 随机参考压力 + 手写算子作用出右端
    std::vector<double> p_ref((size_t)size);
    fill_rand(p_ref, (unsigned)(Nx * 1000 + Ny * 31 + Nz));
    std::vector<double> b((size_t)size);
    apply_operator(Nx, Ny, Nz, p_ref.data(), b.data());

    // 填 fft_data（复数交错，虚部 0）并上传
    for (int t = 0; t < size; t++) { st.h_fft[2 * t] = b[t]; st.h_fft[2 * t + 1] = 0.0; }
    #pragma acc update device(st.fft[0:2*size])

    // (0,0) 列的相容性诊断数组（设备侧长度 1）
    double diag_h[1] = { 0.0 };
    double* diag = diag_h;
    #pragma acc enter data create(diag[0:1])

    solve_pressure(cfg, st.fft, st.tri_w, st.plan_xy, st.p, diag);

    #pragma acc update self(st.p[0:size])
    #pragma acc update self(diag[0:1])
    const double compat = diag_h[0];
    #pragma acc exit data delete(diag[0:1])

    st.finish();

    // d = p - p_ref 应是一个常数（null space = 常向量），去均值后应 ~0
    double mean = 0.0, pmax = 0.0;
    for (int t = 0; t < size; t++) { mean += st.p[t] - p_ref[t]; }
    mean /= (double)size;
    double mx = 0.0;
    for (int t = 0; t < size; t++)
    {
        const double d = (st.p[t] - p_ref[t]) - mean;
        const double a = std::fabs(d);
        if (a > mx) { mx = a; }
        const double pm = std::fabs(p_ref[t]);
        if (pm > pmax) { pmax = pm; }
    }
    const double tol = 1e-10 * (pmax > 0.0 ? pmax : 1.0);
    const bool ok = mx < tol;

    std::printf("  %-56s %s\n", tag, ok ? "PASS" : "FAIL");
    std::printf("      %dx%dx%d   max|(p-p_ref)-mean| = %.3e   阈值 %.3e\n",
                Nx, Ny, Nz, mx, tol);
    std::printf("      相容性残差 |Σb̂| = %.3e   （应 ~机器精度，非 O(1)）\n", compat);
    if (compat > 1e-9) { g_fail++; }   // 相容投影没生效的硬信号
    return ok ? 0 : 1;
}

// ============================================================================
// P6：相容性残差诊断对【不相容】右端应报非零。
// 喂一个任意随机 b（不来自 L p_ref，Σb ≠ 0），壁面路径的 diag 应报出 |Σ_{ijk} b|
// 的量级（(0,0) 模式 = Σ_{i,j} b，再沿 k 求和 = Σ_{ijk} b）。据此验证
// fix_singular_column 确实在算、且诊断有意义，而不是恒 0。
// ============================================================================
static void check_compat_diag()
{
    const int Nx = 6, Ny = 4, Nz = 5;
    const int size = Nx * Ny * Nz;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 999ULL);

    std::vector<double> b((size_t)size);
    fill_rand(b, 31415u);
    double host_sum = 0.0;
    for (int t = 0; t < size; t++) { host_sum += b[t]; }

    for (int t = 0; t < size; t++) { st.h_fft[2 * t] = b[t]; st.h_fft[2 * t + 1] = 0.0; }
    #pragma acc update device(st.fft[0:2*size])

    double diag_h[1] = { 0.0 };
    double* diag = diag_h;
    #pragma acc enter data create(diag[0:1])

    solve_pressure(cfg, st.fft, st.tri_w, st.plan_xy, st.p, diag);

    #pragma acc update self(diag[0:1])
    const double compat = diag_h[0];
    #pragma acc exit data delete(diag[0:1])
    st.finish();

    const bool ok = (std::fabs(host_sum) > 1e-3) && (std::fabs(compat - std::fabs(host_sum)) < 1e-9 * (1.0 + std::fabs(host_sum)));
    check(ok, "P6 相容性残差诊断（不相容右端应报非零且匹配 |Σb|）");
    std::printf("      |Σ_{ijk} b|（主机） = %.3e    diag = %.3e\n", std::fabs(host_sum), compat);
}

// ============================================================================
// 变系数投影判据 V0-V5（需 GPU）
//
// a 场用解析函数（高斯或球），不依赖 Mass.cpp —— 这样「算子/PCG」与「密度场」
// 是两套独立判据，各自能被单独判死。
// ============================================================================
static void build_var_a(std::vector<double>& a, int loc, NS_Config cfg,
                        double amp, double sigma, bool piecewise, double R)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const double cx = 0.5 * Nx, cy = 0.5 * Ny, cz = 0.5 * Nz;
    const double ox = (loc == 0) ? 0.5 : 0.0;
    const double oy = (loc == 1) ? 0.5 : 0.0;
    const double oz = (loc == 2) ? 0.5 : 0.0;
    a.assign((size_t)Nx * Ny * Nz, 1.0);
    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                const double dx = (double)i + ox - cx;
                const double dy = (double)j + oy - cy;
                const double dz = (double)k + oz - cz;
                const double r2 = dx * dx + dy * dy + dz * dz;
                const double phi = piecewise ? ((r2 <= R * R) ? 1.0 : 0.0)
                                             : std::exp(-r2 / (2.0 * sigma * sigma));
                a[(size_t)IDX(i, j, k)] = 1.0 / (1.0 + amp * phi);
            }
}

static double dot_host(const std::vector<double>& x, const std::vector<double>& y)
{
    double s = 0.0;
    for (size_t t = 0; t < x.size(); t++) { s += x[t] * y[t]; }
    return s;
}

// ============================================================================
// V0：对流加速度【不】随 1/rho 缩放。
//
// 合成均匀 a、eta=0、f=0、p=0 的构型，跑一步显式步，比较 v*（tmp_f）。
// 均匀密度下非守恒方程是 ∂t v + (v·∇)v = 0，所以 v* 必须与 a 无关。
// 修正前把整个括号乘 a，v* 会差一个因子 a；这条判据就是抓它的。
// ============================================================================
static void v0_run_step(FpdState& st, double a,
                        const std::vector<double>& vx0, const std::vector<double>& vy0,
                        const std::vector<double>& vz0,
                        std::vector<double>& ox, std::vector<double>& oy, std::vector<double>& oz)
{
    const size_t size = st.size;
    for (size_t t = 0; t < size; t++)
    {
        st.h_vx[t] = vx0[t]; st.h_vy[t] = vy0[t]; st.h_vz[t] = vz0[t];
        st.h_fx[t] = 0.0; st.h_fy[t] = 0.0; st.h_fz[t] = 0.0;
        st.h_inv_rho_x[t] = a; st.h_inv_rho_y[t] = a; st.h_inv_rho_z[t] = a;
    }
    // eta=0：只留下对流项（噪声 W=0 已在 make_ns_config 里关掉）
    std::fill(st.h_eta.begin(), st.h_eta.end(), 0.0);
    std::fill(st.h_etaXY.begin(), st.h_etaXY.end(), 0.0);
    std::fill(st.h_etaYZ.begin(), st.h_etaYZ.end(), 0.0);
    std::fill(st.h_etaZX.begin(), st.h_etaZX.end(), 0.0);
    st.upload(ST_VELOCITY | ST_PHI);

    step_navier_stokes(st.cfg, st.vx, st.vy, st.vz, st.p,
                       st.fx, st.fy, st.fz,
                       st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                       st.cg_r, st.cg_s, st.cg_zq,
                       st.eta, st.etaXY, st.etaYZ, st.etaZX,
                       st.pi_dx, st.pi_dy, st.pi_dz,
                       st.pi_nx, st.pi_ny, st.pi_nz,
                       st.fft, st.plan_xy, st.tri_w, st.diag,
                       st.gen, st.randD, st.randN,
                       st.tmp_fx, st.tmp_fy, st.tmp_fz, 0L);

    #pragma acc update self(st.tmp_fx[0:size])
    #pragma acc update self(st.tmp_fy[0:size])
    #pragma acc update self(st.tmp_fz[0:size])
    ox.assign(st.tmp_fx, st.tmp_fx + size);
    oy.assign(st.tmp_fy, st.tmp_fy + size);
    oz.assign(st.tmp_fz, st.tmp_fz + size);
}

static void v0_convection_independent_of_a(int Nx, int Ny, int Nz)
{
    const size_t size = (size_t)Nx * Ny * Nz;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.002, 0.0, false);
    cfg.density_on = 1;

    FpdState st;
    st.init(cfg, 0, ST_FULL, 20261018ULL);

    std::vector<double> vx0(size), vy0(size), vz0(size);
    fill_rand(vx0, 777u); fill_rand(vy0, 888u); fill_rand(vz0, 999u);
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) { vz0[(size_t)IDX(i, j, Nz - 1)] = 0.0; }

    std::vector<double> ax, ay, az, bx, by, bz;
    v0_run_step(st, 0.5, vx0, vy0, vz0, ax, ay, az);
    v0_run_step(st, 1.0, vx0, vy0, vz0, bx, by, bz);

    double dmax = 0.0, vmax = 0.0;
    for (size_t t = 0; t < size; t++)
    {
        dmax = std::fmax(dmax, std::fabs(ax[t] - bx[t]));
        dmax = std::fmax(dmax, std::fabs(ay[t] - by[t]));
        dmax = std::fmax(dmax, std::fabs(az[t] - bz[t]));
        vmax = std::fmax(vmax, std::fmax(std::fabs(bx[t]),
                             std::fmax(std::fabs(by[t]), std::fabs(bz[t]))));
    }
    const double rel = vmax > 0.0 ? dmax / vmax : dmax;
    std::printf("      V0 均匀 a=0.5 vs a=1.0：max|Δv*|/max|v*| = %.3e\n", rel);
    check(rel == 0.0, "V0 对流加速度与均匀 a 无关（v* 逐位相同）");
    st.finish();
}

// 用给定的 a 场跑一遍「往返 + 对称 + 修正步散度」。
static void var_poisson_battery(int Nx, int Ny, int Nz, const char* tag,
                                double amp, double sigma, bool piecewise, double R)
{
    const size_t size = (size_t)Nx * Ny * Nz;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.002, 0.0, false);
    cfg.density_on = 1;
    PhiParams pp = make_phi_params(3.2, 1.0, 50.0);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 20261017ULL);

    std::vector<double> ax, ay, az;
    build_var_a(ax, 0, cfg, amp, sigma, piecewise, R);
    build_var_a(ay, 1, cfg, amp, sigma, piecewise, R);
    build_var_a(az, 2, cfg, amp, sigma, piecewise, R);
    for (size_t t = 0; t < size; t++)
    { st.h_inv_rho_x[t] = ax[t]; st.h_inv_rho_y[t] = ay[t]; st.h_inv_rho_z[t] = az[t]; }
    #pragma acc update device(st.inv_rho_x[0:size], st.inv_rho_y[0:size], st.inv_rho_z[0:size])

    // --- V1：L * 1 = 0（逐位）---
    {
        std::vector<double> one(size, 1.0), out(size);
        for (size_t t = 0; t < size; t++) { st.h_p[t] = 1.0; }
        #pragma acc update device(st.p[0:size])
        apply_laplacian_var(cfg, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z, st.p, st.cg_r);
        #pragma acc update self(st.cg_r[0:size])
        double mx = 0.0;
        for (size_t t = 0; t < size; t++) { mx = std::fmax(mx, std::fabs(st.cg_r[t])); }
        check(mx == 0.0, (std::string("V1 ") + tag + "：L*1 = 0（逐位）").c_str());
    }

    // --- V2：对称性 (Lp,q) = (p,Lq) ---
    {
        std::vector<double> p(size), q(size);
        fill_rand(p, 111u + (unsigned)amp * 7u);
        fill_rand(q, 222u + (unsigned)amp * 7u);
        double pm = 0.0, qm = 0.0;
        for (size_t t = 0; t < size; t++) { pm += p[t]; qm += q[t]; }
        pm /= (double)size; qm /= (double)size;
        for (size_t t = 0; t < size; t++) { p[t] -= pm; q[t] -= qm; st.h_p[t] = p[t]; st.h_cg_s[t] = q[t]; }
        #pragma acc update device(st.p[0:size], st.cg_s[0:size])
        apply_laplacian_var(cfg, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z, st.p, st.cg_r);   // Lp
        apply_laplacian_var(cfg, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z, st.cg_s, st.cg_zq); // Lq
        #pragma acc update self(st.cg_r[0:size], st.cg_zq[0:size])
        const double lhs = dot_host(std::vector<double>(st.cg_r, st.cg_r + size), q);
        const double rhs = dot_host(p, std::vector<double>(st.cg_zq, st.cg_zq + size));
        const double scale = std::fmax(std::fabs(lhs), std::fabs(rhs));
        const double rel = scale > 0.0 ? std::fabs(lhs - rhs) / scale : 0.0;
        std::printf("      %s 对称性 (Lp,q)=%.6e  (p,Lq)=%.6e  相对差 %.3e\n",
                    tag, lhs, rhs, rel);
        check(rel < 1e-12, (std::string("V2 ") + tag + "：变系数算子对称").c_str());
    }

    // --- V3：往返 ---
    int iters = 0;
    {
        std::vector<double> pref(size);
        fill_rand(pref, 333u + (unsigned)amp * 7u);
        double m = 0.0;
        for (size_t t = 0; t < size; t++) { m += pref[t]; }
        m /= (double)size;
        for (size_t t = 0; t < size; t++) { pref[t] -= m; st.h_p[t] = pref[t]; }
        #pragma acc update device(st.p[0:size])
        apply_laplacian_var(cfg, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z, st.p, st.cg_r); // b = L pref
        #pragma acc update self(st.cg_r[0:size])
        std::vector<double> b(st.cg_r, st.cg_r + size);
        for (size_t t = 0; t < size; t++) { st.h_cg_r[t] = b[t]; st.h_p[t] = 0.0; }
        #pragma acc update device(st.cg_r[0:size], st.p[0:size])

        // 兼容性：b 的和应 ~0
        double bsum = 0.0;
        for (size_t t = 0; t < size; t++) { bsum += b[t]; }

        iters = solve_pressure_var(cfg, st.cg_r, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                                   st.fft, st.tri_w, st.plan_xy, st.p, st.diag,
                                   st.cg_s, st.cg_zq);
        #pragma acc update self(st.p[0:size], st.diag[0:2])

        double pm = 0.0;
        for (size_t t = 0; t < size; t++) { pm += st.p[t] - pref[t]; }
        pm /= (double)size;
        double mx = 0.0, rmax = 0.0;
        for (size_t t = 0; t < size; t++)
        {
            mx = std::fmax(mx, std::fabs((st.p[t] - pref[t]) - pm));
            rmax = std::fmax(rmax, std::fabs(pref[t]));
        }
        const double rel = rmax > 0.0 ? mx / rmax : mx;
        std::printf("      %s 往返：max|(p-p_ref)-mean| = %.3e   PCG 迭代 %d   |Σb|=%.2e   残差 %.2e\n",
                    tag, rel, iters, std::fabs(bsum), st.diag[1]);
        check(rel < 1e-10, (std::string("V3 ") + tag + "：变系数算子往返").c_str());
        check(iters < 200 && st.diag[1] < 1e-11,
              (std::string("V5 ") + tag + "：PCG 收敛（迭代 < 200，相对残差 < 1e-11）").c_str());
    }

    // --- V4：一步完整投影后 div v = 0 ---
    {
        const int size_i = (int)size;
        std::vector<double> vx(size), vy(size), vz(size);
        fill_rand(vx, 444u);
        fill_rand(vy, 555u);
        fill_rand(vz, 666u);
        // 壁面法向面钉死
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++) { vz[(size_t)IDX(i, j, Nz - 1)] = 0.0; }
        for (int t = 0; t < size_i; t++)
        {
            st.h_vx[t] = vx[t]; st.h_vy[t] = vy[t]; st.h_vz[t] = vz[t];
            st.h_fx[t] = 0.0; st.h_fy[t] = 0.0; st.h_fz[t] = 0.0;
        }
        #pragma acc update device(st.vx[0:size], st.vy[0:size], st.vz[0:size], \
                                  st.fx[0:size], st.fy[0:size], st.fz[0:size])

        update_viscosity_fields(cfg, pp, 0, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);

        step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p,
                           st.fx, st.fy, st.fz,
                           st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                           st.cg_r, st.cg_s, st.cg_zq,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz,
                           st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan_xy, st.tri_w, st.diag,
                           st.gen, st.randD, st.randN,
                           st.tmp_fx, st.tmp_fy, st.tmp_fz, 0L);
        st.download(ST_VELOCITY);

        double vmax = 0.0, dmax = 0.0;
        for (int k = 0; k < Nz; k++)
            for (int j = 0; j < Ny; j++)
                for (int i = 0; i < Nx; i++)
                {
                    const int ijk = IDX(i, j, k);
                    const int im = (i - 1 + Nx) % Nx;
                    const int jm = (j - 1 + Ny) % Ny;
                    const int km = (k - 1 + Nz) % Nz;
                    const double div = (st.vx[ijk] - st.vx[IDX(im, j, k)])
                                     + (st.vy[ijk] - st.vy[IDX(i, jm, k)])
                                     + (st.vz[ijk] - st.vz[IDX(i, j, km)]);
                    dmax = std::fmax(dmax, std::fabs(div));
                    vmax = std::fmax(vmax, std::fmax(std::fabs(st.vx[ijk]),
                                        std::fmax(std::fabs(st.vy[ijk]), std::fabs(st.vz[ijk]))));
                }
        const double rel = vmax > 0.0 ? dmax / vmax : dmax;
        std::printf("      %s 一步投影：max|div v|/max|v| = %.3e\n", tag, rel);
        check(rel < 1e-10, (std::string("V4 ") + tag + "：一步投影后 div v = 0").c_str());
    }

    st.finish();
}

// ============================================================================
// --debug-pcg（需 GPU）：测 PCG 到双精度地板的迭代数。
//
// p_ref 是零均值参考解，b = L p_ref 用生产算子算；从 p=0 跑生产
// solve_pressure_var（trace 覆盖 tol=0 / max_iter=300），逐迭代记录
//   res = ||r_k||/||b||，err = ||(p_k-p_ref)-mean||/||p_ref||。
// 只打印标量；p/p_ref/b 不进任何文件。
// ============================================================================
struct PcgStats
{
    int it_res[6] = {0, 0, 0, 0, 0, 0};
    int it_err[6] = {0, 0, 0, 0, 0, 0};
    int min_res_it = 0, min_err_it = 0;
    double min_res = 1e300, min_err = 1e300;
    double final_res = 0.0, final_err = 0.0;
    double direct_err = -1.0;
};

static int first_at(const std::vector<double>& h, int n, double th)
{
    for (int i = 0; i < n; i++) { if (h[i] <= th) { return i + 1; } }
    return 0;
}

static void print_row(const char* tag, const std::vector<double>& h, int n)
{
    static const double th[6] = {1e-8, 1e-10, 1e-12, 1e-14, 1e-15, 1e-16};
    std::printf("  %-24s", tag);
    for (int q = 0; q < 6; q++)
    {
        const int it = first_at(h, n, th[q]);
        if (it > 0) { std::printf(" %9d", it); }
        else        { std::printf(" %9s", "--"); }
    }
    int mi = 0;
    for (int i = 1; i < n; i++) { if (h[i] < h[mi]) { mi = i; } }
    std::printf("   %.3e@%d\n", n > 0 ? h[mi] : 0.0, n > 0 ? mi + 1 : 0);
}

static void pcg_measure(FpdState& st, const char* tag, bool direct, PcgStats& out)
{
    const size_t size = st.size;

    // 零均值参考解 p_ref（确定性 LCG）
    std::vector<double> pref(size);
    fill_rand(pref, 9102u + (unsigned)(size & 0xFFFFu));
    double m = 0.0;
    for (size_t t = 0; t < size; t++) { m += pref[t]; }
    m /= (double)size;
    double pn2 = 0.0;
    for (size_t t = 0; t < size; t++) { pref[t] -= m; pn2 += pref[t] * pref[t]; }
    const double pn = std::sqrt(pn2);

    acc_copyin(pref.data(), size * sizeof(double));

    // b = L p_ref（生产算子；a 由调用方设好）
    apply_laplacian_var(st.cfg, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                        pref.data(), st.cg_r);
    #pragma acc update self(st.cg_r[0:size])
    std::vector<double> b_host(st.cg_r, st.cg_r + size);

    // p0 = 0
    for (size_t t = 0; t < size; t++) { st.h_p[t] = 0.0; }
    #pragma acc update device(st.p[0:size])

    const int max_iter = 300;
    std::vector<double> res(max_iter + 1, 0.0), err(max_iter + 1, 0.0);
    PcgTrace tr;
    tr.p_ref = pref.data();
    tr.res = res.data();
    tr.err = err.data();
    tr.cap = max_iter + 1;
    tr.tol = 0.0;                 // 不提前停：要看双精度地板
    tr.max_iter = max_iter;

    solve_pressure_var(st.cfg, st.cg_r, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                       st.fft, st.tri_w, st.plan_xy, st.p, st.diag,
                       st.cg_s, st.cg_zq, &tr);
    #pragma acc update self(st.diag[0:2])
    const int n = tr.n;
    for (int i = 0; i < n; i++) { err[i] /= pn; }

    // 统计
    static const double th[6] = {1e-8, 1e-10, 1e-12, 1e-14, 1e-15, 1e-16};
    for (int q = 0; q < 6; q++)
    {
        out.it_res[q] = first_at(res, n, th[q]);
        out.it_err[q] = first_at(err, n, th[q]);
    }
    out.min_res = res[0]; out.min_res_it = 1;
    out.min_err = err[0]; out.min_err_it = 1;
    for (int i = 1; i < n; i++)
    {
        if (res[i] < out.min_res) { out.min_res = res[i]; out.min_res_it = i + 1; }
        if (err[i] < out.min_err) { out.min_err = err[i]; out.min_err_it = i + 1; }
    }
    out.final_res = res[n - 1];
    out.final_err = err[n - 1];

    // 直解对照（仅 a≡1 有意义）
    if (direct)
    {
        for (size_t t = 0; t < size; t++)
        {
            st.h_fft[2 * t]     = b_host[t];
            st.h_fft[2 * t + 1] = 0.0;
        }
        #pragma acc update device(st.fft[0:2 * size])
        solve_pressure(st.cfg, st.fft, st.tri_w, st.plan_xy, st.p, st.diag);
        #pragma acc update self(st.p[0:size])

        double em = 0.0, e2 = 0.0;
        for (size_t t = 0; t < size; t++) { em += st.p[t] - pref[t]; }
        em /= (double)size;
        for (size_t t = 0; t < size; t++)
        {
            const double d = (st.p[t] - pref[t]) - em;
            e2 += d * d;
        }
        out.direct_err = std::sqrt(e2) / pn;
    }

    acc_delete(pref.data(), size * sizeof(double));

    std::string rtag = std::string(tag) + " [res]";
    std::string etag = std::string(tag) + " [err]";
    print_row(rtag.c_str(), res, n);
    print_row(etag.c_str(), err, n);
    std::printf("  %-24s final_res=%.3e  final_err=%.3e  n=%d\n",
                (std::string(tag) + " [末]").c_str(), out.final_res, out.final_err, n);
    if (direct)
    {
        std::printf("  %-24s direct FFT: err=%.3e（双精度地板参照）\n",
                    (std::string(tag) + " [直解]").c_str(), out.direct_err);
    }
}

int run_debug_pcg(int Nx, int Ny, int Nz)
{
    g_fail = 0;
    if (Nx <= 0 || Ny <= 0 || Nz <= 0) { Nx = Ny = Nz = 32; }

    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.002, 0.0, false);
    cfg.density_on = 1;
    PhiParams pp = make_phi_params(3.2, 1.0, 50.0);
    const size_t size = (size_t)Nx * Ny * Nz;

    std::printf("=== PCG 精度调试（tol=0, max_iter=300；只输出标量，不导出 p）===\n");
    std::printf("  网格 %d x %d x %d   参考解 p_ref 为零均值确定性随机场\n", Nx, Ny, Nz);
    std::printf("  %-24s %9s %9s %9s %9s %9s %9s   %s\n",
                "case", "1e-8", "1e-10", "1e-12", "1e-14", "1e-15", "1e-16", "min(it)");

    PcgStats s1, sg, ss, sp;

    // 1) a ≡ 1（预条件器精确，应 1 步到位）
    {
        FpdState st; st.init(cfg, 0, ST_FULL, 20261018ULL);
        for (size_t t = 0; t < size; t++)
        { st.h_inv_rho_x[t] = 1.0; st.h_inv_rho_y[t] = 1.0; st.h_inv_rho_z[t] = 1.0; }
        st.upload(ST_PHI);
        pcg_measure(st, "a=1 (PCG)", true, s1);
        st.finish();
    }

    // 2) 高斯 a（amp=1, sigma=4）
    {
        FpdState st; st.init(cfg, 0, ST_FULL, 20261018ULL);
        std::vector<double> ax, ay, az;
        build_var_a(ax, 0, cfg, 1.0, 4.0, false, 0.0);
        build_var_a(ay, 1, cfg, 1.0, 4.0, false, 0.0);
        build_var_a(az, 2, cfg, 1.0, 4.0, false, 0.0);
        for (size_t t = 0; t < size; t++)
        { st.h_inv_rho_x[t] = ax[t]; st.h_inv_rho_y[t] = ay[t]; st.h_inv_rho_z[t] = az[t]; }
        st.upload(ST_PHI);
        pcg_measure(st, "gaussian (amp=1)", false, sg);
        st.finish();
    }

    // 3) 球状强对比 a（amp=3, R=5）
    {
        FpdState st; st.init(cfg, 0, ST_FULL, 20261018ULL);
        std::vector<double> ax, ay, az;
        build_var_a(ax, 0, cfg, 3.0, 0.0, true, 5.0);
        build_var_a(ay, 1, cfg, 3.0, 0.0, true, 5.0);
        build_var_a(az, 2, cfg, 3.0, 0.0, true, 5.0);
        for (size_t t = 0; t < size; t++)
        { st.h_inv_rho_x[t] = ax[t]; st.h_inv_rho_y[t] = ay[t]; st.h_inv_rho_z[t] = az[t]; }
        st.upload(ST_PHI);
        pcg_measure(st, "sphere (amp=3)", false, ss);
        st.finish();
    }

    // 4) 真实粒子场（N=1，ratio_rho=2）；只在支撑域放得下时跑
    if (Nx >= 32 && Ny >= 32 && Nz >= 32)
    {
        FpdState st; st.init(cfg, 1, ST_FULL, 20261018ULL);
        st.h_Rx[0] = 0.5 * (double)Nx;
        st.h_Ry[0] = 0.5 * (double)Ny;
        st.h_Rz[0] = 0.5 * (double)Nz;
        st.upload(ST_PARTICLE);
        update_mass_fields(cfg, pp, 1.0, 1, st.Rx, st.Ry, st.Rz,
                           st.inv_rho_x, st.inv_rho_y, st.inv_rho_z);
        pcg_measure(st, "particle (ratio=2)", false, sp);
        st.finish();
    }
    else
    {
        std::printf("  （网格 < 32³：跳过粒子场用例，支撑域会被盒子截断）\n");
    }

    // a=1 的硬判据：一步到 <1e-15
    check(s1.it_err[4] == 1 && s1.min_err < 1e-15,
          "P0 a=1：PCG 1 步达到 <1e-15（预条件器精确）");
    if (s1.direct_err >= 0.0)
    {
        // 直解地板是几 ulp（32³ 实测 ~1.2e-15），不能拿 1e-15 当阈值。
    check(s1.direct_err < 1e-14, "P0 a=1：FFT 直解误差 <1e-14（双精度地板）");
    }

    std::printf("\n%s\n", g_fail == 0 ? "完成（调试结果仅供判断，不改变生产 tol）" : "有失败项");
    return g_fail;
}

// ============================================================================
// --check-poisson（需 GPU）：壁面算子往返，覆盖偶/奇 Nz。
// ============================================================================
int run_check_poisson(int Nx, int Ny, int Nz)
{
    g_fail = 0;
    std::printf("=== 泊松算子往返判据（需 GPU）===\n");

    // ⚠️ 返回值必须累加进 g_fail。曾经这里 6 个调用点全部丢弃返回值，
    //    于是往返判据 FAIL 时只打印一行 FAIL、退出码仍是 0 —— 所有依赖
    //    --check-poisson 退出码的回归（含 Phase 7-B 的壁面接线）都会变成假的。
    if (Nx <= 0)
    {
        // 默认电池：偶 Nz、奇 Nz、中等盒子各一（奇偶不对称最容易在奇数 Nz 暴露）
        g_fail += roundtrip_grid(8, 4, 6, "P5 壁面算子往返（偶 Nz）");
        g_fail += roundtrip_grid(6, 6, 5, "P5 壁面算子往返（奇 Nz）");
        g_fail += roundtrip_grid(16, 8, 8, "P6 壁面算子往返（中盒子）");
        check_compat_diag();

        std::printf("\n--- 变系数投影（div(a grad p)）---\n");
        v0_convection_independent_of_a(8, 6, 5);
        var_poisson_battery(16, 16, 16, "高斯 a（amp=1, sigma=4）", 1.0, 4.0, false, 0.0);
        var_poisson_battery(16, 16, 16, "球状 a（amp=1, R=6）",        1.0, 0.0, true, 6.0);
        var_poisson_battery(16, 16, 16, "强对比 a（amp=3, R=5）",      3.0, 0.0, true, 5.0);
    }
    else
    {
        g_fail += roundtrip_grid(Nx, Ny, Nz, "壁面算子往返");
        check_compat_diag();
        v0_convection_independent_of_a(8, 6, 5);
        var_poisson_battery(Nx, Ny, Nz, "高斯 a（amp=1, sigma=4）", 1.0, 4.0, false, 0.0);
    }

    std::printf("\n%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}
