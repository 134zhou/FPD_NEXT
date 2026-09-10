#include <cstdio>
#include <cmath>
#include <vector>

#include "./include/Poisson.h"
#include "./include/State.h"
#include "./include/Tests.h"
#include "./include/Check.h"

// ============================================================================
// 压力泊松求解器自检。两个入口：
//   --check-tridiag  纯 CPU，无卡可跑：三对角 Thomas + build_tridiag_coeffs
//                    对照稠密 Gauss，含 (0,0) 奇异列的定规。
//   --check-poisson  需 GPU：算子往返判据（手写 div∘grad 作用出右端，求解器
//                    还原 p，两者只差一个常数）。周期与壁面各跑一遍。
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
        NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false, /*wall_z=*/1);

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
        NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false, 1);
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
// 独立参考 3：手写离散算子 div∘grad（周期与壁面各一份，ghost 表述）。
//   周期：km=(k-1+Nz)%Nz, kp=(k+1)%Nz
//   壁面：km=(k==0?0:k-1), kp=(k==Nz-1?Nz-1:k+1)   （Neumann ghost p[-1]=p[0] 等）
// 两者统一为  p[ip]+p[im]+p[jp]+p[jm]+p[kp]+p[km] - 6 p[ijk]。
// ============================================================================
static void apply_operator(int Nx, int Ny, int Nz, bool wall, const double* p, double* b)
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
                const int kp = wall ? (k == Nz - 1 ? Nz - 1 : k + 1) : (k + 1) % Nz;
                const int km = wall ? (k == 0 ? 0 : k - 1)           : (k - 1 + Nz) % Nz;
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
static int roundtrip_grid(int Nx, int Ny, int Nz, int wall_z, const char* tag)
{
    const int size = Nx * Ny * Nz;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false, wall_z);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 12345ULL);

    // 随机参考压力 + 手写算子作用出右端
    std::vector<double> p_ref((size_t)size);
    fill_rand(p_ref, (unsigned)(Nx * 1000 + Ny * 31 + Nz));
    std::vector<double> b((size_t)size);
    apply_operator(Nx, Ny, Nz, wall_z != 0, p_ref.data(), b.data());

    // 填 fft_data（复数交错，虚部 0）并上传
    for (int t = 0; t < size; t++) { st.h_fft[2 * t] = b[t]; st.h_fft[2 * t + 1] = 0.0; }
    #pragma acc update device(st.fft[0:2*size])

    // 壁面模式的相容性诊断数组（设备侧长度 1）
    double diag_h[1] = { 0.0 };
    double* diag = 0;
    if (wall_z)
    {
        diag = diag_h;
        #pragma acc enter data create(diag[0:1])
    }

    solve_pressure(cfg, st.fft, st.tri_w, st.plan, st.plan_xy, st.p, diag);

    #pragma acc update self(st.p[0:size])
    double compat = 0.0;
    if (wall_z)
    {
        #pragma acc update self(diag[0:1])
        compat = diag_h[0];
        #pragma acc exit data delete(diag[0:1])
    }

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
    std::printf("      %dx%dx%d (wall_z=%d)   max|(p-p_ref)-mean| = %.3e   阈值 %.3e\n",
                Nx, Ny, Nz, wall_z, mx, tol);
    if (wall_z)
    {
        std::printf("      相容性残差 |Σb̂| = %.3e   （应 ~机器精度，非 O(1)）\n", compat);
        if (compat > 1e-9) { g_fail++; }   // 相容投影没生效的硬信号
    }
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
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, false, /*wall_z=*/1);

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

    solve_pressure(cfg, st.fft, st.tri_w, st.plan, st.plan_xy, st.p, diag);

    #pragma acc update self(diag[0:1])
    const double compat = diag_h[0];
    #pragma acc exit data delete(diag[0:1])
    st.finish();

    const bool ok = (std::fabs(host_sum) > 1e-3) && (std::fabs(compat - std::fabs(host_sum)) < 1e-9 * (1.0 + std::fabs(host_sum)));
    check(ok, "P6 相容性残差诊断（不相容右端应报非零且匹配 |Σb|）");
    std::printf("      |Σ_{ijk} b|（主机） = %.3e    diag = %.3e\n", std::fabs(host_sum), compat);
}

// ============================================================================
// --check-poisson（需 GPU）：周期 + 壁面算子往返，覆盖偶/奇 Nz。
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
        g_fail += roundtrip_grid(8, 4, 6, 0, "P5 周期算子往返（偶 Nz）");
        g_fail += roundtrip_grid(8, 4, 6, 1, "P5 壁面算子往返（偶 Nz）");
        g_fail += roundtrip_grid(6, 6, 5, 0, "P5 周期算子往返（奇 Nz）");
        g_fail += roundtrip_grid(6, 6, 5, 1, "P5 壁面算子往返（奇 Nz）");
        g_fail += roundtrip_grid(16, 8, 8, 0, "P6 周期算子往返（中盒子）");
        g_fail += roundtrip_grid(16, 8, 8, 1, "P6 壁面算子往返（中盒子）");
        check_compat_diag();
    }
    else
    {
        g_fail += roundtrip_grid(Nx, Ny, Nz, 0, "周期算子往返");
        g_fail += roundtrip_grid(Nx, Ny, Nz, 1, "壁面算子往返");
        check_compat_diag();
    }

    std::printf("\n%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}
