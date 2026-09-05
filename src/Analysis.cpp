#include "./include/Analysis.h"

#include <cstdio>
#include <cmath>

// ---------------------------------------------------------------------------
// 对单个 Loc 做离散求和。复用 Stencil.h 的 stencil_point —— 不要在这里
// 重写模板盒遍历，那正是缺陷 C1/C2 的成因。
// ---------------------------------------------------------------------------
template<Loc L>
static void sum_phi_powers(PhiParams pp, NS_Config cfg,
                           double Rx, double Ry, double Rz,
                           double& s1, double& s2)
{
    const int n3 = stencil_size(pp);
    s1 = 0.0; s2 = 0.0;

    for (int l = 0; l < n3; l++)
    {
        int ijk; double w;
        if (stencil_point<L>(pp, cfg, Rx, Ry, Rz, l, ijk, w))
        {
            s1 += w;
            s2 += w * w;
        }
    }
}

// 由三方向的 (Sum phi, Sum phi^2) 填出派生量
static void fill_derived(FpdConstants& c)
{
    for (int d = 0; d < 3; d++)
    {
        c.lambda_T[d] = c.int_phi[d] / c.int_phi2[d];
        c.mass[d]     = c.int_phi[d] * c.int_phi[d] / c.int_phi2[d];   // rho = dV = 1
        c.mass_eff[d] = 1.5 * c.mass[d];
    }

    double p_lo = c.int_phi[0],  p_hi = c.int_phi[0];
    double q_lo = c.int_phi2[0], q_hi = c.int_phi2[0];
    double m_lo = c.mass[0],     m_hi = c.mass[0];

    for (int d = 1; d < 3; d++)
    {
        if (c.int_phi[d]  < p_lo) { p_lo = c.int_phi[d];  }
        if (c.int_phi[d]  > p_hi) { p_hi = c.int_phi[d];  }
        if (c.int_phi2[d] < q_lo) { q_lo = c.int_phi2[d]; }
        if (c.int_phi2[d] > q_hi) { q_hi = c.int_phi2[d]; }
        if (c.mass[d]     < m_lo) { m_lo = c.mass[d];     }
        if (c.mass[d]     > m_hi) { m_hi = c.mass[d];     }
    }

    c.spread_int_phi  = (p_hi - p_lo) / p_lo;
    c.spread_int_phi2 = (q_hi - q_lo) / q_lo;
    c.spread_mass     = (m_hi - m_lo) / m_lo;
}

FpdConstants compute_fpd_constants(PhiParams pp, NS_Config cfg,
                                   double Rx, double Ry, double Rz, int n_scan)
{
    FpdConstants c;

    // --- 主值 ---
    sum_phi_powers<FACE_X>(pp, cfg, Rx, Ry, Rz, c.int_phi[0], c.int_phi2[0]);
    sum_phi_powers<FACE_Y>(pp, cfg, Rx, Ry, Rz, c.int_phi[1], c.int_phi2[1]);
    sum_phi_powers<FACE_Z>(pp, cfg, Rx, Ry, Rz, c.int_phi[2], c.int_phi2[2]);
    fill_derived(c);

    // --- 亚格点扫描：格胞内 n_scan^3 个偏移 ---
    // 这个散布就是整个方法的系统误差下限：若它是 0.5%，
    // 就没法验证到优于 0.5% 的精度。
    c.n_scan   = n_scan;
    c.mass_min = 1e300;
    c.mass_max = -1e300;
    double acc = 0.0;
    int    cnt = 0;

    for (int a = 0; a < n_scan; a++)
    {
        for (int b = 0; b < n_scan; b++)
        {
            for (int d = 0; d < n_scan; d++)
            {
                const double ox = (double)a / (double)n_scan;
                const double oy = (double)b / (double)n_scan;
                const double oz = (double)d / (double)n_scan;

                double s1, s2;
                sum_phi_powers<FACE_X>(pp, cfg, Rx + ox, Ry + oy, Rz + oz, s1, s2);
                const double m = s1 * s1 / s2;

                acc += m; cnt++;
                if (m < c.mass_min) { c.mass_min = m; }
                if (m > c.mass_max) { c.mass_max = m; }
            }
        }
    }

    c.mass_mean      = acc / (double)cnt;
    c.spread_subgrid = (c.mass_max - c.mass_min) / c.mass_min;

    return c;
}

// ---------------------------------------------------------------------------
// 连续球坐标积分，作为与代码完全无关的独立对照
// ---------------------------------------------------------------------------
static void continuum_integrals(PhiParams pp, double& i1, double& i2)
{
    const double a   = pp.radius;
    const double xi  = 1.0 / pp.inv_xi;
    const double R   = sqrt(pp.range2);
    const int    n   = 200000;
    const double h   = R / (double)n;

    i1 = 0.0; i2 = 0.0;
    for (int i = 0; i <= n; i++)
    {
        const double r = i * h;
        const double p = 0.5 * (tanh((a - r) / xi) + 1.0);
        const double w = (i == 0 || i == n) ? 0.5 : 1.0;
        const double s = w * 4.0 * M_PI * r * r * h;
        i1 += s * p;
        i2 += s * p * p;
    }
}

void print_fpd_constants(FpdConstants fc, PhiParams pp)
{
    const char* dir[3] = {"FACE_X", "FACE_Y", "FACE_Z"};

    printf("=== FPD 常量表 ===\n");
    printf("  a = %.4f   xi = %.4f   eta_c/eta_l = %.1f   截断半径 = %.1f   模板盒 %d^3\n",
           pp.radius, 1.0 / pp.inv_xi, pp.ratio_eta, sqrt(pp.range2), pp.n_range);
    printf("\n");
    printf("  方向      int(phi)      int(phi^2)    lambda^T      M_i        M_eff\n");
    for (int d = 0; d < 3; d++)
    {
        printf("  %-8s %11.6f  %12.6f  %10.6f  %10.4f  %10.4f\n",
               dir[d], fc.int_phi[d], fc.int_phi2[d],
               fc.lambda_T[d], fc.mass[d], fc.mass_eff[d]);
    }
    printf("\n");
    printf("  三方向散布   int(phi) %.3e   int(phi^2) %.3e   M_i %.3e\n",
           fc.spread_int_phi, fc.spread_int_phi2, fc.spread_mass);
    printf("    -> int(phi^2) 的一致性是交错混叠的直接判据（%s，阈值 1e-4）\n",
           fc.spread_int_phi2 < 1e-4 ? "PASS" : "FAIL");
    printf("\n");
    printf("  亚格点扫描 %d^3 个偏移：M_i = %.4f   [%.4f, %.4f]   散布 %.4f%%\n",
           fc.n_scan, fc.mass_mean, fc.mass_min, fc.mass_max, fc.spread_subgrid * 100.0);
    printf("    -> 这是整个方法的【系统误差下限】（%s，阈值 1%%）\n",
           fc.spread_subgrid < 0.01 ? "PASS" : "FAIL");
    printf("\n");

    // --- 独立对照 ---
    double i1, i2;
    continuum_integrals(pp, i1, i2);
    const double m_cont = i1 * i1 / i2;

    const double a  = pp.radius;
    const double xi = 1.0 / pp.inv_xi;
    const double sphere = 4.0 / 3.0 * M_PI * a * a * a;
    const double curv   = 8.0 * M_PI * a * xi * xi * (M_PI * M_PI / 24.0);

    printf("  独立对照（与代码无关）\n");
    printf("    连续球坐标积分    int(phi) = %.4f   int(phi^2) = %.4f   M_i = %.4f\n",
           i1, i2, m_cont);
    printf("    离散 vs 连续      %.4f%%          %.4f%%           %.4f%%\n",
           fabs(fc.int_phi[0] - i1) / i1 * 100.0,
           fabs(fc.int_phi2[0] - i2) / i2 * 100.0,
           fabs(fc.mass[0] - m_cont) / m_cont * 100.0);
    printf("\n");
    printf("    扩散界面解析展开  int(phi) ~ (4/3)pi a^3 + 8 pi a xi^2 pi^2/24\n");
    printf("                             = %.3f + %.3f = %.3f   (vs 连续 %.3f)\n",
           sphere, curv, sphere + curv, i1);
    printf("    注意：解析球体 (4/3)pi a^3 = %.3f 比真值小 %.1f%% ——\n",
           sphere, (i1 / sphere - 1.0) * 100.0);
    printf("          这是【扩散界面的曲率项】，不是格点离散化误差。\n");
    printf("          格点本身精确到 %.1e —— 用解析球体算 M_i 才会毁掉验证。\n",
           fabs(fc.int_phi[0] - i1) / i1);
}

double equipartition_target(FpdConstants fc, double kT, NS_Config cfg)
{
    const double L3 = (double)cfg.Nx * (double)cfg.Ny * (double)cfg.Nz;
    return 2.0 * kT / fc.mass_mean - 2.0 * kT / L3;    // rho = 1
}

// ============================================================================
// 分块平均
// ============================================================================
int blocking_analysis(const double* x, int n, BlockStat* out, int max_out)
{
    // 全样本均值
    double mean = 0.0;
    for (int i = 0; i < n; i++) { mean += x[i]; }
    mean /= (double)n;

    // 朴素（把相关样本当独立）的标准误，用来反推 tau_int
    double var0 = 0.0;
    for (int i = 0; i < n; i++) { const double d = x[i] - mean; var0 += d * d; }
    var0 /= (double)(n - 1);
    const double stderr_naive = sqrt(var0 / (double)n);

    int cnt = 0;
    for (int b = 1; b <= n / 4 && cnt < max_out; b *= 2)
    {
        const int nb = n / b;
        if (nb < 4) { break; }

        double m = 0.0;
        for (int k = 0; k < nb; k++)
        {
            double s = 0.0;
            for (int i = 0; i < b; i++) { s += x[k * b + i]; }
            m += s / (double)b;
        }
        m /= (double)nb;

        double v = 0.0;
        for (int k = 0; k < nb; k++)
        {
            double s = 0.0;
            for (int i = 0; i < b; i++) { s += x[k * b + i]; }
            const double d = s / (double)b - m;
            v += d * d;
        }
        v /= (double)(nb - 1);

        BlockStat st;
        st.block_len    = b;
        st.n_block      = nb;
        st.mean         = m;
        st.stderr_mean  = sqrt(v / (double)nb);
        st.tau_int_samp = 0.5 * (st.stderr_mean / stderr_naive)
                              * (st.stderr_mean / stderr_naive);
        out[cnt++] = st;
    }
    return cnt;
}

BlockStat pick_plateau(const BlockStat* s, int n_s, int min_block)
{
    for (int i = 0; i + 1 < n_s; i++)
    {
        if (s[i + 1].n_block < min_block) { break; }
        const double r = s[i + 1].stderr_mean / s[i].stderr_mean;
        if (r < 1.05) { return s[i + 1]; }
    }

    // 没有平台 —— 误差棒不可信，必须延长模拟
    BlockStat bad = s[n_s > 0 ? n_s - 1 : 0];
    bad.block_len = -1;
    return bad;
}

void print_blocking(const BlockStat* s, int n_s, BlockStat picked)
{
    printf("  分块平均（平台是否存在要肉眼判，不是公式判）\n");
    printf("     b      n_b        均值          标准误      tau_int(采样数)\n");
    for (int i = 0; i < n_s; i++)
    {
        printf("  %6d %8d  %12.6e  %12.4e  %10.2f%s\n",
               s[i].block_len, s[i].n_block, s[i].mean,
               s[i].stderr_mean, s[i].tau_int_samp,
               (picked.block_len == s[i].block_len) ? "   <== 平台" : "");
    }
    if (picked.block_len < 0)
    {
        printf("  ** 找不到平台：误差棒不可信，必须延长模拟。不许降标准接受此结果 **\n");
    }
}
