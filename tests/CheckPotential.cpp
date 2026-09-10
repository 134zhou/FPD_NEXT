#include <cstdio>
#include <cmath>
#include <vector>

#include "Potential.h"
#include "Tests.h"

// ============================================================================
// 势函数自检（纯 CPU，无卡可跑）
//
// 全部判据只依赖 pair_energy_bare / pair_force_over_r_bare / min_image 这些
// routine seq 的 inline 函数 —— 它们主机侧直接可调（与 Analysis.cpp 在主机侧
// 复用 stencil_point 同一种做法）。势公式的单元判据不需要 GPU。
// ============================================================================

static int g_fail = 0;
static void check(bool ok, const char* what)
{
    std::printf("  %-52s %s\n", what, ok ? "PASS" : "FAIL");
    if (!ok) { g_fail++; }
}

// CPU 参考实现：串行 i<j 半矩阵。用于 ΣF 判据与 N=3 解析对照（与 GPU 的
// compute_particle_forces 是两份独立实现，各自验证装配逻辑）。
static void forces_cpu(const PotentialParams& pot, const ExternalField& ext,
                       int N, double Lx, double Ly, double Lz,
                       const double* Rx, const double* Ry, const double* Rz,
                       double* Fx, double* Fy, double* Fz)
{
    for (int n = 0; n < N; n++) { Fx[n] = ext.gx; Fy[n] = ext.gy; Fz[n] = ext.gz; }
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
        {
            const double dx = min_image(Rx[i] - Rx[j], Lx);
            const double dy = min_image(Ry[i] - Ry[j], Ly);
            const double dz = min_image(Rz[i] - Rz[j], Lz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double g  = pair_force_over_r(pot, r2);
            Fx[i] += g*dx;  Fx[j] -= g*dx;
            Fy[i] += g*dy;  Fy[j] -= g*dy;
            Fz[i] += g*dz;  Fz[j] -= g*dz;
        }
}

// ---------------------------------------------------------------------------
// J1：力 = -dU/dr 的四阶中心差分对照。
// pair_energy_bare 与 pair_force_over_r_bare 是两份独立实现，所以这是
// 独立数值对照，不是自洽性检验。
// ---------------------------------------------------------------------------
static void check_force_derivative(const PotentialParams& pot, double r_lo, double r_hi)
{
    const int npt = 200;
    double max_rel = 0.0;
    bool sign_change = false;
    double prev_uprime = 0.0;
    bool have_prev = false;

    for (int k = 0; k <= npt; k++)
    {
        const double t = (double)k / (double)npt;
        const double r = r_lo * pow(r_hi / r_lo, t);   // 对数扫描，覆盖特征点两侧
        const double h = 1e-4 * r;                      // 相对步长

        // 四阶中心差分 U'(r) ≈ [U(r-2h) - 8U(r-h) + 8U(r+h) - U(r+2h)]/(12h)
        const double um2 = pair_energy_bare(pot, r - 2*h);
        const double um1 = pair_energy_bare(pot, r - h);
        const double up1 = pair_energy_bare(pot, r + h);
        const double up2 = pair_energy_bare(pot, r + 2*h);
        const double uprime_num = (um2 - 8.0*um1 + 8.0*up1 - up2) / (12.0*h);

        // 解析：U'(r) = -g(r)*r
        const double uprime_ana = -pair_force_over_r_bare(pot, r) * r;

        const double denom = std::fmax(std::fabs(uprime_ana), 1.0);
        const double rel = std::fabs(uprime_num - uprime_ana) / denom;
        if (rel > max_rel) { max_rel = rel; }

        if (have_prev && ((prev_uprime > 0) != (uprime_num > 0))) { sign_change = true; }
        prev_uprime = uprime_num;
        have_prev = true;
    }

    // 混合容差：Morse 在 r=r_eq 处力恰好为 0，纯相对容差会爆到 1e+21
    const double tol = 1e-7;
    check(max_rel < tol, "J1 力=-dU/dr 四阶差分对照 (max 相对差 < 1e-7)");
    std::printf("      max 相对差 = %.3e   符号变号 = %s\n",
                max_rel, sign_change ? "是（覆盖到吸引/排斥两侧）" : "否 *** 只测了单调段 ***");
    if (!sign_change) { g_fail++; }   // 覆盖率断言：没变号说明扫描区间不覆盖势阱
}

// ---------------------------------------------------------------------------
// J2：Python 黄金表对照（独立实现，非自洽）。
// 数值来自 tests/spike/spike_potential_ref.py，17 位有效数字。
// ---------------------------------------------------------------------------
struct GoldenPt { double r, U, g; };

static void check_golden(const PotentialParams& pot, const GoldenPt* pts, int n,
                         const char* name)
{
    double max_rel = 0.0;
    for (int i = 0; i < n; i++)
    {
        const double r2 = pts[i].r * pts[i].r;
        const double U = pair_energy(pot, r2);
        const double g = pair_force_over_r(pot, r2);
        const double relU = std::fabs(U - pts[i].U) / std::fmax(std::fabs(pts[i].U), 1.0);
        const double relg = std::fabs(g - pts[i].g) / std::fmax(std::fabs(pts[i].g), 1.0);
        max_rel = std::fmax(max_rel, std::fmax(relU, relg));
    }
    std::printf("  %-52s %s\n", name, max_rel < 1e-12 ? "PASS" : "FAIL");
    std::printf("      max 相对差 = %.3e   (阈值 1e-12)\n", max_rel);
    if (!(max_rel < 1e-12)) { g_fail++; }
}

// ---------------------------------------------------------------------------
// J3：最小镜像 vs 暴力枚举。用完全不同的算法复核同一个量。
// min_image 是一维的，所以暴力枚举一维镜像 q + n*L（n ∈ {-2..2}），
// 取 |p - (q+n*L)| 最小者，与 min_image(p-q, L) 对照。
// 容差而非逐位：min_image 的 d-L 与暴力枚举的 p-(q+n*L) 是两个不同的浮点
// 表达式，差 ~1 ULP，逐位会假红。真正的 C5 式错误（负号错）差 ~L 量级。
// ---------------------------------------------------------------------------
static void check_min_image()
{
    const double L = 32.0;
    bool ok = true;
    for (int trial = 0; trial < 1000 && ok; trial++)
    {
        // 用确定性伪随机（避免引入 <random> 的跨平台差异）
        const double a = (double)((trial * 2654435761u) % 1000000) / 1000000.0;
        const double b = (double)((trial * 40503u + 7u) % 1000000) / 1000000.0;
        const double p = a * L, q = b * L;

        const double d_min = min_image(p - q, L);
        // 暴力枚举 5 个镜像，取最小 |p - (q + n*L)|
        double best = 1e300; double bestd = 0.0;
        for (int n = -2; n <= 2; n++)
        {
            const double qm = q + n * L;
            const double d  = std::fabs(p - qm);
            if (d < best) { best = d; bestd = p - qm; }
        }
        // 容差：差 ~L 的错误（C5 式负号错）会被抓，1 ULP 的重排不误报
        if (std::fabs(d_min - bestd) > 1e-12 * L) { ok = false; }
    }
    check(ok, "J3 最小镜像 vs 5 镜像暴力枚举（1000 对，容差 1e-12*L）");
}

// ---------------------------------------------------------------------------
// J4：N=3 等边三角形解析对照。抓双计数 / 漏算某对 / 符号错 —— N=2 盲的。
// ---------------------------------------------------------------------------
static void check_n3_triangle()
{
    // WCA eps=1 sigma=7.4，边长 6.0，放在 32^3 盒心
    PotentialParams pot = make_potential_params(POT_WCA, SHIFT_ENERGY,
                                                1.0, 7.4, 0, 0, 0, 0);
    ExternalField ext{0, 0, 0};

    const double R = 6.0 / sqrt(3.0);
    const double cx = 16.0, cy = 16.0, cz = 16.0;
    double Rx[3], Ry[3], Rz[3], Fx[3], Fy[3], Fz[3];
    for (int i = 0; i < 3; i++)
    {
        const double ang = (90.0 + 120.0 * i) * M_PI / 180.0;
        Rx[i] = cx + R * cos(ang); Ry[i] = cy + R * sin(ang); Rz[i] = cz;
    }
    forces_cpu(pot, ext, 3, 32.0, 32.0, 32.0, Rx, Ry, Rz, Fx, Fy, Fz);

    const double F_expected = 147.25519213042094;   // |F_pair|*sqrt(3)
    bool ok = true;
    double sumx = 0, sumy = 0, sumz = 0;
    for (int i = 0; i < 3; i++)
    {
        const double Fmag = sqrt(Fx[i]*Fx[i] + Fy[i]*Fy[i] + Fz[i]*Fz[i]);
        if (std::fabs(Fmag - F_expected) > 1e-9 * F_expected) { ok = false; }
        sumx += Fx[i]; sumy += Fy[i]; sumz += Fz[i];
    }
    const double fmax = F_expected;
    const double sum_mag = sqrt(sumx*sumx + sumy*sumy + sumz*sumz);
    if (sum_mag > 1e-11 * 3 * fmax) { ok = false; }

    check(ok, "J4 N=3 等边三角形解析对照（|F_n| 相等且 ΣF=0）");
    std::printf("      |F_1|=%.6e  |F_2|=%.6e  |F_3|=%.6e   期望 %.6e\n",
                sqrt(Fx[0]*Fx[0]+Fy[0]*Fy[0]+Fz[0]*Fz[0]),
                sqrt(Fx[1]*Fx[1]+Fy[1]*Fy[1]+Fz[1]*Fz[1]),
                sqrt(Fx[2]*Fx[2]+Fy[2]*Fy[2]+Fz[2]*Fz[2]), F_expected);
    std::printf("      ΣF = (%.2e, %.2e, %.2e)\n", sumx, sumy, sumz);
}

// ---------------------------------------------------------------------------
// J6：势特征点 + POT_NONE 零判据 + 标号交换对称性。
// ---------------------------------------------------------------------------
static void check_feature_points()
{
    // Morse：U'(r_eq) = 0（力在平衡距离为零）
    {
        PotentialParams pot = make_potential_params(POT_MORSE, SHIFT_ENERGY,
                                                    50.0, 0, 50.0, 1.0, 7.4, 15.0);
        const double g = pair_force_over_r(pot, 7.4 * 7.4);
        check(std::fabs(g) < 1e-14, "J6 Morse: F(r_eq) = 0");
        const double U_eq = pair_energy(pot, 7.4 * 7.4);
        std::printf("      U(r_eq) = %.8f (期望 ≈ -49.95，能量移位后)\n", U_eq);
    }
    // WCA：U(rcut)=0 且 F→0
    {
        PotentialParams pot = make_potential_params(POT_WCA, SHIFT_ENERGY,
                                                    1.0, 7.4, 0, 0, 0, 0);
        const double rcut = pot.rcut;
        check(std::fabs(pair_energy(pot, rcut*rcut)) < 1e-14, "J6 WCA: U(rcut) = 0");
        check(std::fabs(pair_force_over_r(pot, rcut*rcut)) < 1e-14, "J6 WCA: F(rcut) → 0");
    }
    // POT_NONE：F 逐位等于外场，U ≡ 0
    {
        PotentialParams pot = make_potential_params(POT_NONE, SHIFT_ENERGY, 0, 0, 0, 0, 0, 0);
        ExternalField ext{1.0, -2.0, 3.0};
        double Rx[1]{5.0}, Ry[1]{6.0}, Rz[1]{7.0}, Fx[1], Fy[1], Fz[1];
        forces_cpu(pot, ext, 1, 32, 32, 32, Rx, Ry, Rz, Fx, Fy, Fz);
        const bool none_ok = (Fx[0] == 1.0 && Fy[0] == -2.0 && Fz[0] == 3.0);
        check(none_ok, "J6 POT_NONE: F 逐位等于外场，U ≡ 0");
    }
    // 标号交换对称性：N=3 任意构型，交换 0/1 标号后力按标号交换逐位相同
    {
        PotentialParams pot = make_potential_params(POT_LJ126, SHIFT_ENERGY,
                                                    57.1428571429, 7.4, 0, 0, 0, 15.0);
        ExternalField ext{0, 0, 0};
        double Rx[3]{5.0, 20.0, 3.0}, Ry[3]{1.0, 2.0, 25.0}, Rz[3]{10.0, 11.0, 12.0};
        double Fx[3], Fy[3], Fz[3];
        forces_cpu(pot, ext, 3, 32, 32, 32, Rx, Ry, Rz, Fx, Fy, Fz);

        double Rx2[3]{20.0, 5.0, 3.0}, Ry2[3]{2.0, 1.0, 25.0}, Rz2[3]{11.0, 10.0, 12.0};
        double Fx2[3], Fy2[3], Fz2[3];
        forces_cpu(pot, ext, 3, 32, 32, 32, Rx2, Ry2, Rz2, Fx2, Fy2, Fz2);

        const bool sym = (Fx[0] == Fx2[1] && Fy[0] == Fy2[1] && Fz[0] == Fz2[1] &&
                          Fx[1] == Fx2[0] && Fy[1] == Fy2[0] && Fz[1] == Fz2[0] &&
                          Fx[2] == Fx2[2] && Fy[2] == Fy2[2] && Fz[2] == Fz2[2]);
        check(sym, "J6 标号交换对称性（逐位）");
    }
}

// ---------------------------------------------------------------------------
// ΣF = 0：无外场时，CPU 参考实现的合力必须为零。
// 注：CPU 半矩阵版 F[i]+=g*d; F[j]-=g*d 让它在代数上恒等，几乎是盲的 ——
// 它主要抓「外场没加对」和「min_image 符号不对称」。真正的装配判据在 GPU
// 版（compute_particle_forces，全矩阵）下才有判别力，见 --check。
// ---------------------------------------------------------------------------
static void check_sum_zero()
{
    PotentialParams pot = make_potential_params(POT_LJ126, SHIFT_ENERGY,
                                                57.1428571429, 7.4, 0, 0, 0, 15.0);
    ExternalField ext{0, 0, 0};
    const int N = 8;
    double Rx[N], Ry[N], Rz[N], Fx[N], Fy[N], Fz[N];
    for (int i = 0; i < N; i++)
    {
        Rx[i] = (double)((i * 37) % 29) + 0.5;   // 确定性散布
        Ry[i] = (double)((i * 53) % 31) + 0.3;
        Rz[i] = (double)((i * 71) % 27) + 0.7;
    }
    forces_cpu(pot, ext, N, 32, 32, 32, Rx, Ry, Rz, Fx, Fy, Fz);
    double sx = 0, sy = 0, sz = 0, fmax = 0;
    for (int i = 0; i < N; i++)
    {
        sx += Fx[i]; sy += Fy[i]; sz += Fz[i];
        fmax = std::fmax(fmax, std::sqrt(Fx[i]*Fx[i]+Fy[i]*Fy[i]+Fz[i]*Fz[i]));
    }
    const double sum_mag = sqrt(sx*sx + sy*sy + sz*sz);
    const bool ok = sum_mag < 1e-11 * N * fmax;
    check(ok, "ΣF = 0（CPU 参考，N=8）");
    std::printf("      |ΣF| = %.3e   阈值 %.3e\n", sum_mag, 1e-11 * N * fmax);
}

// ============================================================================
int run_potential_check()
{
    std::printf("=== 势函数自检（纯 CPU，无卡可跑）===\n");

    // J1：三种势的力 = -dU/dr 对照
    {
        PotentialParams wca  = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 7.4, 0, 0, 0, 0);
        PotentialParams lj   = make_potential_params(POT_LJ126, SHIFT_ENERGY, 57.1428571429, 7.4, 0, 0, 0, 15.0);
        PotentialParams morse = make_potential_params(POT_MORSE, SHIFT_ENERGY, 50.0, 0, 50.0, 1.0, 7.4, 15.0);
        check_force_derivative(wca,   0.6 * 7.4, 1.4 * wca.rcut);
        check_force_derivative(lj,    0.6 * 7.4, 1.4 * 15.0);
        check_force_derivative(morse, 0.6 * 7.4, 1.4 * 15.0);
    }

    // J2：黄金表（数值来自 tests/spike/potential_ref.py）
    {
        const PotentialParams wca = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 7.4, 0, 0, 0, 0);
        const GoldenPt wca_pts[] = {
            {3.0,  202048.3759802687,  269997.16519427137},
            {6.0,  36.469882479364884, 14.169637469344766},
            {7.4,  1.0,                0.43827611395178956},
            {8.0,  0.063905867237099567, 0.059381618509055871},
            {8.30621915748936, 0.0, 0.0},
        };
        check_golden(wca, wca_pts, 5, "J2 WCA 黄金表对照");

        const PotentialParams lj = make_potential_params(POT_LJ126, SHIFT_ENERGY, 57.1428571429, 7.4, 0, 0, 0, 15.0);
        const GoldenPt lj_pts[] = {
            {8.0,  -50.243524349647913,  3.393235343377166},
            {8.30621915748936, -53.895288191770625, -2.2068694945765761e-15},
            {9.0,  -45.555331507188114, -1.9985929805915754},
            {12.0, -8.6308893152436692, -0.46613417631421977},
            {14.0, -1.6284795178803129, -0.14593895835225618},
            {15.0, 0.0, 0.0},
        };
        check_golden(lj, lj_pts, 6, "J2 LJ126 黄金表对照");

        const PotentialParams morse = make_potential_params(POT_MORSE, SHIFT_ENERGY, 50.0, 0, 50.0, 1.0, 7.4, 15.0);
        const GoldenPt morse_pts[] = {
            {7.0,  -37.855390718741418, 10.481660440731403},
            {7.4,  -49.949967379237798, 0.0},
            {8.0,  -39.771420393030347, -3.0952178022728041},
            {12.0, -0.95009898360930456, -0.082923302856637446},
            {14.0, -0.08591165293271083, -0.0097036959739307989},
            {15.0, 0.0, 0.0},
        };
        check_golden(morse, morse_pts, 6, "J2 Morse 黄金表对照");
    }

    // J3：最小镜像
    check_min_image();

    // J4：N=3 三角形
    check_n3_triangle();

    // J6：特征点 + 零判据 + 对称性
    check_feature_points();

    // ΣF = 0
    check_sum_zero();

    std::printf("\n%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}
