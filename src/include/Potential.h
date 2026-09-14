#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <cmath>
#include "Common.h"

// ============================================================================
// 粒子间势函数 + 外场 + 最小镜像 + 力装配
//
// 势函数用 enum + POD + routine seq 的 switch（决策记录：OpenACC 设备端不可靠，
// 不用函数指针）。WCA 不是独立分支：WCA ≡ LJ12-6 + rcut=2^{1/6}σ（派生）+ 移位，
// 且此时能量移位与力移位【恒等】（因 U'(rc)=0），所以 switch 里只有 3 个真实分支
// （NONE / LJ 形式 / MORSE）。这消除了旧代码 cul_WCA 与 cul_LJ6 两份复制粘贴的
// 漂移温床 —— 正是 C1/C2 的教训。
//
// ⚠️ pair_energy_bare 与 pair_force_over_r_bare 必须【独立实现，互不调用】。
//    一旦让 force 去调 energy（或反之），判据「F = -dU/dr 数值微分对照」就退化成
//    自洽的、盲的 —— 与 CLAUDE.md 记录的「力守恒判据对 stencil_point 内部公式
//    错误是盲的」同一个陷阱。
// ============================================================================

enum PotentialType { POT_NONE = 0, POT_WCA = 1, POT_MORSE = 2, POT_LJ126 = 3 };
enum ShiftMode     { SHIFT_NONE = 0, SHIFT_ENERGY = 1, SHIFT_FORCE = 2 };

// POD，按值传进 GPU kernel（沿用 PhiParams 的约定）
struct PotentialParams
{
    int    type;          // PotentialType
    int    shift;         // ShiftMode
    double eps, sigma;    // WCA / LJ12-6 的 ε 与 σ
    double De, alpha, r_eq;  // Morse 的阱深、宽度、平衡距离（宽度参数叫 alpha 不叫 a
                             // —— a 是粒子半径，见 README.md 的相场公式）
    double rcut, rcut2;
    double u_at_rc;       // 派生：U_bare(rcut)
    double dudr_at_rc;    // 派生：U'_bare(rcut)
};

// 每粒子的恒定外力（不是加速度）。gx/gy/gz 是三分量，为将来的壁面/沉降预留。
struct ExternalField { double gx, gy, gz; };

// 壁面排斥势的参数组。
//
// ⚠️ 半径 a 与势参数【必须打包在一起传】：间隙 h = (z + 1/2) - a 的定义依赖 a，
//    拆成两个形参就留出了「传了势参数、忘了传半径」的空间，而那是个静默错误
//    （力的大小对、平衡间隙错），力守恒判据抓不到。
struct WallParams
{
    PotentialParams pot;     // 复用粒子间势的同一套参数（POT_NONE 即关闭壁面势）
    double          radius;  // 粒子半径 a —— 间隙 h 的定义依赖它
};

// 由原始物理量构造，计算派生量（rcut2/u_at_rc/dudr_at_rc）。
// WCA 的 rcut = 2^{1/6}σ 是派生的，忽略传入的 rcut。
PotentialParams make_potential_params(int type, int shift,
                                      double eps, double sigma,
                                      double De, double alpha, double r_eq,
                                      double rcut);

// 最小镜像的【唯一】实现。前提：位置已折叠到 [0,L)，故单次折叠足够
//（由 Velocity.cpp 的 update_particle_position 保证）。
#pragma acc routine seq
static inline double min_image(double d, double L)
{
    if      (d >  0.5 * L) { return d - L; }
    else if (d < -0.5 * L) { return d + L; }
    return d;
}

// z 向最小镜像的长度。无滑移壁面下【最小镜像不作用于 z】。
//
// 实现上不新增分支到 O(N²) 内层：返回一个足够大的数（1e300），min_image 就自动
// 退化为恒等映射（0.5*1e300 不溢出，且 |d| 永远小于它）。这样「最小镜像只有一份
// 实现」这条不变量（PROGRESS.md 的关键决策）得以保住。
//
// ⚠️ 函数【必须保留】，不许内联成常量、也不许只在一侧删掉对它的调用：
//    CLAUDE.md 约束 6 要求 compute_particle_forces（GPU）与
//    compute_particle_forces_cpu 独立维护。若一侧保留了 z 映射、另一侧删了，
//    CPU/GPU 对照判据【仍然通过】（两边都是恒等映射，判据是盲的），
//    但「min_image 唯一真值源」的不变量没了，下一次编辑就会漂开。
#pragma acc routine seq
static inline double pbc_length_z(NS_Config cfg)
{
    (void)cfg;
    return 1.0e300;
}

// ---------------------------------------------------------------------------
// 裸势 / 裸力（无截断、无移位）。g_bare(r) = -U'_bare(r) / r，
// 调用点 F_i = g * (R_i - R_j)。g > 0 表示排斥。
// ---------------------------------------------------------------------------
#pragma acc routine seq
static inline double pair_energy_bare(PotentialParams pp, double r)
{
    switch (pp.type)
    {
        case POT_LJ126:
        case POT_WCA:
        {
            double s2 = pp.sigma / r; s2 *= s2;      // (σ/r)²
            const double s6 = s2 * s2 * s2;          // (σ/r)⁶
            return 4.0 * pp.eps * (s6 * s6 - s6);
        }
        case POT_MORSE:
        {
            const double E = exp(-pp.alpha * (r - pp.r_eq));
            return pp.De * (E * E - 2.0 * E);
        }
        default: return 0.0;
    }
}

#pragma acc routine seq
static inline double pair_force_over_r_bare(PotentialParams pp, double r)
{
    switch (pp.type)
    {
        case POT_LJ126:
        case POT_WCA:
        {
            double s2 = pp.sigma / r; s2 *= s2;      // (σ/r)²
            const double s6 = s2 * s2 * s2;          // (σ/r)⁶
            return 24.0 * pp.eps * s6 * (2.0 * s6 - 1.0) / (r * r);
        }
        case POT_MORSE:
        {
            const double E = exp(-pp.alpha * (r - pp.r_eq));
            return 2.0 * pp.alpha * pp.De * E * (E - 1.0) / r;
        }
        default: return 0.0;
    }
}

// ---------------------------------------------------------------------------
// 热路径接口（含截断 + 移位），传 r² 避免截断外还要开方。
// ---------------------------------------------------------------------------
#pragma acc routine seq
static inline double pair_energy(PotentialParams pp, double r2)
{
    if (r2 >= pp.rcut2) { return 0.0; }
    const double r = sqrt(r2);
    const double U = pair_energy_bare(pp, r);
    if (pp.shift == SHIFT_NONE)   { return U; }
    if (pp.shift == SHIFT_ENERGY) { return U - pp.u_at_rc; }
    // SHIFT_FORCE：U - U(rc) - (r-rc)U'(rc)
    return U - pp.u_at_rc - (r - pp.rcut) * pp.dudr_at_rc;
}

#pragma acc routine seq
static inline double pair_force_over_r(PotentialParams pp, double r2)
{
    if (r2 >= pp.rcut2 || r2 <= 0.0) { return 0.0; }
    const double r = sqrt(r2);
    double g = pair_force_over_r_bare(pp, r);
    if (pp.shift == SHIFT_FORCE) { g += pp.dudr_at_rc / r; }   // F - F(rc)
    return g;
}

// ---------------------------------------------------------------------------
// 壁面排斥势
//
// 无滑移壁面在 z = -1/2（下）与 z = Nz-1/2（上）。本模块用的是【表面到壁面的
// 间隙】而不是「中心到壁面平面」的距离：
//
//     h_bottom = (z + 1/2) - a
//     h_top    = (Nz - 1/2 - z) - a
//
// 选间隙的理由是【参数独立性】：改 radius 时平衡间隙不变；用中心距则会静默漂移。
//
// 势函数族【逐字复用】粒子间势的 pair_force_over_r（在 h 坐标上），所以
// --check-potential 的 J1「力 vs 能量的四阶数值微分」判据自动覆盖本模块的公式。
// 壁面的【几何装配】（h 的定义、符号、哪面墙）由判据 W8a/W8b 单独判。
//
// ⚠️ 这是【模型】，不是从 FPD 第一性原理推出的壁面自由能。它的职责是「让粒子停在
//    一个可控的间隙上」。判据管的是实现对不对，不管这个模型是不是某个真实壁面势。
//
// ⚠️ 前置条件 h > 0 由【调用方】保证，本函数不兜底、不 clamp —— 与 Velocity.cpp
//    的「z 向从不折叠」同一立场：clamp 会把穿墙变成一个看起来正常的轨迹。
//    这不是洁癖：h < 0 时 h² 仍然是正的，force_over_r(h²)·h 会给出【符号反向】的
//    力（把粒子往壁里推），是一条无诊断的错物理路径。调用方（main.cpp 的
//    check_wall_bounds）必须先中止。
// ---------------------------------------------------------------------------
#pragma acc routine seq
static inline double wall_gap_bottom(double a, double z)
{
    return (z + 0.5) - a;
}

#pragma acc routine seq
static inline double wall_gap_top(NS_Config cfg, double a, double z)
{
    return ((double)cfg.Nz - 0.5 - z) - a;
}

// 粒子受两面壁的 z 向合力（> 0 指向 +z）。
// 两侧各自走 pair_force_over_r 的 r² >= rcut² 守卫 ⇒ 超出截断恒 0。
#pragma acc routine seq
static inline double wall_force_z(WallParams wp, NS_Config cfg, double z)
{
    if (wp.pot.type == POT_NONE) { return 0.0; }
    const double hb = wall_gap_bottom(wp.radius, z);
    const double ht = wall_gap_top(cfg, wp.radius, z);
    // 下壁把粒子往 +z 推，上壁往 -z 推。
    return pair_force_over_r(wp.pot, hb * hb) * hb
         - pair_force_over_r(wp.pot, ht * ht) * ht;
}

// ---------------------------------------------------------------------------
// 关闭壁面势的 WallParams。判据里不关心壁面时显式用它 —— 好过就地构造一个有
// 隐患的默认值（POD 无构造函数，聚合初始化容易漏字段）。
static inline WallParams no_wall()
{
    WallParams w;
    w.pot    = make_potential_params(POT_NONE, SHIFT_ENERGY, 0, 0, 0, 0, 0, 0);
    w.radius = 0.0;
    return w;
}

// ---------------------------------------------------------------------------
// 粒子间力装配。GPU 版：gang-over-i + vector-reduction-over-j 全矩阵（零 atomic）。
// 每个粒子独立算自己的力，牛三靠 g(r) 的对称性 + 全矩阵覆盖保证。
//
// 求和顺序在固定 N 与固定 gang/vector 配置下是确定的 —— 与 Viscosity.cpp:40 /
// Velocity.cpp:21 算 sum_phi 和 V 是同一机制（spike_pair_bench 已实测逐位可复现）。
//
// 壁面力并入同一处装配（不另开 kernel）：这样 --verify-forces / J5 / J7 自动覆盖
// 新的力项，不需要各自改判据。代价是两个实现各自多一行，符合 CLAUDE.md 约束 6
// 对「两份实现独立维护」的要求（势函数体共享，装配各写各的）。
// ---------------------------------------------------------------------------
void compute_particle_forces(NS_Config cfg, PotentialParams pot, ExternalField ext,
                             WallParams wall,
                             int N,
                             const double* Rx, const double* Ry, const double* Rz,
                             double* Fx, double* Fy, double* Fz);

// CPU 参考实现：串行 i<j 半矩阵。与 GPU 版是两份独立实现（各自验证装配逻辑），
// 供纯 CPU 判据（CheckPotential.cpp）与 fpd_tool 的 --verify-forces 使用。
void compute_particle_forces_cpu(NS_Config cfg, PotentialParams pot, ExternalField ext,
                                 WallParams wall,
                                 int N,
                                 const double* Rx, const double* Ry, const double* Rz,
                                 double* Fx, double* Fy, double* Fz);

#endif
