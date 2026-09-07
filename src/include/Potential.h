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
// 粒子间力装配。GPU 版：gang-over-i + vector-reduction-over-j 全矩阵（零 atomic）。
// 每个粒子独立算自己的力，牛三靠 g(r) 的对称性 + 全矩阵覆盖保证。
//
// 求和顺序在固定 N 与固定 gang/vector 配置下是确定的 —— 与 Viscosity.cpp:40 /
// Velocity.cpp:21 算 sum_phi 和 V 是同一机制，逐位可复现（spike_pair_bench 已实测）。
// ---------------------------------------------------------------------------
void compute_particle_forces(NS_Config cfg, PotentialParams pot, ExternalField ext,
                             int N,
                             const double* Rx, const double* Ry, const double* Rz,
                             double* Fx, double* Fy, double* Fz);

#endif
