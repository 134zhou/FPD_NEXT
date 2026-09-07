#ifndef FORCE_H
#define FORCE_H

#include "Common.h"
#include "Stencil.h"

/**
 * @brief 把粒子受力 F_n 按相场权重摊布到流体力密度场
 *
 * f(r) = Σ_n φ_n(r) F_n / ∫φ_n dr  +  bg
 *
 * 分子分母共用 stencil_point<FACE_{X,Y,Z}>，因此 Σ_grid f_α == Σ_n F_α[n]（bg=0 时）
 * 是【恒等式】。sum_phi{x,y,z} 由 update_viscosity_fields 用同一个 Loc 算出。
 *
 * bgx/bgy/bgz 是【标量背景力密度】—— 外场均匀时 Σ_n F = N·g ≠ 0，流体 k=0 会被
 * 线性加速、整盒无界漂移。bg = −ΣF/size 抵消它（gravity_compensate）。加标量参数
 * 不新增数组，不触碰任何 present() 子句。
 */
void update_force_field(
    NS_Config cfg, PhiParams pp,
    int N, const double* Rx, const double* Ry, const double* Rz,
    const double* Fx, const double* Fy, const double* Fz,
    const double* sum_phix, const double* sum_phiy, const double* sum_phiz,
    double bgx, double bgy, double bgz,
    double* fx, double* fy, double* fz
);

#endif
