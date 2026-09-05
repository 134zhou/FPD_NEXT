#ifndef FORCE_H
#define FORCE_H

#include "Common.h"
#include "Stencil.h"

/**
 * @brief 把粒子受力 F_n 按相场权重摊布到流体力密度场
 *
 * f(r) = Σ_n φ_n(r) F_n / ∫φ_n dr
 *
 * 分子分母共用 stencil_point<FACE_{X,Y,Z}>，因此 Σ_grid f_α == F_α[n]
 * 是【恒等式】。sum_phi{x,y,z} 由 update_viscosity_fields 用同一个 Loc 算出。
 */
void update_force_field(
    NS_Config cfg, PhiParams pp,
    int N, const double* Rx, const double* Ry, const double* Rz,
    const double* Fx, const double* Fy, const double* Fz,
    const double* sum_phix, const double* sum_phiy, const double* sum_phiz,
    double* fx, double* fy, double* fz
);

#endif
