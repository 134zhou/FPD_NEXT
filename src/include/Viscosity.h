#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "Common.h"
#include "Stencil.h"

/**
 * @brief 由粒子位置构造粘度场，并计算三套错开的归一化因子
 *
 * 粘度：η(r) = η_ℓ + (η_c - η_ℓ) Σ_i φ_i(r)，代码里 η_ℓ = 1
 * 归一化：sum_phi{x,y,z}[n] = Σ_grid φ_n，分别在 FACE_X/Y/Z 上求和
 *
 * ⚠️ sum_phi{x,y,z} 必须与 Force 的力投影、Velocity 的速度平均用【同一个 Loc】，
 *    否则 Σ_grid f ≠ F_n。三者共用 stencil_point<L> 来保证这一点。
 */
void update_viscosity_fields(
    NS_Config cfg, PhiParams pp,
    int N, const double* Rx, const double* Ry, const double* Rz,
    double* sum_phix, double* sum_phiy, double* sum_phiz,
    double* eta, double* etaXY, double* etaYZ, double* etaZX
);

#endif
