#ifndef VELOCITY_H
#define VELOCITY_H

#include "Common.h"
#include "Stencil.h"

/**
 * @brief 由流体速度场的相场加权平均得到粒子速度
 *
 * V_i = ∫ v φ_i dr / ∫ φ_i dr
 *
 * ⚠️ 权重必须是【粒子 i 自己的】φ_i，不能用全场求和的 Σ_n φ_n（这是缺陷 C2，
 *    N=1 时恰好抵消所以曾长期未被发现）。stencil_point 现场重算 φ_i，
 *    因此不再需要、也不再存在 phi_grid 这个场。
 */
void update_particle_velocity(
    NS_Config cfg, PhiParams pp,
    int N, const double* Rx, const double* Ry, const double* Rz,
    const double* sum_phix, const double* sum_phiy, const double* sum_phiz,
    const double* vx, const double* vy, const double* vz,
    double* Vx, double* Vy, double* Vz
);

/**
 * @brief 推进粒子位置
 *
 * R 施加周期边界；Ru 是【不折叠】的位置，供 MSD 使用，二者用同一个 V
 * 在同一个 kernel 里积分，避免累积漂移。
 */
void update_particle_position(
    NS_Config cfg,
    int N,
    double* Rx, double* Ry, double* Rz,
    double* Rux, double* Ruy, double* Ruz,
    const double* Vx, const double* Vy, const double* Vz
);

#endif
