#ifndef MASS_H
#define MASS_H

#include "Common.h"

// ============================================================================
// 质量密度场（Mass）—— 变密度惯性 + 重力体力
//
// 论文 [SM] 2018 式 (3) 把 rho 放在物质导数前，并【明确假设】液体与粒子密度相同
// （"we assume that the density of the liquid is the same as that of the particles"）。
// 本项目原先把这条假设实现成了 rho == 1。本模块把它换成
//
//     rho(r) = rho_l + (rho_p - rho_l) * sum_i phi_i(r),      rho_l == 1
//     a(r)   = 1 / rho(r)
//
// rho_p / rho_l = ratio_rho；只传对比度 d_rho = ratio_rho - 1。
//
// 作用（两件事，互不重复）：
//   1) update_mass_fields：在三个面心（FACE_X/Y/Z）建 a = 1/rho。
//      位置、权重、截断、回绕全部复用 Stencil.h 的 stencil_point<FACE_*>，
//      与 Viscosity.cpp 算 sum_phi、Force.cpp 摊布力、Velocity.cpp 取平均【同一个 Loc】。
//      多粒子支撑域重叠时写网格必须 atomic（与 eta 场同一机制，也同一逐位重启边界）。
//   2) apply_gravity：把重力作【体力】加进面心力场，f_g = (rho - rho_l) g = (1/a - 1) g。
//      均匀项 rho_l*g 被静压解析减掉（Boussinesq 拆分，代数上精确）。于是净重力
//      自动是 (rho_p - rho_l) V_i g —— 浮力由密度差给出，不需要手填净力。
//
// ⚠️ 重力【不再】进 compute_particle_forces。粒子力只装势与壁面力；
//    重力的进入点只有本文件的 apply_gravity。
//
// 代码里没有「粒子质量」这个量：惯性完全由流体密度场 rho(r) 承担，
// 粒子速度只是流场的相场加权平均。
// ============================================================================

// 由粒子位置重建三个面心的 a = 1/rho。
// d_rho = 0 时三个场均逐位为 1.0（N = 0 也一样）。
void update_mass_fields(NS_Config cfg, PhiParams pp, double d_rho, int N,
                        const double* Rx, const double* Ry, const double* Rz,
                        double* inv_rho_x, double* inv_rho_y, double* inv_rho_z);

// 重力体力：f_alpha += (1/a_alpha - 1) * g_alpha。
// 上壁法向面（FACE_Z 的 k = Nz-1）a 恒 1，贡献恒 0；该面本来就被钉死。
void apply_gravity(NS_Config cfg, ExternalField g,
                   const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
                   double* fx, double* fy, double* fz);

#endif
