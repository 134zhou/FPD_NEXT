#ifndef POISSON_H
#define POISSON_H

#include <cufft.h>
#include "Common.h"

// ============================================================================
// 压力泊松求解器（Phase 7-A）
//
// 两条路径，按 cfg.wall_z 分派：
//   wall_z = 0（周期）：3D FFT 全周期，除以精确离散本征值
//               2(cos kx + cos ky + cos kz - 3)，与 Stokes.cpp 的有限差分离散自洽。
//   wall_z = 1（z 向无滑移壁面）：xy 向批量 2D FFT（cufftPlanMany，batch=Nz）
//               + 对每个 (kx,ky) 模式在 z 向解 Nz×Nz 三对角（Thomas）。
//
// 数学推导全文见 doc/PressurePoisson.md。二者【共同定义】同一个离散算子
// （div∘grad），改任一侧必须同步更新另一侧并重跑 --check-poisson。
// ============================================================================

// 三对角前推系数预算（主机侧，一次性）。含 (0,0) 奇异列的定规 d_0 -= 1。
// tri_w 长度 Nx*Ny*Nz，tri_w[IDX(i,j,k)] = 第 (i,j) 列的第 k 个前推系数 w_k。
// w_k 只依赖 λ_xy 与 k，与右端无关，故只需预算一次。
void build_tridiag_coeffs(NS_Config cfg, double* tri_w);

// 唯一入口：按 cfg.wall_z 分派。就地求解：把 fft_data（复数交错，x 最快）
// 变成压力谱 → 逆变换 → 归一化 → 写入 p。
//
//   plan3d：周期路径的 3D FFT 计划（wall_z=0 时用，否则忽略）。
//   plan_xy：壁面路径的 2D 批量 FFT 计划（wall_z=1 时用，否则忽略）。
//   tri_w：壁面路径的三对角系数（wall_z=1 时用，否则忽略）。
//   diag：壁面模式下指向设备侧长度 ≥1 的数组，写入 (0,0) 列的相容性残差
//         |Σ_k b̂_k|（浮点下 ~1e-13，用作诊断）。周期模式不写。
//
// 归一化因子只在本函数内部出现一次，杜绝「照抄 /size 却该除 Nx*Ny」的漂移。
void solve_pressure(NS_Config cfg, double* fft_data, const double* tri_w,
                    cufftHandle plan3d, cufftHandle plan_xy,
                    double* p, double* diag);

#endif
