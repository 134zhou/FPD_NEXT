#ifndef POISSON_H
#define POISSON_H

#include <cufft.h>
#include "Common.h"

// ============================================================================
// 压力泊松求解器（Phase 7-A）
//
// 【唯一路径】z 向无滑移壁面：xy 向批量 2D FFT（cufftPlanMany，batch=Nz）
// + 对每个 (kx,ky) 模式在 z 向解 Nz×Nz 三对角（Thomas）。
//
// 数学推导全文见 doc/PressurePoisson.md。文档与代码【共同定义】同一个离散算子
// （div∘grad），改任一侧必须同步更新另一侧并重跑 --check-poisson。
// 变系数路径（div(a grad p)）的推导、PCG 与预条件器见同文档 §16。
// ============================================================================

// 三对角前推系数预算（主机侧，一次性）。含 (0,0) 奇异列的定规 d_0 -= 1。
// tri_w 长度 Nx*Ny*Nz，tri_w[IDX(i,j,k)] = 第 (i,j) 列的第 k 个前推系数 w_k。
// w_k 只依赖 λ_xy 与 k，与右端无关，故只需预算一次。
void build_tridiag_coeffs(NS_Config cfg, double* tri_w);

// 唯一入口。就地求解：把 fft_data（复数交错，x 最快）变成压力谱 →
// 逆变换 → 归一化 → 写入 p。
//
//   plan_xy：xy 向 2D 批量 FFT 计划（batch = Nz）。
//   tri_w：z 向三对角前推系数。
//   diag：指向设备侧长度 ≥1 的数组，写入 (0,0) 列的相容性残差
//         |Σ_k b̂_k|（浮点下 ~1e-13，用作诊断）。**必须给，不可为 null**。
//         常系数路径只写 diag[0]；变系数路径另写 diag[1]（见下）。
//
// 归一化因子只在本函数内部出现一次，杜绝「照抄 /size 却该除 Nx*Ny」的漂移。
void solve_pressure(NS_Config cfg, double* fft_data, const double* tri_w,
                    cufftHandle plan_xy, double* p, double* diag);

// ---------------------------------------------------------------------------
// 变系数投影（ratio_rho != 1 时走这条路径）
//
// 动量方程带 rho(r) 后，投影步是
//     div( a grad p ) = (1/dt) div v*,     a = 1/rho（面心）
//     v = v* - dt * a * grad p
// 三处（算子、梯度、散度）必须用同一套面心 a 与同一套离散。这是 --check-poisson
// 的 V0-V5 与 doc/PressurePoisson.md §16 共同定义的不变量。
//
// apply_laplacian_var 算的是 L p = div(a grad p)，L 对称半负定，零空间 = 常向量。
// 离散：x/y 周期回绕；k=0 只保留 +a_z[0]*(p1-p0)，k=Nz-1 只保留
// -a_z[Nz-2]*(p-p_{Nz-2})（Neumann ghost 使壁面通量为 0）。
// ---------------------------------------------------------------------------
void apply_laplacian_var(NS_Config cfg,
                         const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
                         const double* p, double* out);

// PCG 解 L p = b（等价地 A p = -b，A = -L），预条件器 B = -a_bar * M
// （M 是上面的常系数 Laplacian，a_bar 是实际参与算子的面心 a 的均值；
//  a_bar == 1 时退化成上面那个常系数求解器，直接复用 solve_pressure，
//  不复制 FFT/Thomas。a_bar 缩放对稠密体系把迭代数降一个量级）。
//
//   cg_r：入口装 b，出口是残差（就地覆盖）。
//   cg_s / cg_zq：调用方提供的 size 暂存。
//   diag：[0] 写 |Σb|（相容性），[1] 写最终相对残差 ||r||/||b||。
//   返回值：迭代数（0 表示 b == 0）。
//
// 常数零空间：r0 与最终 p 都显式去均值。预条件器输出自带的任意常数【不】逐次投影 ——
// 它被 A 零化、且与零均值的 r 正交，对 CG 递推无影响（只给 p 加常数，最后去掉）。
//
// 调试用 trace（生产传 nullptr，走默认 tol/max_iter，算术路径逐位不变）：
//   p_ref 是设备上的零均值参考解；res 存 ||r||/||b||，err 存去均值后的
//   ||p-p_ref||（未归一化，调用方按 ||p_ref|| 归一化）。仅用于
//   --debug-pcg 测「多少步到双精度地板」，不参与生产求解。
struct PcgTrace
{
    const double* p_ref = 0;   // device，零均值；null 则不记误差
    double* res = 0;           // host，长度 cap
    double* err = 0;           // host，长度 cap
    int     cap = 0;
    int     n   = 0;           // 出口：记录的迭代数
    double  tol = 1e-12;       // 覆盖生产 tol；调试设 0
    int     max_iter = 200;    // 覆盖生产上限；调试设 300
};

int solve_pressure_var(NS_Config cfg,
                       double* cg_r,
                       const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
                       double* fft_data, const double* tri_w, cufftHandle plan_xy,
                       double* p, double* diag,
                       double* cg_s, double* cg_zq,
                       PcgTrace* trace = 0);

#endif
