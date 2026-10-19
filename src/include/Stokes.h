#ifndef STOKES_H
#define STOKES_H

#include <cufft.h>
#include <curand.h>
#include <cmath>

#include "Common.h"
#include "Check.h"
#include "Wall.h"


// ---------------------------------------------------------------------------
// 【逐位重启的关键】随机数流的位置由 step 直接算出，不做任何累计记账。
//
// 实测结论（tests/spike/spike_rng_offset.cpp、spike_rng_slice.cpp、spike_rng_advance.cpp）：
//   1. 只有 CURAND_RNG_PSEUDO_PHILOX4_32_10 支持有意义的 setGeneratorOffset；
//      XORWOW（原 CURAND_RNG_PSEUDO_DEFAULT）和 MRG32K3A 都对不上，无法逐位重启
//   2. 「生成 n 个数后序列前进 n」这个模型是【错的】：offset=1*n 能对上第 2 批，
//      但 2*n、3*n 都对不上连续生成的对应批次。所以累计记账不可行
//   3. 但同一 offset 生成同样数量必定逐位可复现，且不同 offset slot 之间
//      零共同值、互相关 ~1/sqrt(n)（即统计独立）
//
// 于是采用 slot 方案：每步占一段长 slotD + slotN 的随机数流，生成前显式定位：
//   base = step * (slotD + slotN)
//   randD（体心对角噪声）从 base 开始，randN（棱边剪切噪声）从 base + slotD 开始
// slotD = 3*Nx*Ny*Nz、slotN = 3*Nx*Ny*(Nz+1)，均向上取偶（见 Wall.h 的 wz_rand_slot_*）。
// 重启时按同一公式算出同一 offset，逐位相等是构造上保证的，与跑了多少步无关。
// ---------------------------------------------------------------------------
// 返回变系数路径的 PCG 迭代数；常系数路径（density_on=0）返回 0。
//
// 密度参数：inv_rho_x/y/z 是面心 1/rho（Mass.cpp 建；density_on=0 时恒 1.0）。
// ⚠️ 1/rho 只乘【力密度】（粘性、噪声、外力、压力梯度），不乘对流加速度：
//    显式步拆成 d_adv（不乘 a）与 d_visc = d_pix - d_adv（乘 a），修正步乘 a。
//    density_on=1 时压力方程解 div(a grad p)=b；判据是 --check-poisson 的 V0。
// 三个 cg_* 是 density_on=1 才分配的 PCG 暂存（=0 时可为 null，不被读取）。
int step_navier_stokes(
    NS_Config cfg,
    double* vx, double* vy, double* vz,
    double* p,
    double* fx, double* fy, double* fz,
    const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
    double* cg_r, double* cg_s, double* cg_zq,
    double* eta,
    double* etaXY, double* etaYZ, double* etaZX,
    double* pi_dx, double* pi_dy, double* pi_dz,
    double* pi_nx, double* pi_ny, double* pi_nz,
    double* fft_data,
    cufftHandle plan_xy,                   // xy 批量 2D FFT 计划（batch = Nz）
    const double* tri_w,                   // z 向三对角前推系数
    double* diag,                          // [0]=|Σb̂|，[1]=PCG 最终相对残差（长度 2）
    curandGenerator_t gen,
    double* randD, double* randN,
    double* tmp_fx, double* tmp_fy, double* tmp_fz,
    long step                              // 决定随机数流的 slot 位置
);

#endif