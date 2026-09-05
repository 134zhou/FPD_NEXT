#ifndef STOKES_H
#define STOKES_H

#include <cufft.h>
#include <curand.h>
#include <cmath>

#include "Common.h"
#include "Check.h"


// ---------------------------------------------------------------------------
// 【逐位重启的关键】随机数流的位置由 step 直接算出，不做任何累计记账。
//
// 实测结论（spike/spike_rng_offset.cpp、spike_rng_slice.cpp、spike_rng_advance.cpp）：
//   1. 只有 CURAND_RNG_PSEUDO_PHILOX4_32_10 支持有意义的 setGeneratorOffset；
//      XORWOW（原 CURAND_RNG_PSEUDO_DEFAULT）和 MRG32K3A 都对不上，无法逐位重启
//   2. 「生成 n 个数后序列前进 n」这个模型是【错的】：offset=1*n 能对上第 2 批，
//      但 2*n、3*n 都对不上连续生成的对应批次。所以累计记账不可行
//   3. 但同一 offset 生成同样数量必定逐位可复现，且不同 offset slot 之间
//      零共同值、互相关 ~1/sqrt(n)（即统计独立）
//
// 于是采用 slot 方案：每次生成前显式 setGeneratorOffset 到
//   offset(step, c) = (2*step + c) * (size*3),   c ∈ {0,1}
// 重启时按同一公式算出同一 offset，逐位相等是构造上保证的，与跑了多少步无关。
// ---------------------------------------------------------------------------
void step_navier_stokes(
    NS_Config cfg,
    double* vx, double* vy, double* vz,
    double* p,
    double* fx, double* fy, double* fz,
    double* eta,
    double* etaXY, double* etaYZ, double* etaZX,
    double* pi_dx, double* pi_dy, double* pi_dz,
    double* pi_nx, double* pi_ny, double* pi_nz,
    double* fft_data,
    cufftHandle plan,
    curandGenerator_t gen,
    double* randD, double* randN,
    double* tmp_fx, double* tmp_fy, double* tmp_fz,
    long step                              // 决定随机数流的 slot 位置
);

#endif