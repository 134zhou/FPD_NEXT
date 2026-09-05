#ifndef STOKES_H
#define STOKES_H

#include <cufft.h>
#include <curand.h>
#include <cmath>

#include "Common.h"


// 函数声明
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
    double* tmp_fx, double* tmp_fy, double* tmp_fz
);

#endif