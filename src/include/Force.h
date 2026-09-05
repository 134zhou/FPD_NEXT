#ifndef FORCE_H
#define FORCE_H

#include <cufft.h>
#include <curand.h>
#include <cmath>
#include "Common.h"

void update_force_field(
    NS_Config cfg,
    double* sum_phi,
    int N, double* Rx, double* Ry, double* Rz,
    double* Fx, double* Fy, double* Fz,
    double* fx, double* fy, double* fz
);

#endif