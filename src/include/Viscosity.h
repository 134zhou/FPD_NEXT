#ifndef VISCOSITY_H
#define VISCOSITY_H
#include <cufft.h>
#include <curand.h>
#include <cmath>

#include "Common.h"

void update_viscosity_fields(
    NS_Config cfg,
    double* sum_phi,
    int N, double* Rx, double* Ry, double* Rz,
    double* phi_grid, double* eta, double* etaXY, double* etaYZ, double* etaZX
);
#endif