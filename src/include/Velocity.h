#ifndef VELOCITY_H
#define VELOCITY_H

#include <cufft.h>
#include <curand.h>
#include <cmath>

#include "Common.h"

void update_particale_VandR
(
    NS_Config cfg,
    double* phi_grid, double* sum_phi,
    int N, double* Rx, double* Ry, double* Rz,
    double* Vx, double* Vy, double* Vz,
    double* vx, double* vy, double* vz
);

#endif