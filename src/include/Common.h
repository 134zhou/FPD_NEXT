#ifndef COMMON_H
#define COMMON_H

#include <cmath>
// #define IDX(i, j, k) (((i) * Ny + (j)) * Nz + (k))
#define IDX(i, j, k) ((i) + (j)*Nx + (k)*Nx*Ny)

const double XI = 1.;
const double INV_XI = 1.;
const double RADIUS = 3.2;
const double RATIO_ETA = 50.0;

const int range = int(2.*(RADIUS+1./INV_XI));//相场函数的作用范围
const int range2 = range * range;
const int range_m1 = range-1;
const int N_range = 2*range;//range的两倍

struct NS_Config
{
    int Nx, Ny, Nz;
    double dt, inv_dt;
    double W; // 噪声强度系数
};

#pragma acc routine seq
static inline double order(double dx, double dy, double dz)
{
	return 0.5*(tanh((RADIUS - sqrt(dx*dx + dy*dy + dz*dz))*INV_XI) + 1.);
}

#endif