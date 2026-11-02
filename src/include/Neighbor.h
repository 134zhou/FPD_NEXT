#ifndef NEIGHBOR_H
#define NEIGHBOR_H
#include <cstddef>
#include "Potential.h"
// 设备工作区由 FpdState 持有，排序只改变索引。
struct NeighborWorkspace
{
    int nx=0, ny=0, nz=0, capacity=0;
    int *begin=0, *end=0;
    unsigned long long *keys=0, *sorted=0;
    void* scratch=0;
    size_t scratch_bytes=0;
};
void neighbor_allocate(NeighborWorkspace& w, NS_Config cfg, int N, double rcut);
void neighbor_release(NeighborWorkspace& w);
void neighbor_build(NeighborWorkspace w, NS_Config cfg, int N,
                    const double* Rx, const double* Ry, const double* Rz);
void compute_particle_forces_cell(NS_Config cfg, PotentialParams pot, WallParams wall,
                                 int N, const double* Rx, const double* Ry, const double* Rz,
                                 double* Fx, double* Fy, double* Fz, NeighborWorkspace w);
bool use_cell_search(int mode, int N, NS_Config cfg, PotentialParams pot);
#endif
