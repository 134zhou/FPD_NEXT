#ifndef NEIGHBOR_H
#define NEIGHBOR_H
#include <cstddef>
#include "Potential.h"
// 每步重建的空间分箱索引，不是跨步缓存的 Verlet 邻居表。
// 设备工作区由 FpdState 持有；这里的指针均是 CUDA 设备地址，不是主机数组。
// 排序只改变索引，原粒子编号与物理数组保持不变；工作区不写入检查点。
// 内存 O(N+B)，B 为箱数；建表含 O(B) 清零和固定宽度整数基数排序。
// 力计算 O(候选对数)：固定密度/截断下通常近似 O(N)，聚集时可退化 O(N²)。
struct NeighborWorkspace
{
    // 三方向箱数及已分配的粒子容量。
    int nx=0, ny=0, nz=0, capacity=0;
    // 每个箱在 sorted 中的半开区间；空箱 begin=end=0。
    int *begin=0, *end=0;
    // 输入/输出键：(箱号 << 32) | 原粒子编号。
    unsigned long long *keys=0, *sorted=0;
    void* scratch=0;
    size_t scratch_bytes=0;
};
void neighbor_allocate(NeighborWorkspace& w, NS_Config cfg, int N, double rcut);
void neighbor_release(NeighborWorkspace& w);
// x/y/z 必须是设备地址，调用方先等待 OpenACC；此函数返回时 CUDA 已同步。
void neighbor_build(NeighborWorkspace w, NS_Config cfg, int N,
                    const double* Rx, const double* Ry, const double* Rz);
void compute_particle_forces_cell(NS_Config cfg, PotentialParams pot, WallParams wall,
                                 int N, const double* Rx, const double* Ry, const double* Rz,
                                 double* Fx, double* Fy, double* Fz, NeighborWorkspace w);
bool use_cell_search(int mode, int N, NS_Config cfg, PotentialParams pot);
#endif
