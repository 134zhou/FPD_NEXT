#include "Neighbor.h"
#include <cub/device/device_radix_sort.cuh>
#include <cuda_runtime.h>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#define CUDA_CHECK(x) do { cudaError_t e=(x); if(e!=cudaSuccess) { fprintf(stderr,"CUDA: %s (%s:%d)\n",cudaGetErrorString(e),__FILE__,__LINE__); std::abort(); } } while(0)

// FpdState 首次启用分箱时调用；N、盒长和截断在该工作区生命周期内不变。
// 工作区只分配一次，后续每步只重建内容；释放也由 FpdState 调用。
void neighbor_allocate(NeighborWorkspace& w, NS_Config cfg, int N, double rc)
{
    // n=floor(L/rc)，有多个箱时箱宽 L/n >= rc，因此只需本箱及两侧邻箱。
    // rc>L 时退化为单箱，虽然箱宽小于 rc，但所有粒子都已包含在这个箱中。
    w.nx = std::max(1, int(cfg.Nx/rc));
    w.ny = std::max(1, int(cfg.Ny/rc));
    w.nz = std::max(1, int(cfg.Nz/rc));
    w.capacity = N;
    size_t boxes = size_t(w.nx)*w.ny*w.nz;
    CUDA_CHECK(cudaMalloc(&w.begin, boxes*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&w.end, boxes*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&w.keys, size_t(N)*sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&w.sorted, size_t(N)*sizeof(unsigned long long)));
    // 查询并一次性分配 CUB 临时空间，避免每步排序触发设备内存分配。
    CUDA_CHECK(cub::DeviceRadixSort::SortKeys(nullptr,w.scratch_bytes,w.keys,w.sorted,N));
    CUDA_CHECK(cudaMalloc(&w.scratch,w.scratch_bytes));
}
void neighbor_release(NeighborWorkspace& w)
{
    CUDA_CHECK(cudaFree(w.scratch));
    CUDA_CHECK(cudaFree(w.sorted));
    CUDA_CHECK(cudaFree(w.keys));
    CUDA_CHECK(cudaFree(w.end));
    CUDA_CHECK(cudaFree(w.begin));
    w = NeighborWorkspace{};
}
// 一个线程生成一个粒子的排序键，不交换或复制 R/Ru/V 等物理量。
__global__ void make_keys(NeighborWorkspace w, NS_Config cfg, int N,
                          const double* x,const double* y,const double* z)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=N) return;
    int ix=min(w.nx-1,int(x[i]*w.nx/cfg.Nx));
    int iy=min(w.ny-1,int(y[i]*w.ny/cfg.Ny));
    // z 的物理范围是 [-0.5, Nz-0.5]；加 0.5 后转换到分箱坐标。
    int iz=min(w.nz-1,max(0,int((z[i]+0.5)*w.nz/cfg.Nz)));
    unsigned long long cell=ix+size_t(w.nx)*(iy+size_t(w.ny)*iz);
    // 高 32 位是箱号，低 32 位是原粒子号；要求两者均能用 32 位表示。
    // 完整键唯一，按键排序后同箱连续且箱内编号升序，不依赖原子插入顺序。
    w.keys[i]=(cell<<32)|static_cast<unsigned int>(i);
}
// 排序后只识别箱号变化处；每箱起点/终点各有唯一写者，无需 atomic。
__global__ void make_ranges(NeighborWorkspace w,int N)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=N) return;
    unsigned int c=w.sorted[i]>>32;
    if(i==0 || (w.sorted[i-1]>>32)!=c) w.begin[c]=i;
    if(i==N-1 || (w.sorted[i+1]>>32)!=c) w.end[c]=i+1;
}
void neighbor_build(NeighborWorkspace w, NS_Config cfg, int N,
                    const double* x,const double* y,const double* z)
{
    // 每步从当前 R 重建，不保留上一时刻邻居，也不需要 skin 或位移阈值。
    // 先把所有箱设为空区间，make_ranges 只覆盖非空箱的半开区间。
    CUDA_CHECK(cudaMemset(w.begin,0,size_t(w.nx)*w.ny*w.nz*sizeof(int)));
    CUDA_CHECK(cudaMemset(w.end,0,size_t(w.nx)*w.ny*w.nz*sizeof(int)));
    make_keys<<<(N+255)/256,256>>>(w,cfg,N,x,y,z);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cub::DeviceRadixSort::SortKeys(w.scratch,w.scratch_bytes,w.keys,w.sorted,N));
    make_ranges<<<(N+255)/256,256>>>(w,N);
    CUDA_CHECK(cudaGetLastError());
    // CUDA 默认流上的建键、排序、建区间依次完成后，才交回 OpenACC 力 kernel。
    CUDA_CHECK(cudaDeviceSynchronize());
}
