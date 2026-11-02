#include "Neighbor.h"
#include <cub/device/device_radix_sort.cuh>
#include <cuda_runtime.h>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#define CUDA_CHECK(x) do { cudaError_t e=(x); if(e!=cudaSuccess) { fprintf(stderr,"CUDA: %s (%s:%d)\n",cudaGetErrorString(e),__FILE__,__LINE__); std::abort(); } } while(0)

void neighbor_allocate(NeighborWorkspace& w, NS_Config cfg, int N, double rc)
{
    w.nx = std::max(1, int(cfg.Nx/rc));
    w.ny = std::max(1, int(cfg.Ny/rc));
    w.nz = std::max(1, int(cfg.Nz/rc));
    w.capacity = N;
    size_t boxes = size_t(w.nx)*w.ny*w.nz;
    CUDA_CHECK(cudaMalloc(&w.begin, boxes*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&w.end, boxes*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&w.keys, size_t(N)*sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&w.sorted, size_t(N)*sizeof(unsigned long long)));
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
__global__ void make_keys(NeighborWorkspace w, NS_Config cfg, int N,
                          const double* x,const double* y,const double* z)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=N) return;
    int ix=min(w.nx-1,int(x[i]*w.nx/cfg.Nx));
    int iy=min(w.ny-1,int(y[i]*w.ny/cfg.Ny));
    int iz=min(w.nz-1,max(0,int((z[i]+0.5)*w.nz/cfg.Nz)));
    unsigned long long cell=ix+size_t(w.nx)*(iy+size_t(w.ny)*iz);
    w.keys[i]=(cell<<32)|static_cast<unsigned int>(i);
}
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
    CUDA_CHECK(cudaMemset(w.begin,0,size_t(w.nx)*w.ny*w.nz*sizeof(int)));
    CUDA_CHECK(cudaMemset(w.end,0,size_t(w.nx)*w.ny*w.nz*sizeof(int)));
    make_keys<<<(N+255)/256,256>>>(w,cfg,N,x,y,z);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cub::DeviceRadixSort::SortKeys(w.scratch,w.scratch_bytes,w.keys,w.sorted,N));
    make_ranges<<<(N+255)/256,256>>>(w,N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}
