#include "Neighbor.h"
#include <algorithm>
#include <openacc.h>

// mode：0 自动、1 全遍历、2 强制分箱。只有粒子对势使用邻箱，壁面力不参与。
// auto 的保守规则：小体系或各方向都不超过三箱时，使用全遍历以省去建箱开销。
bool use_cell_search(int mode, int N, NS_Config cfg, PotentialParams pot)
{
    if (pot.type == POT_NONE || N < 2 || mode == 1) { return false; }
    if (mode == 2) { return true; }
    return N >= 10000 && (cfg.Nx/pot.rcut >= 4 || cfg.Ny/pot.rcut >= 4 || cfg.Nz/pot.rcut >= 4);
}

void compute_particle_forces_cell(NS_Config cfg, PotentialParams pot, WallParams wall,
                                 int N, const double* Rx, const double* Ry, const double* Rz,
                                 double* Fx, double* Fy, double* Fz, NeighborWorkspace w)
{
    // R 由 OpenACC 更新；先等待其完成，再把映射后的设备地址交给 CUDA。
    // neighbor_build 返回前同步 CUDA，保证下面的 OpenACC kernel 能读到完整箱表。
    acc_wait_all();
    neighbor_build(w,cfg,N,(const double*)acc_deviceptr((void*)Rx),
                   (const double*)acc_deviceptr((void*)Ry),(const double*)acc_deviceptr((void*)Rz));
    const int nx=w.nx, ny=w.ny, nz=w.nz;
    const int *begin=w.begin, *end=w.end;
    const unsigned long long* sorted=w.sorted;
    // R/F 是 OpenACC 映射数组，用 present；箱表由 cudaMalloc 分配，用 deviceptr。
    // 每个粒子一个 gang，只写自己的力：同一对分别由 i 和 j 计算，不用力 atomic。
    // 每个箱只访问一次；箱内按原编号排序，归约顺序不依赖排序的调度。
    #pragma acc parallel loop gang present(Rx,Ry,Rz,Fx,Fy,Fz) deviceptr(begin,end,sorted)
    for (int i=0; i<N; ++i)
    {
        // 与 make_keys 使用同一分箱坐标：x/y 从 0 开始，z 从下壁 -0.5 开始。
        // 边缘索引限制只用于选箱，不移动粒子，也不代替实际穿墙检查。
        int cx=int(Rx[i]*nx/cfg.Nx), cy=int(Ry[i]*ny/cfg.Ny);
        int cz=int((Rz[i]+0.5)*nz/cfg.Nz);
        cx=cx<nx?cx:nx-1; cy=cy<ny?cy:ny-1;
        cz=cz<0?0:(cz<nz?cz:nz-1);
        double sx=0, sy=0, sz=0;
        // 周期方向最多访问三个不同箱；一箱/两箱时减少循环次数，避免回绕后重计。
        // 例如 nx=2 时偏移只取 -1、0，+1 与 -1 是同一个箱。
        int ax=nx<3?nx:3, ay=ny<3?ny:3;
        for (int dz=-1; dz<=1; ++dz)
        {
            // z 不周期；越过上下壁的邻箱直接跳过。
            int z=cz+dz;
            if(z<0 || z>=nz) { continue; }
            for(int dy=0; dy<ay; ++dy)
            {
                int y=(cy+dy-1+ny)%ny;
                for(int dx=0; dx<ax; ++dx)
                {
                    int x=(cx+dx-1+nx)%nx;
                    int c=x+nx*(y+ny*z);
                    // [begin[c], end[c]) 是箱 c 在排序数组中的连续区间；空箱区间为空。
                    double tx=0,ty=0,tz=0;
                    #pragma acc loop vector reduction(+:tx,ty,tz)
                    for(int q=begin[c]; q<end[c]; ++q)
                    {
                        // 排序键低 32 位仍是原粒子编号，物理数组从未重排。
                        int j=(unsigned int)sorted[q];
                        if(j==i) { continue; }
                        double rx=min_image(Rx[i]-Rx[j],double(cfg.Nx));
                        double ry=min_image(Ry[i]-Ry[j],double(cfg.Ny));
                        double rz=Rz[i]-Rz[j];
                        // 邻箱只是候选集；球形截断仍由原势函数按 r² 精确判断。
                        double g=pair_force_over_r(pot,rx*rx+ry*ry+rz*rz);
                        tx+=g*rx; ty+=g*ry; tz+=g*rz;
                    }
                    // 每箱先做 vector 归约，再按固定邻箱顺序累加。
                    // 顺序区别于 direct，跨路径仅要求容差一致。
                    sx+=tx; sy+=ty; sz+=tz;
                }
            }
        }
        // 壁面力直接按表面间隙计算，与粒子对搜索分开。
        Fx[i]=sx; Fy[i]=sy; Fz[i]=sz+wall_force_z(wall,cfg,Rz[i]);
    }
}
