#include "Neighbor.h"
#include <algorithm>
#include <openacc.h>

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
    acc_wait_all();
    neighbor_build(w,cfg,N,(const double*)acc_deviceptr((void*)Rx),
                   (const double*)acc_deviceptr((void*)Ry),(const double*)acc_deviceptr((void*)Rz));
    const int nx=w.nx, ny=w.ny, nz=w.nz;
    const int *begin=w.begin, *end=w.end;
    const unsigned long long* sorted=w.sorted;
    // 每个箱只访问一次；箱内按原编号排序，归约顺序不依赖排序的调度。
    #pragma acc parallel loop gang present(Rx,Ry,Rz,Fx,Fy,Fz) deviceptr(begin,end,sorted)
    for (int i=0; i<N; ++i)
    {
        int cx=int(Rx[i]*nx/cfg.Nx), cy=int(Ry[i]*ny/cfg.Ny);
        int cz=int((Rz[i]+0.5)*nz/cfg.Nz);
        cx=cx<nx?cx:nx-1; cy=cy<ny?cy:ny-1;
        cz=cz<0?0:(cz<nz?cz:nz-1);
        double sx=0, sy=0, sz=0;
        int ax=nx<3?nx:3, ay=ny<3?ny:3;
        for (int dz=-1; dz<=1; ++dz)
        {
            int z=cz+dz;
            if(z<0 || z>=nz) { continue; }
            for(int dy=0; dy<ay; ++dy)
            {
                int y=(cy+dy-1+ny)%ny;
                for(int dx=0; dx<ax; ++dx)
                {
                    int x=(cx+dx-1+nx)%nx;
                    int c=x+nx*(y+ny*z);
                    double tx=0,ty=0,tz=0;
                    #pragma acc loop vector reduction(+:tx,ty,tz)
                    for(int q=begin[c]; q<end[c]; ++q)
                    {
                        int j=(unsigned int)sorted[q];
                        if(j==i) { continue; }
                        double rx=min_image(Rx[i]-Rx[j],double(cfg.Nx));
                        double ry=min_image(Ry[i]-Ry[j],double(cfg.Ny));
                        double rz=Rz[i]-Rz[j];
                        double g=pair_force_over_r(pot,rx*rx+ry*ry+rz*rz);
                        tx+=g*rx; ty+=g*ry; tz+=g*rz;
                    }
                    sx+=tx; sy+=ty; sz+=tz;
                }
            }
        }
        Fx[i]=sx; Fy[i]=sy; Fz[i]=sz+wall_force_z(wall,cfg,Rz[i]);
    }
}
