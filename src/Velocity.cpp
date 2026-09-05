#include "./include/Velocity.h"

void update_particale_VandR
(
    NS_Config cfg,
    double* phi_grid, double* sum_phi,
    int N, double* Rx, double* Ry, double* Rz,
    double* Vx, double* Vy, double* Vz,
    double* vx, double* vy, double* vz
)
{
    int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    double DT = cfg.dt;


    #pragma acc parallel loop collapse(4) present(phi_grid, Rx, Ry, Rz, Vx, Vy, Vz, vx, vy, vz)
    for(int n=0;n< N;n++)
    {
        for (int li = 0; li < N_range; li++)
        {
            for (int lj = 0; lj < N_range; lj++)
            {
                for (int lk = 0; lk < N_range; lk++)
                {
                    double Rnx = Rx[n]; double Rny = Ry[n]; double Rnz = Rz[n];
                    int in = (int)Rnx; int jn = (int)Rny; int kn = (int)Rnz;

                    // li, lj, lk 是局部范围索引，映射到全局格点
                    // 映射到全局周期性网格
                    int ir = li + in - range_m1;
                    int jr = lj + jn - range_m1;
                    int kr = lk + kn - range_m1;

                    int irP = (ir + Nx) % Nx;
                    int jrP = (jr + Ny) % Ny;
                    int krP = (kr + Nz) % Nz;
                    int ijkP = IDX(irP, jrP, krP);

                    double temp_phi = phi_grid[ijkP];

                    #pragma acc atomic update
                    Vx[n] += vx[ijkP] * temp_phi;

                    #pragma acc atomic update
                    Vy[n] += vy[ijkP] * temp_phi;

                    #pragma acc atomic update
                    Vz[n] += vz[ijkP] * temp_phi;
                }
            }
        }
    }


    #pragma acc parallel loop collapse(1) present(Vx, Vy, Vz, sum_phi)
    for (int n = 0; n < N; n++)
    {
        Vx[n] = Vx[n] / sum_phi[n];
        Vy[n] = Vy[n] / sum_phi[n];
        Vz[n] = Vz[n] / sum_phi[n];
    }


    // time evolution of particles' positions --			
    #pragma acc parallel loop collapse(1)
    for(int n=0;n<N;n++)
    {
        Rx[n] += DT*Vx[n]; 			
        Ry[n] += DT*Vy[n];
        Rz[n] += DT*Vz[n];
        

        //PBC	
        Rx[n] = fmod(Rx[n], (double)Nx);		
        Ry[n] = fmod(Ry[n], (double)Ny);
        Rz[n] = fmod(Rz[n], (double)Nz);
    }

}