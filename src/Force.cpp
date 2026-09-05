#include "include/Force.h"

/**
 * @brief 将粒子的外部力映射到流体网格受力场
 * @param Fx, Fy, Fz 粒子受力数组 (大小为 N)
 * @param fx, fy, fz 流体网格受力场 (大小为 Nx*Ny*Nz)
 */
void update_force_field(
    NS_Config cfg,
    double* sum_phi,
    int N, double* Rx, double* Ry, double* Rz,
    double* Fx, double* Fy, double* Fz,
    double* fx, double* fy, double* fz
)
{
    int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;

    // --- 1. 重置受力场 (每一帧都需要重新计算) ---
    #pragma acc parallel loop collapse(3) present(fx, fy, fz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                fx[ijk] = 0.0; fy[ijk] = 0.0; fz[ijk] = 0.0;
            }
        }
    }

    // --- 2. 粒子受力投影 ---
    #pragma acc parallel loop collapse(4) present(Rx, Ry, Rz, Fx, Fy, Fz, fx, fy, fz)
    for (int n = 0; n < N; n++)
    {
        for (int li = 0; li < N_range; li++)
        {
            for (int lj = 0; lj < N_range; lj++)
            {
                for (int lk = 0; lk < N_range; lk++)
                {

                    double Rnx = Rx[n], Rny = Ry[n], Rnz = Rz[n];
                    int in = (int)Rnx, jn = (int)Rny, kn = (int)Rnz;

                    int irP = (li + in - range_m1 + Nx) % Nx;
                    int jrP = (lj + jn - range_m1 + Ny) % Ny;
                    int krP = (lk + kn - range_m1 + Nz) % Nz;
                    int ijkP = IDX(irP, jrP, krP);

                    double dx = (li + in - range_m1) - Rnx;
                    double dy = (lj + jn - range_m1) - Rny;
                    double dz = (lk + kn - range_m1) - Rnz;

                    if (dx*dx + dy*dy + dz*dz <= range2)
                    {
                        // 关键点：力场必须投射到交错网格的面中心
                        // 使用原子加法处理多个粒子对同一格点的贡献
                        
                        // x方向分量：投射到 (i+0.5, j, k)
                        double weightX = order(dx + 0.5, dy, dz);
                        #pragma acc atomic update
                        fx[ijkP] += Fx[n] * weightX / sum_phi[n];

                        // y方向分量：投射到 (i, j+0.5, k)
                        double weightY = order(dx, dy + 0.5, dz);
                        #pragma acc atomic update
                        fy[ijkP] += Fy[n] * weightY / sum_phi[n];

                        // z方向分量：投射到 (i, j, k+0.5)
                        double weightZ = order(dx, dy, dz + 0.5);
                        #pragma acc atomic update
                        fz[ijkP] += Fz[n] * weightZ / sum_phi[n];

                    }
                }
            }
        }
        



    }
}