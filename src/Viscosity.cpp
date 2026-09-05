#include "./include/Viscosity.h"

/**
 * @brief 更新粘度场：将粒子位置映射到网格中心及棱边的粘度
 * * 逻辑：
 * 1. 初始化所有 eta 为背景粘度 (1.0)
 * 2. 遍历粒子，对每个粒子影响的局部范围进行累加
 * 3. 使用 atomic 操作处理粒子重叠导致的并行写入冲突
 */
void update_viscosity_fields(
    NS_Config cfg,
    double* sum_phi,
    int N, double* Rx, double* Ry, double* Rz,
    double* phi_grid, double* eta, double* etaXY, double* etaYZ, double* etaZX
)
{
    int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    int size = Nx * Ny * Nz;

    #pragma acc parallel loop collapse(1) present(sum_phi)
    for (int n = 0; n < N; n++)
    {
        sum_phi[n] = 0.0;
    }

    // --- 1. 初始化背景粘度与 phi (1.0) ---
    #pragma acc parallel loop collapse(3) present(phi_grid, eta, etaXY, etaYZ, etaZX)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                eta[ijk]   = 1.0;
                phi_grid[ijk]   = 0.0;
                etaXY[ijk] = 1.0;
                etaYZ[ijk] = 1.0;
                etaZX[ijk] = 1.0;
            }
        }
    }

    // --- 2. 粒子映射到粘度场 ---
    // 并行化策略：外层循环粒子，内层循环局部范围
    // 注意：多个粒子可能影响同一个格点，必须使用 atomic
    #pragma acc parallel loop collapse(4) present(Rx, Ry, Rz, phi_grid, eta, etaXY, etaYZ, etaZX)
    for (int n = 0; n < N; n++)
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

                    // 计算相对于粒子中心的距离
                    double dx = ir - Rnx;
                    double dy = jr - Rny;
                    double dz = kr - Rnz;

                    // 只有在影响半径内的点才进行计算
                    if (dx*dx + dy*dy + dz*dz <= range2)
                    {
                        // 正应力位置 (Cell Center)
                        double val_e  = RATIO_ETA * order(dx, dy, dz);
                        // 切应力位置 (Edges - 带有 0.5 的偏移)
                        double val_xy = RATIO_ETA * order(dx+0.5, dy+0.5, dz);
                        double val_yz = RATIO_ETA * order(dx, dy+0.5, dz+0.5);
                        double val_zx = RATIO_ETA * order(dx+0.5, dy, dz+0.5);

                        #pragma acc atomic update
                        phi_grid[ijkP] += val_e;

                        #pragma acc atomic update
                        sum_phi[n] += val_e;

                        #pragma acc atomic update
                        eta[ijkP] += val_e;

                        #pragma acc atomic update
                        etaXY[ijkP] += val_xy;

                        #pragma acc atomic update
                        etaYZ[ijkP] += val_yz;

                        #pragma acc atomic update
                        etaZX[ijkP] += val_zx;
                    }
                }
            }
        }




    }
}