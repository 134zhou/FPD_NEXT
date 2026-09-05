#include "./include/Force.h"

void update_force_field(
    NS_Config cfg, PhiParams pp,
    int N, const double* Rx, const double* Ry, const double* Rz,
    const double* Fx, const double* Fy, const double* Fz,
    const double* sum_phix, const double* sum_phiy, const double* sum_phiz,
    double* fx, double* fy, double* fz
)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const int n3 = stencil_size(pp);

    // --- 1. 清零 ---
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

    // --- 2. 投影到交错网格的面心 ---
    // Loc 必须与 update_viscosity_fields 里算 sum_phi{x,y,z} 时一致（修 C1）
    #pragma acc parallel loop collapse(2) \
        present(Rx, Ry, Rz, Fx, Fy, Fz, sum_phix, sum_phiy, sum_phiz, fx, fy, fz)
    for (int n = 0; n < N; n++)
    {
        for (int l = 0; l < n3; l++)
        {
            int ijk; double w;

            if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                fx[ijk] += Fx[n] * w / sum_phix[n];
            }
            if (stencil_point<FACE_Y>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                fy[ijk] += Fy[n] * w / sum_phiy[n];
            }
            if (stencil_point<FACE_Z>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                fz[ijk] += Fz[n] * w / sum_phiz[n];
            }
        }
    }
}
