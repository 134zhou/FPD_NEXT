#include "./include/Viscosity.h"

void update_viscosity_fields(
    NS_Config cfg, PhiParams pp,
    int N, const double* Rx, const double* Ry, const double* Rz,
    double* sum_phix, double* sum_phiy, double* sum_phiz,
    double* eta, double* etaXY, double* etaYZ, double* etaZX
)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const int n3 = stencil_size(pp);

    // η = η_ℓ + (η_c - η_ℓ)Σφ，η_ℓ = 1 → 系数是 ratio_eta - 1（修 C3）
    const double d_eta = pp.ratio_eta - 1.0;

    // --- 1. 背景粘度 ---
    #pragma acc parallel loop collapse(3) present(eta, etaXY, etaYZ, etaZX)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                eta[ijk]   = 1.0;
                etaXY[ijk] = 1.0;
                etaYZ[ijk] = 1.0;
                etaZX[ijk] = 1.0;
            }
        }
    }

    // --- 2. 三套错开的归一化因子 ---
    // gang over 粒子 + 内层 vector reduction（修 C7：不再用 atomic 打同一地址）
    #pragma acc parallel loop gang present(Rx, Ry, Rz, sum_phix, sum_phiy, sum_phiz)
    for (int n = 0; n < N; n++)
    {
        double sx = 0.0, sy = 0.0, sz = 0.0;

        #pragma acc loop vector reduction(+:sx, sy, sz)
        for (int l = 0; l < n3; l++)
        {
            int ijk; double w;
            if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sx += w; }
            if (stencil_point<FACE_Y>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sy += w; }
            if (stencil_point<FACE_Z>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sz += w; }
        }

        sum_phix[n] = sx;
        sum_phiy[n] = sy;
        sum_phiz[n] = sz;
    }

    // --- 3. 粒子对粘度场的贡献 ---
    // 写到【网格】上，多粒子支撑域会重叠，必须 atomic；但这是分散写，无竞争热点
    #pragma acc parallel loop collapse(2) present(Rx, Ry, Rz, eta, etaXY, etaYZ, etaZX)
    for (int n = 0; n < N; n++)
    {
        for (int l = 0; l < n3; l++)
        {
            int ijk; double w;

            // 对角应力位置：体心
            if (stencil_point<CELL>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                eta[ijk] += d_eta * w;
            }

            // 剪切应力位置：三个棱心
            if (stencil_point<EDGE_XY>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                etaXY[ijk] += d_eta * w;
            }
            if (stencil_point<EDGE_YZ>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                etaYZ[ijk] += d_eta * w;
            }
            if (stencil_point<EDGE_ZX>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                etaZX[ijk] += d_eta * w;
            }
        }
    }
}
