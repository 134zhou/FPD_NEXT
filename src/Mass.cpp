#include "Mass.h"
#include "Stencil.h"

// ============================================================================
// 三个面心的 a = 1/rho。位置约定全部来自 Stencil.h，见 Mass.h 的说明。
//
// 三段式（与 Viscosity.cpp 的 eta 场同构）：
//   1. 全填 1.0（= rho_l）
//   2. 逐粒子 atomic 累加 d_rho * w（= rho 场）
//   3. 就地取倒数（= a 场）
// 不要在 2 和 3 之间插入任何「顺手优化」：先建 rho 再取倒数是唯一正确的顺序，
// 直接对 a 累加 1/rho 会把重叠支撑域的插值变成非线性的错误平均。
// ============================================================================
void update_mass_fields(NS_Config cfg, PhiParams pp, double d_rho, int N,
                        const double* Rx, const double* Ry, const double* Rz,
                        double* inv_rho_x, double* inv_rho_y, double* inv_rho_z)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const int n3 = stencil_size(pp);

    // --- 1. 背景密度 rho = 1 ---
    #pragma acc parallel loop collapse(3) present(inv_rho_x, inv_rho_y, inv_rho_z)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                const int ijk = IDX(i, j, k);
                inv_rho_x[ijk] = 1.0;
                inv_rho_y[ijk] = 1.0;
                inv_rho_z[ijk] = 1.0;
            }
        }
    }

    // --- 2. 粒子贡献：rho += d_rho * phi ---
    // FACE_Z 的 k = Nz-1（上壁法向面）由 LocZKind<FACE_Z> = ZFACE 自动排除，
    // 所以那里保持 1.0；与 sum_phiz 的求和域严格一致。
    #pragma acc parallel loop collapse(2) present(Rx, Ry, Rz, inv_rho_x, inv_rho_y, inv_rho_z)
    for (int n = 0; n < N; n++)
    {
        for (int l = 0; l < n3; l++)
        {
            int ijk; double w;
            if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                inv_rho_x[ijk] += d_rho * w;
            }
            if (stencil_point<FACE_Y>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                inv_rho_y[ijk] += d_rho * w;
            }
            if (stencil_point<FACE_Z>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                inv_rho_z[ijk] += d_rho * w;
            }
        }
    }

    // --- 3. rho -> 1/rho ---
    #pragma acc parallel loop collapse(3) present(inv_rho_x, inv_rho_y, inv_rho_z)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                const int ijk = IDX(i, j, k);
                inv_rho_x[ijk] = 1.0 / inv_rho_x[ijk];
                inv_rho_y[ijk] = 1.0 / inv_rho_y[ijk];
                inv_rho_z[ijk] = 1.0 / inv_rho_z[ijk];
            }
        }
    }
}

// ============================================================================
// 重力体力。f_g = (rho - rho_l) g = (1/a - 1) g。
//
// 必须在 update_force_field 之后调用：它把体力【加】到已经摊好粒子力的面心力场上。
// 除以 rho 的那一步在 Stokes.cpp 的显式步里（面心 a 用同一个数组），所以这里
// 只写体密度，不写加速度。
// ============================================================================
void apply_gravity(NS_Config cfg, ExternalField g,
                   const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
                   double* fx, double* fy, double* fz)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;

    #pragma acc parallel loop collapse(3) \
        present(inv_rho_x, inv_rho_y, inv_rho_z, fx, fy, fz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                const int ijk = IDX(i, j, k);
                // (1/a - 1) 就是 (rho - 1)；除一次比重建 rho 场更省。
                fx[ijk] += (1.0 / inv_rho_x[ijk] - 1.0) * g.gx;
                fy[ijk] += (1.0 / inv_rho_y[ijk] - 1.0) * g.gy;
                fz[ijk] += (1.0 / inv_rho_z[ijk] - 1.0) * g.gz;
            }
        }
    }
}
