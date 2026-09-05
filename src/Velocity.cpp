#include "./include/Velocity.h"

void update_particle_velocity(
    NS_Config cfg, PhiParams pp,
    int N, const double* Rx, const double* Ry, const double* Rz,
    const double* sum_phix, const double* sum_phiy, const double* sum_phiz,
    const double* vx, const double* vy, const double* vz,
    double* Vx, double* Vy, double* Vz
)
{
    const int n3 = stencil_size(pp);

    // gang over 粒子 + 内层 vector reduction。
    // 局部累加变量天然清零，顺带修掉 C4（旧代码 Vx[n] += 前未清零）。
    #pragma acc parallel loop gang \
        present(Rx, Ry, Rz, sum_phix, sum_phiy, sum_phiz, vx, vy, vz, Vx, Vy, Vz)
    for (int n = 0; n < N; n++)
    {
        double sx = 0.0, sy = 0.0, sz = 0.0;

        #pragma acc loop vector reduction(+:sx, sy, sz)
        for (int l = 0; l < n3; l++)
        {
            int ijk; double w;

            // 权重与力投影、归一化因子同源：同一个 stencil_point<FACE_*>
            if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sx += vx[ijk] * w; }
            if (stencil_point<FACE_Y>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sy += vy[ijk] * w; }
            if (stencil_point<FACE_Z>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sz += vz[ijk] * w; }
        }

        Vx[n] = sx / sum_phix[n];
        Vy[n] = sy / sum_phiy[n];
        Vz[n] = sz / sum_phiz[n];
    }
}

void update_particle_position(
    NS_Config cfg,
    int N,
    double* Rx, double* Ry, double* Rz,
    double* Rux, double* Ruy, double* Ruz,
    const double* Vx, const double* Vy, const double* Vz
)
{
    const double DT = cfg.dt;
    const double Lx = (double)cfg.Nx, Ly = (double)cfg.Ny, Lz = (double)cfg.Nz;

    #pragma acc parallel loop present(Rx, Ry, Rz, Rux, Ruy, Ruz, Vx, Vy, Vz)
    for (int n = 0; n < N; n++)
    {
        const double dx = DT * Vx[n];
        const double dy = DT * Vy[n];
        const double dz = DT * Vz[n];

        // 不折叠的位置，供 MSD 使用
        Rux[n] += dx;
        Ruy[n] += dy;
        Ruz[n] += dz;

        Rx[n] += dx;
        Ry[n] += dy;
        Rz[n] += dz;

        // 周期边界（修 C5：fmod 对负数返回负值，边界会破）
        if      (Rx[n] <  0.0) { Rx[n] += Lx; }
        else if (Rx[n] >= Lx ) { Rx[n] -= Lx; }
        if      (Ry[n] <  0.0) { Ry[n] += Ly; }
        else if (Ry[n] >= Ly ) { Ry[n] -= Ly; }
        if      (Rz[n] <  0.0) { Rz[n] += Lz; }
        else if (Rz[n] >= Lz ) { Rz[n] -= Lz; }
    }
}
