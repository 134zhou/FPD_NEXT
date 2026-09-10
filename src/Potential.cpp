#include "Potential.h"

PotentialParams make_potential_params(int type, int shift,
                                      double eps, double sigma,
                                      double De, double alpha, double r_eq,
                                      double rcut)
{
    PotentialParams pp;
    pp.type   = type;
    pp.shift  = shift;
    pp.eps    = eps;
    pp.sigma  = sigma;
    pp.De     = De;
    pp.alpha  = alpha;
    pp.r_eq   = r_eq;

    // WCA 的截断是派生的：r = 2^{1/6}σ（LJ 势最小值处，力天然连续为零）
    if (type == POT_WCA) { rcut = pow(2.0, 1.0 / 6.0) * sigma; }
    pp.rcut  = rcut;
    pp.rcut2 = rcut * rcut;

    // 派生量：U_bare(rcut) 与 U'_bare(rcut)。复用 pair_*_bare，不重写公式
    //（避免 C1/C2 式的公式漂移）。U' = -g*r。
    if (type == POT_NONE)
    {
        pp.u_at_rc    = 0.0;
        pp.dudr_at_rc = 0.0;
    }
    else
    {
        pp.u_at_rc    = pair_energy_bare(pp, rcut);
        pp.dudr_at_rc = -pair_force_over_r_bare(pp, rcut) * rcut;
    }
    return pp;
}

void compute_particle_forces(NS_Config cfg, PotentialParams pot, ExternalField ext,
                             int N,
                             const double* Rx, const double* Ry, const double* Rz,
                             double* Fx, double* Fy, double* Fz)
{
    const double Lx = (double)cfg.Nx, Ly = (double)cfg.Ny;
    const double Lz = pbc_length_z(cfg);   // 壁面模式下退化为恒等映射（z 不折叠）

    // 全矩阵 + gang/vector reduction（零 atomic）：每个粒子 i 由一个 gang 独立
    // 归约它对所有 j 的力。算两倍的对（j<i 和 j>i 各一次），在 GPU 上不值一提，
    // 却换掉了 atomic —— 正是 C7 的教训要求的写法。
    #pragma acc parallel loop gang present(Rx, Ry, Rz, Fx, Fy, Fz)
    for (int i = 0; i < N; i++)
    {
        double sx = 0.0, sy = 0.0, sz = 0.0;
        #pragma acc loop vector reduction(+:sx, sy, sz)
        for (int j = 0; j < N; j++)
        {
            if (j == i) { continue; }
            const double dx = min_image(Rx[i] - Rx[j], Lx);
            const double dy = min_image(Ry[i] - Ry[j], Ly);
            const double dz = min_image(Rz[i] - Rz[j], Lz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double g  = pair_force_over_r(pot, r2);
            sx += g * dx;
            sy += g * dy;
            sz += g * dz;
        }
        Fx[i] = ext.gx + sx;
        Fy[i] = ext.gy + sy;
        Fz[i] = ext.gz + sz;
    }
}

void compute_particle_forces_cpu(NS_Config cfg, PotentialParams pot, ExternalField ext,
                                 int N,
                                 const double* Rx, const double* Ry, const double* Rz,
                                 double* Fx, double* Fy, double* Fz)
{
    const double Lx = (double)cfg.Nx, Ly = (double)cfg.Ny;
    const double Lz = pbc_length_z(cfg);   // 壁面模式下退化为恒等映射（z 不折叠）

    for (int n = 0; n < N; n++) { Fx[n] = ext.gx; Fy[n] = ext.gy; Fz[n] = ext.gz; }
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
        {
            const double dx = min_image(Rx[i] - Rx[j], Lx);
            const double dy = min_image(Ry[i] - Ry[j], Ly);
            const double dz = min_image(Rz[i] - Rz[j], Lz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double g  = pair_force_over_r(pot, r2);
            Fx[i] += g*dx;  Fx[j] -= g*dx;
            Fy[i] += g*dy;  Fy[j] -= g*dy;
            Fz[i] += g*dz;  Fz[j] -= g*dz;
        }
}
