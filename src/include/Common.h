#ifndef COMMON_H
#define COMMON_H

#include <cmath>

// 索引约定：x 变化最快。与 main.cpp 的 cufftPlan3d(Nz, Ny, Nx) 匹配
#define IDX(i, j, k) ((i) + (j)*Nx + (k)*Nx*Ny)

// --- 流体求解器参数（POD，按值进 GPU kernel）---
struct NS_Config
{
    int    Nx, Ny, Nz;
    double dt, inv_dt;
    double W;           // 噪声强度系数，应由 kT 派生：W = sqrt(2*kT/dt)
    int    wall_z;      // 0 = z 周期（默认）；1 = z 上下无滑移壁面（Phase 7）
};

// --- 相场参数（POD，按值进 GPU kernel）---
// 从编译期常量改为运行时成员，为 Phase 2 的配置文件让路。
// 循环上界用运行时值对 OpenACC 无影响。
struct PhiParams
{
    double radius;      // 粒子半径 a
    double inv_xi;      // 1/ξ，界面宽度的倒数
    double ratio_eta;   // η_c / η_ℓ
    int    range;       // 相场作用范围
    int    range_m1;
    int    n_range;     // 2*range，局部模板盒的边长
    double range2;      // 球形截断判据（比较的是距离平方）
};

// 相场函数 φ(r) = ½[tanh((a - |r|)/ξ) + 1]
#pragma acc routine seq
static inline double order(double dx, double dy, double dz,
                           double radius, double inv_xi)
{
	return 0.5*(tanh((radius - sqrt(dx*dx + dy*dy + dz*dz))*inv_xi) + 1.);
}

// 由物理量构造 PhiParams。Phase 2 后改由 FpdConfig 驱动。
static inline PhiParams make_phi_params(double radius, double xi, double ratio_eta)
{
    PhiParams pp;
    pp.radius    = radius;
    pp.inv_xi    = 1.0 / xi;
    pp.ratio_eta = ratio_eta;
    pp.range     = (int)(2.*(radius + xi));
    pp.range_m1  = pp.range - 1;
    pp.n_range   = 2 * pp.range;
    pp.range2    = (double)(pp.range * pp.range);
    return pp;
}

// 由 dt 构造 NS_Config，保证 inv_dt 与 W 永远自洽（修 C10 / C6）
static inline NS_Config make_ns_config(int Nx, int Ny, int Nz,
                                       double dt, double kT, bool noise_on,
                                       int wall_z = 0)
{
    NS_Config cfg;
    cfg.Nx = Nx; cfg.Ny = Ny; cfg.Nz = Nz;
    cfg.dt     = dt;
    cfg.inv_dt = 1.0 / dt;
    cfg.W      = noise_on ? sqrt(2.0 * kT / dt) : 0.0;
    cfg.wall_z = wall_z;
    return cfg;
}

#endif
