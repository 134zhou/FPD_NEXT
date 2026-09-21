#ifndef COMMON_H
#define COMMON_H

#include <cmath>

// 索引约定：x 变化最快。与 cuFFT 的 row-major 约定匹配：2D 批量变换
// （cufftPlanMany，n = {Ny, Nx}）把 x 当最后一维，z-slab 在 IDX 下天然连续。
#define IDX(i, j, k) ((i) + (j)*Nx + (k)*Nx*Ny)

// --- 流体求解器参数（POD，按值进 GPU kernel）---
struct NS_Config
{
    // z 边界【只有一种】：x/y 周期，z 上下无滑移硬壁。曾有的 wall_z 字段与
    // 三维 FFT 全周期路径已于 Phase 8-A 整条删除 —— 单一语义 = 单一漂移面。
    // 棱边存储与法向面钉死由 Wall.h 定义；切向壁面通量在 Stokes.cpp 显式展开。
    int    Nx, Ny, Nz;
    double dt, inv_dt;
    double W;           // 噪声强度系数，应由 kT 派生：W = sqrt(2*kT/dt)

    // --- 判据专用开关（生产路径恒为默认值）---
    // noise_gamma1 = 1 时把壁面棱边的噪声因子强制成 1（关掉 √2）。**只为判据对照**，
    // 用来证明 √2 是被数据选中的、而不是被假设的。生产路径恒为 0。
    int    noise_gamma1;
};

// --- 相场参数（POD，按值进 GPU kernel）---
// 从编译期常量改为运行时成员，为 Phase 2 的配置文件让路。
// 循环上界用运行时值对 OpenACC 无影响。
struct PhiParams
{
    double radius;      // 粒子半径 a
    double inv_xi;      // 1/ξ，界面宽度的倒数
    double ratio_eta;   // η_c / η_ℓ
    int    range;       // 相场作用范围（球形截断半径）
    int    n_range;     // 2*range + 1，局部模板盒的边长（见下）
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
    // ⚠️ n_range 必须是 2*range + 1，【不能】是 2*range。
    //    盒取 [kn-range, kn+range]（见 Stencil.h）。半径为 range 的球心落在格胞
    //    任意位置时，能触及的整数层最多有 2*range+1 个 —— 宽 2*range 的盒会
    //    在负方向【静默漏掉一层】。漏掉的权重随亚格点位置和 Loc 变化，属于
    //    C1/C2 那一类「方向相关的伪偏差」，而且力守恒判据抓不到（分子分母
    //    一起漏，比值自洽）。
    pp.n_range   = 2 * pp.range + 1;
    pp.range2    = (double)(pp.range * pp.range);
    return pp;
}

// 由 dt 构造 NS_Config，保证 inv_dt 与 W 永远自洽（修 C10 / C6）
static inline NS_Config make_ns_config(int Nx, int Ny, int Nz,
                                       double dt, double kT, bool noise_on)
{
    NS_Config cfg;
    cfg.Nx = Nx; cfg.Ny = Ny; cfg.Nz = Nz;
    cfg.dt     = dt;
    cfg.inv_dt = 1.0 / dt;
    cfg.W      = noise_on ? sqrt(2.0 * kT / dt) : 0.0;
    cfg.noise_gamma1 = 0;
    return cfg;
}

#endif
