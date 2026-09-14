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

    // --- 判据专用开关（生产路径恒为默认值）---
    // adv_on = 0 关掉对流项：系统退化成【严格线性高斯】，于是
    //   「平稳协方差 == kT·I」变成一条精确的矩阵恒等式，可以用逐个注入单位向量的
    //   方式秒级判定，不需要跑几小时做统计。这是能造出来的最强单元测试。
    // noise_skip = 1 时不重新生成随机数，直接用 randD/randN 的现值 —— 供判据注入。
    // 两者取默认值时，算术表达式与没有它们时【逐位一致】（乘 1.0 是精确的）。
    int    adv_on;
    int    noise_skip;
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
//
// ⚠️ wall_z 【没有默认值】，这是故意的：带默认值时漏写的调用点会静默取周期语义，
//    在只支持壁面的构建里变成一个不报错的错误结果。去掉默认值让编译器穷举调用点。
static inline NS_Config make_ns_config(int Nx, int Ny, int Nz,
                                       double dt, double kT, bool noise_on,
                                       int wall_z)
{
    NS_Config cfg;
    cfg.Nx = Nx; cfg.Ny = Ny; cfg.Nz = Nz;
    cfg.dt     = dt;
    cfg.inv_dt = 1.0 / dt;
    cfg.W      = noise_on ? sqrt(2.0 * kT / dt) : 0.0;
    cfg.wall_z = wall_z;
    cfg.adv_on       = 1;
    cfg.noise_skip   = 0;
    cfg.noise_gamma1 = 0;
    return cfg;
}

#endif
