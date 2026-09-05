#ifndef STENCIL_H
#define STENCIL_H

#include "Common.h"

// ============================================================================
// 交错网格模板盒遍历的【唯一真值源】
//
// 半格偏移约定、球形截断判据、周期回绕 —— 这三样东西只在本文件出现一次。
//
// 背景：Viscosity / Force / Velocity 三个模块曾各自复制粘贴同一段遍历循环，
// 然后漂移出两个致命缺陷（见 README 的 C1/C2）：力投影用 x-face 的权重、
// 归一化却用 cell-center 的求和，导致 Σ_grid f ≠ F_n，产生位置相关的伪力。
//
// 现在三个模块共用 stencil_point<L>()，只要 Loc 一致，
// 「Σ_grid fx == Fx[n]」就是恒等式而非巧合。
//
// ⚠️ 代价：本文件成了单点故障。它错了三个模块一起错，
//    而「力守恒」判据对它的内部公式错误是【盲的】（投影和归一化会一起错）。
//    所以必须同时有独立数值对照，见 spike/spike_template_routine.cpp。
// ============================================================================

// 交错网格上的位置。MAC 布局：
//   CELL              体心   —— p, eta, 对角应力 Π_xx/yy/zz
//   FACE_X/Y/Z        面心   —— vx/vy/vz, fx/fy/fz
//   EDGE_XY/YZ/ZX     棱心   —— etaXY/YZ/ZX, 非对角应力 Π_xy/yz/zx
enum Loc { CELL, FACE_X, FACE_Y, FACE_Z, EDGE_XY, EDGE_YZ, EDGE_ZX };

// 每个位置相对体心的半格偏移
template<Loc L> struct LocOffset;
template<> struct LocOffset<CELL>    { static constexpr double dx=0.0, dy=0.0, dz=0.0; };
template<> struct LocOffset<FACE_X>  { static constexpr double dx=0.5, dy=0.0, dz=0.0; };
template<> struct LocOffset<FACE_Y>  { static constexpr double dx=0.0, dy=0.5, dz=0.0; };
template<> struct LocOffset<FACE_Z>  { static constexpr double dx=0.0, dy=0.0, dz=0.5; };
template<> struct LocOffset<EDGE_XY> { static constexpr double dx=0.5, dy=0.5, dz=0.0; };
template<> struct LocOffset<EDGE_YZ> { static constexpr double dx=0.0, dy=0.5, dz=0.5; };
template<> struct LocOffset<EDGE_ZX> { static constexpr double dx=0.5, dy=0.0, dz=0.5; };

// 局部模板盒的平铺点数。调用方的内层循环上界。
static inline int stencil_size(PhiParams pp)
{
    return pp.n_range * pp.n_range * pp.n_range;
}

// 把局部平铺索引 l ∈ [0, n_range³) 映射到全局周期索引 ijk 和相场权重 w。
//
// 返回 false 表示该点在球形截断之外，调用方应跳过。
//
// 注意 pp / cfg 必须【按值】传：acc routine seq 里的引用参数
// 会要求被引用对象位于设备内存。
#pragma acc routine seq
template<Loc L>
inline bool stencil_point(PhiParams pp, NS_Config cfg,
                          double Rnx, double Rny, double Rnz,
                          int l, int& ijk, double& w)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const int nr = pp.n_range;

    // 平铺索引 -> 模板盒内的三维局部坐标
    const int li = l / (nr*nr);
    const int lj = (l / nr) % nr;
    const int lk = l % nr;

    // 模板盒以粒子所在格胞为中心
    const int in = (int)Rnx, jn = (int)Rny, kn = (int)Rnz;
    const int ir = li + in - pp.range_m1;
    const int jr = lj + jn - pp.range_m1;
    const int kr = lk + kn - pp.range_m1;

    // 到粒子中心的距离，带该 Loc 的半格偏移
    const double dx = ir - Rnx + LocOffset<L>::dx;
    const double dy = jr - Rny + LocOffset<L>::dy;
    const double dz = kr - Rnz + LocOffset<L>::dz;

    if (dx*dx + dy*dy + dz*dz > pp.range2) { return false; }   // 唯一的截断判据

    // 唯一的周期回绕
    ijk = IDX((ir + Nx) % Nx, (jr + Ny) % Ny, (kr + Nz) % Nz);
    w   = order(dx, dy, dz, pp.radius, pp.inv_xi);
    return true;
}

#endif
