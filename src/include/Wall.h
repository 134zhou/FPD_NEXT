#ifndef WALL_H
#define WALL_H

#include "Common.h"

// ============================================================================
// z 向棱边布局与法向边界的【唯一真值源】（Phase 7-B）
//
// z 边界只有一种：z = -1/2 与 z = Nz-1/2 两片无滑移硬壁，x/y 周期。
// （曾有的 wall_z=0 全周期分支在 Phase 8-A 随三维 FFT 路径一起删除。）
//
// 核心设计：**逻辑 k 的物理含义是 z = k+1/2（棱边/面心）或 z = k（体心）**，
// 与存储映射解耦。Stokes.cpp 把内部棱边与两片壁面拆成独立 kernel；本文件只负责
// 棱边存储映射、法向面钉死与 RNG 长度，避免把物理逻辑下标误当成存储层号。
//
// 量按 z 半格偏移分三类：
//   体心层  p, eta, Π_xx/yy/zz, vx/fx, vy/fy, etaXY, Π_xy   逻辑 k ∈ [0, Nz-1]   存 IDX
//   z 面层  vz, fz                                          自由 k ∈ [0, Nz-2]   存 IDX
//                                                            k=Nz-1 是上壁，恒 0；k=-1 不存在
//   棱边层  etaYZ/etaZX, Π_yz/Π_zx                          逻辑 k ∈ [-1, Nz-1]  存 Nz+1 层
//
// ⚠️ 为什么棱边层必须多一层：下壁棱边（逻辑 k=-1）与上壁棱边（逻辑 k=Nz-1）承载的
//    vx 完全不同（一个是 -vx[0]，一个是 -vx[Nz-1]）。Nz 层布局下 k=-1 会回绕到
//    Nz-1，两者【撞同一个槽位】。所以 Nz+1 不是图方便，是必需。
//
// ⚠️ 切向无滑移的 ghost 关系是 vx[-1]=-vx[0]、vx[Nz]=-vx[Nz-1]（vy 同理）。
//    Stokes.cpp 的壁面 kernel 已把它代入通量：对流为 0，单侧应变率为 ±2vx/±2vy。
//
// ⚠️ 壁面棱边的剪切噪声必须 ×√2。推导见
//    doc/PressurePoisson.md §14：噪声力幅度 ∝ √(η·γ)，离散涨落耗散要求 γ·w² = w，
//    壁面棱边的控制体只有内部的一半 ⇒ w = 1/2 ⇒ γ = 2。漏掉它不会报错，只会让
//    壁面附近的 kT 系统性偏低。
//
// 本文件全部 inline + `acc routine seq`，cfg 按值传（照 Stencil.h 的约定：
// acc routine seq 里的引用参数会要求被引用对象位于设备内存）。
// ============================================================================

// 棱边数组的 z 层数（逻辑 k ∈ [-1, Nz-1] 共 Nz+1 个）
#pragma acc routine seq
static inline int wz_edge_nz(NS_Config cfg)
{
    return cfg.Nz + 1;
}

// 棱边数组元素总数（也是 randN 每个分量的跨步）
#pragma acc routine seq
static inline int wz_edge_size(NS_Config cfg)
{
    return cfg.Nx * cfg.Ny * wz_edge_nz(cfg);
}

// 棱边存储下标。逻辑 k 的物理位置是 z = k + 1/2，k ∈ [-1, Nz-1]。
// 平移 +1 存进 Nz+1 层，下壁棱边落在槽位 0。
//
// ⚠️ 越界【没有被兜住】：这里没有任何回绕，k ∉ [-1, Nz-1] 就是一次无诊断的越界写。
//    这份安全是【调用方的前置条件】而非访问器的性质 —— 删除周期分支时丢掉的
//    正是原先 `%Nz` 提供的静默吸收。调用方必须自己保证范围（Stencil.h 先做
//    球形+ZKind 判据，Stokes.cpp / Viscosity.cpp 的循环都在 [-1, Nz-1]）。
//    `acc routine seq` 里做不了运行时检查。
#pragma acc routine seq
static inline int wz_edge_idx(NS_Config cfg, int i, int j, int k)
{
    // 注意 IDX 是宏，内部引用裸的 Nx/Ny，必须靠同名局部变量代入
    const int Nx = cfg.Nx, Ny = cfg.Ny;
    return i + j * Nx + (k + 1) * Nx * Ny;
}

// 棱边的 z 向下邻居的【逻辑】下标。z 向不再回绕（k=0 时给 -1，那正是下壁棱边，
// 合法且必须被访问到）。
//
// 保留逻辑下标形式，让调用点显式表达「当前棱边的下邻居」，不改写成依赖存储布局的
// `eb - Nx*Ny` 指针算术。
#pragma acc routine seq
static inline int wz_edge_km(NS_Config cfg, int k)
{
    (void)cfg;
    return k - 1;
}

// z 面（vz/fz 的位置）是否被钉死。上壁那个元素存在但恒 0；
// 下壁（k=-1）根本不在数组里，Stokes.cpp 的壁面公式直接使用 0。
#pragma acc routine seq
static inline int wz_pinned_zface(NS_Config cfg, int k)
{
    return k >= cfg.Nz - 1 ? 1 : 0;
}

// --- RNG slot 的两个数组长度 ----------------------------------------------
// ⚠️ curandGenerateNormalDouble 要求生成个数为【偶数】（Box–Muller 成对）。
//    3*Nx*Ny*(Nz+1) 在奇数盒子上会是奇数（例如 3x3x3 的壁面盒 = 243），
//    不取偶会在运行时报 CURAND_STATUS_LENGTH_NOT_MULTIPLE。
//    取偶不改变前 n 个值，所以 slot 的定位语义不受影响。
static inline unsigned long long wz_even(unsigned long long n)
{
    return (n + 1ULL) & ~1ULL;
}

#pragma acc routine seq
static inline unsigned long long wz_rand_slot_d(NS_Config cfg)
{
    return wz_even(3ULL * (unsigned long long)cfg.Nx * (unsigned long long)cfg.Ny
                        * (unsigned long long)cfg.Nz);
}

#pragma acc routine seq
static inline unsigned long long wz_rand_slot_n(NS_Config cfg)
{
    return wz_even(3ULL * (unsigned long long)wz_edge_size(cfg));
}

#endif
