#ifndef TESTS_H
#define TESTS_H

#include "Common.h"

// ============================================================================
// fpd_check 的各自检/验证路径的声明。
//
// 约定：返回值 = 失败项数，0 表示全绿。不返回 bool —— 需要知道「错了几条」
// 来判断是单点故障还是整片崩。与 main 的退出码直接对应。
// ============================================================================

// Phase 1 自检：力守恒 + 亚格点不变性 + 交错混叠体检（需要 GPU）
int run_force_conservation_check(NS_Config cfg, PhiParams pp);

// Phase 1 自检：多粒子重叠（N=2，C2 的判决性判据）（需要 GPU）
int run_overlap_check(NS_Config cfg, PhiParams pp);

// 测试 3A：纯流体噪声谱（无粒子，需要 GPU）
int run_noise_check(int L, double dt, double kT, long n_steps);

// 测试 3B：粒子能量均分（需要 GPU）
enum EquipartMode { EQ_GHOST, EQ_FROZEN, EQ_MOVING };
int run_equipartition_check(EquipartMode mode, int L, double dt, double kT,
                            long n_steps, unsigned long long seed);

// 势函数自检（纯 CPU，无卡可跑）：力=-dU/dr、黄金表、最小镜像、N=3、特征点
int run_potential_check();

#endif
