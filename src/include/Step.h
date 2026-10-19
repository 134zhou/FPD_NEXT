#ifndef STEP_H
#define STEP_H

#include "Common.h"
#include "Config.h"
#include "Potential.h"
#include "State.h"

// ============================================================================
// 生产主循环一步的各段【唯一实现】。main.cpp 与 fpd_calib（容差标定）共用，
// 标定程序因此测到的就是生产路径本身，不是一份手抄的副本。
//
// 一步的顺序（main.cpp）：
//   step_forces → [检查点 / 日志] → step_fields → step_fluid → step_particles
// 检查点夹在力与场之间，是为了让文件里的 F 对应文件里的 R。
// ============================================================================

// 一步所需的全部 POD 参数，由 FpdConfig（尺寸已取自 init_file）构造一次。
struct StepParams
{
    NS_Config       cfg;
    PhiParams       pp;
    PotentialParams pot;
    WallParams      wall;
    ExternalField   ext;
    double          d_rho;     // ratio_rho - 1
};

StepParams make_step_params(const FpdConfig& effective);

// F = 粒子对势 + 壁面势（重力不在这里）
void step_forces(const StepParams& sp, FpdState& st);

// η、S_i、面心 1/ρ（density_on）、力密度 f（背景 0）、重力体力（density_on 且 g≠0）
void step_fields(const StepParams& sp, FpdState& st);

// 噪声 → Π → v* → 压力 → 投影。cfg 单独传，标定程序要用 W=0 的副本。返回 PCG 迭代数。
int step_fluid(NS_Config cfg, FpdState& st, long step);

// V ← v，R += dt·V
void step_particles(const StepParams& sp, FpdState& st);

#endif
