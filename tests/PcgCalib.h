#ifndef PCG_CALIB_H
#define PCG_CALIB_H

#include <vector>
#include "Step.h"

// ============================================================================
// PCG 容差标定的核心：在生产状态上，量「残差导致的粒子速度误差」与
// 「同一步热噪声导致的粒子速度增量」之比。fpd_calib 与 fpd_check --pcg-noise 共用。
//
// 一个测量点（calib_measure）：
//   1. step_forces + step_fields（生产路径），记下面心 a = 1/ρ
//   2. 同一 step、同一 RNG slot 跑两次显式步：W 正常 → v*_A，W = 0 → v*_B；
//      v* 对 W 线性，Δv* = v*_A - v*_B 就是本步噪声增量
//   3. 噪声增量投影：PΔv* = Δv* - dt·a·∇p_N，L p_N = div Δv*/dt（参考精度求解）；
//      粒子平均（生产 update_particle_velocity）得 ΔV_n，分母 Vn_rms = sqrt(mean_n |ΔV_n|²)
//   4. 参考压力 p_ref：L p_ref = div v*_A/dt（参考精度）；对每个候选容差 τ 用生产判停
//      （tol = τ，上限 PCG_MAX_ITER）求 p_τ，δv = -dt·a·∇(p_τ - p_ref)：
//        ε_V(τ) = max_n |δV_n| / Vn_rms
//        ε_E(τ) = sqrt( Σ_f ρ_f |δv_f|² / Σ_f ρ_f |PΔv*_f|² )
//   5. 把 v、V 恢复成测量前的值；调用方接着走一个正常的生产步
//
// 右端、修正步、能量都在主机侧手写（独立实现），生产代码只作为被测对象（PCG）与
// 物理装配（Step.h、update_particle_velocity）被调用。只产生标量，不写场数据。
// ============================================================================

struct CalibLevel
{
    double tol;        // 候选容差 τ
    int    iters;      // 生产判停下的迭代数
    double res;        // 实际相对残差 ||r||/||b||
    double eV, eE;     // ε_V、ε_E
};

struct CalibPoint
{
    long   step;
    double noise_frac;     // ||b_noise|| / ||b||
    double Vn_rms;         // rms_n |ΔV_n|（分母）
    double ref_res;        // 参考解的相对残差
    std::vector<CalibLevel> lv;
};

struct CalibOpts
{
    long   warm   = 200;   // 首个测量点前的生产步数
    long   gap    = 50;    // 测量点间隔
    int    meas   = 4;     // 测量点个数
    double target = 1e-3;  // ε 门槛
    double margin = 10.0;  // 推荐值 = 最松合格档 / margin
};

// 候选容差：1e-2, 3e-3, 1e-3, …, 1e-12（每半个量级一档）
std::vector<double> calib_levels();

// 从 start_step 起推进生产轨迹，在 warm、warm+gap、… 处测量。返回 0 表示正常。
int calib_run(const StepParams& sp, FpdState& st, long start_step,
              const CalibOpts& o, std::vector<CalibPoint>& out);

// 打印逐点与汇总表；返回推荐容差（<= 0 表示连 1e-12 也不满足）。
// tau_star 返回满足门槛的最松档，cur_tol 是配置当前值（报告它是否满足）。
double calib_report(const std::vector<CalibPoint>& pts, const CalibOpts& o,
                    double cur_tol, double& tau_star);

#endif
