#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Common.h"
#include "Stencil.h"

// ============================================================================
// FPD 常量表与平衡态统计分析
// ============================================================================

// ---------------------------------------------------------------------------
// 【为什么三个方向要分开列】
// 交错网格上 phi_x/phi_y/phi_z 是同一个连续函数在三组错开半格的格点上的
// 【重新采样】，不是彼此的平移 —— 平移只改相位不改模，重采样会改混叠分支的
// 组合方式，所以 |phi_hat_x(k)|^2 != |phi_hat_y(k)|^2 严格成立。
//
// 由 Parseval，Sum_k |phi_hat_alpha(k)|^2 = N_g * Sum_r phi_alpha(r)^2，
// 所以三方向 Sum(phi^2) 的一致性【本身】就是交错混叠的直接数值判据。
//
// 注意 main 的 --check 目前只检了 Sum(phi_alpha)（那只是 k=0 分量），
// Sum(phi_alpha^2) 是漏掉的那一半。
// ---------------------------------------------------------------------------
struct FpdConstants
{
    // 与模拟完全相同的离散求和（stencil_point<FACE_X/Y/Z>）
    double int_phi[3];        // Sum phi_alpha
    double int_phi2[3];       // Sum phi_alpha^2
    double lambda_T[3];       // lambda^T = int(phi) / int(phi^2)
    double mass[3];           // M_i = rho (int phi)^2 / int(phi^2)   (rho = dV = 1)
    double mass_eff[3];       // M_eff = 1.5 M_i（附加质量）

    // 三方向一致性
    double spread_int_phi;    // (max-min)/min
    double spread_int_phi2;   // 交错混叠的一阶体检
    double spread_mass;

    // 亚格点扫描（格胞内 n_scan^3 个偏移）
    int    n_scan;
    double mass_mean, mass_min, mass_max;
    double spread_subgrid;    // 【整个方法的系统误差下限】
};

// 用与模拟【同一个】 stencil_point 做离散求和。纯 CPU，一次性调用。
// 在 (Rx,Ry,Rz) 求主值，另在格胞内扫 n_scan^3 个偏移给误差带。
FpdConstants compute_fpd_constants(PhiParams pp, NS_Config cfg,
                                   double Rx, double Ry, double Rz, int n_scan);

// 打印常量表，并与连续球坐标积分、扩散界面解析展开对照
void print_fpd_constants(FpdConstants fc, PhiParams pp);

// ⚠️ 这里曾有 equipartition_target()（粒子能量均分的目标值
//    <|V|^2> = 2kT/M_i − 2kT/(rho L^3)）。它是【周期专属】的：第二项是周期盒
//    k=0 模式被冻结在零的修正，无滑移壁面下总动量本就不守恒（壁面是真实的
//    动量汇），没有可平移的 k=0 模式，该修正不成立。随 z 周期路径一起删除。
//    替代判据需要新推导，见 PROGRESS.md 的覆盖损失第 3 条。

// ============================================================================
// 分块平均（block averaging）
//
// 把相关样本当独立样本会把误差棒低估约 sqrt(2 tau_int / dt_samp) 倍
// —— 本项目参数下约 7.5 倍。不分块的误差棒是假的。
// ============================================================================
struct BlockStat
{
    int    block_len;         // 每块的样本数 b
    int    n_block;           // n_b = n / b
    double mean;
    double stderr_mean;       // sqrt( Var(块均值) / (n_b - 1) )
    double tau_int_samp;      // 0.5 * (stderr_b / stderr_naive)^2，以采样间隔为单位
};

// 分块长度扫描 b = 1,2,4,8,...，结果写进 out[0..返回值)
int blocking_analysis(const double* x, int n, BlockStat* out, int max_out);

// 平台判据：最小的 b 使 stderr(2b)/stderr(b) < 1.05 且 n_block(2b) >= min_block。
// min_block 默认 20 —— 块方差自身的相对误差是 1/sqrt(2(n_b-1))，n_b=20 时已有 16%。
// 找不到平台时返回的 block_len = -1：此时【误差棒不可信，必须延长模拟】。
BlockStat pick_plateau(const BlockStat* s, int n_s, int min_block);

// 打印整条 stderr(b) 曲线。平台是否存在是肉眼一秒能判、公式判不了的事。
void print_blocking(const BlockStat* s, int n_s, BlockStat picked);

#endif
