#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include "Common.h"
#include "Potential.h"

// ============================================================================
// 主机侧完整参数集。含 std::string，【永不进 GPU】。
//
// 与 NS_Config / PhiParams 的关系：
//   FpdConfig（原始物理量） --make_ns_config--> NS_Config（POD，按值进 kernel）
//                          --make_phi_params-> PhiParams（POD，按值进 kernel）
//
// 配置文件里只出现【原始】物理量。派生量 W / inv_dt / range / n_range 一律由
// Common.h 的两个工厂函数计算 —— 那是它们的唯一真值源（修 C10/C6 的保证）。
// ============================================================================
struct FpdConfig
{
    // --- 网格与时间（必填）---
    int    Nx = 0, Ny = 0, Nz = 0;
    double dt = 0.0;
    long   n_steps = 0;          // 【绝对】目标步数，不是「再跑多少步」

    // --- 物理（必填）---
    double kT = 0.0;
    double radius = 0.0;
    double xi = 0.0;
    double ratio_eta = 0.0;

    // --- 物理（可选）---
    double gravity_x = 0.0;
    double gravity_y = 0.0;
    double gravity_z = 0.0;
    int    gravity_compensate = 1;   // 外场叠加背景力密度，抵消 k=0 漂移（过渡；壁面 Phase 后转 0）
    int    noise_on = 1;

    // --- 粒子间势（可选；potential=none 时以下全部不用）---
    std::string potential = "none";   // none | wca | morse | lj126
    double pot_eps    = 0.0;          // WCA/LJ 的 ε
    double pot_sigma  = 0.0;          // WCA/LJ 的 σ
    double pot_De     = 0.0;          // Morse 阱深
    double pot_alpha  = 0.0;          // Morse 宽度参数（叫 alpha 不叫 a —— a 是粒子半径）
    double pot_r_eq   = 0.0;          // Morse 平衡距离
    double pot_rcut   = 0.0;          // LJ/Morse 截断（WCA 派生，用户不应给）
    std::string pot_shift = "energy"; // none | energy | force

    // --- IO（可选）---
    std::string init_file;       // .fpd；留空则用内置默认（单粒子放盒心）
    std::string out_dir  = "out";
    std::string run_name = "run";
    long interval_ckpt = 10000;
    long interval_log  = 10000;
    int  save_pressure = 1;      // 存 p 便于可视化；重启时忽略

    unsigned long long seed = 1234ULL;

    // 显式设置过的 key 名（防「用不到的参数被静默忽略」）
    std::vector<std::string> keys_given;
};

// ---------------------------------------------------------------------------
// 字段表：解析、覆盖、dump 三者【共用同一张表】。
// 分开维护三份是 C1/C2 式漂移的经典温床。
// ---------------------------------------------------------------------------
enum FieldType { F_INT, F_LONG, F_DOUBLE, F_STRING, F_ULL };

struct FieldDesc
{
    const char* key;
    FieldType   type;
    void*       ptr;         // 指向某个 FpdConfig 实例的成员
    bool        required;
    const char* help;
};

std::vector<FieldDesc> config_fields(FpdConfig& c);

// 全部返回 true 表示成功；失败时 err 里是可直接给用户看的中文说明。
bool load_config(const char* path, FpdConfig& c, std::string& err);
bool apply_override(FpdConfig& c, const char* kv, std::string& err);   // "kT=1.0"
bool validate_config(const FpdConfig& c, std::string& err);

// 归档最终配置。附带打印派生量（标注为 derived，不可作为输入）。
bool dump_config(const FpdConfig& c, const char* path, std::string& err);
void print_config_help();

NS_Config make_ns_config(const FpdConfig& c);
PhiParams make_phi_params(const FpdConfig& c);
PotentialParams make_potential_params(const FpdConfig& c);
ExternalField  make_external_field(const FpdConfig& c);

#endif
