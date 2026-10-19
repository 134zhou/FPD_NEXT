#ifndef IOBIN_H
#define IOBIN_H

#include <string>
#include <cstdint>
#include "Common.h"

// ============================================================================
// .fpd —— C++ 侧【唯一】的输入输出格式
//
// 初始构型、检查点、重启文件都是这一个格式：文件里的 step 决定从哪继续。
// 因此 C++ 除配置文件外没有任何文本解析。
//
// 布局（小端，裸二进制）：
//   magic     char[8]   "FPDCKPT\0"
//   endian    uint32    0x01020304        读端校验
//   version   uint32    4（写）；读端接受 3 与 4
//   flags     uint32    仅 v4；bit0 = 含压力历史（FPD_FLAG_PRESSURE）
//   Nx,Ny,Nz,N          int32 x4
//   step                int64
//   --- 数组，f64，按线性 IDX = i + j*Nx + k*Nx*Ny 顺序【连续】存放 ---
//   vx[size], vy[size], vz[size]
//   Rx[N],Ry[N],Rz[N], Rux[N],Ruy[N],Ruz[N]
//   Vx[N],Vy[N],Vz[N], Fx[N],Fy[N],Fz[N]
//   [p[size], p_prev[size]]                   仅 flags bit0：压力分裂法的压力历史
//   --- 尾部 ---
//   checksum  uint64    FNV-1a，覆盖此前所有字节
//
// ⚠️ 绝不存派生量 W / inv_dt / range / n_range。重启时重新过
//    make_ns_config() / make_phi_params()，否则 C10/C6 的唯一真值源保证被绕过。
//
// ⚠️ 也绝不存 z 边界类型：v2 起它不再是变量（恒为无滑移壁面）。
//
// ⚠️ 密度场 a=1/rho 同样是派生量：每次由 R 经 Mass.cpp 重算，不入档。
//    所以 ratio_rho 不在 .fpd 里，必须与初态/续跑用同一个配置（.config.used）。
//
// ⚠️ 压力默认【不】入档：PCG 路径的 p 是求解器内部状态（p0=0 冷启动），不是续跑状态。
//    唯一的例外是压力分裂法（pressure_solver=split）：外推压力 p̂ = 2pⁿ - pⁿ⁻¹ 依赖历史，
//    要逐位重启就必须保存 pⁿ、pⁿ⁻¹，此时置 FPD_FLAG_PRESSURE 并追加两个整场数组。
//    合法文件大小恒为 v3：48 + 24*size + 96*N；v4：52 + 24*size + 96*N [+ 16*size]。
//    fpd_tool --dump-ckpt 会校验这条，防止整场数组悄悄回流。
//
// ⚠️ 数组是【一次 fread/fwrite】的连续块，不做三重循环。
//    Python 侧 reshape 成 (Nz, Ny, Nx) —— x 是最后一维。
// ============================================================================

// v3 删除物理参数、RNG 记录、压力标志及压力数据，只保留续跑状态。
// v4 在 version 之后加 flags，可选追加压力历史；不带压力时数据区与 v3 逐字节相同。
// v2 及更早版本由版本检查明确拒绝。
// 历史：version 1 → 2（Phase 8-A）时 z 周期存档语义被删除。
//    v1 的存档可能是 z 周期构型（vz[...,Nz-1] 非 0）。那样的初态会直接破坏
//    泊松的相容性条件 Σ_k b̂_k = 0，表现为压力的线性漂移，且【看不出错】。
#define FPD_MAGIC       "FPDCKPT"
#define FPD_VERSION     4u
#define FPD_VERSION_MIN 3u
#define FPD_ENDIAN_TAG  0x01020304u
#define FPD_FLAG_PRESSURE 1u

struct CkptHeader
{
    uint32_t version   = FPD_VERSION;
    uint32_t flags     = 0;          // v3 文件读入时为 0
    int32_t  Nx = 0, Ny = 0, Nz = 0, N = 0;
    int64_t  step = 0;

    size_t size() const { return (size_t)Nx * (size_t)Ny * (size_t)Nz; }
    bool has_pressure() const { return (flags & FPD_FLAG_PRESSURE) != 0; }
    // 合法文件的精确字节数（--dump-ckpt 的尺寸守卫）
    long long file_bytes() const
    {
        return (version >= 4u ? 52LL : 48LL) + 24LL * (long long)size() + 96LL * (long long)N
             + (has_pressure() ? 16LL * (long long)size() : 0LL);
    }
};

// 数组视图。IOBin 不拥有内存，只负责搬运。
struct CkptArrays
{
    double *vx = 0, *vy = 0, *vz = 0;
    double *Rx = 0, *Ry = 0, *Rz = 0;
    double *Rux = 0, *Ruy = 0, *Ruz = 0;
    double *Vx = 0, *Vy = 0, *Vz = 0;
    double *Fx = 0, *Fy = 0, *Fz = 0;
    // 压力历史：写时 h.flags 含 FPD_FLAG_PRESSURE 才使用；读时为 null 则读入后丢弃（仍参与校验和）
    double *p = 0, *p_prev = 0;
};

// 生产路径可先读头、分配数组，再从同一句柄读取数据与校验和。
struct CkptReader
{
    void* f = 0;
    uint64_t hash = 0;
};
bool open_checkpoint(const char* path, CkptHeader& h, CkptReader& r, std::string& err);
bool read_checkpoint(CkptReader& r, const CkptHeader& h,
                     const CkptArrays& a, std::string& err);
void close_checkpoint(CkptReader& r);

// 只读 header，供离线工具按尺寸分配数组
bool read_ckpt_header(const char* path, CkptHeader& h, std::string& err);

// 数组必须已由调用方按 header 的尺寸分配好
bool save_checkpoint(const char* path, const CkptHeader& h,
                     const CkptArrays& a, std::string& err);
bool load_checkpoint(const char* path, CkptHeader& h,
                     const CkptArrays& a, std::string& err);

// <out_dir>/<run_name>_<step 补零到 7 位>.fpd
std::string ckpt_path(const std::string& out_dir, const std::string& run_name, long step);

// 目录不存在就创建（不检查会静默丢数据）
bool ensure_dir(const std::string& path, std::string& err);

#endif
