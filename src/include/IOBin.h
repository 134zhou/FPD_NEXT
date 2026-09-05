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
//   version   uint32    1
//   endian    uint32    0x01020304        读端校验
//   flags     uint32    bit0 = has_pressure
//   Nx,Ny,Nz,N          int32 x4
//   step                int64
//   dt,kT,radius,xi,ratio_eta   double x5     只存【原始】物理量
//   noise_on            int32
//   seed                uint64
//   rng_draws           uint64
//   --- 数组，f64，按线性 IDX = i + j*Nx + k*Nx*Ny 顺序【连续】存放 ---
//   vx[size], vy[size], vz[size]
//   [p[size]]                              仅当 has_pressure
//   Rx[N],Ry[N],Rz[N], Rux[N],Ruy[N],Ruz[N]
//   Vx[N],Vy[N],Vz[N], Fx[N],Fy[N],Fz[N]
//   --- 尾部 ---
//   checksum  uint64    FNV-1a，覆盖此前所有字节
//
// ⚠️ 绝不存派生量 W / inv_dt / range / n_range。重启时重新过
//    make_ns_config() / make_phi_params()，否则 C10/C6 的唯一真值源保证被绕过。
//
// ⚠️ 数组是【一次 fread/fwrite】的连续块，不做三重循环。
//    Python 侧 reshape 成 (Nz, Ny, Nx) —— x 是最后一维。
// ============================================================================

#define FPD_MAGIC       "FPDCKPT"
#define FPD_VERSION     1u
#define FPD_ENDIAN_TAG  0x01020304u
#define FPD_FLAG_PRESSURE 0x1u

struct CkptHeader
{
    uint32_t version   = FPD_VERSION;
    uint32_t flags     = 0;
    int32_t  Nx = 0, Ny = 0, Nz = 0, N = 0;
    int64_t  step = 0;
    double   dt = 0, kT = 0, radius = 0, xi = 0, ratio_eta = 0;
    int32_t  noise_on = 1;
    uint64_t seed = 0;
    uint64_t rng_draws = 0;

    bool   has_pressure() const { return (flags & FPD_FLAG_PRESSURE) != 0; }
    size_t size() const { return (size_t)Nx * (size_t)Ny * (size_t)Nz; }
};

// 数组视图。IOBin 不拥有内存，只负责搬运。
// p 可为 nullptr（此时写出时不置 has_pressure）。
struct CkptArrays
{
    double *vx = 0, *vy = 0, *vz = 0, *p = 0;
    double *Rx = 0, *Ry = 0, *Rz = 0;
    double *Rux = 0, *Ruy = 0, *Ruz = 0;
    double *Vx = 0, *Vy = 0, *Vz = 0;
    double *Fx = 0, *Fy = 0, *Fz = 0;
};

// 只读 header（用于在分配数组前拿到 Nx/Ny/Nz/N）
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
