#ifndef STATE_H
#define STATE_H

#include <vector>
#include <cufft.h>
#include <curand.h>
#include "Common.h"
#include "IOBin.h"
#include "Wall.h"

// ============================================================================
// FpdState —— 主机侧存储 + 设备映射 + FFT/RNG 句柄的唯一持有者
//
// 背景：四个自检路径曾各自复制一份「std::vector 分配 → #pragma acc enter data
// → exit data」样板，共 4 份 —— 这是 C1/C2 式漂移的温床。FpdState 把它收成一份。
//
// 映射机制：用 OpenACC【运行时 API】（acc_copyin/acc_create/acc_update_self/
// acc_update_device/acc_delete），不用 pragma。实测（tests/spike/spike_state_map.cpp）：
// 运行时 API 建的映射能被【另一个 TU】里的 present() 子句完全认账，
// 且 acc_is_present 可用于运行时断言。这是它相对 pragma 的最大优势 ——
// 见 require()。
//
// 与 tool_main.cpp 的 Holder 同源：Holder 是 CPU 版先例（vector 持有 + 装配
// 裸指针视图），FpdState 把它扩展到设备映射。不是新抽象。
// ============================================================================

// 数组分组。用 flag 位而非分层结构：构造顺序与生命周期问题最小，
// 且符合「带 vector 的 C」的既定风格。
enum StatePart : unsigned
{
    ST_PARTICLE = 1u,   // Rx..sum_phiz（15 个 N 数组）
    ST_PHI      = 2u,   // fx,fy,fz,eta,etaXY,etaYZ,etaZX（7 个 size 数组）
    ST_VELOCITY = 4u,   // vx,vy,vz,p（4 个 size 数组）
    ST_SOLVER   = 8u,   // pi_d*/pi_n*/tmp_f*/fft/randD/randN + plan + gen
};

constexpr unsigned ST_STENCIL   = ST_PARTICLE | ST_PHI;
constexpr unsigned ST_KINEMATIC = ST_STENCIL  | ST_VELOCITY;
constexpr unsigned ST_FULL      = ST_KINEMATIC| ST_SOLVER;

struct FpdState
{
    NS_Config cfg;
    int    N       = 0;      // 粒子数（可为 0，见 init）
    int    nalloc  = 1;      // = max(N,1)；acc_* 对 0 字节危险，clamp 到 1
    size_t size    = 0;      // Nx*Ny*Nz
    size_t esize   = 0;      // 棱边数组元素数：Nx*Ny*(Nz+1)
    size_t nbEdge  = 0;      // = esize * sizeof(double)
    size_t slotD   = 0;      // randD 元素数（向上取偶，见 Wall.h）
    size_t slotN   = 0;      // randN 元素数（向上取偶）
    unsigned parts = 0;      // 已映射的位集

    // --- 主机侧存储（RAII）---
    std::vector<double> h_Rx, h_Ry, h_Rz, h_Rux, h_Ruy, h_Ruz;
    std::vector<double> h_Vx, h_Vy, h_Vz, h_Fx, h_Fy, h_Fz;
    std::vector<double> h_sum_phix, h_sum_phiy, h_sum_phiz;
    std::vector<double> h_fx, h_fy, h_fz;
    std::vector<double> h_eta, h_etaXY, h_etaYZ, h_etaZX;
    std::vector<double> h_vx, h_vy, h_vz, h_p;
    std::vector<double> h_pi_dx, h_pi_dy, h_pi_dz, h_pi_nx, h_pi_ny, h_pi_nz;
    std::vector<double> h_tmp_fx, h_tmp_fy, h_tmp_fz;
    std::vector<double> h_fft, h_randD, h_randN;
    std::vector<double> h_tri_w;                    // z 向三对角前推系数
    std::vector<double> h_diag;                     // 长度 1：(0,0) 列相容性残差 |Σb̂|

    // --- 裸指针视图（= vector::data()，kernel 参数用）---
    double *Rx=0,*Ry=0,*Rz=0, *Rux=0,*Ruy=0,*Ruz=0;
    double *Vx=0,*Vy=0,*Vz=0, *Fx=0,*Fy=0,*Fz=0;
    double *sum_phix=0,*sum_phiy=0,*sum_phiz=0;
    double *fx=0,*fy=0,*fz=0;
    double *eta=0,*etaXY=0,*etaYZ=0,*etaZX=0;
    double *vx=0,*vy=0,*vz=0,*p=0;
    double *pi_dx=0,*pi_dy=0,*pi_dz=0, *pi_nx=0,*pi_ny=0,*pi_nz=0;
    double *tmp_fx=0,*tmp_fy=0,*tmp_fz=0;
    double *fft=0,*randD=0,*randN=0;
    double *tri_w=0;                                // z 向三对角系数
    double *diag=0;                                 // 长度 1 的诊断量（见 Poisson.h）

    cufftHandle       plan_xy = 0;                  // xy 向 2D 批量 FFT 计划（batch=Nz）
    curandGenerator_t gen  = 0;

    // 分配 + 装配指针 + 设备映射 + 建 plan/gen。
    // N_ 可为 0（纯流体测试），内部 nalloc = max(N_,1)。
    // want 是 ST_* 的位集；依赖关系在 init 里断言（SOLVER ⇒ VELOCITY|PHI）。
    void init(NS_Config cfg_, int N_, unsigned want, unsigned long long seed);
    void finish();                                  // 逆序 delete + destroy，幂等

    // 未映射就 abort 并报「是谁要的」—— 把 present() 遗漏从「跑出垃圾数」
    // 变成「启动即报错」。
    void require(unsigned bits, const char* who) const;

    // 给 IOBin 的数组视图
    CkptArrays ckpt_arrays() const;

    // 拉回 / 推回某组数组（按 StatePart 位）。字节数在这里统一算，调用方不碰。
    void download(unsigned bits, bool pressure = true); // 生产检查点无需回传 p；判据可按需读取
    void upload(unsigned bits);                     // = update device（host→device）
};

#endif
