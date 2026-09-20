#include "State.h"

#include <cstdio>
#include <cstdlib>
#include <openacc.h>
#include "Check.h"
#include "Poisson.h"

// 运行时 API 的映射建立/删除（S0-A spike 已验证这些映射能被 present() 认账）。
// 注意字节数一律用 (size_t) 显式转换，size*sizeof(double) 是 size_t 运算。

static void api_create(double* p, size_t nbytes)  { acc_create(p, nbytes); }
static void api_copyin(double* p, size_t nbytes)  { acc_copyin(p, nbytes); }
static void api_delete(double* p, size_t nbytes)  { acc_delete(p, nbytes); }
static void api_down(double* p, size_t nbytes)    { acc_update_self(p, nbytes); }
static void api_up(double* p, size_t nbytes)      { acc_update_device(p, nbytes); }

void FpdState::init(NS_Config cfg_, int N_, unsigned want, unsigned long long seed)
{
    cfg  = cfg_;
    N    = N_;
    nalloc = N_ > 1 ? N_ : 1;          // acc_* 对 0 字节危险
    size = (size_t)cfg.Nx * (size_t)cfg.Ny * (size_t)cfg.Nz;
    parts = want;

    const size_t nbN    = (size_t)nalloc * sizeof(double);
    const size_t nbSize = size * sizeof(double);

    // 棱边数组（etaYZ/etaZX/pi_nx/pi_ny/randN 的两个分量）多一层：
    // 逻辑 k ∈ [-1, Nz-1] 共 Nz+1 层（含 z=∓1/2 两片壁面棱边）。见 Wall.h。
    esize = (size_t)wz_edge_size(cfg);
    nbEdge = esize * sizeof(double);
    slotD  = (size_t)wz_rand_slot_d(cfg);   // randD 长度（偶数）
    slotN  = (size_t)wz_rand_slot_n(cfg);   // randN 长度（偶数）

    // 依赖关系：SOLVER 需要 VELOCITY 与 PHI 的场（step_navier_stokes 的参数）
    if (want & ST_SOLVER)
    {
        if (!(want & ST_VELOCITY) || !(want & ST_PHI))
        {
            fprintf(stderr, "FpdState::init 内部错误：ST_SOLVER 需要 ST_VELOCITY|ST_PHI\n");
            std::abort();
        }
    }

    // --- 分配主机存储 + 装配裸指针 ---
    if (want & ST_PARTICLE)
    {
        h_Rx.assign(nalloc, 0); h_Ry.assign(nalloc, 0); h_Rz.assign(nalloc, 0);
        h_Rux.assign(nalloc, 0); h_Ruy.assign(nalloc, 0); h_Ruz.assign(nalloc, 0);
        h_Vx.assign(nalloc, 0); h_Vy.assign(nalloc, 0); h_Vz.assign(nalloc, 0);
        h_Fx.assign(nalloc, 0); h_Fy.assign(nalloc, 0); h_Fz.assign(nalloc, 0);
        h_sum_phix.assign(nalloc, 0); h_sum_phiy.assign(nalloc, 0); h_sum_phiz.assign(nalloc, 0);
        Rx = h_Rx.data(); Ry = h_Ry.data(); Rz = h_Rz.data();
        Rux = h_Rux.data(); Ruy = h_Ruy.data(); Ruz = h_Ruz.data();
        Vx = h_Vx.data(); Vy = h_Vy.data(); Vz = h_Vz.data();
        Fx = h_Fx.data(); Fy = h_Fy.data(); Fz = h_Fz.data();
        sum_phix = h_sum_phix.data(); sum_phiy = h_sum_phiy.data(); sum_phiz = h_sum_phiz.data();
    }
    if (want & ST_PHI)
    {
        h_fx.assign(size, 0); h_fy.assign(size, 0); h_fz.assign(size, 0);
        h_eta.assign(size, 0); h_etaXY.assign(size, 0);
        h_etaYZ.assign(esize, 0); h_etaZX.assign(esize, 0);   // 棱边：Nz+1 层
        fx = h_fx.data(); fy = h_fy.data(); fz = h_fz.data();
        eta = h_eta.data(); etaXY = h_etaXY.data(); etaYZ = h_etaYZ.data(); etaZX = h_etaZX.data();
    }
    if (want & ST_VELOCITY)
    {
        h_vx.assign(size, 0); h_vy.assign(size, 0); h_vz.assign(size, 0); h_p.assign(size, 0);
        vx = h_vx.data(); vy = h_vy.data(); vz = h_vz.data(); p = h_p.data();
    }
    if (want & ST_SOLVER)
    {
        h_pi_dx.assign(size, 0); h_pi_dy.assign(size, 0); h_pi_dz.assign(size, 0);
        h_pi_nx.assign(esize, 0); h_pi_ny.assign(esize, 0); h_pi_nz.assign(size, 0);
        h_tmp_fx.assign(size, 0); h_tmp_fy.assign(size, 0); h_tmp_fz.assign(size, 0);
        h_fft.assign(size * 2, 0);
        h_randD.assign(slotD, 0); h_randN.assign(slotN, 0);
        h_diag.assign(1, 0);
        pi_dx = h_pi_dx.data(); pi_dy = h_pi_dy.data(); pi_dz = h_pi_dz.data();
        pi_nx = h_pi_nx.data(); pi_ny = h_pi_ny.data(); pi_nz = h_pi_nz.data();
        tmp_fx = h_tmp_fx.data(); tmp_fy = h_tmp_fy.data(); tmp_fz = h_tmp_fz.data();
        fft = h_fft.data(); randD = h_randD.data(); randN = h_randN.data();
        diag = h_diag.data();

        // z 向三对角前推系数（主机侧一次性预算，与右端无关）
        h_tri_w.assign(size, 0);
        tri_w = h_tri_w.data();
        build_tridiag_coeffs(cfg, tri_w);
    }

    // --- 设备映射 ---
    if (want & ST_PARTICLE)
    {
        api_copyin(Rx, nbN); api_copyin(Ry, nbN); api_copyin(Rz, nbN);
        api_copyin(Rux, nbN); api_copyin(Ruy, nbN); api_copyin(Ruz, nbN);
        api_copyin(Vx, nbN); api_copyin(Vy, nbN); api_copyin(Vz, nbN);
        api_copyin(Fx, nbN); api_copyin(Fy, nbN); api_copyin(Fz, nbN);
        api_create(sum_phix, nbN); api_create(sum_phiy, nbN); api_create(sum_phiz, nbN);
    }
    if (want & ST_PHI)
    {
        api_copyin(fx, nbSize); api_copyin(fy, nbSize); api_copyin(fz, nbSize);
        api_create(eta, nbSize); api_create(etaXY, nbSize);
        api_create(etaYZ, nbEdge); api_create(etaZX, nbEdge);
    }
    if (want & ST_VELOCITY)
    {
        api_copyin(vx, nbSize); api_copyin(vy, nbSize); api_copyin(vz, nbSize); api_copyin(p, nbSize);
    }
    if (want & ST_SOLVER)
    {
        api_create(pi_dx, nbSize); api_create(pi_dy, nbSize); api_create(pi_dz, nbSize);
        api_create(pi_nx, nbEdge); api_create(pi_ny, nbEdge); api_create(pi_nz, nbSize);
        api_create(tmp_fx, nbSize); api_create(tmp_fy, nbSize); api_create(tmp_fz, nbSize);
        api_create(fft, size * 2 * sizeof(double));
        api_create(randD, slotD * sizeof(double));
        api_create(randN, slotN * sizeof(double));
        api_create(diag, sizeof(double));       // 单元素诊断量，create 即可（求解器整体覆写）
        api_copyin(tri_w, nbSize);
    }

    // --- FFT plan + RNG ---
    if (want & ST_SOLVER)
    {
        // xy 向批量 2D FFT（batch=Nz）。cuFFT row-major，n[] 最后一维最快，
        // 我们 x 最快 ⇒ n={Ny,Nx}。z-slab 在 IDX 布局下天然连续（idist=Nx*Ny），无需重排。
        int n[2] = { cfg.Ny, cfg.Nx };
        CUFFT_CHECK(cufftPlanMany(&plan_xy, 2, n, NULL, 1, cfg.Nx * cfg.Ny,
                                  NULL, 1, cfg.Nx * cfg.Ny, CUFFT_Z2Z, cfg.Nz));
        // 【统一 Philox】3A/3B 原来用 XORWOW（CURAND_RNG_PSEUDO_DEFAULT），
        // 生产用 Philox。CMakeLists 要求验证与生产走同一代码路径，且 XORWOW
        // 无法逐位重启（见 PROGRESS.md），故统一为 Philox。
        CURAND_CHECK(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_PHILOX4_32_10));
        CURAND_CHECK(curandSetPseudoRandomGeneratorSeed(gen, seed));
    }
}

void FpdState::finish()
{
    const size_t nbN    = (size_t)nalloc * sizeof(double);
    const size_t nbSize = size * sizeof(double);

    // 逆序 delete；acc_delete 幂等（不在 present 表时是 no-op）
    if (parts & ST_SOLVER)
    {
        api_delete(diag, sizeof(double));
        api_delete(randN, slotN * sizeof(double));
        api_delete(randD, slotD * sizeof(double));
        api_delete(fft, size * 2 * sizeof(double));
        api_delete(tri_w, nbSize); tri_w = 0;
        api_delete(tmp_fz, nbSize); api_delete(tmp_fy, nbSize); api_delete(tmp_fx, nbSize);
        api_delete(pi_nz, nbSize); api_delete(pi_ny, nbEdge); api_delete(pi_nx, nbEdge);
        api_delete(pi_dz, nbSize); api_delete(pi_dy, nbSize); api_delete(pi_dx, nbSize);
    }
    if (parts & ST_VELOCITY)
    {
        api_delete(p, nbSize); api_delete(vz, nbSize); api_delete(vy, nbSize); api_delete(vx, nbSize);
    }
    if (parts & ST_PHI)
    {
        api_delete(etaZX, nbEdge); api_delete(etaYZ, nbEdge); api_delete(etaXY, nbSize); api_delete(eta, nbSize);
        api_delete(fz, nbSize); api_delete(fy, nbSize); api_delete(fx, nbSize);
    }
    if (parts & ST_PARTICLE)
    {
        api_delete(sum_phiz, nbN); api_delete(sum_phiy, nbN); api_delete(sum_phix, nbN);
        api_delete(Fz, nbN); api_delete(Fy, nbN); api_delete(Fx, nbN);
        api_delete(Vz, nbN); api_delete(Vy, nbN); api_delete(Vx, nbN);
        api_delete(Ruz, nbN); api_delete(Ruy, nbN); api_delete(Rux, nbN);
        api_delete(Rz, nbN); api_delete(Ry, nbN); api_delete(Rx, nbN);
    }
    if (plan_xy) { CUFFT_CHECK(cufftDestroy(plan_xy)); plan_xy = 0; }
    if (gen)  { CURAND_CHECK(curandDestroyGenerator(gen)); gen = 0; }
    parts = 0;
}

void FpdState::require(unsigned bits, const char* who) const
{
    if ((parts & bits) != bits)
    {
        fprintf(stderr, "FpdState::require 失败：%s 需要 0x%x，但只映射了 0x%x\n",
                who, bits, parts);
        std::abort();
    }
}

CkptArrays FpdState::ckpt_arrays() const
{
    CkptArrays a;
    a.vx = vx; a.vy = vy; a.vz = vz;
    a.Rx = Rx; a.Ry = Ry; a.Rz = Rz;
    a.Rux = Rux; a.Ruy = Ruy; a.Ruz = Ruz;
    a.Vx = Vx; a.Vy = Vy; a.Vz = Vz;
    a.Fx = Fx; a.Fy = Fy; a.Fz = Fz;
    return a;
}

void FpdState::download(unsigned bits, bool pressure)
{
    const size_t nbN    = (size_t)nalloc * sizeof(double);
    const size_t nbSize = size * sizeof(double);

    if (bits & ST_PARTICLE)
    {
        api_down(Rx, nbN); api_down(Ry, nbN); api_down(Rz, nbN);
        api_down(Rux, nbN); api_down(Ruy, nbN); api_down(Ruz, nbN);
        api_down(Vx, nbN); api_down(Vy, nbN); api_down(Vz, nbN);
        api_down(Fx, nbN); api_down(Fy, nbN); api_down(Fz, nbN);
        api_down(sum_phix, nbN); api_down(sum_phiy, nbN); api_down(sum_phiz, nbN);
    }
    if (bits & ST_PHI)
    {
        api_down(fx, nbSize); api_down(fy, nbSize); api_down(fz, nbSize);
        api_down(eta, nbSize); api_down(etaXY, nbSize);
        api_down(etaYZ, nbEdge); api_down(etaZX, nbEdge);
    }
    if (bits & ST_VELOCITY)
    {
        api_down(vx, nbSize); api_down(vy, nbSize); api_down(vz, nbSize);
        if (pressure) { api_down(p, nbSize); }
    }
    // ST_SOLVER 的临时场不下载（仅供 step_navier_stokes 内部使用）；
    // 例外是 diag：它是壁面相容性诊断量，生产运行要定期读回监控（判据 W8）。
    if (bits & ST_SOLVER) { api_down(diag, sizeof(double)); }
}

void FpdState::upload(unsigned bits)
{
    const size_t nbN    = (size_t)nalloc * sizeof(double);
    const size_t nbSize = size * sizeof(double);

    if (bits & ST_PARTICLE)
    {
        api_up(Rx, nbN); api_up(Ry, nbN); api_up(Rz, nbN);
        api_up(Rux, nbN); api_up(Ruy, nbN); api_up(Ruz, nbN);
        api_up(Vx, nbN); api_up(Vy, nbN); api_up(Vz, nbN);
        api_up(Fx, nbN); api_up(Fy, nbN); api_up(Fz, nbN);
        api_up(sum_phix, nbN); api_up(sum_phiy, nbN); api_up(sum_phiz, nbN);
    }
    if (bits & ST_PHI)
    {
        api_up(fx, nbSize); api_up(fy, nbSize); api_up(fz, nbSize);
        api_up(eta, nbSize); api_up(etaXY, nbSize);
        api_up(etaYZ, nbEdge); api_up(etaZX, nbEdge);
    }
    if (bits & ST_VELOCITY)
    {
        api_up(vx, nbSize); api_up(vy, nbSize); api_up(vz, nbSize); api_up(p, nbSize);
    }
}
