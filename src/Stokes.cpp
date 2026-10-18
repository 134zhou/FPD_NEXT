#include "Stokes.h"
#include "Poisson.h"

// ============================================================================
// 对流通量 v⊗v 的【唯一实现】。§2a/§2b 用它写 Π，§3a 用它把对流散度从
// 总通量散度里拆出来（变密度时需要单独处理对流项，见 §3a 的说明）。
// 两处必须逐字同源：先算的面心速度（半格平均）再相乘，不能一处用 (a+b)*0.5、
// 另一处用 0.5*(a+b) 之外的任何重排 —— 拆出来的 d_visc 才与 Π 的粘性部分严格一致。
// ============================================================================
#pragma acc routine seq
static inline double conv_xx(double v0, double v1)
{
    const double v = (v0 + v1) * 0.5;
    return v * v;
}

#pragma acc routine seq
static inline double conv_yy(double v0, double v1)
{
    const double v = (v0 + v1) * 0.5;
    return v * v;
}

#pragma acc routine seq
static inline double conv_zz(double v0, double v1)
{
    const double v = (v0 + v1) * 0.5;
    return v * v;
}

#pragma acc routine seq
static inline double conv_xy(double vx0, double vx1, double vy0, double vy1)
{
    return ((vx0 + vx1) * 0.5) * ((vy0 + vy1) * 0.5);
}

#pragma acc routine seq
static inline double conv_xz(double vz0, double vz1, double vx0, double vx1)
{
    return ((vz0 + vz1) * 0.5) * ((vx0 + vx1) * 0.5);
}

#pragma acc routine seq
static inline double conv_yz(double vy0, double vy1, double vz0, double vz1)
{
    return ((vy0 + vy1) * 0.5) * ((vz0 + vz1) * 0.5);
}

int step_navier_stokes
(
    NS_Config cfg,
    double* vx, double* vy, double* vz,
    double* p,
    double* fx, double* fy, double* fz,
    const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
    double* cg_r, double* cg_s, double* cg_zq,
    double* eta,                            // 正应力位置粘度 (Cell Center)
    double* etaXY, double* etaYZ, double* etaZX, // 切应力位置粘度 (Edges)
    double* pi_dx, double* pi_dy, double* pi_dz,
    double* pi_nx, double* pi_ny, double* pi_nz,
    double* fft_data,
    cufftHandle plan_xy,
    const double* tri_w,
    double* diag,
    curandGenerator_t gen,
    double* randD, double* randN,
    double* tmp_fx, double* tmp_fy, double* tmp_fz,
    long step
)
{
    int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    int size = Nx * Ny * Nz;
    double DT = cfg.dt;
    const int esize = wz_edge_size(cfg);   // 棱边数组元素数（Nz+1 层，见 Wall.h）

    // --- 1. 生成随机噪声 (Σ 项) ---
    // randD: 0,1,2 对应 xx, yy, zz 方向（体心，size 个）；randN: 0,1,2 对应 xy, yz, zx
    // 方向（棱边，esize 个）。每次生成前显式定位到本步专属的 slot，不依赖「上次生成后
    // 序列前进了多少」—— 实测那个前进量并非 n，累计记账会失配（详见 Stokes.h）。
    //
    // slot 按两个数组各自的实际长度算。
    const unsigned long long slotD = wz_rand_slot_d(cfg);
    const unsigned long long slotN = wz_rand_slot_n(cfg);
    const unsigned long long base  = (unsigned long long)step * (slotD + slotN);

    // 每一步按绝对 step 重新定位，保证连续运行与重启使用同一段随机数流。
    {
        #pragma acc host_data use_device(randD, randN)
        {
            CURAND_CHECK(curandSetGeneratorOffset(gen, base));
            CURAND_CHECK(curandGenerateNormalDouble(gen, randD, (int)slotD, 0.0, 1.0));
            CURAND_CHECK(curandSetGeneratorOffset(gen, base + slotD));
            CURAND_CHECK(curandGenerateNormalDouble(gen, randN, (int)slotN, 0.0, 1.0));
        }
    }

    // --- 2a. 体心与 XY 棱：k ∈ [0, Nz-1] ---
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, eta, etaXY, randD, randN, pi_dx, pi_dy, pi_dz, pi_nz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ip = (i + 1) % Nx; int im = (i - 1 + Nx) % Nx;
                int jp = (j + 1) % Ny; int jm = (j - 1 + Ny) % Ny;

                // ⚠️ z 向邻居用 (k-1+Nz)%Nz 是【承重的壁面代码，不是周期残留】：
                //    k=0 时读到 vz[Nz-1]，而那个元素被下面 3a 段钉死为 0 —— 正是
                //    下壁 vz[-1]=0 所需的 ghost。所以这两处【一行都不用改】，
                //    也不许「修」成 k-1（会越界）。
                const int km = (k - 1 + Nz) % Nz;
                const int ijk = IDX(i, j, k);

                // 正应力分量：体心 (i,j,k)。对流项用 conv_*（与 §3a 的拆分同源）
                pi_dx[ijk] = conv_xx(vx[ijk], vx[IDX(im, j, k)]);
                pi_dx[ijk] -= eta[ijk] * 2.0 * (vx[ijk] - vx[IDX(im, j, k)]);
                pi_dx[ijk] -= sqrt(2.0 * eta[ijk]) * cfg.W * randD[0 * size + ijk];

                pi_dy[ijk] = conv_yy(vy[ijk], vy[IDX(i, jm, k)]);
                pi_dy[ijk] -= eta[ijk] * 2.0 * (vy[ijk] - vy[IDX(i, jm, k)]);
                pi_dy[ijk] -= sqrt(2.0 * eta[ijk]) * cfg.W * randD[1 * size + ijk];

                pi_dz[ijk] = conv_zz(vz[ijk], vz[IDX(i, j, km)]);
                pi_dz[ijk] -= eta[ijk] * 2.0 * (vz[ijk] - vz[IDX(i, j, km)]);
                pi_dz[ijk] -= sqrt(2.0 * eta[ijk]) * cfg.W * randD[2 * size + ijk];

                // Πxy (pi_nz)：XY 棱 (i+1/2,j+1/2,k)
                pi_nz[ijk] = conv_xy(vx[ijk], vx[IDX(i, jp, k)],
                                     vy[ijk], vy[IDX(ip, j, k)]);
                pi_nz[ijk] -= etaXY[ijk] * ((vx[IDX(i, jp, k)] - vx[ijk])
                                           + (vy[IDX(ip, j, k)] - vy[ijk]));
                pi_nz[ijk] -= sqrt(etaXY[ijk]) * cfg.W * randN[0 * esize + ijk];
            }
        }
    }

    // --- 2b. 内部 YZ/ZX 棱：k ∈ [0, Nz-2] ---
    // 这一段不触及 ghost，所有速度都能用普通 IDX 直接访问，噪声控制体完整（γ=1）。
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, etaYZ, etaZX, randN, pi_nx, pi_ny)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz - 1; k++)
            {
                int ip = (i + 1) % Nx;
                int jp = (j + 1) % Ny;
                const int ijk = IDX(i, j, k);
                const int eb  = wz_edge_idx(cfg, i,  j,  k);

                // Πyz (pi_nx)：YZ 棱 (i,j+1/2,k+1/2)
                pi_nx[eb] = conv_yz(vy[ijk], vy[IDX(i, j, k + 1)],
                                    vz[ijk], vz[IDX(i, jp, k)]);
                pi_nx[eb] -= etaYZ[eb] * ((vy[IDX(i, j, k + 1)] - vy[ijk])
                                        + (vz[IDX(i, jp, k)] - vz[ijk]));
                pi_nx[eb] -= sqrt(etaYZ[eb]) * cfg.W * randN[1 * esize + eb];

                // Πzx (pi_ny)：ZX 棱 (i+1/2,j,k+1/2)
                pi_ny[eb] = conv_xz(vz[ijk], vz[IDX(ip, j, k)],
                                    vx[ijk], vx[IDX(i, j, k + 1)]);
                pi_ny[eb] -= etaZX[eb] * ((vz[IDX(ip, j, k)] - vz[ijk])
                                        + (vx[IDX(i, j, k + 1)] - vx[ijk]));
                pi_ny[eb] -= sqrt(etaZX[eb]) * cfg.W * randN[2 * esize + eb];
            }
        }
    }

    // --- 2c. 两片壁面的 YZ/ZX 棱：k=-1 与 k=Nz-1 ---
    // 无滑移 ghost：vα[-1]=-vα[0]，vα[Nz]=-vα[Nz-1]（α=x,y），法向 vz=0。
    // 代入一般式后，两侧对流与横向 vz 导数都为 0；半格单侧导数分别为
    // +2vα[0] 与 -2vα[Nz-1]。壁面控制体只有一半，故 γ=2（幅度 ×√2）。
    const double wall_gamma = cfg.noise_gamma1 ? 1.0 : 2.0;
    #pragma acc parallel loop collapse(2) present(vx, vy, etaYZ, etaZX, randN, pi_nx, pi_ny)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            const int bottom = IDX(i, j, 0);
            const int top = IDX(i, j, Nz - 1);
            const int eb_bottom = wz_edge_idx(cfg, i, j, -1);
            const int eb_top = wz_edge_idx(cfg, i, j, Nz - 1);

            pi_nx[eb_bottom] = -etaYZ[eb_bottom] * (2.0 * vy[bottom]);
            pi_nx[eb_bottom] -= sqrt(etaYZ[eb_bottom] * wall_gamma) * cfg.W
                              * randN[1 * esize + eb_bottom];
            pi_nx[eb_top] = -etaYZ[eb_top] * (-2.0 * vy[top]);
            pi_nx[eb_top] -= sqrt(etaYZ[eb_top] * wall_gamma) * cfg.W
                           * randN[1 * esize + eb_top];

            pi_ny[eb_bottom] = -etaZX[eb_bottom] * (2.0 * vx[bottom]);
            pi_ny[eb_bottom] -= sqrt(etaZX[eb_bottom] * wall_gamma) * cfg.W
                              * randN[2 * esize + eb_bottom];
            pi_ny[eb_top] = -etaZX[eb_top] * (-2.0 * vx[top]);
            pi_ny[eb_top] -= sqrt(etaZX[eb_top] * wall_gamma) * cfg.W
                           * randN[2 * esize + eb_top];
        }
    }

    // --- 3a. 显式一步：v* = v + dt·[ -(v·∇)v + a·(-∇·τ_vis - ∇·noise + f) ] ---
    // 每个面心只算【一次】(∇·Π)_α，直接写 tmp_fα = v*（tmp_f* 复用为 v* 的存储，
    // 零额外显存）。原实现把每个面算了两次（d_pix_ijk 与左邻的 d_pix_im 是同一个量），
    // 这里消掉重复，Π 的读取减半 —— 正对 README 记录的缓存瓶颈。
    //
    // ⚠️ 1/rho 只乘【力密度】，不乘对流加速度。ρDv/Dt 两边除以 ρ 后对流项在左边、
    //    系数是 1；把整个 ∇·Π 都乘 a 会得到 ρ∂t v + ∇·(vv) = ...，既不是非守恒形式
    //    也不是守恒形式（tmp.txt 的 V0 反例：均匀 a 下对流加速度会被错误地缩放）。
    //    d_pix 仍是总通量散度（W6 的动量收支用它的 Π_zz），这里把它拆成
    //    d_adv（对流）与 d_visc = d_pix - d_adv（粘性+噪声）。
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, fx, fy, fz, inv_rho_x, inv_rho_y, inv_rho_z, pi_dx, pi_dy, pi_dz, pi_nx, pi_ny, pi_nz, tmp_fx, tmp_fy, tmp_fz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                int ip = (i + 1) % Nx; int im = (i - 1 + Nx) % Nx;
                int jp = (j + 1) % Ny; int jm = (j - 1 + Ny) % Ny;
                int kp = (k + 1) % Nz;

                // 棱边量的 z 向邻居走逻辑索引：k=0 的下面那层是【下壁】（逻辑 k=-1），
                // 不是回绕过来的上壁。搞错会静默用错壁面剪切。
                const int eb   = wz_edge_idx(cfg, i,  j,  k);
                const int ebim = wz_edge_idx(cfg, im, j,  k);
                const int ebjm = wz_edge_idx(cfg, i,  jm, k);
                const int ebkm = wz_edge_idx(cfg, i,  j,  wz_edge_km(cfg, k));

                // (∇·Π)_x，物理位置 (i+½, j, k)：d(Πxx)/dx + d(Πxy)/dy + d(Πxz)/dz
                double d_pix = (pi_dx[IDX(ip, j, k)] - pi_dx[ijk])
                             + (pi_nz[ijk] - pi_nz[IDX(i, jm, k)])
                             + (pi_ny[eb] - pi_ny[ebkm]);

                // (∇·Π)_y，物理位置 (i, j+½, k)
                double d_piy = (pi_nz[ijk] - pi_nz[IDX(im, j, k)])
                             + (pi_dy[IDX(i, jp, k)] - pi_dy[ijk])
                             + (pi_nx[eb] - pi_nx[ebkm]);

                // (∇·Π)_z，物理位置 (i, j, k+½)
                double d_piz = (pi_ny[eb] - pi_ny[ebim])
                             + (pi_nx[eb] - pi_nx[ebjm])
                             + (pi_dz[IDX(i, j, kp)] - pi_dz[ijk]);

                // --- 变密度：对流项【不】乘 1/rho ---
                // 非守恒形式 rho(∂t v + (v·∇)v) = -∇p + ∇·τ + f 除以 rho 后是
                //   ∂t v = -(v·∇)v + a(∇·τ + ∇·noise + f)
                // 对流项是加速度，不带 a；只有粘性/噪声/外力这些力密度带 a。
                // 代码里 d_pi = ∇·Π = d_adv + d_visc（d_visc = ∇·(-τ - noise)），
                // 所以 v* = v + dt[ -d_adv + a(-d_visc + f) ]。
                // density_on=0 时保留原来的合并写法（与 HEAD 逐位相同）。
                if (cfg.density_on)
                {
                    const int km = (k - 1 + Nz) % Nz;
                    const double cxx_i  = conv_xx(vx[ijk], vx[IDX(im, j, k)]);
                    const double cxx_ip = conv_xx(vx[IDX(ip, j, k)], vx[ijk]);
                    const double cyy_i  = conv_yy(vy[ijk], vy[IDX(i, jm, k)]);
                    const double cyy_jp = conv_yy(vy[IDX(i, jp, k)], vy[ijk]);
                    const double czz_i  = conv_zz(vz[ijk], vz[IDX(i, j, km)]);
                    const double czz_kp = conv_zz(vz[IDX(i, j, kp)], vz[ijk]);

                    const double cxy_i  = conv_xy(vx[ijk], vx[IDX(i, jp, k)],
                                                  vy[ijk], vy[IDX(ip, j, k)]);
                    const double cxy_jm = conv_xy(vx[IDX(i, jm, k)], vx[ijk],
                                                  vy[IDX(i, jm, k)], vy[IDX(ip, jm, k)]);
                    const double cxy_im = conv_xy(vx[IDX(im, j, k)], vx[IDX(im, jp, k)],
                                                  vy[IDX(im, j, k)], vy[ijk]);

                    // ZX 棱（pi_ny）：eb 在 k、ebkm 在 k-1、ebim 在 (i-1,j,k)。
                    // 壁面棱边（逻辑 k=-1 / Nz-1）的对流项恒为 0（§2c 只写粘性+噪声）。
                    const double cxz_eb   = (k == Nz - 1) ? 0.0
                        : conv_xz(vz[ijk], vz[IDX(ip, j, k)], vx[ijk], vx[IDX(i, j, k + 1)]);
                    const double cxz_ebkm = (k == 0) ? 0.0
                        : conv_xz(vz[IDX(i, j, k - 1)], vz[IDX(ip, j, k - 1)],
                                  vx[IDX(i, j, k - 1)], vx[ijk]);
                    const double cxz_ebim = (k == Nz - 1) ? 0.0
                        : conv_xz(vz[IDX(im, j, k)], vz[ijk],
                                  vx[IDX(im, j, k)], vx[IDX(im, j, k + 1)]);

                    // YZ 棱（pi_nx）：eb 在 k、ebkm 在 k-1、ebjm 在 (i,j-1,k)
                    const double cyz_eb   = (k == Nz - 1) ? 0.0
                        : conv_yz(vy[ijk], vy[IDX(i, j, k + 1)], vz[ijk], vz[IDX(i, jp, k)]);
                    const double cyz_ebkm = (k == 0) ? 0.0
                        : conv_yz(vy[IDX(i, j, k - 1)], vy[ijk],
                                  vz[IDX(i, j, k - 1)], vz[IDX(i, jp, k - 1)]);
                    const double cyz_ebjm = (k == Nz - 1) ? 0.0
                        : conv_yz(vy[IDX(i, jm, k)], vy[IDX(i, jm, k + 1)],
                                  vz[IDX(i, jm, k)], vz[ijk]);

                    const double adv_x = (cxx_ip - cxx_i) + (cxy_i - cxy_jm) + (cxz_eb - cxz_ebkm);
                    const double adv_y = (cxy_i - cxy_im) + (cyy_jp - cyy_i) + (cyz_eb - cyz_ebkm);
                    const double adv_z = (cxz_eb - cxz_ebim) + (cyz_eb - cyz_ebjm) + (czz_kp - czz_i);

                    const double dvx = d_pix - adv_x;
                    const double dvy = d_piy - adv_y;
                    const double dvz = d_piz - adv_z;

                    tmp_fx[ijk] = vx[ijk] + DT * (-adv_x + inv_rho_x[ijk] * (-dvx + fx[ijk]));
                    tmp_fy[ijk] = vy[ijk] + DT * (-adv_y + inv_rho_y[ijk] * (-dvy + fy[ijk]));
                    tmp_fz[ijk] = vz[ijk] + DT * (-adv_z + inv_rho_z[ijk] * (-dvz + fz[ijk]));
                }
                else
                {
                    tmp_fx[ijk] = vx[ijk] + DT * (-d_pix + fx[ijk]);
                    tmp_fy[ijk] = vy[ijk] + DT * (-d_piy + fy[ijk]);
                    tmp_fz[ijk] = vz[ijk] + DT * (-d_piz + fz[ijk]);
                }

                // 壁面法向面【永不修正】：v*z 在两壁钉死为 0。
                // ⚠️ 这个置 0 是【承重的】：上面体心段的 km=(k-1+Nz)%Nz 与下面
                //    散度段的 km 都靠它把「绕回来的上一层」变成正确的 ghost
                //    v*z[-1]=0。这正是 doc/PressurePoisson.md §8a 相容性条件
                //    Σ_k b̂_k = 0 的前提，动了它 W1/W2/W6 一起破。
                if (wz_pinned_zface(cfg, k)) { tmp_fz[ijk] = 0.0; }
            }
        }
    }

    // --- 3b. 泊松右端 b = (1/dt)·∇·v* ---
    // 显式 v* 的意义：壁面条件就落在「v*z 在两壁被钉死」这一句话上，泊松的相容性
    // （Σ_k b̂_k = 0）随之在精确算术下自动成立（见 doc/PressurePoisson.md §8a）。
    // 从面心差分到体心：v*x[ijk] 是右面，tmp_fx[im] 是左面。
    #pragma acc parallel loop collapse(3) present(tmp_fx, tmp_fy, tmp_fz, fft_data)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                int im = (i - 1 + Nx) % Nx;
                int jm = (j - 1 + Ny) % Ny;
                int km = (k - 1 + Nz) % Nz;

                double div_vs = (tmp_fx[ijk] - tmp_fx[IDX(im, j, k)])
                              + (tmp_fy[ijk] - tmp_fy[IDX(i, jm, k)])
                              + (tmp_fz[ijk] - tmp_fz[IDX(i, j, km)]);

                fft_data[ijk * 2]     = cfg.inv_dt * div_vs;
                fft_data[ijk * 2 + 1] = 0.0;
            }
        }
    }

    // --- 4. 求解压力泊松方程 ---
    // 常密度：div(grad p) = b，xy 2D 批量 FFT + z 向 Thomas，直接解。
    // 变密度：div(a grad p) = b，用同一个 FFT 求解器作预条件器跑 PCG。
    // 数学推导与公式↔代码对照见 doc/PressurePoisson.md（变系数见 §16）。
    int cg_iter = 0;
    if (cfg.density_on)
    {
        // 把右端（fft_data 的实部）拷进 cg_r；solve_pressure_var 入口约定 = b。
        #pragma acc parallel loop collapse(3) present(fft_data, cg_r)
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    const int ijk = IDX(i, j, k);
                    cg_r[ijk] = fft_data[2 * ijk];
                }
            }
        }
        cg_iter = solve_pressure_var(cfg, cg_r, inv_rho_x, inv_rho_y, inv_rho_z,
                                     fft_data, tri_w, plan_xy, p, diag, cg_s, cg_zq);
    }
    else
    {
        solve_pressure(cfg, fft_data, tri_w, plan_xy, p, diag);
    }

    // --- 5. 修正步 (Correction / projection) ---
    // v^{n+1} = v* - dt·a·∇p，其中 v* 就是 tmp_f*，a = 1/rho 在面心。
    // 必须与上面的求解分成两个 kernel：p 的梯度要读相邻线程写的 p。
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, p, inv_rho_x, inv_rho_y, inv_rho_z, tmp_fx, tmp_fy, tmp_fz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                int ip = (i + 1) % Nx;
                int jp = (j + 1) % Ny;
                int kp = (k + 1) % Nz;

                vx[ijk] = tmp_fx[ijk] - DT * inv_rho_x[ijk] * (p[IDX(ip, j, k)] - p[ijk]);
                vy[ijk] = tmp_fy[ijk] - DT * inv_rho_y[ijk] * (p[IDX(i, jp, k)] - p[ijk]);
                // 壁面法向面永不修正：vz 在那两处恒 0。k=Nz-1 时 kp 会回绕到 0，
                // 那是错的，所以这一层干脆不写（vz[Nz-1] 保持初值 0）。
                if (!wz_pinned_zface(cfg, k))
                {
                    vz[ijk] = tmp_fz[ijk] - DT * inv_rho_z[ijk] * (p[IDX(i, j, kp)] - p[ijk]);
                }
            }
        }
    }

    return cg_iter;
}
