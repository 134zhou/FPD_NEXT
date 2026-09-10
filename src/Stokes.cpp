#include "./include/Stokes.h"
#include "./include/Poisson.h"

void step_navier_stokes
(
    NS_Config cfg,
    double* vx, double* vy, double* vz,
    double* p,
    double* fx, double* fy, double* fz,
    double* eta,                            // 正应力位置粘度 (Cell Center)
    double* etaXY, double* etaYZ, double* etaZX, // 切应力位置粘度 (Edges)
    double* pi_dx, double* pi_dy, double* pi_dz,
    double* pi_nx, double* pi_ny, double* pi_nz,
    double* fft_data,
    cufftHandle plan,
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
    const int wz   = cfg.wall_z;      // 0 = 周期，1 = z 向无滑移壁面
    const int esize = wz_edge_size(cfg);   // 棱边数组元素数（壁面模式是 Nz+1 层）

    // --- 1. 生成随机噪声 (Σ 项) ---
    // randD: 0,1,2 对应 xx, yy, zz 方向（体心，size 个）；randN: 0,1,2 对应 xy, yz, zx
    // 方向（棱边，esize 个）。每次生成前显式定位到本步专属的 slot，不依赖「上次生成后
    // 序列前进了多少」—— 实测那个前进量并非 n，累计记账会失配（详见 Stokes.h）。
    //
    // slot 按两个数组各自的实际长度算：周期模式下 esize == size ⇒ 与
    // (2*step+c)*3*size 逐位相同（Phase 7-B 之前的行为原样保住）。
    const unsigned long long slotD = wz_rand_slot_d(cfg);
    const unsigned long long slotN = wz_rand_slot_n(cfg);
    const unsigned long long base  = (unsigned long long)step * (slotD + slotN);

    #pragma acc host_data use_device(randD, randN)
    {
        CURAND_CHECK(curandSetGeneratorOffset(gen, base));
        CURAND_CHECK(curandGenerateNormalDouble(gen, randD, (int)slotD, 0.0, 1.0));
        CURAND_CHECK(curandSetGeneratorOffset(gen, base + slotD));
        CURAND_CHECK(curandGenerateNormalDouble(gen, randN, (int)slotN, 0.0, 1.0));
    }

    // --- 2. 计算总动量通量 Π = Advection - Viscosity - Noise ---
    // 逻辑 k 统一遍历：体心量在 k ∈ [0,Nz-1]；棱边量 Π_yz/Π_zx 在 k ∈ [-1,Nz-1]
    // （壁面模式下多出 k=-1 这一层 = 下壁棱边，上壁棱边在 k=Nz-1）。
    // 周期模式下 wz=0 ⇒ kk≡k、范围不变，与 Phase 7-B 之前逐位一致。
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, eta, etaXY, etaYZ, etaZX, randD, randN, pi_dx, pi_dy, pi_dz, pi_nx, pi_ny, pi_nz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int kk = 0; kk < Nz + wz; kk++)
            {
                const int k = kk - wz;
                int ip = (i + 1) % Nx; int im = (i - 1 + Nx) % Nx;
                int jp = (j + 1) % Ny; int jm = (j - 1 + Ny) % Ny;

                // ---------------- 体心量：仅 k ∈ [0, Nz-1] ----------------
                // 注意 z 向邻居仍用 (k-1+Nz)%Nz：壁面模式下 k=0 时 vaz 与 Π_zz 的
                // 粘性项会读到 vz[Nz-1]，而那个元素被钉死为 0 —— 正是下壁
                // vz[-1]=0 所需的 ghost。所以这两处【一行都不用改】。
                if (k >= 0)
                {
                    const int ijk = IDX(i, j, k);
                    const int kp = (k + 1) % Nz; const int km = (k - 1 + Nz) % Nz;

                    // --- 正应力分量 (Diagonal terms: Πxx, Πyy, Πzz) ---
                    // 物理位置：Cell Center (i,j,k)
                    double vax = (vx[ijk] + vx[IDX(im, j, k)]) * 0.5;
                    double vay = (vy[ijk] + vy[IDX(i, jm, k)]) * 0.5;
                    double vaz = (vz[ijk] + vz[IDX(i, j, km)]) * 0.5;

                    pi_dx[ijk] = vax * vax; // Advection: v_i * v_j
                    pi_dx[ijk] -= eta[ijk] * 2.0 * (vx[ijk] - vx[IDX(im, j, k)]); // Viscosity: 2*eta*Exx
                    pi_dx[ijk] -= sqrt(2.0 * eta[ijk]) * cfg.W * randD[0 * size + ijk]; // Noise: Σxx

                    pi_dy[ijk] = vay * vay;
                    pi_dy[ijk] -= eta[ijk] * 2.0 * (vy[ijk] - vy[IDX(i, jm, k)]);
                    pi_dy[ijk] -= sqrt(2.0 * eta[ijk]) * cfg.W * randD[1 * size + ijk];

                    pi_dz[ijk] = vaz * vaz;
                    pi_dz[ijk] -= eta[ijk] * 2.0 * (vz[ijk] - vz[IDX(i, j, km)]);
                    pi_dz[ijk] -= sqrt(2.0 * eta[ijk]) * cfg.W * randD[2 * size + ijk];

                    // Πxy (pi_nz): 位于 XY 面的棱边（体心层，不含 z 导数，控制体完整）
                    double vsx_y = (vx[ijk] + vx[IDX(i, jp, k)]) * 0.5;
                    double vsy_x = (vy[ijk] + vy[IDX(ip, j, k)]) * 0.5;
                    pi_nz[ijk] = vsx_y * vsy_x; // 对流
                    pi_nz[ijk] -= etaXY[ijk] * ((vx[IDX(i, jp, k)] - vx[ijk]) + (vy[IDX(ip, j, k)] - vy[ijk])); // 粘性
                    pi_nz[ijk] -= sqrt(etaXY[ijk]) * cfg.W * randN[0 * esize + ijk]; // 随机
                }

                // ---------------- 棱边量：k ∈ [-1, Nz-1] ----------------
                // 这两个量在 z 方向【没有】周期回绕（壁面模式下不环绕到另一片壁），
                // 所以要用 wz_edge_idx 的逻辑索引 + wz_vx/wz_vy/wz_vz 的 ghost 取值。
                // 代入 ghost 后表达式不变就自动给出正确结果：
                //   k=-1（下壁）对流两个因子都是 0；粘性 (vx[0]-vx[-1]) = 2vx[0]。
                const int eb  = wz_edge_idx(cfg, i,  j,  k);
                const int ebp = wz_edge_idx(cfg, ip, j,  k);
                const int ebj = wz_edge_idx(cfg, i,  jp, k);
                const double gam = wz_noise_gamma(cfg, k);   // 壁面棱边 = 2（幅度 ×√2）

                // Πyz (pi_nx): 位于 YZ 面的棱边
                double vsy_z = (wz_vy(cfg, vy, i, j, k) + wz_vy(cfg, vy, i, j, k + 1)) * 0.5;
                double vsz_y = (wz_vz(cfg, vz, i, j, k) + wz_vz(cfg, vz, i, jp, k)) * 0.5;
                pi_nx[eb] = vsy_z * vsz_y; // 对流
                pi_nx[eb] -= etaYZ[eb] * ((wz_vy(cfg, vy, i, j, k + 1) - wz_vy(cfg, vy, i, j, k))
                                        + (wz_vz(cfg, vz, i, jp, k) - wz_vz(cfg, vz, i, j, k))); // 粘性
                pi_nx[eb] -= sqrt(etaYZ[eb] * gam) * cfg.W * randN[1 * esize + eb]; // 随机

                // Πzx (pi_ny): 位于 ZX 面的棱边
                double vsz_x = (wz_vz(cfg, vz, i, j, k) + wz_vz(cfg, vz, ip, j, k)) * 0.5;
                double vsx_z = (wz_vx(cfg, vx, i, j, k) + wz_vx(cfg, vx, i, j, k + 1)) * 0.5;
                pi_ny[eb] = vsz_x * vsx_z; // 对流
                pi_ny[eb] -= etaZX[eb] * ((wz_vz(cfg, vz, ip, j, k) - wz_vz(cfg, vz, i, j, k))
                                        + (wz_vx(cfg, vx, i, j, k + 1) - wz_vx(cfg, vx, i, j, k))); // 粘性
                pi_ny[eb] -= sqrt(etaZX[eb] * gam) * cfg.W * randN[2 * esize + eb]; // 随机
            }
        }
    }

    // --- 3a. 显式一步：v* = v + dt·(-∇·Π + f) ---
    // 每个面心只算【一次】(∇·Π)_α，直接写 tmp_fα = v*（tmp_f* 复用为 v* 的存储，
    // 零额外显存）。原实现把每个面算了两次（d_pix_ijk 与左邻的 d_pix_im 是同一个量），
    // 这里消掉重复，Π 的读取减半 —— 正对 README 记录的缓存瓶颈。
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, fx, fy, fz, pi_dx, pi_dy, pi_dz, pi_nx, pi_ny, pi_nz, tmp_fx, tmp_fy, tmp_fz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                int ip = (i + 1) % Nx; int im = (i - 1 + Nx) % Nx;
                int jp = (j + 1) % Ny; int jm = (j - 1 + Ny) % Ny;
                int kp = (k + 1) % Nz; int km = (k - 1 + Nz) % Nz;

                // 棱边量的 z 向邻居走逻辑索引：壁面模式下 k=0 的下面那层是【下壁】
                // （逻辑 k=-1），不是周期回绕过来的上壁。搞错会静默用错壁面剪切。
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

                tmp_fx[ijk] = vx[ijk] + DT * (-d_pix + fx[ijk]);   // v*x
                tmp_fy[ijk] = vy[ijk] + DT * (-d_piy + fy[ijk]);   // v*y
                tmp_fz[ijk] = vz[ijk] + DT * (-d_piz + fz[ijk]);   // v*z

                // 壁面法向面【永不修正】：v*z 在两壁钉死为 0。这正是
                // doc/PressurePoisson.md §8a 相容性条件 Σ_k b̂_k = 0 的前提。
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
    // 按 cfg.wall_z 分派：周期走 3D FFT + 精确离散本征值；壁面走 xy 2D 批量 FFT
    // + z 向 Thomas。归一化因子（size / Nx·Ny）只在 solve_pressure 内部出现一次。
    // 数学推导与公式↔代码对照见 doc/PressurePoisson.md。
    solve_pressure(cfg, fft_data, tri_w, plan, plan_xy, p, diag);

    // --- 5. 修正步 (Correction / projection) ---
    // v^{n+1} = v* - dt·∇p，其中 v* 就是 tmp_f*。
    // 必须与上面的求解分成两个 kernel：p 的梯度要读相邻线程写的 p。
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, p, tmp_fx, tmp_fy, tmp_fz)
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

                vx[ijk] = tmp_fx[ijk] - DT * (p[IDX(ip, j, k)] - p[ijk]);
                vy[ijk] = tmp_fy[ijk] - DT * (p[IDX(i, jp, k)] - p[ijk]);
                // 壁面法向面永不修正：vz 在那两处恒 0。k=Nz-1 时 kp 会回绕到 0，
                // 那在壁面模式下是错的，所以这一层干脆不写（vz[Nz-1] 保持初值 0）。
                if (!wz_pinned_zface(cfg, k))
                {
                    vz[ijk] = tmp_fz[ijk] - DT * (p[IDX(i, j, kp)] - p[ijk]);
                }
            }
        }
    }

}