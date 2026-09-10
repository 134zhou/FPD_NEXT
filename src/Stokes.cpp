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

    // --- 1. 生成随机噪声 (Σ 项) ---
    // randD: 0,1,2 对应 xx, yy, zz 方向；randN: 0,1,2 对应 xy, yz, zx 方向
    // 每次生成前显式定位到本步专属的 slot，不依赖「上次生成后序列前进了多少」
    // —— 实测那个前进量并非 n，累计记账会失配（详见 Stokes.h 的说明）。
    const unsigned long long slot = 3ULL * (unsigned long long)size;
    const unsigned long long base = 2ULL * (unsigned long long)step * slot;

    #pragma acc host_data use_device(randD, randN)
    {
        CURAND_CHECK(curandSetGeneratorOffset(gen, base));
        CURAND_CHECK(curandGenerateNormalDouble(gen, randD, size * 3, 0.0, 1.0));
        CURAND_CHECK(curandSetGeneratorOffset(gen, base + slot));
        CURAND_CHECK(curandGenerateNormalDouble(gen, randN, size * 3, 0.0, 1.0));
    }
    
    // --- 2. 计算总动量通量 Π = Advection - Viscosity - Noise ---
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, eta, etaXY, etaYZ, etaZX, randD, randN, pi_dx, pi_dy, pi_dz, pi_nx, pi_ny, pi_nz)
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

                // --- 切应力分量 (Off-diagonal: Πxy, Πyz, Πzx) ---
                
                // Πxy (pi_nz): 位于 XY 面的棱边
                double vsx_y = (vx[ijk] + vx[IDX(i, jp, k)]) * 0.5;
                double vsy_x = (vy[ijk] + vy[IDX(ip, j, k)]) * 0.5;
                pi_nz[ijk] = vsx_y * vsy_x; // 对流
                pi_nz[ijk] -= etaXY[ijk] * ((vx[IDX(i, jp, k)] - vx[ijk]) + (vy[IDX(ip, j, k)] - vy[ijk])); // 粘性
                pi_nz[ijk] -= sqrt(etaXY[ijk]) * cfg.W * randN[0 * size + ijk]; // 随机

                // Πyz (pi_nx): 位于 YZ 面的棱边
                double vsy_z = (vy[ijk] + vy[IDX(i, j, kp)]) * 0.5;
                double vsz_y = (vz[ijk] + vz[IDX(i, jp, k)]) * 0.5;
                pi_nx[ijk] = vsy_z * vsz_y; // 对流
                pi_nx[ijk] -= etaYZ[ijk] * ((vy[IDX(i, j, kp)] - vy[ijk]) + (vz[IDX(i, jp, k)] - vz[ijk])); // 粘性
                pi_nx[ijk] -= sqrt(etaYZ[ijk]) * cfg.W * randN[1 * size + ijk]; // 随机

                // Πzx (pi_ny): 位于 ZX 面的棱边
                double vsz_x = (vz[ijk] + vz[IDX(ip, j, k)]) * 0.5;
                double vsx_z = (vx[ijk] + vx[IDX(i, j, kp)]) * 0.5;
                pi_ny[ijk] = vsz_x * vsx_z; // 对流
                pi_ny[ijk] -= etaZX[ijk] * ((vz[IDX(ip, j, k)] - vz[ijk]) + (vx[IDX(i, j, kp)] - vx[ijk])); // 粘性
                pi_ny[ijk] -= sqrt(etaZX[ijk]) * cfg.W * randN[2 * size + ijk]; // 随机
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

                // (∇·Π)_x，物理位置 (i+½, j, k)：d(Πxx)/dx + d(Πxy)/dy + d(Πxz)/dz
                double d_pix = (pi_dx[IDX(ip, j, k)] - pi_dx[ijk])
                             + (pi_nz[ijk] - pi_nz[IDX(i, jm, k)])
                             + (pi_ny[ijk] - pi_ny[IDX(i, j, km)]);

                // (∇·Π)_y，物理位置 (i, j+½, k)
                double d_piy = (pi_nz[ijk] - pi_nz[IDX(im, j, k)])
                             + (pi_dy[IDX(i, jp, k)] - pi_dy[ijk])
                             + (pi_nx[ijk] - pi_nx[IDX(i, j, km)]);

                // (∇·Π)_z，物理位置 (i, j, k+½)
                double d_piz = (pi_ny[ijk] - pi_ny[IDX(im, j, k)])
                             + (pi_nx[ijk] - pi_nx[IDX(i, jm, k)])
                             + (pi_dz[IDX(i, j, kp)] - pi_dz[ijk]);

                tmp_fx[ijk] = vx[ijk] + DT * (-d_pix + fx[ijk]);   // v*x
                tmp_fy[ijk] = vy[ijk] + DT * (-d_piy + fy[ijk]);   // v*y
                tmp_fz[ijk] = vz[ijk] + DT * (-d_piz + fz[ijk]);   // v*z
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
                vz[ijk] = tmp_fz[ijk] - DT * (p[IDX(i, j, kp)] - p[ijk]);
            }
        }
    }

}