#include "./include/Stokes.h"

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
    curandGenerator_t gen,
    double* randD, double* randN,
    double* tmp_fx, double* tmp_fy, double* tmp_fz
)
{
    int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    int size = Nx * Ny * Nz;
    double DT = cfg.dt;

    // --- 1. 生成随机噪声 (Σ 项) ---
    // randD: 0,1,2 对应 xx, yy, zz 方向；randN: 0,1,2 对应 xy, yz, zx 方向
    #pragma acc host_data use_device(randD, randN)
    {
        CURAND_CHECK(curandGenerateNormalDouble(gen, randD, size * 3, 0.0, 1.0));
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

    #pragma acc parallel loop collapse(3) present(vx, vy, vz, fx, fy, fz, pi_dx, pi_dy, pi_dz, pi_nx, pi_ny, pi_nz, fft_data, tmp_fx, tmp_fy, tmp_fz)
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

                // --- A. 速度散度项 (∇·v / dt) ---
                // 物理位置：从面中心差分到体中心。vx[ijk]是右面，vx[im]是左面。
                double div_v = (vx[ijk] - vx[IDX(im, j, k)]) + 
                            (vy[ijk] - vy[IDX(i, jm, k)]) + 
                            (vz[ijk] - vz[IDX(i, j, km)]);

                // --- B. 动量通量张量的散度 (∇·Π) ---
                // 我们需要计算当前格子的三个分量力，它们分别位于 x, y, z 面中心。
                
                // 1. x-面上的总通量散度 (f_star_x)
                // 包含：d(Pixx)/dx + d(Pixy)/dy + d(Pixz)/dz
                double d_pix_ijk = (pi_dx[IDX(ip, j, k)] - pi_dx[ijk]) +      // 正应力梯度
                                (pi_nz[ijk] - pi_nz[IDX(i, jm, k)]) +     // 切应力梯度 (xy)
                                (pi_ny[ijk] - pi_ny[IDX(i, j, km)]);      // 切应力梯度 (xz)

                // 2. x-左面上的总通量散度 (f_star_x_im)
                // 物理位置：(i-0.5, j, k)
                double d_pix_im = (pi_dx[ijk] - pi_dx[IDX(im, j, k)]) + 
                                (pi_nz[IDX(im, j, k)] - pi_nz[IDX(im, jm, k)]) + 
                                (pi_ny[IDX(im, j, k)] - pi_ny[IDX(im, j, km)]);

                // 3. y-面上的总通量散度 (f_star_y)
                double d_piy_ijk = (pi_nz[ijk] - pi_nz[IDX(im, j, k)]) + 
                                (pi_dy[IDX(i, jp, k)] - pi_dy[ijk]) + 
                                (pi_nx[ijk] - pi_nx[IDX(i, j, km)]);

                // 4. y-后面上的总通量散度 (f_star_y_jm)
                double d_piy_jm = (pi_nz[IDX(i, jm, k)] - pi_nz[IDX(im, jm, k)]) + 
                                (pi_dy[ijk] - pi_dy[IDX(i, jm, k)]) + 
                                (pi_nx[IDX(i, jm, k)] - pi_nx[IDX(i, jm, km)]);

                // 5. z-面上的总通量散度 (f_star_z)
                double d_piz_ijk = (pi_ny[ijk] - pi_ny[IDX(im, j, k)]) + 
                                (pi_nx[ijk] - pi_nx[IDX(i, jm, k)]) + 
                                (pi_dz[IDX(i, j, kp)] - pi_dz[ijk]);

                // 6. z-底面上的总通量散度 (f_star_z_km)
                double d_piz_km = (pi_ny[IDX(i, j, km)] - pi_ny[IDX(im, j, km)]) + 
                                (pi_nx[IDX(i, j, km)] - pi_nx[IDX(i, jm, km)]) + 
                                (pi_dz[ijk] - pi_dz[IDX(i, j, km)]);

                // --- C. 外力的散度 (∇·f) ---
                double div_f = (fx[ijk] - fx[IDX(im, j, k)]) + 
                            (fy[ijk] - fy[IDX(i, jm, k)]) + 
                            (fz[ijk] - fz[IDX(i, j, km)]);

                // --- D. 最终源项合成 ---
                // 根据泊松方程：∇²p = (1/dt)∇·v - ∇·(∇·Π) + ∇·f
                // 注意：((d_pix_ijk - d_pix_im) + ...) 这一步完成了第二次散度运算
                fft_data[ijk * 2] = cfg.inv_dt * div_v - 
                                    ((d_pix_ijk - d_pix_im) + (d_piy_ijk - d_piy_jm) + (d_piz_ijk - d_piz_km)) + 
                                    div_f;
                fft_data[ijk * 2 + 1] = 0.0;


                // 计算 x 方向的面心合力 (∇·Π)_x 并存储
                tmp_fx[ijk] = d_pix_ijk;

                // 计算 y 方向的面心合力 (∇·Π)_y 并存储
                tmp_fy[ijk] = d_piy_ijk;

                // 计算 z 方向的面心合力 (∇·Π)_z 并存储
                tmp_fz[ijk] = d_piz_ijk;

                // 构造泊松源项时使用这些计算好的值计算“散度的散度”
            }
        }
    }

    // --- 4. FFT 求解泊松方程 (Pressure Solver) ---
    #pragma acc host_data use_device(fft_data)
    {
        CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex*)fft_data,
                                 (cufftDoubleComplex*)fft_data, CUFFT_FORWARD));
    }

    // 除以 7 点差分格式的【精确离散】拉普拉斯本征值 2(cos kx + cos ky + cos kz - 3)，
    // 而非连续谱的 -k²。这样与第 3 步的有限差分离散严格自洽。
    #pragma acc parallel loop collapse(3) present(fft_data)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                if (i == 0 && j == 0 && k == 0) { fft_data[ijk*2] = 0; fft_data[ijk*2+1] = 0; continue; }
                double nrm = 0.5 / (cos(2.*M_PI*i/Nx) + cos(2.*M_PI*j/Ny) + cos(2.*M_PI*k/Nz) - 3.0);
                fft_data[ijk * 2] *= nrm;
                fft_data[ijk * 2 + 1] *= nrm;
            }
        }
    }



    #pragma acc host_data use_device(fft_data)
    {
        CUFFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex*)fft_data,
                                 (cufftDoubleComplex*)fft_data, CUFFT_INVERSE));
    }

    // --- 5. 最终速度更新 (Correction step) ---
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, p, fft_data, tmp_fx, tmp_fy, tmp_fz, fx, fy, fz)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);

                // 1. 归一化并存储压力
                double p_ijk = fft_data[ijk * 2] / (double)size;
                p[ijk] = p_ijk;

                // 此时必须等待所有 p[ijk] 写入完成，或使用索引访问 fft_data 计算压力梯度
                // 在 GPU 上，由于 p[ip] 属于相邻线程，需确保同步或直接计算
            }
        }
    }

    // 分离出更新步骤以确保压力场完全算出 (或者使用 fft_data 归一化后的值)
    #pragma acc parallel loop collapse(3) present(vx, vy, vz, p, tmp_fx, tmp_fy, tmp_fz, fx, fy, fz)
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

                // 2. 利用缓存的散度项 tmp_f 直接更新速度
                // v(n+1) = v(n) + dt * [ -∇p - (∇·Π) + f ]
                vx[ijk] += DT * ( -(p[IDX(ip, j, k)] - p[ijk]) - tmp_fx[ijk] + fx[ijk] );
                vy[ijk] += DT * ( -(p[IDX(i, jp, k)] - p[ijk]) - tmp_fy[ijk] + fy[ijk] );
                vz[ijk] += DT * ( -(p[IDX(i, j, kp)] - p[ijk]) - tmp_fz[ijk] + fz[ijk] );
            }
        }
    }

}