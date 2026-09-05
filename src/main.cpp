#include <iostream>
#include <vector>
#include <openacc.h>

#include "./include/Force.h"
#include "./include/Viscosity.h"
#include "./include/Stokes.h"
#include "./include/IO.h"
#include "./include/Velocity.h"

int main()
{
    NS_Config cfg = {128, 64, 32, 0.001, 1000.0, 10};
    int size = cfg.Nx * cfg.Ny * cfg.Nz;
    int N = 1;
    
    std::vector<double> vx(size, 0.0), vy(size, 0.0), vz(size, 0.0), p(size, 0.0);// 流体速度场
    std::vector<double> fx(size, 0.0), fy(size, 0.0), fz(size, 0.0);// 流体受力场

    std::vector<double> Rx(N, 0.0), Ry(N, 0.0), Rz(N, 0.0);// 粒子位置
    std::vector<double> Vx(N, 0.0), Vy(N, 0.0), Vz(N, 0.0);// 粒子速度
    std::vector<double> Fx(N, 0.0), Fy(N, 0.0), Fz(N, 0.0);// 粒子受力

    std::vector<double> phi(size), sum_phi(N);
    std::vector<double> eta(size), etaXY(size), etaYZ(size), etaZX(size);
    std::vector<double> fft_data(size * 2);


    std::vector<double> randD(size * 3), randN(size * 3);

    std::vector<double> h_tmp_fx(size), h_tmp_fy(size), h_tmp_fz(size);// 中间量
    std::vector<double> pi_dx(size), pi_dy(size), pi_dz(size);// 矩阵对角元
    std::vector<double> pi_nx(size), pi_ny(size), pi_nz(size);// 矩阵非对角元

    Rx[0] = 64.0; Ry[0] = 32.0; Rz[0] = 16.0;
    // Fx[0] = 0.0;

    // 获取原始指针
    double* d_vx = vx.data(); double* d_vy = vy.data(); double* d_vz = vz.data();
    double* d_p = p.data();
    double* d_fx = fx.data(); double* d_fy = fy.data(); double* d_fz = fz.data();
    

    double* d_Rx = Rx.data(); double* d_Ry = Ry.data(); double* d_Rz = Rz.data();
    double* d_Vx = Vx.data(); double* d_Vy = Vy.data(); double* d_Vz = Vz.data();
    double* d_Fx = Fx.data(); double* d_Fy = Fy.data(); double* d_Fz = Fz.data();

    // 流体phi和粘度场
    double* d_phi = phi.data(); 
    double* d_sum_phi = sum_phi.data();
    double* d_eta = eta.data(); 
    double* d_etaXY = etaXY.data(); 
    double* d_etaYZ = etaYZ.data(); 
    double* d_etaZX = etaZX.data();
    double* d_fft = fft_data.data();

    double* tmp_fx = h_tmp_fx.data();
    double* tmp_fy = h_tmp_fy.data();
    double* tmp_fz = h_tmp_fz.data();

    double* d_pi_dx = pi_dx.data(); double* d_pi_dy = pi_dy.data(); double* d_pi_dz = pi_dz.data();
    double* d_pi_nx = pi_nx.data(); double* d_pi_ny = pi_ny.data(); double* d_pi_nz = pi_nz.data();

    double* d_randD = randD.data(); double* d_randN = randN.data();
    // 初始化库
    cufftHandle plan;

    cufftPlan3d(&plan, cfg.Nz, cfg.Ny, cfg.Nx, CUFFT_Z2Z);
    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);

    // 使用原始指针进行 OpenACC 数据映射
    #pragma acc enter data copyin(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size], d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc enter data copyin(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Vx[0:N], d_Vy[0:N], d_Vz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])

    #pragma acc enter data create(d_phi[0:size], d_sum_phi[0:N])
    #pragma acc enter data create(d_eta[0:size], d_etaXY[0:size], d_etaYZ[0:size], d_etaZX[0:size])
    #pragma acc enter data create(d_fft[0:size*2])

    #pragma acc enter data create(tmp_fx[0:size], tmp_fy[0:size], tmp_fz[0:size])
    #pragma acc enter data create(d_pi_dx[0:size], d_pi_dy[0:size], d_pi_dz[0:size], d_pi_nx[0:size], d_pi_ny[0:size], d_pi_nz[0:size], d_randD[0:size*3], d_randN[0:size*3])

    std::cout << "Starting Simulation..." << std::endl;

    for(int step = 0; step < 200000; step++)
    {
        update_viscosity_fields
        (
            cfg,
            d_sum_phi,
            N, d_Rx, d_Ry, d_Rz,
            d_phi, d_eta, d_etaXY, d_etaYZ, d_etaZX
        );

        update_force_field
        (
            cfg,
            d_sum_phi,
            N, d_Rx, d_Ry, d_Rz,
            d_Fx, d_Fy, d_Fz,
            d_fx, d_fy, d_fz
        );

        step_navier_stokes
        (
            cfg, d_vx, d_vy, d_vz, d_p, d_fx, d_fy, d_fz, 
            d_eta, d_etaXY, d_etaYZ, d_etaZX, 
            d_pi_dx, d_pi_dy, d_pi_dz, d_pi_nx, d_pi_ny, d_pi_nz, 
            d_fft, plan, gen, d_randD, d_randN,
            tmp_fx, tmp_fy, tmp_fz
        );

        update_particale_VandR
        (
            cfg,
            d_phi, d_sum_phi,
            N, d_Rx, d_Ry, d_Rz,
            d_Vx, d_Vy, d_Vz,
            d_vx, d_vy, d_vz
        );
        std::cout << "\rStep " << step << " completed." << std::flush;
        // 每 20 步将数据拷回 CPU 并保存一次
        if (step % 10000 == 0)
        {
            #pragma acc update host(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size])
            //#pragma acc update host(d_phi[0:size])
            //#pragma acc update host(d_sum_phi[0:N])
            #pragma acc update host(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N])
            #pragma acc update host(d_Vx[0:N], d_Vy[0:N], d_Vz[0:N])
            save_fluid_vtk("result", step, cfg, vx, vy, vz, p);
            save_particles_vtk("result", step, N, Rx, Ry, Rz, Vx, Vy, Vz, Fx, Fy, Fz);
        }
    }

    // 退出并拷贝结果

    #pragma acc exit data delete(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size], d_fx[0:size], d_fy[0:size], d_fz[0:size])
    #pragma acc exit data delete(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Vx[0:N], d_Vy[0:N], d_Vz[0:N], d_Fx[0:N], d_Fy[0:N], d_Fz[0:N])

    #pragma acc exit data delete(d_phi[0:size], d_sum_phi[0:N])
    #pragma acc exit data delete(d_eta[0:size], d_etaXY[0:size], d_etaYZ[0:size], d_etaZX[0:size])
    #pragma acc exit data delete(d_fft[0:size*2])

    #pragma acc exit data delete(tmp_fx[0:size], tmp_fy[0:size], tmp_fz[0:size])
    #pragma acc exit data delete(d_pi_dx[0:size], d_pi_dy[0:size], d_pi_dz[0:size], d_pi_nx[0:size], d_pi_ny[0:size], d_pi_nz[0:size], d_randD[0:size*3], d_randN[0:size*3])

    
    cufftDestroy(plan);
    curandDestroyGenerator(gen);
    
    std::cout << "Done. Final vx[0]: " << vx[0] << std::endl;
    return 0;
}