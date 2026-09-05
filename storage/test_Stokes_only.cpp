#include <iostream>
#include <vector>
#include <openacc.h>

#include "./include/Stokes.h"
#include "./include/IO.h"
int main()
{
    NS_Config cfg = {64, 64, 64, 0.01, 100.0, 0.1};
    int size = cfg.Nx * cfg.Ny * cfg.Nz;

    // 申请容器
    std::vector<double> vx(size, 0.0), vy(size, 0.0), vz(size, 0.0), p(size, 0.0);
    std::vector<double> fx(size, 0.0), fy(size, 0.0), fz(size, 0.0);
    std::vector<double> eta(size, 1.0), etaXY(size, 1.0), etaYZ(size, 1.0), etaZX(size, 1.0);
    std::vector<double> pi_dx(size), pi_dy(size), pi_dz(size);
    std::vector<double> pi_nx(size), pi_ny(size), pi_nz(size);
    std::vector<double> fft_data(size * 2, 0.0);
    std::vector<double> randD(size * 3), randN(size * 3);

    // 1. 在 main 中为临时变量分配真实的 Host 内存，或者使用 vector
    std::vector<double> h_tmp_fx(size), h_tmp_fy(size), h_tmp_fz(size);

    // Y 的上半部分向右 (vx=1), 下半部分向左 (vx=-1)
    for (int i = 0; i < cfg.Nx; ++i)
    {
        for (int j = 0; j < cfg.Ny; ++j)
        {
            for (int k = 0; k < cfg.Nz; ++k)
            {
                int ijk = ((i * cfg.Ny + j) * cfg.Nz + k);
                if(16<i && i<48 && 16<j && j<48 && 16<k && k<48)
                {
                    vx[ijk] = 1;
                }
            }
        }
    }

    // 获取原始指针
    double* d_vx = vx.data(); double* d_vy = vy.data(); double* d_vz = vz.data();
    double* d_p = p.data();
    double* d_fx = fx.data(); double* d_fy = fy.data(); double* d_fz = fz.data();
    double* d_eta = eta.data(); double* d_etaXY = etaXY.data(); double* d_etaYZ = etaYZ.data(); double* d_etaZX = etaZX.data();
    double* d_pi_dx = pi_dx.data(); double* d_pi_dy = pi_dy.data(); double* d_pi_dz = pi_dz.data();
    double* d_pi_nx = pi_nx.data(); double* d_pi_ny = pi_ny.data(); double* d_pi_nz = pi_nz.data();
    double* d_fft = fft_data.data();
    double* d_randD = randD.data(); double* d_randN = randN.data();

    double* tmp_fx = h_tmp_fx.data();
    double* tmp_fy = h_tmp_fy.data();
    double* tmp_fz = h_tmp_fz.data();

    // 初始化库
    cufftHandle plan;
    cufftPlan3d(&plan, cfg.Nx, cfg.Ny, cfg.Nz, CUFFT_Z2Z);
    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);

    // 使用原始指针进行 OpenACC 数据映射
    #pragma acc enter data copyin(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size], d_eta[0:size], d_etaXY[0:size], d_etaYZ[0:size], d_etaZX[0:size], d_fx[0:size], d_fy[0:size], d_fz[0:size]) \
                         create(d_pi_dx[0:size], d_pi_dy[0:size], d_pi_dz[0:size], d_pi_nx[0:size], d_pi_ny[0:size], d_pi_nz[0:size], d_fft[0:size*2], d_randD[0:size*3], d_randN[0:size*3])

    #pragma acc enter data create(tmp_fx[0:size], tmp_fy[0:size], tmp_fz[0:size])

    std::cout << "Starting Simulation..." << std::endl;

    for(int step = 0; step < 1000; step++)
    {
        step_navier_stokes
                        (
                            cfg, d_vx, d_vy, d_vz, d_p, d_fx, d_fy, d_fz, 
                            d_eta, d_etaXY, d_etaYZ, d_etaZX, 
                            d_pi_dx, d_pi_dy, d_pi_dz, d_pi_nx, d_pi_ny, d_pi_nz, 
                            d_fft, plan, gen, d_randD, d_randN,
                            tmp_fx, tmp_fy, tmp_fz
                        );
        std::cout << "Step " << step << " completed." << std::endl;
        // 每 20 步将数据拷回 CPU 并保存一次
        if (step % 100 == 0)
        {
            #pragma acc update host(d_vx[0:size], d_vy[0:size], d_vz[0:size], d_p[0:size])
            save_vtk(step, cfg, vx, vy, vz, p);
        }
    }

    // 退出并拷贝结果
    #pragma acc exit data copyout(d_vx[0:size]) delete(d_vy[0:size], d_vz[0:size], d_p[0:size], d_eta[0:size], d_etaXY[0:size], d_etaYZ[0:size], d_etaZX[0:size], d_fx[0:size], d_fy[0:size], d_fz[0:size], d_pi_dx[0:size], d_pi_dy[0:size], d_pi_dz[0:size], d_pi_nx[0:size], d_pi_ny[0:size], d_pi_nz[0:size], d_fft[0:size*2], d_randD[0:size*3], d_randN[0:size*3])
    
    cufftDestroy(plan);
    curandDestroyGenerator(gen);
    
    std::cout << "Done. Final vx[0]: " << vx[0] << std::endl;
    return 0;
}