// fpd_check —— 自检与验证路径的独立可执行入口
//
//   fpd_check --check
//   fpd_check --check-potential
//   fpd_check --check-tridiag
//   fpd_check --check-poisson [Nx Ny Nz]
//   fpd_check --check-wall [Nx Ny Nz]
//   fpd_check --lambda [L]
//
// 从 fpd 拆出：main.cpp 里的自检路径搬到这里，生产主循环留在 fpd。
// 拆分后 fpd 的 main 只做配置解析 + run_production。
//
// ⚠️ 原 --noise（测试 3A 纯流体噪声谱）与 --equipart（测试 3B 粒子能量均分）
//    是【周期专属】判据，随 z 周期路径一起删除。见 PROGRESS.md 的覆盖损失记录。
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "Tests.h"
#include "Analysis.h"

int main(int argc, char** argv)
{
    // 自检与验证路径用内置参数，不读配置文件
    const int    Nx = 128, Ny = 64, Nz = 32;
    const double dt        = 0.001;
    const double kT        = 0.05;
    const bool   noise_on  = true;
    const double radius    = 3.2;
    const double xi        = 1.0;
    const double ratio_eta = 50.0;

    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, kT, noise_on);
    PhiParams pp  = make_phi_params(radius, xi, ratio_eta);

    if (argc > 1 && std::strcmp(argv[1], "--check") == 0)
    {
        int fails = run_force_conservation_check(cfg, pp);
        std::cout << "\n";
        fails += run_overlap_check(cfg, pp);
        std::cout << "\n";
        fails += run_force_pipeline_check(cfg, pp);
        return fails;
    }

    if (argc > 1 && std::strcmp(argv[1], "--check-config") == 0)
    {
        return run_config_check();
    }

    // 势函数自检：纯 CPU，无卡可跑
    if (argc > 1 && std::strcmp(argv[1], "--check-potential") == 0)
    {
        return run_potential_check();
    }

    // 压力泊松三对角自检：纯 CPU，无卡可跑
    if (argc > 1 && std::strcmp(argv[1], "--check-tridiag") == 0)
    {
        return run_check_tridiag();
    }

    // 压力泊松算子往返判据：需 GPU（走 cuFFT）
    if (argc > 1 && std::strcmp(argv[1], "--check-poisson") == 0)
    {
        const int Px = (argc > 2) ? std::atoi(argv[2]) : 0;
        const int Py = (argc > 3) ? std::atoi(argv[3]) : 0;
        const int Pz = (argc > 4) ? std::atoi(argv[4]) : 0;
        return run_check_poisson(Px, Py, Pz);
    }

    // Phase 7-B 壁面判据：需 GPU（走 cuFFT），但不开噪声
    if (argc > 1 && std::strcmp(argv[1], "--check-wall") == 0)
    {
        const int Px = (argc > 2) ? std::atoi(argv[2]) : 0;
        const int Py = (argc > 3) ? std::atoi(argv[3]) : 0;
        const int Pz = (argc > 4) ? std::atoi(argv[4]) : 0;
        return run_check_wall(Px, Py, Pz);
    }

    // 常量表：--lambda [L]   （纯 CPU，秒级）
    if (argc > 1 && std::strcmp(argv[1], "--lambda") == 0)
    {
        const int L_ = (argc > 2) ? std::atoi(argv[2]) : 32;
        // ⚠️ L 必须大到测试粒子离两壁 >= range（相场模板盒的 z 向截断），
        //    否则量到的是「被壁面截断的支撑域」而不是体相值。默认 32 足够。
        NS_Config c2 = make_ns_config(L_, L_, L_, dt, kT, noise_on);
        FpdConstants fc = compute_fpd_constants(pp, c2,
                                                L_/2.0 + 0.3, L_/2.0 + 0.7, L_/2.0 + 0.5, 4);
        print_fpd_constants(fc, pp);
        return 0;
    }

    std::cerr << "用法: fpd_check --check\n"
              << "      fpd_check --check-config\n"
              << "      fpd_check --check-potential\n"
              << "      fpd_check --check-tridiag\n"
              << "      fpd_check --check-poisson [Nx Ny Nz]\n"
              << "      fpd_check --check-wall [Nx Ny Nz]\n"
              << "      fpd_check --lambda [L]\n";
    return 2;
}
