// fpd_check —— 自检与验证路径的独立可执行入口
//
//   fpd_check --check
//   fpd_check --lambda [L]
//   fpd_check --noise [L] [dt] [kT] [steps]
//   fpd_check --equipart <ghost|frozen|moving> [L] [dt] [kT] [steps] [seed]
//
// 从 fpd 拆出：main.cpp 里的四个自检路径 + --lambda 搬到这里，生产主循环
// 留在 fpd。拆分后 fpd 的 main 只做配置解析 + run_production。
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "./include/Tests.h"
#include "./include/Analysis.h"

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

    // 常量表：--lambda [L]   （纯 CPU，秒级）
    if (argc > 1 && std::strcmp(argv[1], "--lambda") == 0)
    {
        const int L_ = (argc > 2) ? std::atoi(argv[2]) : 32;
        NS_Config c2 = make_ns_config(L_, L_, L_, dt, kT, noise_on);
        FpdConstants fc = compute_fpd_constants(pp, c2,
                                                L_/2.0 + 0.3, L_/2.0 + 0.7, L_/2.0 + 0.5, 4);
        print_fpd_constants(fc, pp);

        std::cout << "\n  能量均分目标（L=" << L_ << "）\n";
        const double L3 = (double)L_ * L_ * L_;
        std::printf("    2kT/M_i              = %.6e * kT\n", 2.0 / fc.mass_mean);
        std::printf("    k=0 冻结修正 -2kT/L^3 = %.6e * kT   (%.3f%%)\n",
                    -2.0 / L3, -fc.mass_mean / L3 * 100.0);
        std::printf("    目标 <|V|^2>         = %.6e * kT\n",
                    equipartition_target(fc, 1.0, c2));
        return 0;
    }

    // 测试 3B：--equipart <ghost|frozen|moving> [L] [dt] [kT] [steps] [seed]
    if (argc > 2 && std::strcmp(argv[1], "--equipart") == 0)
    {
        EquipartMode m;
        if      (std::strcmp(argv[2], "ghost")  == 0) { m = EQ_GHOST;  }
        else if (std::strcmp(argv[2], "frozen") == 0) { m = EQ_FROZEN; }
        else if (std::strcmp(argv[2], "moving") == 0) { m = EQ_MOVING; }
        else { std::cerr << "mode 必须是 ghost | frozen | moving\n"; return 2; }

        const int    L_  = (argc > 3) ? std::atoi(argv[3]) : 32;
        const double dt_ = (argc > 4) ? std::atof(argv[4]) : 0.002;
        const double kT_ = (argc > 5) ? std::atof(argv[5]) : 0.25;
        const long   ns_ = (argc > 6) ? std::atol(argv[6]) : 300000;
        const unsigned long long sd_ = (argc > 7) ? std::strtoull(argv[7], 0, 10) : 1ULL;
        return run_equipartition_check(m, L_, dt_, kT_, ns_, sd_);
    }

    // 测试 3A：--noise [L] [dt] [kT] [steps]
    if (argc > 1 && std::strcmp(argv[1], "--noise") == 0)
    {
        const int    L_    = (argc > 2) ? std::atoi(argv[2]) : 32;
        const double dt_   = (argc > 3) ? std::atof(argv[3]) : 0.01;
        const double kT_   = (argc > 4) ? std::atof(argv[4]) : 1.0;
        const long   nst_  = (argc > 5) ? std::atol(argv[5]) : 200000;
        return run_noise_check(L_, dt_, kT_, nst_);
    }

    std::cerr << "用法: fpd_check --check\n"
              << "      fpd_check --check-potential\n"
              << "      fpd_check --check-tridiag\n"
              << "      fpd_check --check-poisson [Nx Ny Nz]\n"
              << "      fpd_check --lambda [L]\n"
              << "      fpd_check --noise [L] [dt] [kT] [steps]\n"
              << "      fpd_check --equipart <ghost|frozen|moving> [L] [dt] [kT] [steps] [seed]\n";
    return 2;
}
