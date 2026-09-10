// spike：实测决定「粒子间力」算在 CPU 还是 GPU。
//
// 背景：PROGRESS.md 决策记录写「粒子间力在 CPU 侧算」，论据是「传输 17KB/步、
// 延迟主导 20-40µs，NS 求解 1-3ms，占比 <3%」。复核发现该论据漏算了 O(N²) 计算
// 本身，且「CPU 才能逐位重启」不成立（Viscosity/Velocity 已用 GPU reduction 并
// 随 Phase 2 的逐位重启判据验证过）。
//
// 本 spike 用【同一份 routine seq 的 pair_force_over_r】分别驱动：
//   方案 A：CPU 串行 i<j 半矩阵
//   方案 B：GPU gang-over-i + vector-reduction-over-j 全矩阵（零 atomic）
// 测四件事：
//   (1) 两版结果一致（相对差 < 1e-14，浮点归约顺序不同，不要求逐位）
//   (2) GPU 版重复运行【逐位可复现】← 决定方案 B 可不可用
//   (3) 耗时：CPU 纯计算 / GPU 纯计算 / 传输往返，合成方案 A 与方案 B 的总开销
//   (4) N 标度：CPU 应见 O(N²) 上翘，GPU 应基本平坦
//
// 判决规则：GPU 版逐位可复现 且 不慢于 CPU 版 → 选 GPU（方案 B）；
//           GPU 版逐位不可复现 → 必须选 CPU（方案 A）。
//
// 编译：
//   nvc++ -acc -O3 -Isrc/include -Minfo=accel tests/spike/spike_pair_bench.cpp -o build/spike_pair
//   ./build/spike_pair

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <openacc.h>

// ---- 势函数（最小版，与将来 src/include/Potential.h 同源）----
// 这里只写 Morse 一支（含 exp()，是三种势里最贵的分支），bench 结论对它最保守。
// 力公式：g(r) = -(1/r)dU/dr = 2*alpha*De*E*(E-1)/r，E = exp(-alpha*(r-r_eq))。
// F_i = g(r) * (R_i - R_j)。返回标量 g，调用点乘位移向量。
#pragma acc routine seq
static inline double pair_force_over_r(double r2, double De, double alpha, double r_eq,
                                       double rcut2)
{
    if (r2 >= rcut2 || r2 <= 0.0) { return 0.0; }
    const double r = sqrt(r2);
    const double E = exp(-alpha * (r - r_eq));
    return 2.0 * alpha * De * E * (E - 1.0) / r;
}

// 最小镜像（与将来 Potential.cpp 同源；前提：位置已折叠到 [0,L)）
#pragma acc routine seq
static inline double min_image(double d, double L)
{
    if      (d >  0.5 * L) { return d - L; }
    else if (d < -0.5 * L) { return d + L; }
    return d;
}

// ---- 方案 A：CPU 串行 i<j 半矩阵 ----
static void forces_cpu(int N, double Lx, double Ly, double Lz,
                       double De, double alpha, double r_eq, double rcut2,
                       const double* Rx, const double* Ry, const double* Rz,
                       double* Fx, double* Fy, double* Fz)
{
    for (int n = 0; n < N; n++) { Fx[n] = 0; Fy[n] = 0; Fz[n] = 0; }
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
        {
            const double dx = min_image(Rx[i] - Rx[j], Lx);
            const double dy = min_image(Ry[i] - Ry[j], Ly);
            const double dz = min_image(Rz[i] - Rz[j], Lz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double g = pair_force_over_r(r2, De, alpha, r_eq, rcut2);
            Fx[i] += g * dx;  Fx[j] -= g * dx;
            Fy[i] += g * dy;  Fy[j] -= g * dy;
            Fz[i] += g * dz;  Fz[j] -= g * dz;
        }
}

// ---- 方案 B：GPU 全矩阵 + gang/vector reduction（零 atomic）----
static void forces_gpu(int N, double Lx, double Ly, double Lz,
                       double De, double alpha, double r_eq, double rcut2,
                       const double* Rx, const double* Ry, const double* Rz,
                       double* Fx, double* Fy, double* Fz)
{
    #pragma acc parallel loop gang present(Rx, Ry, Rz, Fx, Fy, Fz)
    for (int i = 0; i < N; i++)
    {
        double sx = 0.0, sy = 0.0, sz = 0.0;
        #pragma acc loop vector reduction(+:sx, sy, sz)
        for (int j = 0; j < N; j++)
        {
            if (j == i) { continue; }
            const double dx = min_image(Rx[i] - Rx[j], Lx);
            const double dy = min_image(Ry[i] - Ry[j], Ly);
            const double dz = min_image(Rz[i] - Rz[j], Lz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double g = pair_force_over_r(r2, De, alpha, r_eq, rcut2);
            sx += g * dx;  sy += g * dy;  sz += g * dz;
        }
        Fx[i] = sx;  Fy[i] = sy;  Fz[i] = sz;
    }
}

static double now_s()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

int main()
{
    const double Lx = 128.0, Ly = 64.0, Lz = 32.0;
    const double De = 50.0, alpha = 1.0, r_eq = 7.4;
    const double rcut = 15.0, rcut2 = rcut * rcut;   // 最小镜像上限 min(N)/2 = 16

    printf("=== spike_pair_bench：粒子间力 CPU vs GPU ===\n");
    printf("盒子 128x64x32   Morse(De=%.0f, alpha=%.1f, r_eq=%.1f, rcut=%.0f)\n\n",
           De, alpha, r_eq, rcut);

    const int Ns[] = {363, 1000, 4000};

    for (int it = 0; it < 3; it++)
    {
        const int N = Ns[it];

        // 固定 seed 的随机位置（可复现），折叠在 [0,L) 内
        std::mt19937 rng(12345u);
        std::uniform_real_distribution<double> ux(0, Lx), uy(0, Ly), uz(0, Lz);
        std::vector<double> Rx(N), Ry(N), Rz(N), Fx(N), Fy(N), Fz(N), Gx(N), Gy(N), Gz(N);
        for (int n = 0; n < N; n++)
        {
            Rx[n] = ux(rng); Ry[n] = uy(rng); Rz[n] = uz(rng);
        }

        // ---- 方案 A：CPU ----
        const double t_cpu0 = now_s();
        forces_cpu(N, Lx, Ly, Lz, De, alpha, r_eq, rcut2,
                   Rx.data(), Ry.data(), Rz.data(), Fx.data(), Fy.data(), Fz.data());
        const double t_cpu = now_s() - t_cpu0;

        // ---- 方案 B：GPU ----
        double *d_Rx = Rx.data(), *d_Ry = Ry.data(), *d_Rz = Rz.data();
        double *d_Gx = Gx.data(), *d_Gy = Gy.data(), *d_Gz = Gz.data();
        #pragma acc enter data copyin(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N]) create(d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])

        // 预热一次
        forces_gpu(N, Lx, Ly, Lz, De, alpha, r_eq, rcut2, d_Rx, d_Ry, d_Rz, d_Gx, d_Gy, d_Gz);

        // 计时（纯 kernel，不含传输）
        const int REPS = 20;
        const double t_gpu0 = now_s();
        for (int r = 0; r < REPS; r++)
        {
            forces_gpu(N, Lx, Ly, Lz, De, alpha, r_eq, rcut2, d_Rx, d_Ry, d_Rz, d_Gx, d_Gy, d_Gz);
        }
        const double t_gpu = (now_s() - t_gpu0) / REPS;

        #pragma acc update host(d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])

        // (1) 两版结果一致
        double relmax = 0.0, fmax = 0.0;
        for (int n = 0; n < N; n++)
        {
            fmax = fmax > fabs(Fx[n]) ? fmax : fabs(Fx[n]);
            fmax = fmax > fabs(Fy[n]) ? fmax : fabs(Fy[n]);
            fmax = fmax > fabs(Fz[n]) ? fmax : fabs(Fz[n]);
        }
        for (int n = 0; n < N; n++)
        {
            double d1 = fabs(Fx[n] - Gx[n]), d2 = fabs(Fy[n] - Gy[n]), d3 = fabs(Fz[n] - Gz[n]);
            double d  = d1 > d2 ? d1 : d2; d = d > d3 ? d : d3;
            relmax = relmax > d ? relmax : d;
        }
        relmax = fmax > 0 ? relmax / fmax : 0.0;

        // (2) GPU 版逐位可复现：再跑 10 次，每次拉回与第一次逐位比较
        bool bitwise = true;
        for (int r = 0; r < 10 && bitwise; r++)
        {
            forces_gpu(N, Lx, Ly, Lz, De, alpha, r_eq, rcut2, d_Rx, d_Ry, d_Rz, d_Gx, d_Gy, d_Gz);
            #pragma acc update host(d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])
            // 与 Fx 逐位比较 —— 注意 Fx 是 CPU 结果，不是 GPU 第一次结果！
            // 这里改为：把第一次 GPU 结果存下来再比较。为简洁，直接与第一次 Gx 的副本比。
        }
        // 上面的逐位比较逻辑有误，重做：存第一份 GPU 结果，跑 9 次逐位比
        {
            std::vector<double> g0x(N), g0y(N), g0z(N);
            for (int n = 0; n < N; n++) { g0x[n] = Gx[n]; g0y[n] = Gy[n]; g0z[n] = Gz[n]; }
            bitwise = true;
            for (int r = 0; r < 9 && bitwise; r++)
            {
                forces_gpu(N, Lx, Ly, Lz, De, alpha, r_eq, rcut2, d_Rx, d_Ry, d_Rz, d_Gx, d_Gy, d_Gz);
                #pragma acc update host(d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])
                for (int n = 0; n < N && bitwise; n++)
                {
                    if (Gx[n] != g0x[n] || Gy[n] != g0y[n] || Gz[n] != g0z[n]) { bitwise = false; }
                }
            }
        }

        #pragma acc exit data delete(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])

        // (3) 传输往返计时：方案 A 需要 update host(R) + update device(F)
        {
            #pragma acc enter data copyin(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N]) create(d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])
            const int TREPS = 50;
            const double t_x0 = now_s();
            for (int r = 0; r < TREPS; r++)
            {
                #pragma acc update host(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N])
                #pragma acc update device(d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])
            }
            const double t_xfer = (now_s() - t_x0) / TREPS;
            #pragma acc exit data delete(d_Rx[0:N], d_Ry[0:N], d_Rz[0:N], d_Gx[0:N], d_Gy[0:N], d_Gz[0:N])

            const double tA = t_cpu + t_xfer;   // 方案 A 每步总开销
            const double tB = t_gpu;            // 方案 B 每步总开销（无往返）

            printf("N=%4d   对 %9d\n", N, N * (N - 1) / 2);
            printf("   CPU 计算      %8.1f us\n", t_cpu * 1e6);
            printf("   GPU 计算      %8.1f us\n", t_gpu * 1e6);
            printf("   传输往返      %8.1f us\n", t_xfer * 1e6);
            printf("   ---- 合成 ----\n");
            printf("   方案 A (CPU+传输)  %8.1f us   占 6.5ms/步的 %.2f%%\n", tA * 1e6, tA / 6.5e-3 * 100.0);
            printf("   方案 B (GPU)       %8.1f us   占 6.5ms/步的 %.2f%%\n", tB * 1e6, tB / 6.5e-3 * 100.0);
            printf("   两版结果 相对差 = %.2e   %s\n", relmax, relmax < 1e-14 ? "PASS" : "FAIL");
            printf("   GPU 重复 10 次逐位 = %s\n\n", bitwise ? "可复现" : "*** 不可复现 ***");
        }
    }

    printf("判决：GPU 版逐位可复现且不慢于 CPU 版 → 选 GPU（方案 B）；\n");
    printf("      否则必须选 CPU（方案 A），并把结论写进 PROGRESS.md。\n");
    return 0;
}
