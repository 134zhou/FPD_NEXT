// spike：验证 nvc++ 下 #pragma acc routine seq 作用于模板函数是否可行
// 编译：nvc++ -acc -O2 -Minfo=accel spike_template_routine.cpp -o spike
#include <cstdio>
#include <cmath>

#define IDX(i, j, k) ((i) + (j)*Nx + (k)*Nx*Ny)

enum Loc { CELL, FACE_X, FACE_Y, FACE_Z, EDGE_XY, EDGE_YZ, EDGE_ZX };

// 每个位置的半格偏移，唯一真值源
template<Loc L> struct LocOffset;
template<> struct LocOffset<CELL>    { static constexpr double dx=0.0, dy=0.0, dz=0.0; };
template<> struct LocOffset<FACE_X>  { static constexpr double dx=0.5, dy=0.0, dz=0.0; };
template<> struct LocOffset<FACE_Y>  { static constexpr double dx=0.0, dy=0.5, dz=0.0; };
template<> struct LocOffset<FACE_Z>  { static constexpr double dx=0.0, dy=0.0, dz=0.5; };
template<> struct LocOffset<EDGE_XY> { static constexpr double dx=0.5, dy=0.5, dz=0.0; };
template<> struct LocOffset<EDGE_YZ> { static constexpr double dx=0.0, dy=0.5, dz=0.5; };
template<> struct LocOffset<EDGE_ZX> { static constexpr double dx=0.5, dy=0.0, dz=0.5; };

struct PhiParams
{
    double radius, inv_xi;
    int    range, range_m1, n_range;
    double range2;
};

struct NS_Config { int Nx, Ny, Nz; double dt, inv_dt, W; };

#pragma acc routine seq
inline double order(double dx, double dy, double dz, double radius, double inv_xi)
{
    return 0.5*(tanh((radius - sqrt(dx*dx + dy*dy + dz*dz))*inv_xi) + 1.);
}

// 【被测对象】模板函数 + routine seq
#pragma acc routine seq
template<Loc L>
inline bool stencil_point(PhiParams pp, NS_Config cfg,
                          double Rnx, double Rny, double Rnz,
                          int l, int& ijk, double& w)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const int nr = pp.n_range;
    int li = l / (nr*nr), lj = (l / nr) % nr, lk = l % nr;
    int in = (int)Rnx, jn = (int)Rny, kn = (int)Rnz;
    int ir = li + in - pp.range_m1;
    int jr = lj + jn - pp.range_m1;
    int kr = lk + kn - pp.range_m1;

    double dx = ir - Rnx + LocOffset<L>::dx;
    double dy = jr - Rny + LocOffset<L>::dy;
    double dz = kr - Rnz + LocOffset<L>::dz;

    if (dx*dx + dy*dy + dz*dz > pp.range2) { return false; }

    ijk = IDX((ir + Nx) % Nx, (jr + Ny) % Ny, (kr + Nz) % Nz);
    w   = order(dx, dy, dz, pp.radius, pp.inv_xi);
    return true;
}

int main()
{
    NS_Config cfg = {128, 64, 32, 0.001, 1000.0, 10.0};
    PhiParams pp;
    pp.radius = 3.2; pp.inv_xi = 1.0;
    pp.range = (int)(2.*(3.2 + 1.0));
    pp.range_m1 = pp.range - 1;
    pp.n_range  = 2*pp.range;
    pp.range2   = (double)(pp.range*pp.range);

    const int size = cfg.Nx*cfg.Ny*cfg.Nz;
    const int n3   = pp.n_range*pp.n_range*pp.n_range;
    const int N    = 1;

    double* Rx = new double[N]; double* Ry = new double[N]; double* Rz = new double[N];
    double* Fx = new double[N];
    double* sum_phix = new double[N];
    double* fx = new double[size];
    Rx[0] = 64.3; Ry[0] = 32.7; Rz[0] = 16.5;   // 故意用非整数亚格点位置
    Fx[0] = 1.0;

    #pragma acc enter data copyin(Rx[0:N], Ry[0:N], Rz[0:N], Fx[0:N]) create(sum_phix[0:N], fx[0:size])

    // --- 归约：gang over n + vector reduction ---
    #pragma acc parallel loop gang present(Rx, Ry, Rz, sum_phix)
    for (int n = 0; n < N; n++)
    {
        double sx = 0.0;
        #pragma acc loop vector reduction(+:sx)
        for (int l = 0; l < n3; l++)
        {
            int ijk; double w;
            if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sx += w; }
        }
        sum_phix[n] = sx;
    }

    #pragma acc parallel loop present(fx)
    for (int i = 0; i < size; i++) { fx[i] = 0.0; }

    // --- 投影：与上面同一个 Loc ---
    #pragma acc parallel loop collapse(2) present(Rx, Ry, Rz, Fx, sum_phix, fx)
    for (int n = 0; n < N; n++)
        for (int l = 0; l < n3; l++)
        {
            int ijk; double w;
            if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
            {
                #pragma acc atomic update
                fx[ijk] += Fx[n] * w / sum_phix[n];
            }
        }

    #pragma acc update host(fx[0:size], sum_phix[0:N])

    // 力守恒判据：Σ_grid fx 必须精确等于 Fx[0]
    double tot = 0.0;
    for (int i = 0; i < size; i++) { tot += fx[i]; }

    printf("sum_phix[0] = %.10f\n", sum_phix[0]);
    printf("Sigma fx    = %.16f\n", tot);
    printf("|Sigma fx - 1.0| = %.3e   -> %s\n",
           fabs(tot - 1.0), (fabs(tot - 1.0) < 1e-12) ? "PASS" : "FAIL");

    #pragma acc exit data delete(Rx[0:N], Ry[0:N], Rz[0:N], Fx[0:N], sum_phix[0:N], fx[0:size])
    return 0;
}
