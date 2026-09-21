#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <openacc.h>

#include "State.h"
#include "Stokes.h"
#include "Viscosity.h"
#include "Force.h"
#include "Velocity.h"
#include "Potential.h"
#include "Tests.h"
#include "Analysis.h"

// ============================================================================
// z 向无滑移壁面的端到端判据（Phase 7-B）。
//
//   W1  max|div v| 含两个壁面层        —— doc/PressurePoisson.md §13 盲区 1 的正解
//   W2  vz 在两壁被钉死                —— 法向边界条件
//   W3  平面 Poiseuille 的【精确离散】闭式 —— 切向 ghost 与散度端系数的判决性判据
//   W3b 上下壁 Π_yz/Π_zx 展开式            —— 两片壁面、两个切向分量的直接判据
//   W6  z 向动量收支恒等式（逐步精确）  —— fz 截断 / vz 钉死 / ∂p/∂z / Π_zz 四件事一次全查
//
// 这些判据都【不开噪声】。把「算子对不对」与「噪声幅度对不对」彻底分离 ——
// 噪声幅度（√2 因子）由 W5b（能量均分 + γ≡1 对照）负责。
//
// 独立实现约定：本文件的散度、Poiseuille 闭式、动量收支都是【手写】的，
// 不调用 Stokes.cpp 的内部量（只能从 FpdState 的公开数组取）。
// ============================================================================

static int g_fail = 0;
static void check(bool ok, const char* what)
{
    std::printf("  %-58s %s\n", what, ok ? "PASS" : "FAIL");
    if (!ok) { g_fail++; }
}

// 壁面感知的离散散度（独立手写，不调用生产代码）。
// 体心 (i,j,k)：z 向是 vz[k] - vz[k-1]，两壁的 ghost 都是 0。
static double max_abs_div(const FpdState& st)
{
    const int Nx = st.cfg.Nx, Ny = st.cfg.Ny, Nz = st.cfg.Nz;
    double mx = 0.0;
    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                const int ijk = IDX(i, j, k);
                const int im = (i - 1 + Nx) % Nx;
                const int jm = (j - 1 + Ny) % Ny;
                const double vzb = (k > 0) ? st.vz[IDX(i, j, k - 1)] : 0.0;  // k=0 下面就是下壁
                const double d = (st.vx[ijk] - st.vx[IDX(im, j, k)])
                               + (st.vy[ijk] - st.vy[IDX(i, jm, k)])
                               + (st.vz[ijk] - vzb);
                mx = std::max(mx, std::fabs(d));
            }
        }
    }
    return mx;
}

static double max_abs_v(const FpdState& st)
{
    const int size = (int)st.size;
    double mx = 0.0;
    for (int t = 0; t < size; t++)
    {
        mx = std::max(mx, std::fabs(st.vx[t]));
        mx = std::max(mx, std::fabs(st.vy[t]));
        mx = std::max(mx, std::fabs(st.vz[t]));
    }
    return mx;
}

// 跑 n 步（无粒子、无噪声）。bgx 是均匀体力。
static void run_steps(FpdState& st, PhiParams pp, int N, double bgx,
                      double bgy, double bgz, int n)
{
    for (int s = 0; s < n; s++)
    {
        update_viscosity_fields(st.cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(st.cfg, pp, N, st.Rx, st.Ry, st.Rz,
                           st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz,
                           bgx, bgy, bgz, st.fx, st.fy, st.fz);
        step_navier_stokes(st.cfg, st.vx, st.vy, st.vz, st.p,
                           st.fx, st.fy, st.fz,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz,
                           st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan_xy, st.tri_w, st.diag,
                           st.gen, st.randD, st.randN,
                           st.tmp_fx, st.tmp_fy, st.tmp_fz, s);
    }
}

// ============================================================================
// W1 + W2：随机速度场跑若干步后，∇·v 必须在【所有】层上到机器精度。
// 这是 Phase 7-A 盲区 1 的正解 —— 往返判据对「算子与修正步不匹配」是盲的。
// ============================================================================
static void w1_w2_divergence(int Nx, int Ny, int Nz)
{
    const int size = Nx * Ny * Nz;
    const double dt = 0.002, kT = 0.0;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, kT, /*noise_on=*/false);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 4242ULL);

    // 随机初速度，并把上壁 vz 钉死（壁面构型的合法初态）
    unsigned s = 12345u;
    for (int t = 0; t < size; t++)
    {
        s = s * 1664525u + 1013904223u; st.vx[t] = ((double)(s & 0xFFFF) / 32768.0) - 1.0;
        s = s * 1664525u + 1013904223u; st.vy[t] = ((double)(s & 0xFFFF) / 32768.0) - 1.0;
        s = s * 1664525u + 1013904223u; st.vz[t] = ((double)(s & 0xFFFF) / 32768.0) - 1.0;
    }
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++) { st.vz[IDX(i, j, Nz - 1)] = 0.0; }
    }
    #pragma acc update device(st.vx[0:size], st.vy[0:size], st.vz[0:size])

    const int n_steps = 20;
    run_steps(st, make_phi_params(3.2, 1.0, 50.0), 0, 0.0, 0.0, 0.0, n_steps);
    st.download(ST_VELOCITY);

    const double divm = max_abs_div(st);
    const double vm   = max_abs_v(st);

    std::printf("      %dx%dx%d  %d 步后 max|div v| = %.4e   max|v| = %.4e   比值 = %.3e\n",
                Nx, Ny, Nz, n_steps, divm, vm, divm / vm);
    check(divm / vm < 1e-12, "W1 壁面：max|div v|/max|v| < 1e-12（含 k=0 与 k=Nz-1 两层）");

    // W2：上壁那一层必须【逐位】为 0
    double vztop = 0.0;
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            vztop = std::max(vztop, std::fabs(st.vz[IDX(i, j, Nz - 1)]));
        }
    }
    std::printf("      上壁层的 max|vz| = %.4e\n", vztop);
    check(vztop == 0.0, "W2 壁面：vz[:,:,Nz-1] 逐位为 0");

    st.finish();
}

// ============================================================================
// W3：平面 Poiseuille 的精确离散闭式。
//
// 均匀体力 fx = g、η ≡ 1（无粒子）、无噪声下的稳态。稳态时 vx 与 x,y 无关 ⇒
// 对流项恒 0、div v = 0 ⇒ p ≡ 0，所以这纯粹是粘性-外力平衡。
//
// 离散方程（见 doc/PressurePoisson.md §14）：
//   内部 k      vx[k+1] - 2vx[k] + vx[k-1] = -g/eta
//   壁面 k=0    vx[1]   - 3vx[0]           = -g/eta     （ghost vx[-1] = -vx[0]）
// 解得闭式
//   vx[k] = (g/2eta) * [ (Nz^2+1)/4 - (k - (Nz-1)/2)^2 ]
//
// ⚠️ 判据必须用【这个离散闭式】，不能用连续抛物线 (g/2)(Nz^2/4 - z^2)：
//    后者带 O(1/Nz^2) 的假偏差（Nz=32 时 0.1%），足以让人误判 ghost 写错了。
//
// 做法：把初值【直接设成闭式】，跑若干步，验证它是一个不动点。
// 这比「跑到稳态」快几个数量级，且检验的正是离散算子的平衡关系。
// ============================================================================
static void w3_poiseuille(int Nz, double g, bool verbose, int& fails_out,
                          double& err_out)
{
    const int Nx = 8, Ny = 8;
    const double dt = 0.01;             // eta=1 时稳定上界 dt < 1/6，留足余量
    const double eta = 1.0;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, 0.0, /*noise_on=*/false);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 777ULL);

    const int size = Nx * Ny * Nz;
    for (int t = 0; t < size; t++) { st.vx[t] = 0.0; st.vy[t] = 0.0; st.vz[t] = 0.0; }

    const double c = 0.5 * (Nz - 1);
    for (int k = 0; k < Nz; k++)
    {
        const double u = (g / (2.0 * eta)) * ((double)(Nz * Nz + 1) / 4.0 - (k - c) * (k - c));
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++) { st.vx[IDX(i, j, k)] = u; }
        }
    }
    #pragma acc update device(st.vx[0:size], st.vy[0:size], st.vz[0:size])

    // 闭式是精确不动点 ⇒ 每步的变化应为机器精度量级。
    // 20 步是刻意的：若闭式【不是】不动点，每步的 O(dt) 偏差会累积出来。
    run_steps(st, make_phi_params(3.2, 1.0, 50.0), 0, g, 0.0, 0.0, 20);
    st.download(ST_VELOCITY);

    const double c2 = 0.5 * (Nz - 1);
    double mx = 0.0, umax = 0.0;
    for (int k = 0; k < Nz; k++)
    {
        const double exact = (g / (2.0 * eta)) * ((double)(Nz * Nz + 1) / 4.0 - (k - c2) * (k - c2));
        umax = std::max(umax, std::fabs(exact));
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                mx = std::max(mx, std::fabs(st.vx[IDX(i, j, k)] - exact));
            }
        }
    }
    const double rel = umax > 0.0 ? mx / umax : mx;
    err_out = rel;

    if (verbose)
    {
        std::printf("      Nz=%2d  与离散闭式的 max 相对差 = %.3e\n", Nz, rel);
    }
    else if (rel >= 1e-12)
    {
        fails_out++;
    }
    st.finish();
}

// ============================================================================
// W3b：直接核对两片壁面上的剪切通量展开式。
//
// kT=0、vz=0 时壁面对流和横向 vz 导数都为 0，只剩
//   bottom: Π_yz=-2η_yz vy[0]，Π_zx=-2η_zx vx[0]
//   top:    Π_yz=+2η_yz vy[Nz-1]，Π_zx=+2η_zx vx[Nz-1]
// 这里直接读取生产 kernel 写出的 pi_nx/pi_ny；期望值独立地在 CPU 上计算。
// ============================================================================
static void w3b_wall_shear_flux()
{
    const int Nx = 4, Ny = 3, Nz = 4;
    const int size = Nx * Ny * Nz;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.001, 0.0, /*noise_on=*/false);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 314159ULL);

    for (int t = 0; t < size; t++)
    {
        st.vx[t] = 0.15 + 0.003 * t;
        st.vy[t] = -0.20 + 0.002 * t;
        st.vz[t] = 0.0;
        st.fx[t] = st.fy[t] = st.fz[t] = 0.0;
        st.eta[t] = st.etaXY[t] = 1.0;
    }
    for (size_t t = 0; t < st.esize; t++)
    {
        st.etaYZ[t] = 1.0 + 0.001 * (double)t;
        st.etaZX[t] = 1.3 + 0.001 * (double)t;
    }
    st.upload(ST_VELOCITY | ST_PHI);

    step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p,
                       st.fx, st.fy, st.fz,
                       st.eta, st.etaXY, st.etaYZ, st.etaZX,
                       st.pi_dx, st.pi_dy, st.pi_dz,
                       st.pi_nx, st.pi_ny, st.pi_nz,
                       st.fft, st.plan_xy, st.tri_w, st.diag,
                       st.gen, st.randD, st.randN,
                       st.tmp_fx, st.tmp_fy, st.tmp_fz, 0L);

    acc_update_self(st.pi_nx, st.nbEdge);
    acc_update_self(st.pi_ny, st.nbEdge);

    double worst = 0.0;
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            const int bottom = IDX(i, j, 0);
            const int top = IDX(i, j, Nz - 1);
            const int eb_bottom = wz_edge_idx(cfg, i, j, -1);
            const int eb_top = wz_edge_idx(cfg, i, j, Nz - 1);

            worst = std::max(worst, std::fabs(st.pi_nx[eb_bottom]
                                           - (-2.0 * st.etaYZ[eb_bottom] * st.vy[bottom])));
            worst = std::max(worst, std::fabs(st.pi_nx[eb_top]
                                           - ( 2.0 * st.etaYZ[eb_top] * st.vy[top])));
            worst = std::max(worst, std::fabs(st.pi_ny[eb_bottom]
                                           - (-2.0 * st.etaZX[eb_bottom] * st.vx[bottom])));
            worst = std::max(worst, std::fabs(st.pi_ny[eb_top]
                                           - ( 2.0 * st.etaZX[eb_top] * st.vx[top])));
        }
    }

    std::printf("      上下壁 Π_yz/Π_zx 展开式 max|误差| = %.3e\n", worst);
    check(worst < 1e-13, "W3b 两片壁面的 Π_yz/Π_zx 展开式与独立 CPU 公式一致");
    st.finish();
}

// ============================================================================
// W6：z 向动量收支恒等式，【每一步】精确成立。
//
//   Σ_{ij, k=0..Nz-2} (vz^{n+1} - vz^n)/dt = F_z - (P̄_{Nz-1} - P̄_0) - (Z̄_{Nz-1} - Z̄_0)
//
// P̄_k = Σ_ij p[i,j,k]，Z̄_k = Σ_ij Π_zz[i,j,k]，F_z = Σ_grid fz。
// x/y 剪切项在周期方向望远镜到 0；z 向项在自由面 k=0..Nz-2 上望远镜。
//
// 一次检验四件事：fz 在壁面被截断、vz 在两壁被钉死、压力梯度、Π_zz。
// 稳态时左边为 0，右边就是「壁面通过静压 + 法向应力支撑粒子重量」。
// ============================================================================
static void w6_momentum_budget(int Nx, int Ny, int Nz)
{
    const int size = Nx * Ny * Nz;
    const double dt = 0.002;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, 0.0, /*noise_on=*/false);

    FpdState st;
    st.init(cfg, 0, ST_FULL, 999ULL);

    unsigned s = 31415u;
    for (int t = 0; t < size; t++)
    {
        s = s * 1664525u + 1013904223u; st.vx[t] = ((double)(s & 0xFFFF) / 32768.0) - 1.0;
        s = s * 1664525u + 1013904223u; st.vy[t] = ((double)(s & 0xFFFF) / 32768.0) - 1.0;
        s = s * 1664525u + 1013904223u; st.vz[t] = ((double)(s & 0xFFFF) / 32768.0) - 1.0;
    }
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++) { st.vz[IDX(i, j, Nz - 1)] = 0.0; }
    }
    #pragma acc update device(st.vx[0:size], st.vy[0:size], st.vz[0:size])

    PhiParams pp = make_phi_params(3.2, 1.0, 50.0);
    const double gx = 0.3, gz = -0.7;   // 非零外力，让 F_z 那一项真的有内容
    double worst = 0.0, scale = 0.0;

    std::vector<double> vz_old(size), vz_new(size);
    for (int step = 0; step < 10; step++)
    {
        st.download(ST_VELOCITY);
        for (int t = 0; t < size; t++) { vz_old[t] = st.vz[t]; }

        update_viscosity_fields(cfg, pp, 0, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, 0, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz,
                           gx, 0.0, gz, st.fx, st.fy, st.fz);
        step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p,
                           st.fx, st.fy, st.fz,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz,
                           st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan_xy, st.tri_w, st.diag,
                           st.gen, st.randD, st.randN,
                           st.tmp_fx, st.tmp_fy, st.tmp_fz, step);

        st.download(ST_VELOCITY | ST_PHI);
        for (int t = 0; t < size; t++) { vz_new[t] = st.vz[t]; }

        // 左边：所有【自由】z 面（k = 0 .. Nz-2）的 vz 变化之和 / dt
        double lhs = 0.0;
        for (int k = 0; k < Nz - 1; k++)
        {
            for (int j = 0; j < Ny; j++)
            {
                for (int i = 0; i < Nx; i++)
                {
                    const int ijk = IDX(i, j, k);
                    lhs += (vz_new[ijk] - vz_old[ijk]) / dt;
                }
            }
        }

        // F_z：只累加【自由】面。k=Nz-1 是上壁，那里的 fz 不进入任何一项
        // （K2a 算完之后 tmp_fz 在该层被直接钉 0）。
        double fz_sum = 0.0;
        for (int k = 0; k < Nz - 1; k++)
        {
            for (int j = 0; j < Ny; j++)
            {
                for (int i = 0; i < Nx; i++) { fz_sum += st.fz[IDX(i, j, k)]; }
            }
        }

        // P̄_k 从主机侧 p 取（ST_VELOCITY 已下载）；Z̄_k 必须把 pi_dz 拉下来 ——
        // 它是每步重算的中间量，download() 按设计不搬运它。
        // ⚠️ 用【运行时 API】而不是 `#pragma acc parallel loop reduction`：
        //    后者在本文件里读到了主机侧的陈旧零值（present 子句与运行时映射
        //    在此 TU 的配合不可靠），诊断为「Z̄ 恒等于 0」—— 正是 CLAUDE.md 点名的
        //    「present 遗漏不报错，只给垃圾数据」那一类。
        acc_update_self(st.pi_dz, st.size * sizeof(double));
        double Pb[2] = {0.0, 0.0}, Zb[2] = {0.0, 0.0};
        const int kk[2] = { 0, Nz - 1 };
        for (int q = 0; q < 2; q++)
        {
            double ps = 0.0, zs = 0.0;
            const int base = kk[q] * Nx * Ny;
            for (int t = 0; t < Nx * Ny; t++)
            {
                ps += st.p[base + t];
                zs += st.pi_dz[base + t];
            }
            Pb[q] = ps; Zb[q] = zs;
        }
        const double rhs = fz_sum - (Pb[1] - Pb[0]) - (Zb[1] - Zb[0]);

        // ⚠️ 归一化必须用【各项的量级】，不能用 max(|lhs|,|rhs|)：
        //    本算例第一步之后就到达稳态，两项都掉到 ~1e-11，
        //    再用残差的相对量会放大成 O(1) 的假失败。
        const double err  = std::fabs(lhs - rhs);
        const double dP   = Pb[1] - Pb[0];
        const double dZ   = Zb[1] - Zb[0];
        const double scal = std::fabs(fz_sum) + std::fabs(dP) + std::fabs(dZ) + 1e-30;
        worst = std::max(worst, err / scal);
        scale = std::max(scale, scal);
        if (step == 0)
        {
            std::printf("      第 1 步：lhs = %.6e  F_z = %.6e  ΔP̄ = %.6e  ΔZ̄ = %.6e"
                        "  rhs = %.6e  |lhs-rhs| = %.3e\n",
                        lhs, fz_sum, Pb[1] - Pb[0], Zb[1] - Zb[0], rhs,
                        std::fabs(lhs - rhs));
        }
    }

    std::printf("      10 步中 max||lhs-rhs| / (|F_z|+|ΔP̄|+|ΔZ̄|)| = %.4e   （各项量级 %.3e）\n",
                worst, scale);
    check(worst < 1e-11, "W6 壁面：z 向动量收支恒等式逐步成立");
    st.finish();
}

// ============================================================================
// W4：壁面截断下的力守恒 + 均匀流场的粒子速度。
//
// 力守恒 Σ_grid f_α == F_α 在截断下【仍然精确成立】，理由是分子（力投影）与分母
// （sum_phi_α）是同一个 stencil_point<L> 在同一组点上的求和 —— 截断只是让那组点
// 变少。这正是 CLAUDE.md 第 1 条约束保证的性质。
//
// ⚠️ 但这条判据对 stencil_point 的【内部公式错误】仍然是盲的（投影和归一化会一起
//    错）。所以这里同时跑「均匀流场 v≡1 ⇒ V_i == 1」—— 那条对内部错误不盲。
// ============================================================================
static void w4_truncation()
{
    const int Nx = 32, Ny = 32, Nz = 32;
    const int size = Nx * Ny * Nz;
    const double TOL = 1e-12;

    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.002, 0.25, true);
    PhiParams pp  = make_phi_params(3.2, 1.0, 50.0);

    // 四种构型：远离壁 / 贴【下】壁 / 贴【上】壁 / N=2 重叠且靠近下壁。
    // 上壁与下壁走的是【不同】代码路径（下壁棱边是逻辑 k=-1 的额外一层，
    // 上壁是数组里的 k=Nz-1），所以必须两片都测 —— 只有一侧的判据对
    // 「上下不对称」是盲的。
    // ⚠️ 贴上壁那一例的 Rz 必须取【精确离散镜像】Nz-1-Rz：格胞在 z=0..Nz-1、
    //    棱边在 z=k+1/2，关于 z=(Nz-1)/2 的反射 z -> Nz-1-z 把格胞映到格胞、
    //    壁面棱边映到壁面棱边。取 Rz=27.0 而不是"看起来对称"的 28.0，
    //    这样两例的 z/x 必须逐个数字相等，构成真正的上下对称判据。
    const int    NC[4] = {1, 1, 1, 2};
    const char*  TG[4] = {"远离壁 Rz=16.5", "贴下壁 Rz=4.0（支撑域被截断）",
                          "贴上壁 Rz=27.0（= 31-4，下壁的离散镜像）",
                          "N=2 重叠（间距 6）+ 靠近下壁"};
    const double RZ[4][2] = {{16.5, 0.0}, {4.0, 0.0}, {27.0, 0.0}, {5.0, 11.0}};
    double zx_ratio[4] = {0.0, 0.0, 0.0, 0.0};

    for (int c = 0; c < 4; c++)
    {
        const int N = NC[c];
        FpdState st;
        st.init(cfg, N, ST_KINEMATIC, 0);

        for (int n = 0; n < N; n++)
        {
            st.Rx[n] = 16.3; st.Ry[n] = 16.7; st.Rz[n] = RZ[c][n];
            st.Fx[n] = 1.0;  st.Fy[n] = 2.0;  st.Fz[n] = 3.0;
        }
        st.upload(ST_PARTICLE);

        update_viscosity_fields(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                           st.fx, st.fy, st.fz);
        st.download(ST_PHI | ST_PARTICLE);

        double sx = 0.0, sy = 0.0, sz = 0.0;
        for (int t = 0; t < size; t++) { sx += st.fx[t]; sy += st.fy[t]; sz += st.fz[t]; }
        double FX = 0.0, FY = 0.0, FZ = 0.0;
        for (int n = 0; n < N; n++) { FX += st.Fx[n]; FY += st.Fy[n]; FZ += st.Fz[n]; }

        const double ex = std::fabs(sx - FX), ey = std::fabs(sy - FY), ez = std::fabs(sz - FZ);
        const bool ok = (ex < TOL && ey < TOL && ez < TOL);
        std::printf("      %-34s |err| = %.3e / %.3e / %.3e  %s\n",
                    TG[c], ex, ey, ez, ok ? "PASS" : "FAIL");
        if (!ok) { g_fail++; }

        // sum_phiz < sum_phix 证明截断【真的在发生】——否则本判据是空的。
        // 贴壁那一例的支撑域被下壁切掉一块，z 向归一化必然小于 x 向。
        zx_ratio[c] = st.sum_phiz[0] / st.sum_phix[0];
        std::printf("      %-34s sum_phi = %.4f / %.4f / %.4f   z/x = %.6f\n",
                    "", st.sum_phix[0], st.sum_phiy[0], st.sum_phiz[0], zx_ratio[c]);

        // 均匀流场：vx = vy = 1，vz = 0（壁面法向必须为 0）。
        // V_i = Σ v w / sum_phi 必须【精确】等于 1（x/y）与 0（z），与截断、与重叠无关。
        for (int t = 0; t < size; t++) { st.vx[t] = 1.0; st.vy[t] = 1.0; st.vz[t] = 0.0; }
        st.upload(ST_VELOCITY);
        update_particle_velocity(cfg, pp, N, st.Rx, st.Ry, st.Rz,
                                 st.sum_phix, st.sum_phiy, st.sum_phiz,
                                 st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);
        st.download(ST_PARTICLE);

        double vmax = 0.0;
        for (int n = 0; n < N; n++)
        {
            vmax = std::max(vmax, std::fabs(st.Vx[n] - 1.0));
            vmax = std::max(vmax, std::fabs(st.Vy[n] - 1.0));
            vmax = std::max(vmax, std::fabs(st.Vz[n] - 0.0));
        }
        std::printf("      %-34s 均匀流场 max|V_i - 期望| = %.3e  %s\n",
                    "", vmax, (vmax < TOL) ? "PASS" : "FAIL");
        if (!(vmax < TOL)) { g_fail++; }

        st.finish();
    }

    // 上下对称性：离散镜像位置上的 z 向归一化必须【逐个数字】一致。
    // 这是对「上下必须对称排除」那条要求的判决性检验 —— 只排一侧会让
    // 力守恒仍然通过（分子分母一起错），但这里会立刻显形。
    const double asym = std::fabs(zx_ratio[1] - zx_ratio[2]);
    std::printf("      上下对称：|z/x(下) - z/x(上)| = %.3e   (%.6f vs %.6f)  %s\n",
                asym, zx_ratio[1], zx_ratio[2], (asym < 1e-12) ? "PASS" : "FAIL");
    if (!(asym < 1e-12)) { g_fail++; }
}

// ---------------------------------------------------------------------------
// W5：壁面能量均分，带 γ≡1 对照 —— √2 因子的判决性判据。
//
//   <Σ_r |v|²> = kT · dim ，  dim = Nx·Ny·(2Nz−1) + 1
//
// 采样只累加【自由】分量（壁面模式下 vz 的上壁那一层不是自由度）。
// 误差棒必须来自分块平均（复用 Analysis.h 的 blocking_analysis / pick_plateau）：
// 把相关样本当独立样本会低估误差约 7.5 倍（PROGRESS.md 的实测结论）。
//
// 盒子取小（4×4×4）：最慢模式 τ = (2Nz/π)²/ν ≈ 1.6，几十秒就能跑几百个 τ。
// 大盒子要 41 万步（L=128），这儿只要几千步。
//
// 判别力：γ 取错的效应集中在壁面剪切模式上，其占的总自由度比例随 Nz 减小而增大
// （Nz=4 时壁面棱边对应的自由度约 2/(2Nz−1) = 29%）。所以小 Nz 不但便宜，
// 而且判别力最强 —— 与 3A/3B 为了统计量而选大盒子的直觉正好相反。
// ---------------------------------------------------------------------------
static double last_ratio = 0.0, last_sigma = 0.0;
static void w5_wall_equipart(int Nx, int Ny, int Nz, double dt, double kT,
                             bool gamma1, long n_steps, unsigned long long seed,
                             bool verbose = true)
{
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, kT, /*noise_on=*/true);
    cfg.noise_gamma1 = gamma1 ? 1 : 0;

    FpdState st;
    st.init(cfg, 0, ST_FULL, seed);
    PhiParams pp = make_phi_params(3.2, 1.0, 50.0);

    const int size = Nx * Ny * Nz;
    // dim = n_v − rank(D) = n_v − size + 1（见 W5' 的推导）。
    // 壁面 n_v = Nx·Ny·(3Nz−1)（上壁那一层 vz 不是自由度）。
    // ⚠️ 这个 dim 与 doc/PressurePoisson.md §14.4 的自由度计数推导必须一致。
    const double dim = (double)(Nx * Ny * (2 * Nz - 1) + 1);

    // 平衡时间：z 向 Neumann 的最长波长是 2H ⇒ k = π/H
    const double tau    = 1.0 / ((M_PI / (double)Nz) * (M_PI / (double)Nz));
    const long   n_equil = (long)(8.0 * tau / dt);
    // 采样间隔固定在【物理时间】上（0.2），否则 dt 越小样本越相关
    const long samp_every = std::max(1L, (long)(0.2 / dt));
    const unsigned long long n_slot = (unsigned long long)st.slotD + st.slotN;

    std::vector<double> series;
    for (long step = 0; step < n_steps; step++)
    {
        update_viscosity_fields(cfg, pp, 0, st.Rx, st.Ry, st.Rz,
                                st.sum_phix, st.sum_phiy, st.sum_phiz,
                                st.eta, st.etaXY, st.etaYZ, st.etaZX);
        update_force_field(cfg, pp, 0, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                           st.sum_phix, st.sum_phiy, st.sum_phiz,
                           0.0, 0.0, 0.0, st.fx, st.fy, st.fz);
        step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p,
                           st.fx, st.fy, st.fz,
                           st.eta, st.etaXY, st.etaYZ, st.etaZX,
                           st.pi_dx, st.pi_dy, st.pi_dz,
                           st.pi_nx, st.pi_ny, st.pi_nz,
                           st.fft, st.plan_xy, st.tri_w, st.diag,
                           st.gen, st.randD, st.randN,
                           st.tmp_fx, st.tmp_fy, st.tmp_fz, step);

        if (step >= n_equil && step % samp_every == 0)
        {
            acc_update_self(st.vx, (size_t)size * sizeof(double));
            acc_update_self(st.vy, (size_t)size * sizeof(double));
            acc_update_self(st.vz, (size_t)size * sizeof(double));
            double s = 0.0;
            for (int t = 0; t < size; t++)
            {
                s += st.vx[t] * st.vx[t] + st.vy[t] * st.vy[t];
                // 上壁那一层 vz 不是自由度（恒 0，永不更新）
                if (t / (Nx * Ny) != Nz - 1) { s += st.vz[t] * st.vz[t]; }
            }
            series.push_back(s);
        }
    }
    (void)n_slot;
    st.finish();

    const int n = (int)series.size();
    BlockStat bs[64];
    const int ns = blocking_analysis(series.data(), n, bs, 64);
    const BlockStat pk = pick_plateau(bs, ns, 20);

    const double target = kT * dim;
    const char* lbl = gamma1 ? "壁面 γ≡1" : "壁面 γ=2";
    if (pk.block_len < 0)
    {
        std::printf("      %-10s 未找到分块平台（误差棒不可信）→ 按设计拒绝背书 FAIL\n", lbl);
        g_fail++;
        return;
    }
    const double ratio = pk.mean / target;
    const double sigma = pk.stderr_mean / target;
    const double dev   = sigma > 0.0 ? std::fabs(ratio - 1.0) / sigma : 0.0;
    std::printf("      %-10s dt=%-7g n=%d  b=%-4d  <Σ|v|²>/(kT·dim) = %.4f ± %.4f  (%.1fσ)\n",
                lbl, dt, n, pk.block_len, ratio, sigma, dev);
    last_ratio = ratio; last_sigma = sigma;
    // ⚠️ 这里【不做】绝对判据。绝对量测（比值是否精确为 1）被这台小盒子量测机器
    //    的 ~1.7% 系统偏差污染（周期控制也落在 0.98），判据由调用方以【差分】形式
    //    给出。这个函数只负责产生数值。
    (void)verbose;
}

// ============================================================================
// W8：粒子-壁面排斥势
//
//   W8a  F_z(z) = -dU_wall/dz 的四阶数值微分对照。能量用本文件【独立手写】的
//        间隙表达式（不调用 wall_gap_*，那正是被检验的对象），且对 z 求导而不是
//        对 h —— 链式法则会自动产生符号结构，手写的符号错了这里立刻显形。
//   W8b  上下壁镜像：F_z(Nz-1-z) == -F_z(z)，且壁面势在 x/y 上恒为 0
//   W8c  体相逐位零力；且 wallpotential=none 时【所有】z 上逐位为 0（回归锚）
//   W8d  静态力平衡：独立二分反解 z_rest，断言重力+壁面力为零，且两侧符号
//        相反（吸引子）、且偏离时残差非零（判别力断言）
//   W8e  GPU 与 CPU 两份装配在壁面开启时一致
//
// ⚠️ W8d 是【力平衡】判据，不是动力学判据。完整的「粒子真的沉到 z_rest」由生产
//    算例端到端确认（见 PROGRESS 的 Phase 6），不放进单元判据 —— 流体驰豫时间
//    τ ≈ (L/π)²/ν 是几万步量级，放进来会让 --check-wall 从 2 分钟变成几十分钟，
//    而判据的判别力并不增加（力平衡是动力学的必要条件，且是这里唯一可精确判的量）。
// ============================================================================

// 本文件【独立手写】的壁面势能量，供 W8a 做数值微分对照。
// ⚠️ 刻意不调用 Potential.h 的 wall_gap_bottom/top：那正是被检验的对象。
//    这里用「中心到壁面的距离」写成 (z - z_wall) - a 的另一种展开式，
//    能抓到「间隙定义差一个常数」这类错误。
//    扫描区间保证 hb, ht > 0（h<=0 的物理无意义，且不在本判据的覆盖范围内）。
static double u_wall_ref(const WallParams& w, int Nz, double z)
{
    const double zb = -0.5, zt = (double)Nz - 0.5;
    const double hb = (z - zb) - w.radius;
    const double ht = (zt - z) - w.radius;
    return pair_energy(w.pot, hb * hb) + pair_energy(w.pot, ht * ht);
}

static void w8a_force_vs_energy()
{
    const int Nz = 32;
    NS_Config cfg = make_ns_config(32, 32, Nz, 0.002, 0.0, false);
    const double a = 3.2;

    struct Case { const char* name; int type; int shift;
                  double eps, sigma, De, alpha, r_eq, rcut; double z_lo, z_hi; };
    // z_lo/z_hi 是扫描区间：必须整段落在 hb < rcut（也就是壁面势真的非零）之内，
    // 否则整段都是恒 0，判据退化成空的。
    const Case CS[3] = {
        {"WCA shift=energy",  POT_WCA,   SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0, 3.0, 4.85},
        {"Morse shift=force", POT_MORSE, SHIFT_FORCE,  0, 0, 5.0, 2.0, 1.5, 2.0, 3.0, 4.85},
        {"LJ126 shift=force", POT_LJ126, SHIFT_FORCE,  1.0, 2.0, 0, 0, 0, 3.0, 3.0, 4.85},
    };

    for (int c = 0; c < 3; c++)
    {
        const Case& K = CS[c];
        WallParams w;
        w.pot    = make_potential_params(K.type, K.shift, K.eps, K.sigma,
                                         K.De, K.alpha, K.r_eq, K.rcut);
        w.radius = a;

        const int npt = 400;
        double max_err = 0.0, fmax = 0.0, fmin = 1e300;
        for (int k = 0; k <= npt; k++)
        {
            const double z = K.z_lo + (K.z_hi - K.z_lo) * (double)k / (double)npt;
            const double h = 1e-6 * z;

            // dU/dz 的四阶中心差分：[-U(z+2h) + 8U(z+h) - 8U(z-h) + U(z-2h)]/(12h)
            const double dUdz = ( -u_wall_ref(w, Nz, z + 2 * h)
                                + 8.0 * u_wall_ref(w, Nz, z + h)
                                - 8.0 * u_wall_ref(w, Nz, z - h)
                                + u_wall_ref(w, Nz, z - 2 * h) ) / (12.0 * h);
            // 保守力 F_z = -dU/dz。⚠️ 这个负号必须【从微分定义】来，不能照抄
            //    wall_force_z 的写法 —— 抄了就把判据变成了自洽的、盲的。
            const double num = -dUdz;
            const double ana = wall_force_z(w, cfg, z);

            // 混合容差：|F| 跨十几个量级，纯相对容差在小力端会被舍入淹没
            const double denom = std::fmax(std::fabs(ana), 1.0);
            const double err = std::fabs(num - ana) / denom;
            if (err > max_err) { max_err = err; }
            fmax = std::fmax(fmax, std::fabs(ana));
            if (std::fabs(ana) > 1e-12) { fmin = std::fmin(fmin, std::fabs(ana)); }
        }
        std::printf("      %-18s 扫描 z∈[%g, %g]  max 相对差 = %.3e   |F| 跨 %.1e → %.1e\n",
                    K.name, K.z_lo, K.z_hi, max_err, fmin, fmax);
        // 覆盖率断言：力必须跨至少 2 个量级，否则「扫描区全在势外」会让判据变空
        const bool covered = (fmax / std::fmax(fmin, 1e-300)) > 100.0;
        check(max_err < 1e-6 && covered,
              "W8a 壁面力 = -dU_wall/dz 四阶差分对照（力跨 ≥2 量级）");
    }

    // POT_NONE：壁面势关闭时 F_z 必须【逐位】为 0 —— 这是 S1 回归锚的单元形式
    {
        WallParams w = no_wall();
        w.radius = a;
        bool allzero = true;
        for (int k = 0; k <= 400; k++)
        {
            const double z = 0.0 + 32.0 * (double)k / 400.0;
            if (wall_force_z(w, cfg, z) != 0.0) { allzero = false; }
        }
        check(allzero, "W8a wallpotential=none 时 F_z 逐位为 0（所有 z）");
    }
}

static void w8b_mirror_symmetry()
{
    const int Nz = 32;
    NS_Config cfg = make_ns_config(32, 32, Nz, 0.002, 0.0, false);

    WallParams w;
    w.pot    = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0);
    w.radius = 3.2;

    // 离散镜像 z -> Nz-1-z 把下壁映到上壁。整数/半整数位置上 (z+1/2) 与
    // ((Nz-0.5)-z) 的计算【逐位相同】，所以这里可以断言逐位相等。
    const double ZS[4] = {4.0, 4.5, 5.0, 3.75};
    double worst_bit = 0.0, worst_sub = 0.0;
    for (int k = 0; k < 4; k++)
    {
        const double z  = ZS[k];
        const double zm = (double)(Nz - 1) - z;
        const double f1 = wall_force_z(w, cfg, z);
        const double f2 = wall_force_z(w, cfg, zm);
        const double d  = std::fabs(f1 + f2);
        if (k < 3) { worst_bit = std::fmax(worst_bit, d); }         // 逐位
        else       { worst_sub = std::fmax(worst_sub, d / std::fmax(std::fabs(f1), 1.0)); }
    }
    std::printf("      整数/半整数镜像：|F_z(z) + F_z(Nz-1-z)| = %.3e（逐位）\n", worst_bit);
    std::printf("      亚格点镜像 z=3.75：相对差 = %.3e（(z+0.5) 与 ((Nz-0.5)-z) 差 1 ulp）\n",
                worst_sub);
    check(worst_bit == 0.0, "W8b 上下壁镜像逐位反对称（整数/半整数位置）");
    check(worst_sub < 1e-14, "W8b 上下壁镜像反对称（亚格点，容 1e-14）");
}

static void w8c_bulk_zero()
{
    const int Nz = 32;
    NS_Config cfg = make_ns_config(32, 32, Nz, 0.002, 0.0, false);

    WallParams w;
    w.pot    = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0);
    w.radius = 3.2;
    const double rcut = w.pot.rcut;

    // 体相：min(hb, ht) > rcut ⇒ 两面壁都超出截断 ⇒ 逐位 0
    const double z_lo = -0.5 + 3.2 + rcut + 1e-9;
    const double z_hi = (double)Nz - 0.5 - 3.2 - rcut - 1e-9;
    bool allzero = true;
    for (int k = 0; k <= 200; k++)
    {
        const double z = z_lo + (z_hi - z_lo) * (double)k / 200.0;
        if (wall_force_z(w, cfg, z) != 0.0) { allzero = false; }
    }
    std::printf("      体相 z∈[%.4f, %.4f]（rcut=%.4f）\n", z_lo, z_hi, rcut);
    check(allzero, "W8c 体相内壁面力逐位为 0（不污染远离壁面的物理）");

    // 边界刚过截断的【外侧】也必须为 0（守卫写反了会在这一侧漏）
    const double z_out = -0.5 + 3.2 + rcut + 1e-6;
    check(wall_force_z(w, cfg, z_out) == 0.0,
          "W8c 恰好超出截断的外侧逐位为 0");
}

// W8d：静态力平衡。
// z_rest 用【本文件自己的】二分反解（Config.cpp 的 wall_rest_gap 是 static，
// 本文件根本调不到 —— 独立性是结构保证的，不是靠自觉）。
static double w8d_solve_rest(const WallParams& w, NS_Config cfg, double gz)
{
    auto F = [&](double z) { return gz + wall_force_z(w, cfg, z); };
    // 下界取「间隙 = 0.2*rcut」（力很大，F>0），上界取「间隙 = rcut」（F<0 由重力给）
    double lo = -0.5 + w.radius + 0.2 * w.pot.rcut;
    double hi = -0.5 + w.radius + 0.999 * w.pot.rcut;
    for (int it = 0; it < 200; it++)
    {
        const double mid = 0.5 * (lo + hi);
        if (F(mid) > 0.0) { lo = mid; } else { hi = mid; }
    }
    return 0.5 * (lo + hi);
}

static void w8d_static_balance()
{
    const int Nz = 32;
    NS_Config cfg = make_ns_config(32, 32, Nz, 0.002, 0.0, false);
    const double gz = -10.0;

    WallParams w;
    w.pot    = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0);
    w.radius = 3.2;

    const double z_rest = w8d_solve_rest(w, cfg, gz);
    const double f0 = gz + wall_force_z(w, cfg, z_rest);
    std::printf("      z_rest = %.6f（间隙 %.6f）   F_total(z_rest) = %.3e\n",
                z_rest, z_rest + 0.5 - w.radius, f0);
    check(std::fabs(f0) < 1e-10 * std::fabs(gz),
          "W8d 平衡点处 重力 + 壁面力 = 0（独立二分解出的 z_rest）");

    // 吸引子：偏下必须被推上去（F_total > 0），偏上必须被推下来（F_total < 0）
    const double d = 0.05;
    const double fdn = gz + wall_force_z(w, cfg, z_rest - d);
    const double fup = gz + wall_force_z(w, cfg, z_rest + d);
    std::printf("      z_rest-%.2f: F_total = %+.4f   z_rest+%.2f: F_total = %+.4f\n",
                d, fdn, d, fup);
    check(fdn > 0.0 && fup < 0.0,
          "W8d 平衡点是【吸引子】（下方被推回、上方被压下）");

    // 判别力：残差随偏离线性增长，且不能处处为 0（否则上面的判据是空的）
    const double d2 = 0.20;
    const double f2 = gz + wall_force_z(w, cfg, z_rest - d2);
    check(std::fabs(f2) > 100.0 * std::fabs(f0) && f2 > 0.0,
          "W8d 判据非空：偏离 0.2 的残差比平衡点残差大 >100 倍");

    // 刚度余量：这个势必须能在半个平衡间隙处提供远超重力的排斥
    //（否则热涨落/一步的力脉冲就足以把粒子压穿，check_wall_bounds 会中止）
    const double f_half = wall_force_z(w, cfg, -0.5 + w.radius + 0.5 * (z_rest + 0.5 - w.radius));
    std::printf("      半间隙处排斥力 = %.3e   重力 = %.1f   倍数 = %.1e\n",
                f_half, std::fabs(gz), f_half / std::fabs(gz));
    check(f_half > 100.0 * std::fabs(gz),
          "W8d 刚度余量：半间隙处排斥 > 100×重力（不会被压穿）");
}

// W8e：GPU 与 CPU 两份装配在壁面【开启】时一致，且壁面力把重力【精确】抵掉。
//
// 刻意把粒子间势关掉（pot = none）：粒子间力与壁面力混在一起时，前者可以大几个
// 量级，会把后者的装配错误淹没（相对差仍小到「看起来通过」）。本判据要的就是
// 单独照见壁面项。
//
// ⚠️ 诚实标注：GPU 与 CPU 共用 wall_force_z，所以本判据对「wall_force_z 内部公式
//    错误」是【盲的】—— 那由 W8a 的四阶微分对照负责。这里能抓的是装配侧的遗漏：
//    某一侧忘了加壁面项、加了两遍、或错误地放进了 j 归约。
static void w8e_gpu_cpu_wall()
{
    const int Nx = 32, Ny = 32, Nz = 32;
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.002, 0.0, false);

    const int N = 3;
    WallParams w;
    w.pot    = make_potential_params(POT_WCA, SHIFT_ENERGY, 1.0, 2.0, 0, 0, 0, 0);
    w.radius = 3.2;
    const PotentialParams pot = no_wall().pot;      // 关掉粒子间势，隔离壁面项
    const ExternalField   ext{0.0, 0.0, -10.0};

    const double z_rest = w8d_solve_rest(w, cfg, ext.gz);
    FpdState st;
    st.init(cfg, N, ST_STENCIL, 0);
    // 下壁平衡点 / 体相 / 上壁平衡点（离散镜像）
    const double RZ[3] = {z_rest, 15.5, (double)(Nz - 1) - z_rest};
    for (int n = 0; n < N; n++)
    {
        st.Rx[n] = 16.3 + 5.0 * n; st.Ry[n] = 16.7; st.Rz[n] = RZ[n];
    }
    st.upload(ST_PARTICLE);

    compute_particle_forces(cfg, pot, ext, w, N, st.Rx, st.Ry, st.Rz,
                            st.Fx, st.Fy, st.Fz);
    st.download(ST_PARTICLE);

    double cFx[3], cFy[3], cFz[3];
    compute_particle_forces_cpu(cfg, pot, ext, w, N,
                                st.Rx, st.Ry, st.Rz, cFx, cFy, cFz);

    double fmax = 0.0, dmax = 0.0;
    for (int n = 0; n < N; n++)
    {
        fmax = std::fmax(fmax, std::fmax(std::fabs(st.Fx[n]),
                                std::fmax(std::fabs(st.Fy[n]), std::fabs(st.Fz[n]))));
        dmax = std::fmax(dmax, std::fabs(st.Fx[n] - cFx[n]));
        dmax = std::fmax(dmax, std::fabs(st.Fy[n] - cFy[n]));
        dmax = std::fmax(dmax, std::fabs(st.Fz[n] - cFz[n]));
    }
    const double rel = fmax > 0.0 ? dmax / fmax : 0.0;

    std::printf("      下壁平衡点 z=%.4f  F_z = %+.3e\n", st.Rz[0], st.Fz[0]);
    std::printf("      体相       z=%.4f  F_z = %+.3e\n", st.Rz[1], st.Fz[1]);
    std::printf("      上壁平衡点 z=%.4f  F_z = %+.3e\n", st.Rz[2], st.Fz[2]);
    std::printf("      GPU vs CPU 装配 max 相对差 = %.3e\n", rel);

    check(rel < 1e-12, "W8e GPU 与 CPU 装配在壁面开启时一致");

    // 定量断言①：下壁平衡点上【总力恰好为 0】—— 这就是沉降的力平衡，
    // 跑在 GPU 装配路径上。壁面项被漏掉时这里会读到 -10，一目了然。
    check(std::fabs(st.Fz[0]) < 1e-10,
          "W8e 下壁平衡点上的总力 = 0（壁面力精确抵消重力）");

    // 定量断言②：镜像位置上的壁面力必须【逐位反对称】⇒ 总力 = 2·g_z。
    //
    // ⚠️ 上壁【不存在】平衡点，这不是缺陷：上壁的排斥方向朝下（远离壁面），
    //    与重力同向，两者叠加永远托不住粒子。只有下壁能承接沉降堆积。
    //    把这里写成「上壁也应该平衡」会得到一个永远 FAIL 的假判据。
    check(std::fabs(st.Fz[2] - 2.0 * ext.gz) < 1e-10,
          "W8e 上壁镜像处壁面力反对称 ⇒ 总力 = 2·g_z（上壁无平衡点，符合物理）");

    // 定量断言③：体相只受外场
    check(std::fabs(st.Fz[1] - ext.gz) < 1e-12 && st.Fz[1] == cFz[1],
          "W8e 体相粒子只受外场（壁面力逐位为 0）");

    st.finish();
}

// ============================================================================
// 入口。Nx<=0 时跑内置电池。
// ============================================================================
int run_check_wall(int Nx, int Ny, int Nz)
{
    g_fail = 0;
    std::printf("=== z 向无滑移壁面判据（需 GPU）===\n");

    if (Nx <= 0) { Nx = 16; Ny = 16; Nz = 16; }

    std::printf("\n[W1/W2] 壁面模式下的 ∇·v 与法向速度钉死\n");
    w1_w2_divergence(Nx, Ny, Nz);

    std::printf("\n[W3] 平面 Poiseuille：与【精确离散】闭式的对照\n");
    {
        const double g = 0.4;
        double e = 0.0;
        w3_poiseuille(Nz, g, true, g_fail, e);
        check(e < 1e-12, "W3 壁面：稳态与离散闭式 vx[k] = (g/2η)[(Nz²+1)/4 − (k−(Nz−1)/2)²] 一致");

        // 顺带报告与【连续】抛物线的偏差 —— 它应该是 O(1/Nz²)，不是判据
        const double c = 0.5 * (Nz - 1);
        double dmax = 0.0, umax = 0.0;
        for (int k = 0; k < Nz; k++)
        {
            const double zd = (double)k - c;
            const double cont  = (g / 2.0) * ((double)(Nz * Nz) / 4.0 - zd * zd);
            const double disc  = (g / 2.0) * ((double)(Nz * Nz + 1) / 4.0 - zd * zd);
            umax = std::max(umax, std::fabs(disc));
            dmax = std::max(dmax, std::fabs(disc - cont));
        }
        std::printf("      （参考）离散闭式与连续抛物线的相对差 = %.3e ≈ 1/(2Nz²) = %.3e\n",
                    dmax / umax, 1.0 / (2.0 * (double)Nz * Nz));
    }

    std::printf("\n[W3b] 两片壁面的剪切通量展开式\n");
    w3b_wall_shear_flux();

    std::printf("\n[W4] 壁面截断下的力守恒与粒子速度（含贴壁与 N=2 重叠）\n");
    w4_truncation();

    std::printf("\n[W6] z 向动量收支恒等式\n");
    w6_momentum_budget(Nx, Ny, Nz);

    std::printf("\n[W8] 粒子-壁面排斥势\n");
    w8a_force_vs_energy();
    w8b_mirror_symmetry();
    w8c_bulk_zero();
    w8d_static_balance();
    w8e_gpu_cpu_wall();

    std::printf("\n[W5b] 壁面能量均分 + γ≡1 对照（√2 的判决性判据）\n");
    double r_g2 = 0.0, s_g2 = 0.0, r_g1 = 0.0, s_g1 = 0.0;
    w5_wall_equipart(4, 4, 4, 0.002, 1.0, false, 300000, 1ULL);
    r_g2 = last_ratio; s_g2 = last_sigma;
    w5_wall_equipart(4, 4, 4, 0.002, 1.0, true,  300000, 1ULL);
    r_g1 = last_ratio; s_g1 = last_sigma;

    // ⚠️ 判据的形态：**差分**，不是绝对。
    //    4³ 小盒子里这台量测机器本身带 ~1.7% 的系统偏差，来源是 dim 的约定
    //    （冻结模式数）+ 显式 Euler 的 O(dt) 偏差，与壁面无关。所以能可靠交付的
    //    判据是「γ=2 是否被数据选中」，而不是「比值是否精确为 1」。
    //    绝对一致需要更大的盒子 + dt 外推，见 PROGRESS.md 的遗留问题。
    {
        const double diff = r_g2 - r_g1;
        const double sd   = std::sqrt(s_g2 * s_g2 + s_g1 * s_g1);
        const double nsig = sd > 0.0 ? std::fabs(diff) / sd : 0.0;
        std::printf("      γ=2 与 γ=1 之差 = %+.4f ± %.4f  (%.1fσ)；"
                    "γ=2 更接近 1：%s\n", diff, sd, nsig, (r_g2 > r_g1) ? "是" : "否");

        // 两者必须显著可分（否则本判据没有判别力）
        check(nsig > 5.0 && r_g2 > r_g1,
              "W5 对照：γ=2 与 γ=1 显著可分且 γ=2 更接近 1（√2 被数据选中）");
    }

    // ⚠️ 这里曾有子判据②：「壁面 γ=2 与【周期控制】的系统偏差同量级」。
    //    它靠同一台量测机器上的周期运行提供系统偏差基线。z 周期路径删除后
    //    基线不存在了 —— 用硬编码常数替代会是【假判据】，所以直接删掉。
    //    更强的替代（留待后续）：把周期控制换成 dt 控制，跑 dt 与 dt/2，
    //    断言 |ratio−1| 大致减半 —— 它识别偏差来源，而不是只比幅度。

    std::printf("\n%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}
