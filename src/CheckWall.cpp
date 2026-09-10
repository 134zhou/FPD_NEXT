#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <openacc.h>

#include "./include/State.h"
#include "./include/Stokes.h"
#include "./include/Viscosity.h"
#include "./include/Force.h"
#include "./include/Velocity.h"
#include "./include/Tests.h"

// ============================================================================
// z 向无滑移壁面的端到端判据（Phase 7-B）。
//
//   W1  max|div v| 含两个壁面层        —— doc/PressurePoisson.md §13 盲区 1 的正解
//   W2  vz 在两壁被钉死                —— 法向边界条件
//   W3  平面 Poiseuille 的【精确离散】闭式 —— 切向 ghost 与散度端系数的判决性判据
//   W6  z 向动量收支恒等式（逐步精确）  —— fz 截断 / vz 钉死 / ∂p/∂z / Π_zz 四件事一次全查
//
// 这些判据都【不开噪声】。把「算子对不对」与「噪声幅度对不对」彻底分离 ——
// 噪声幅度（√2 因子）由 --noise-wall 的 W5 单独负责。
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
                           st.fft, st.plan, st.plan_xy, st.tri_w, st.diag,
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
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, kT, /*noise_on=*/false, /*wall_z=*/1);

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
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, 0.0, /*noise_on=*/false, /*wall_z=*/1);

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
    NS_Config cfg = make_ns_config(Nx, Ny, Nz, dt, 0.0, /*noise_on=*/false, /*wall_z=*/1);

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
                           st.fft, st.plan, st.plan_xy, st.tri_w, st.diag,
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

    NS_Config cfg = make_ns_config(Nx, Ny, Nz, 0.002, 0.25, true, /*wall_z=*/1);
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

    std::printf("\n[W4] 壁面截断下的力守恒与粒子速度（含贴壁与 N=2 重叠）\n");
    w4_truncation();

    std::printf("\n[W6] z 向动量收支恒等式\n");
    w6_momentum_budget(Nx, Ny, Nz);

    std::printf("\n%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}
