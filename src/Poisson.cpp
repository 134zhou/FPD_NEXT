#include "Poisson.h"

#include "Check.h"

// ============================================================================
// 压力泊松求解器：xy 2D 批量 FFT + z 向三对角（Thomas）【唯一路径】。
//
// 数学推导见 doc/PressurePoisson.md，本文件是它的代码对照。核心不变量：
// 泊松算子必须是 div∘grad 的精确复合 —— 修正步用的梯度、散度检验用的散度、
// 这里解算的算子，三者是同一套离散。三对角（Neumann）是「壁面法向面不修正」
// 直接导出的结论，不是假设。
//
// ⚠️ 曾有第二条路径：全周期 3D FFT + 精确离散本征值
//    0.5/(cos kx + cos ky + cos kz - 3)。它随 z 周期边界一起在 Phase 8-A 删除。
// ============================================================================

// ---------------------------------------------------------------------------
// (0,0) 奇异列的相容性投影 + 残差记录。
// λ_xy(0,0)=0 时该列矩阵 A0 有零空间 = 常向量，可解 ⟺ Σ_k b̂_k = 0（精确算术下
// 自动成立）。浮点下残差 ~1e-13 会积成 p 的线性漂移 → 恒定伪竖直加速度。
// 这里显式减掉均值（相容投影）并把 |Σ_k b̂_k| 记入 diag[0]。
// 注意：这是减【均值】而非把整列清零 —— 该列的 z 结构是支撑粒子重量的静压。
// ---------------------------------------------------------------------------
static void fix_singular_column(int Nx, int Ny, int Nz, double* fft_data, double* diag)
{
    #pragma acc serial present(fft_data, diag)
    {
        double sr = 0.0, si = 0.0;
        for (int k = 0; k < Nz; k++)
        {
            sr += fft_data[2 * IDX(0, 0, k)];
            si += fft_data[2 * IDX(0, 0, k) + 1];
        }
        diag[0] = sqrt(sr * sr + si * si);
        const double mr = sr / (double)Nz;
        const double mi = si / (double)Nz;
        for (int k = 0; k < Nz; k++)
        {
            fft_data[2 * IDX(0, 0, k)]     -= mr;
            fft_data[2 * IDX(0, 0, k) + 1] -= mi;
        }
    }
}

// ---------------------------------------------------------------------------
// Thomas 算法，就地解 Nx·Ny 个独立的三对角组（每个 (i,j) 一组，z 向 Nz 个未知数）。
//   次/超对角 a_k = c_k = 1
//   对角   d_k = λ_xy - 1（k=0 或 Nz-1），λ_xy - 2（其余）
// 前推系数 w_k 已由 build_tridiag_coeffs 预算进 tri_w（与右端无关）。
// 实部虚部用同一组 w_k，在同一遍 k 循环里同时推进（不是两遍扫描）。
//
// 并行/合并：collapse(2) over (i,j)，内层 k 串行。固定 k 时相邻线程（相邻 i）
// 访问相邻地址（fft_data 跨 2 个 double、tri_w 跨 1 个 double）→ 访存合并。
// ---------------------------------------------------------------------------
static void thomas_sweep(int Nx, int Ny, int Nz, double* fft_data, const double* tri_w)
{
    const int stride = Nx * Ny;   // z 向相邻层跨 Nx*Ny 个元素（IDX 的 k 最快…实际 x 最快）
    #pragma acc parallel loop collapse(2) present(fft_data, tri_w)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            const int base = IDX(i, j, 0);   // = i + j*Nx

            // --- 前代（y_0 = b_0·w_0,  y_k = (b_k - y_{k-1})·w_k）---
            double prev_re = fft_data[2 * base]     * tri_w[base];
            double prev_im = fft_data[2 * base + 1] * tri_w[base];
            fft_data[2 * base]     = prev_re;
            fft_data[2 * base + 1] = prev_im;
            for (int k = 1; k < Nz; k++)
            {
                const int idx = base + k * stride;
                const double w  = tri_w[idx];
                const double re = (fft_data[2 * idx]     - prev_re) * w;
                const double im = (fft_data[2 * idx + 1] - prev_im) * w;
                fft_data[2 * idx]     = re;
                fft_data[2 * idx + 1] = im;
                prev_re = re; prev_im = im;
            }

            // --- 回代（x_{Nz-1} = y_{Nz-1},  x_k = y_k - w_k·x_{k+1}）---
            for (int k = Nz - 2; k >= 0; k--)
            {
                const int idx = base + k * stride;
                const double w = tri_w[idx];
                fft_data[2 * idx]     -= w * fft_data[2 * (idx + stride)];
                fft_data[2 * idx + 1] -= w * fft_data[2 * (idx + stride) + 1];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2D 批量 FFT → 奇异列相容投影 → Thomas → 逆 2D FFT → p = Re/(Nx*Ny)。
// 归一化因子是 Nx*Ny 而非 size（2D 变换点数），见 doc/PressurePoisson.md §9。
// ---------------------------------------------------------------------------
void solve_pressure(NS_Config cfg, double* fft_data, const double* tri_w,
                       cufftHandle plan_xy, double* p, double* diag)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;

    #pragma acc host_data use_device(fft_data)
    {
        CUFFT_CHECK(cufftExecZ2Z(plan_xy, (cufftDoubleComplex*)fft_data,
                                 (cufftDoubleComplex*)fft_data, CUFFT_FORWARD));
    }

    fix_singular_column(Nx, Ny, Nz, fft_data, diag);
    thomas_sweep(Nx, Ny, Nz, fft_data, tri_w);

    #pragma acc host_data use_device(fft_data)
    {
        CUFFT_CHECK(cufftExecZ2Z(plan_xy, (cufftDoubleComplex*)fft_data,
                                 (cufftDoubleComplex*)fft_data, CUFFT_INVERSE));
    }

    const double nrm = 1.0 / (double)(Nx * Ny);
    #pragma acc parallel loop collapse(3) present(fft_data, p)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                int ijk = IDX(i, j, k);
                p[ijk] = fft_data[ijk * 2] * nrm;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 三对角前推系数（主机侧一次性预算）。含 (0,0) 列的定规 d_0 -= 1。
// ---------------------------------------------------------------------------
void build_tridiag_coeffs(NS_Config cfg, double* tri_w)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            const double lam = 2.0 * (cos(2.0 * M_PI * i / Nx) - 1.0)
                             + 2.0 * (cos(2.0 * M_PI * j / Ny) - 1.0);
            // d_0 = lam - 1；奇异列 (0,0) 加定规 d_0 -= 1（见 doc §8）
            double d0 = lam - 1.0;
            if (i == 0 && j == 0) { d0 -= 1.0; }
            double w = 1.0 / d0;
            tri_w[IDX(i, j, 0)] = w;
            for (int k = 1; k < Nz; k++)
            {
                const double dk = (k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
                w = 1.0 / (dk - w);
                tri_w[IDX(i, j, k)] = w;
            }
        }
    }
}

// ============================================================================
// 变系数投影（ratio_rho != 1）：算子 + PCG
//
// 与常系数路径的【唯一】区别是面心权重 a = 1/rho。离散、BC、零空间完全同构：
// a == 1 时本算子逐位退化成 solve_pressure 解的那个 div∘grad，所以 PCG 一步收敛。
// 文档对照见 doc/PressurePoisson.md §16。
// ============================================================================

void apply_laplacian_var(NS_Config cfg,
                         const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
                         const double* p, double* out)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;

    #pragma acc parallel loop collapse(3) present(inv_rho_x, inv_rho_y, inv_rho_z, p, out)
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                const int ijk = IDX(i, j, k);
                const int ip = (i + 1) % Nx;
                const int im = (i - 1 + Nx) % Nx;
                const int jp = (j + 1) % Ny;
                const int jm = (j - 1 + Ny) % Ny;

                double lap = inv_rho_x[ijk] * (p[IDX(ip, j, k)] - p[ijk])
                           - inv_rho_x[IDX(im, j, k)] * (p[ijk] - p[IDX(im, j, k)])
                           + inv_rho_y[ijk] * (p[IDX(i, jp, k)] - p[ijk])
                           - inv_rho_y[IDX(i, jm, k)] * (p[ijk] - p[IDX(i, jm, k)]);

                // 壁面通量：Neumann ghost 使壁面那一面的梯度为 0，所以直接不写该项。
                // 与常系数算子的 d_0 = lam-1（k=0/Nz-1 少一个 -1）是同一件事。
                if (k < Nz - 1) { lap += inv_rho_z[ijk] * (p[IDX(i, j, k + 1)] - p[ijk]); }
                if (k > 0)      { lap -= inv_rho_z[IDX(i, j, k - 1)] * (p[ijk] - p[IDX(i, j, k - 1)]); }

                out[ijk] = lap;
            }
        }
    }
}

int solve_pressure_var(NS_Config cfg,
                       double* cg_r,
                       const double* inv_rho_x, const double* inv_rho_y, const double* inv_rho_z,
                       double* fft_data, const double* tri_w, cufftHandle plan_xy,
                       double* p, double* diag,
                       double* cg_s, double* cg_zq)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny, Nz = cfg.Nz;
    const size_t size = (size_t)Nx * (size_t)Ny * (size_t)Nz;
    const double tol = 1e-12;
    const int max_iter = 200;

    // --- 1. 右端统计：sum(b) 用于相容性诊断与去均值，||b|| 用于相对判据 ---
    double sb = 0.0, nb = 0.0;
    #pragma acc parallel loop collapse(3) reduction(+:sb,nb) present(cg_r)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
            {
                const double v = cg_r[IDX(i, j, k)];
                sb += v;
                nb += v * v;
            }

    const double norm_b = sqrt(nb);

    // b == 0：p 取 0（差分为 0，修正步不动），不进入迭代。
    if (norm_b <= 0.0)
    {
        #pragma acc parallel loop collapse(3) present(p)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) { p[IDX(i, j, k)] = 0.0; }
        diag[0] = fabs(sb);
        diag[1] = 0.0;
        #pragma acc update device(diag[0:2])
        return 0;
    }

    // 预条件器的标量缩放：B = -a_bar * M，a_bar 是【实际参与算子的面】的平均 a。
    // 这不是新机制，只是把 M 的系数从 1 换成 a_bar：FFT 求解后再乘 1/a_bar。
    // 壁面 z 面（k = Nz-1）不参与算子，不计入平均。
    double sa = 0.0;
    #pragma acc parallel loop collapse(3) reduction(+:sa) present(inv_rho_x, inv_rho_y, inv_rho_z)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
            {
                const int ijk = IDX(i, j, k);
                sa += inv_rho_x[ijk] + inv_rho_y[ijk];
                if (k < Nz - 1) { sa += inv_rho_z[ijk]; }
            }
    const double nfaces = 2.0 * (double)size + (double)Nx * (double)Ny * (double)(Nz - 1);
    const double a_bar = (nfaces > 0.0) ? (sa / nfaces) : 1.0;
    const double inv_abar = (a_bar > 0.0) ? (1.0 / a_bar) : 1.0;

    // --- 2. r0 = -b 去均值 -------------------------------------------------
    // 【不用 warm start】：p 不入存档，若用上一步的 p 作初值，重启后的第一个 CG
    // 会走另一条（虽然收敛到同一解的）迭代路径，逐位重启就破了。这里固定 p0=0，
    // 只要 a、b 由 R 确定性地重算，CG 的每一轮算术都可逐位复现。
    const double bmean = sb / (double)size;
    #pragma acc parallel loop collapse(3) present(cg_r)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
            {
                const int ijk = IDX(i, j, k);
                cg_r[ijk] = -cg_r[ijk] + bmean;      // r0 = -b - mean(-b)
            }
    #pragma acc parallel loop collapse(3) present(p)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) { p[IDX(i, j, k)] = 0.0; }

    double rr = norm_b;          // 当前 ||r||；初始值只用于诊断

    // --- 3. 首次预条件：z = B^-1 r = -(1/a_bar) M^-1 r -----------------------
    // 载入 fft_data = -r，solve_pressure 解出 u = M^-1(-r) = -M^-1 r；
    // 再乘 1/a_bar 得 B^-1 r = -(1/a_bar) M^-1 r。
    // 它的任意常数不逐次去掉（见 Poisson.h 的说明）。
    #pragma acc parallel loop collapse(3) present(cg_r, fft_data)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
            {
                const int ijk = IDX(i, j, k);
                fft_data[2 * ijk]     = -cg_r[ijk];
                fft_data[2 * ijk + 1] = 0.0;
            }
    solve_pressure(cfg, fft_data, tri_w, plan_xy, cg_zq, diag);

    #pragma acc parallel loop collapse(3) present(cg_zq)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) { cg_zq[IDX(i, j, k)] *= inv_abar; }

    #pragma acc parallel loop collapse(3) present(cg_s, cg_zq)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) { cg_s[IDX(i, j, k)] = cg_zq[IDX(i, j, k)]; }

    double rz = 0.0;
    #pragma acc parallel loop collapse(3) reduction(+:rz) present(cg_r, cg_zq)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) { rz += cg_r[IDX(i, j, k)] * cg_zq[IDX(i, j, k)]; }

    int niter = 0;
    if (rz <= 0.0)
    {
        diag[0] = fabs(sb);
        diag[1] = rr / norm_b;
        #pragma acc update device(diag[0:2])
        return 0;
    }

    // --- 4. PCG 主循环 ------------------------------------------------------
    // A = -L。apply_laplacian_var 给的是 L s，所以 (s, A s) = -(s, L s)，
    // r 的更新是 r += alpha*(L s)（不是 -=）。
    for (int it = 0; it < max_iter; it++)
    {
        niter = it + 1;

        apply_laplacian_var(cfg, inv_rho_x, inv_rho_y, inv_rho_z, cg_s, cg_zq);

        double sq = 0.0;
        #pragma acc parallel loop collapse(3) reduction(+:sq) present(cg_s, cg_zq)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) { sq += cg_s[IDX(i, j, k)] * cg_zq[IDX(i, j, k)]; }

        if (sq >= 0.0) { break; }            // A 正定 => sq < 0；数值异常就停

        const double alpha = rz / (-sq);
        #pragma acc parallel loop collapse(3) present(p, cg_r, cg_s, cg_zq)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++)
                {
                    const int ijk = IDX(i, j, k);
                    p[ijk]     += alpha * cg_s[ijk];
                    cg_r[ijk]  += alpha * cg_zq[ijk];
                }

        double rr2 = 0.0;
        #pragma acc parallel loop collapse(3) reduction(+:rr2) present(cg_r)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) { const double v = cg_r[IDX(i, j, k)]; rr2 += v * v; }
        rr = sqrt(rr2);

        if (rr <= tol * norm_b) { break; }

        // 下一次预条件：z = B^-1 r
        #pragma acc parallel loop collapse(3) present(cg_r, fft_data)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++)
                {
                    const int ijk = IDX(i, j, k);
                    fft_data[2 * ijk]     = -cg_r[ijk];
                    fft_data[2 * ijk + 1] = 0.0;
                }
        solve_pressure(cfg, fft_data, tri_w, plan_xy, cg_zq, diag);

        #pragma acc parallel loop collapse(3) present(cg_zq)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) { cg_zq[IDX(i, j, k)] *= inv_abar; }

        double rz_new = 0.0;
        #pragma acc parallel loop collapse(3) reduction(+:rz_new) present(cg_r, cg_zq)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) { rz_new += cg_r[IDX(i, j, k)] * cg_zq[IDX(i, j, k)]; }

        if (rz_new <= 0.0) { break; }
        const double beta = rz_new / rz;
        rz = rz_new;

        #pragma acc parallel loop collapse(3) present(cg_s, cg_zq)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++)
                {
                    const int ijk = IDX(i, j, k);
                    cg_s[ijk] = cg_zq[ijk] + beta * cg_s[ijk];
                }
    }

    // --- 5. p 去均值（压力只以差分进入修正步；去掉累计的任意常数）---
    double pm = 0.0;
    #pragma acc parallel loop collapse(3) reduction(+:pm) present(p)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) { pm += p[IDX(i, j, k)]; }
    const double pmean = pm / (double)size;
    #pragma acc parallel loop collapse(3) present(p)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) { p[IDX(i, j, k)] -= pmean; }

    // 诊断量：最后一次预条件调用把 diag[0] 覆写成了预条件器自己的 |Σb̂|，
    // 这里恢复成真实右端的 |Σb|，并写入最终相对残差。
    diag[0] = fabs(sb);
    diag[1] = rr / norm_b;
    #pragma acc update device(diag[0:2])

    return niter;
}