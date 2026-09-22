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