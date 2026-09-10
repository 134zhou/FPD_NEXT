// spike：实测「xy 向批量 2D FFT + 手写 z 向 1D DFT」是否等价于 3D FFT。
//
// 背景：Phase 7-A 的壁面泊松求解器把 3D 全周期 FFT 换成
//   (1) 沿每个 z-slab 做 2D FFT（cufftPlanMany，batch=Nz）
//   (2) 对每个 (kx,ky) 模式在 z 向解三对角（Thomas）
// 在动手写 Poisson.cpp 之前，必须先证掉两件「不报错、只给错结果」的事：
//   (a) 轴序：cuFFT 是 row-major，n[] 最后一维最快。我们 x 最快 ⇒ n={Ny,Nx}。
//       搞反会得到一个转置的场（与 .fpd 的 (Nz,Ny,Nx) 是同一类陷阱）。
//   (b) 归一化：cuFFT 不归一化，forward+inverse 放大【变换点数】。
//       3D 版是 Nx*Ny*Nz，2D 批量版是 Nx*Ny。照抄 /size 会让 p 小 Nz 倍。
//
// 判据（对固定随机数据）：
//   P1  批量 2D FFT + 手写 z-DFT  vs  3D FFT，逐点相对差 < 1e-13
//   P2  forward→inverse 放大因子 == Nx*Ny（精确），而非 Nx*Ny*Nz
//   P3  inembed=NULL 与 inembed=n 给出相同结果（确认基本布局假设）
//
// 编译：
//   nvc++ -acc -O3 -I<CUDA_INC> -Isrc/include tests/spike/spike_fft2d_batch.cpp -lcufft -o /tmp/sf
//   /tmp/sf

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cufft.h>
#include <openacc.h>

static inline int IDX3(int i, int j, int k, int Nx, int Ny) { return i + j*Nx + k*Nx*Ny; }

// 确定性伪随机复数填充（避免 <random> 跨平台差异），直接填交错数组
static void fill_interleaved(std::vector<double>& buf, size_t n, unsigned seed)
{
    buf.assign(n * 2, 0.0);
    unsigned s = seed;
    for (size_t t = 0; t < n; t++)
    {
        s = s * 1664525u + 1013904223u;
        buf[t*2]   = ((double)(s & 0xFFFFFF) / 16777216.0) * 2.0 - 1.0;
        s = s * 1664525u + 1013904223u;
        buf[t*2+1] = ((double)(s & 0xFFFFFF) / 16777216.0) * 2.0 - 1.0;
    }
}

// 在设备上执行 forward FFT（原地）。buf 须已在设备上。
static void fft_forward(cufftHandle plan, double* d, size_t n)
{
    #pragma acc host_data use_device(d)
    {
        cufftExecZ2Z(plan, (cufftDoubleComplex*)d, (cufftDoubleComplex*)d, CUFFT_FORWARD);
    }
}

static double max_rel_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    const size_t n = a.size();   // 交错数组，n 为元素数（= 2*size）
    double mx = 0.0;
    for (size_t t = 0; t < n; t += 2)
    {
        const double ar = a[t], ai = a[t+1];
        const double br = b[t], bi = b[t+1];
        const double dre = fabs(ar - br), dim = fabs(ai - bi);
        const double amag = sqrt(ar*ar + ai*ai);
        const double bmag = sqrt(br*br + bi*bi);
        const double mag = amag > bmag ? amag : bmag;
        if (mag > 0.0)
        {
            const double rel = sqrt(dre*dre + dim*dim) / mag;
            if (rel > mx) { mx = rel; }
        }
    }
    return mx;
}

int main(int argc, char** argv)
{
    const int Nx = (argc > 1) ? atoi(argv[1]) : 8;
    const int Ny = (argc > 2) ? atoi(argv[2]) : 6;
    const int Nz = (argc > 3) ? atoi(argv[3]) : 5;
    const size_t size = (size_t)Nx * Ny * Nz;

    printf("=== spike_fft2d_batch：2D 批量 FFT vs 3D FFT ===\n");
    printf("盒子 %d x %d x %d   (size = %zu)\n\n", Nx, Ny, Nz, size);

    int n2[2] = { Ny, Nx };
    const int idist = Nx * Ny;

    int fails = 0;

    // ---------- P1：批量 2D FFT + 手写 z-DFT  vs  3D FFT ----------
    {
        std::vector<double> ref, bat;
        fill_interleaved(ref, size, 12345u);
        bat = ref;

        // 3D FFT
        cufftHandle p3 = 0;
        cufftPlan3d(&p3, Nz, Ny, Nx, CUFFT_Z2Z);
        double* d_ref = ref.data();
        #pragma acc enter data copyin(d_ref[0:2*size])
        fft_forward(p3, d_ref, size);
        #pragma acc exit data copyout(d_ref[0:2*size])
        cufftDestroy(p3);

        // 2D 批量 FFT
        cufftHandle p2 = 0;
        cufftPlanMany(&p2, 2, n2, NULL, 1, idist, NULL, 1, idist, CUFFT_Z2Z, Nz);
        double* d_bat = bat.data();
        #pragma acc enter data copyin(d_bat[0:2*size])
        fft_forward(p2, d_bat, size);
        #pragma acc exit data copyout(d_bat[0:2*size])
        cufftDestroy(p2);

        // 手写 z 向 1D DFT：对每个 (i,j) 模式，沿 k 做 forward DFT
        // X[i,j,l] = sum_k bat[i,j,k] * exp(-i 2pi l k / Nz)
        std::vector<double> z(size * 2, 0.0);
        for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
        for (int l = 0; l < Nz; l++)
        {
            double sr = 0.0, si = 0.0;
            for (int k = 0; k < Nz; k++)
            {
                const double th = 2.0 * M_PI * (double)l * (double)k / (double)Nz;
                const double c = cos(th), s = sin(th);
                const size_t idx = (size_t)IDX3(i, j, k, Nx, Ny);
                const double ar = bat[idx*2], ai = bat[idx*2+1];
                sr += ar * c + ai * s;   // (a+ib)(cos - i sin) 的实部
                si += ai * c - ar * s;   // 虚部
            }
            const size_t o = (size_t)IDX3(i, j, l, Nx, Ny);
            z[o*2] = sr; z[o*2+1] = si;
        }

        const double rel1 = max_rel_diff(z, ref);
        printf("  P1  2D批量+zDFT vs 3DFFT    max 相对差 = %.3e   %s\n",
               rel1, rel1 < 1e-13 ? "PASS" : "FAIL");
        if (!(rel1 < 1e-13)) { fails++; }
    }

    // ---------- P2：forward→inverse 放大因子 == Nx*Ny ----------
    {
        std::vector<double> buf, orig;
        fill_interleaved(orig, size, 4242u);
        buf = orig;
        cufftHandle p2 = 0;
        cufftPlanMany(&p2, 2, n2, NULL, 1, idist, NULL, 1, idist, CUFFT_Z2Z, Nz);
        double* d_buf = buf.data();
        #pragma acc enter data copyin(d_buf[0:2*size])
        #pragma acc host_data use_device(d_buf)
        {
            cufftExecZ2Z(p2, (cufftDoubleComplex*)d_buf, (cufftDoubleComplex*)d_buf, CUFFT_FORWARD);
            cufftExecZ2Z(p2, (cufftDoubleComplex*)d_buf, (cufftDoubleComplex*)d_buf, CUFFT_INVERSE);
        }
        #pragma acc exit data copyout(d_buf[0:2*size])
        cufftDestroy(p2);

        // buf 应 == 原始 * (Nx*Ny)
        const double scale = (double)(Nx * Ny);
        double mxabs = 0.0;
        for (size_t t = 0; t < 2*size; t++)
        {
            const double d = fabs(buf[t] - orig[t]*scale);
            if (d > mxabs) { mxabs = d; }
        }
        const bool ok2 = mxabs < 1e-12;
        printf("  P2  2D 批量 forward+inverse 放大 == Nx*Ny=%d   max|err| = %.3e   %s\n",
               Nx*Ny, mxabs, ok2 ? "PASS" : "FAIL");
        printf("      (若误用 size=Nx*Ny*Nz=%d，p 会小 Nz=%d 倍)\n", Nx*Ny*Nz, Nz);
        if (!ok2) { fails++; }
    }

    // ---------- P3：inembed=NULL 与 inembed=n 结果相同 ----------
    {
        std::vector<double> a, b;
        fill_interleaved(a, size, 7777u);
        b = a;

        cufftHandle pa = 0, pb = 0;
        cufftPlanMany(&pa, 2, n2, NULL, 1, idist, NULL, 1, idist, CUFFT_Z2Z, Nz);
        cufftPlanMany(&pb, 2, n2, n2,   1, idist, n2,   1, idist, CUFFT_Z2Z, Nz);
        double* da = a.data(); double* db = b.data();
        #pragma acc enter data copyin(da[0:2*size], db[0:2*size])
        fft_forward(pa, da, size);
        fft_forward(pb, db, size);
        #pragma acc exit data copyout(da[0:2*size], db[0:2*size])
        cufftDestroy(pa); cufftDestroy(pb);

        double mxabs = 0.0;
        for (size_t t = 0; t < 2*size; t++)
        {
            const double d = fabs(a[t] - b[t]);
            if (d > mxabs) { mxabs = d; }
        }
        const bool ok3 = mxabs == 0.0;   // 期望逐位相同
        printf("  P3  inembed=NULL vs inembed=n   max 绝对差 = %.3e   %s\n",
               mxabs, ok3 ? "PASS" : "FAIL");
        if (!ok3) { fails++; }
    }

    printf("\n%s\n", fails == 0 ? "全部通过" : "有失败项");
    return fails;
}
