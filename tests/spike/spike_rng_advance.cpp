// spike：测出 cuRAND 生成 n 个正态 double 后，内部序列实际前进了多少
//
// 背景：spike_rng_offset 发现「小批次多次调用能对齐、大批次对不上」，
// 说明生成 n 个数后内部位置不是前进 n，而是按并行网格容量向上取整。
// 这直接推翻了 rng_draws += 6*size 这种记账方式。
//
// 本 spike 用二分/扫描测出真实前进量 adv(n)，看它是否是 n 向上取整到某个粒度。
//
// 编译：nvc++ -acc -cuda -cudalib=curand -O2 tests/spike/spike_rng_advance.cpp -o build/spike_adv
#include <cstdio>
#include <cstring>
#include <vector>
#include <curand.h>
#include <cuda_runtime.h>

#define CK(x) do { if ((x) != CURAND_STATUS_SUCCESS) { \
    printf("cuRAND 错误 @ %d: %s\n", __LINE__, #x); return -1; } } while (0)

static const unsigned long long SEED = 20260905ULL;
static const curandRngType_t    TYPE = CURAND_RNG_PSEUDO_PHILOX4_32_10;

static int gen_at(unsigned long long off, double* host, size_t n)
{
    curandGenerator_t g;
    double* d = 0;
    if (cudaMalloc(&d, n * sizeof(double)) != cudaSuccess) { return -1; }
    CK(curandCreateGenerator(&g, TYPE));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    if (off > 0) { CK(curandSetGeneratorOffset(g, off)); }
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));
    cudaMemcpy(host, d, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d);
    CK(curandDestroyGenerator(g));
    return 0;
}

// 连续生成两批 n，返回第二批
static int gen_second(size_t n, double* host)
{
    curandGenerator_t g;
    double* d = 0;
    if (cudaMalloc(&d, n * sizeof(double)) != cudaSuccess) { return -1; }
    CK(curandCreateGenerator(&g, TYPE));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));   // 第一批，丢弃
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));   // 第二批
    cudaMemcpy(host, d, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d);
    CK(curandDestroyGenerator(g));
    return 0;
}

// 找出使 gen_at(off, n) == 连续生成的第二批 的 off
static long long find_advance(size_t n)
{
    const size_t CMP = n < 64 ? n : 64;      // 比对前 CMP 个就够判定
    std::vector<double> want(n), got(n);
    if (gen_second(n, want.data()) != 0) { return -1; }

    // 候选：n 向上取整到 2 的幂次粒度
    for (unsigned long long gran = 1; gran <= (1ULL << 22); gran <<= 1)
    {
        const unsigned long long adv = ((n + gran - 1) / gran) * gran;
        if (gen_at(adv, got.data(), n) != 0) { return -1; }
        if (memcmp(got.data(), want.data(), CMP * sizeof(double)) == 0)
        {
            return (long long)adv;
        }
    }
    return -1;
}

int main()
{
    printf("测出 cuRAND(Philox) 生成 n 个正态 double 后内部序列的实际前进量\n\n");
    printf("       n        实际前进 adv      adv/n     adv-n\n");

    const size_t ns[] = {8, 1024, 6144, 24576, 98304, 196608, 1572864};
    for (size_t t = 0; t < sizeof(ns)/sizeof(ns[0]); t++)
    {
        const size_t n = ns[t];
        const long long adv = find_advance(n);
        if (adv < 0)
        {
            printf("  %8zu        未找到（不是 2 的幂次取整）\n", n);
            continue;
        }
        printf("  %8zu   %14lld   %8.4f  %8lld\n",
               n, adv, (double)adv / (double)n, adv - (long long)n);
    }

    printf("\n若 adv 恒等于 n，则累计记账可用；否则必须在每次生成前显式 setOffset。\n");
    return 0;
}
