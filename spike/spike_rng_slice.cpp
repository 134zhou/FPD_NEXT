// spike：cuRAND 正态生成的结果是否只依赖 (offset, n)？分批与整批是否等价？
//
// 这是逐位重启的症结。三个问题：
//   Q1 同样 (offset, n) 生成两次，是否逐位相同？（确定性）
//   Q2 一次生成 2n，与在同一发生器上连续生成两次 n，是否拼接相等？（分批等价性）
//   Q3 setOffset(n) 后生成 n，是否等于「连续两次生成 n」的第二次？（可跳转性）
//
// 编译：nvc++ -acc -cuda -cudalib=curand -O2 spike/spike_rng_slice.cpp -o build/spike_slice
#include <cstdio>
#include <cstring>
#include <vector>
#include <curand.h>
#include <cuda_runtime.h>

#define CK(x) do { if ((x) != CURAND_STATUS_SUCCESS) { \
    printf("cuRAND 错误 @ %d\n", __LINE__); return 1; } } while (0)

static const unsigned long long SEED = 20260905ULL;
static const curandRngType_t    TYPE = CURAND_RNG_PSEUDO_PHILOX4_32_10;

// 返回首个不同元素的下标，全同返回 -1
static long first_diff(const double* a, const double* b, size_t n)
{
    for (size_t i = 0; i < n; i++)
    {
        if (memcmp(&a[i], &b[i], sizeof(double)) != 0) { return (long)i; }
    }
    return -1;
}

static int probe(size_t n)
{
    printf("\n=== n = %zu ===\n", n);

    double* d = 0;
    if (cudaMalloc(&d, 2 * n * sizeof(double)) != cudaSuccess) { return 1; }
    curandGenerator_t g;

    std::vector<double> one(2 * n), two(2 * n), rep(n), jump(n);

    // 整批：一次生成 2n
    CK(curandCreateGenerator(&g, TYPE));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    CK(curandGenerateNormalDouble(g, d, 2 * n, 0.0, 1.0));
    cudaMemcpy(one.data(), d, 2 * n * sizeof(double), cudaMemcpyDeviceToHost);
    CK(curandDestroyGenerator(g));

    // 分批：同一发生器连续两次 n
    CK(curandCreateGenerator(&g, TYPE));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));
    cudaMemcpy(two.data(), d, n * sizeof(double), cudaMemcpyDeviceToHost);
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));
    cudaMemcpy(two.data() + n, d, n * sizeof(double), cudaMemcpyDeviceToHost);
    CK(curandDestroyGenerator(g));

    // 重复：再来一次「分批」的第一批，验证确定性
    CK(curandCreateGenerator(&g, TYPE));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));
    cudaMemcpy(rep.data(), d, n * sizeof(double), cudaMemcpyDeviceToHost);
    CK(curandDestroyGenerator(g));

    // 跳转：setOffset(n) 后生成 n
    CK(curandCreateGenerator(&g, TYPE));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    CK(curandSetGeneratorOffset(g, (unsigned long long)n));
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));
    cudaMemcpy(jump.data(), d, n * sizeof(double), cudaMemcpyDeviceToHost);
    CK(curandDestroyGenerator(g));

    cudaFree(d);

    const long q1 = first_diff(two.data(), rep.data(), n);
    printf("  Q1 确定性（同 offset 同 n 生成两次）      : %s\n",
           q1 < 0 ? "逐位相同 ✓" : "不同 ✗");

    const long q2 = first_diff(one.data(), two.data(), 2 * n);
    if (q2 < 0) { printf("  Q2 整批 2n == 分两批 n+n                  : 逐位相同 ✓\n"); }
    else        { printf("  Q2 整批 2n == 分两批 n+n                  : 第 %ld 个起不同 ✗\n", q2); }

    const long q3 = first_diff(jump.data(), two.data() + n, n);
    if (q3 < 0) { printf("  Q3 setOffset(n)+gen(n) == 分批的第二批    : 逐位相同 ✓\n"); }
    else        { printf("  Q3 setOffset(n)+gen(n) == 分批的第二批    : 第 %ld 个起不同 ✗ (共 %zu)\n",
                         q3, n); }
    return 0;
}

int main()
{
    printf("探查 cuRAND(Philox) 正态生成的分批/跳转语义\n");
    probe(8);
    probe(1024);
    probe(98304);       // 32^3 的 size*3，生产实际用的批量
    printf("\nQ3 是逐位重启能否成立的关键。\n");
    return 0;
}
