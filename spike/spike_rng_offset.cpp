// spike：验证 curandSetGeneratorOffset 的计数单位
//
// 逐位重启依赖于「跳到序列第 n 个位置」，但 curandGenerateNormalDouble 走
// Box-Muller，offset 的计数单位与产出的 double 个数不一定 1:1，且没有稳定文档保证。
//
// 本 spike 直接测出换算系数：连续生成 M 个 vs 从 offset 处生成后半段，比对是否逐位相同。
//
// 编译：nvc++ -acc -cuda -cudalib=curand -O2 spike/spike_rng_offset.cpp -o build/spike_rng
#include <cstdio>
#include <cstring>
#include <vector>
#include <curand.h>
#include <cuda_runtime.h>

#define CK(x) do { if ((x) != CURAND_STATUS_SUCCESS) { \
    printf("cuRAND 错误 @ %d: %s\n", __LINE__, #x); return 1; } } while (0)

static const unsigned long long SEED = 20260905ULL;

// 生成 n 个正态 double；skip>0 时先 setOffset(skip)
static int gen(curandGenerator_t g, unsigned long long skip,
               double* host, size_t n)
{
    double* d = 0;
    if (cudaMalloc(&d, n * sizeof(double)) != cudaSuccess) { return 1; }
    if (skip > 0) { CK(curandSetGeneratorOffset(g, skip)); }
    CK(curandGenerateNormalDouble(g, d, n, 0.0, 1.0));
    cudaMemcpy(host, d, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d);
    return 0;
}

static int run_for(curandRngType_t type, const char* name)
{
    printf("\n=== %s ===\n", name);

    const size_t HALF = 8;
    const size_t FULL = 2 * HALF;

    std::vector<double> ref(FULL), tst(HALF);

    // 参考：从头连续生成 FULL 个
    curandGenerator_t g;
    CK(curandCreateGenerator(&g, type));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    if (gen(g, 0, ref.data(), FULL)) { return 1; }
    CK(curandDestroyGenerator(g));

    printf("  参考序列前 4 个: %.6f %.6f %.6f %.6f\n", ref[0], ref[1], ref[2], ref[3]);
    printf("  参考序列后半段前 2 个: %.6f %.6f\n", ref[HALF], ref[HALF + 1]);

    // 试各种换算系数，找出哪一个能让 setOffset 跳到后半段
    const unsigned long long factors[] = {1, 2, 4};
    int found = -1;
    for (int fi = 0; fi < 3; fi++)
    {
        const unsigned long long off = HALF * factors[fi];
        CK(curandCreateGenerator(&g, type));
        CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
        if (gen(g, off, tst.data(), HALF)) { return 1; }
        CK(curandDestroyGenerator(g));

        const bool same = (memcmp(tst.data(), ref.data() + HALF,
                                  HALF * sizeof(double)) == 0);
        printf("  offset = %zu * %llu = %-4llu  -> 首值 %.6f   %s\n",
               HALF, factors[fi], off, tst[0], same ? "逐位相同 ✓" : "不同");
        if (same && found < 0) { found = fi; }
    }

    if (found >= 0)
    {
        printf("  => 换算系数 = %llu（offset 单位 = double 个数 x %llu）\n",
               factors[found], factors[found]);
    }
    else
    {
        printf("  => 【没有一个系数对得上】此发生器不能用于逐位重启\n");
    }
    return found >= 0 ? 0 : 1;
}

// 生产场景：setOffset 之后【连续多次】调用 generate，看是否仍与连续生成对齐。
// spike 的第一部分只测了「setOffset 后生成一次」，两者未必等价。
static int run_multicall(curandRngType_t type, const char* name, size_t batch, int nbatch)
{
    printf("\n=== %s：多批次连续生成（每批 %zu 个，共 %d 批）===\n", name, batch, nbatch);

    const size_t total = batch * (size_t)nbatch;
    std::vector<double> ref(total), tst(total);

    // 参考：一个发生器连续跑 nbatch 批
    curandGenerator_t g;
    CK(curandCreateGenerator(&g, type));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    for (int b = 0; b < nbatch; b++)
    {
        if (gen(g, 0, ref.data() + b * batch, batch)) { return 1; }
    }
    CK(curandDestroyGenerator(g));

    // 测试：新发生器，跳过前 half 批，再连续跑剩下的
    const int half = nbatch / 2;
    CK(curandCreateGenerator(&g, type));
    CK(curandSetPseudoRandomGeneratorSeed(g, SEED));
    CK(curandSetGeneratorOffset(g, (unsigned long long)(batch * half)));
    for (int b = half; b < nbatch; b++)
    {
        if (gen(g, 0, tst.data() + b * batch, batch)) { return 1; }
    }
    CK(curandDestroyGenerator(g));

    int firstbad = -1;
    for (int b = half; b < nbatch; b++)
    {
        if (memcmp(tst.data() + b * batch, ref.data() + b * batch,
                   batch * sizeof(double)) != 0)
        { firstbad = b; break; }
    }
    if (firstbad < 0)
    {
        printf("  跳过前 %d 批后，其余各批全部逐位相同 ✓\n", half);
        return 0;
    }
    printf("  第 %d 批开始不同（前 %d 批跳过）\n", firstbad, half);
    printf("    参考 %.6f %.6f   实测 %.6f %.6f\n",
           ref[firstbad * batch], ref[firstbad * batch + 1],
           tst[firstbad * batch], tst[firstbad * batch + 1]);
    return 1;
}

int main()
{
    printf("测试 curandSetGeneratorOffset 对 curandGenerateNormalDouble 的计数单位\n");

    int bad = 0;
    bad |= run_for(CURAND_RNG_PSEUDO_PHILOX4_32_10, "PHILOX4_32_10");
    bad |= run_for(CURAND_RNG_PSEUDO_XORWOW,        "XORWOW (当前在用)");
    bad |= run_for(CURAND_RNG_PSEUDO_MRG32K3A,      "MRG32K3A");

    bad |= run_multicall(CURAND_RNG_PSEUDO_PHILOX4_32_10, "PHILOX4_32_10", 8, 6);
    bad |= run_multicall(CURAND_RNG_PSEUDO_PHILOX4_32_10, "PHILOX4_32_10 (大批次)", 98304, 6);

    printf("\n%s\n", bad ? "有情形不支持精确 offset" : "全部可用");
    return 0;
}
