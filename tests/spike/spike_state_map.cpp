// spike：验证 FpdState 的映射机制 —— 这是 FpdState 设计的唯一硬阻塞，必须先测。
//
// 要回答的问题：
//   Q1  运行时 API（acc_create/acc_copyin/acc_update_device/acc_update_self/acc_delete）
//       建的映射，能被【另一个 TU】里的 present() 子句认账吗？
//   Q2  成员函数里的 #pragma acc enter data / update host 是否和 lambda 一样失效？
//   Q3  自由函数 + 指针参数上的 pragma（退路 B）是否可靠？
//   Q4  acc_create 与 acc_delete 的引用计数是否对称？
//
// 编译：
//   nvc++ -acc -O2 -Minfo=accel tests/spike/spike_state_map.cpp tests/spike/spike_state_map_kernel.cpp -o build/spike_state
//   ./build/spike_state
//
// 判据：kernel 里 p[i]*2 写到 out[i]，update host 拉回后 out[i] 必须【逐位】等于
//       2*p[i]；acc_is_present 在 map 后为真、unmap 后为假。

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <openacc.h>

extern "C" void kernel_times2(double* p, double* out, int n);

static int g_fail = 0;

static void check(bool ok, const char* what)
{
    printf("  %-58s %s\n", what, ok ? "PASS" : "FAIL");
    if (!ok) { g_fail++; }
}

// 逐位核对 out[i] == 2*p[i]
static bool verify_doubled(const std::vector<double>& p,
                           const std::vector<double>& out, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (out[i] != 2.0 * p[i]) { return false; }   // 用 == 要求逐位
    }
    return true;
}

// ---------------------------------------------------------------------------
// 方式 1：运行时 API，放进成员函数
// ---------------------------------------------------------------------------
struct ViaRuntimeAPI
{
    std::vector<double> p, out;
    double* d_p = 0;      // 裸指针 = vector::data()
    double* d_out = 0;
    int n = 0;

    void init(int n_)
    {
        n = n_;
        p.assign(n, 0); out.assign(n, 0);
        for (int i = 0; i < n; i++) { p[i] = (double)(i + 1) * 0.5; }
        d_p = p.data(); d_out = out.data();
    }

    // 关键：全部用运行时 API，不用任何 pragma
    void map()
    {
        acc_copyin(d_p, n * sizeof(double));      // = enter data copyin
        acc_create(d_out, n * sizeof(double));    // = enter data create
    }
    void download() { acc_update_self(d_out, n * sizeof(double)); }
    void unmap()
    {
        acc_delete(d_p, n * sizeof(double));
        acc_delete(d_out, n * sizeof(double));
    }
    bool is_present() const
    {
        return acc_is_present(d_p, n * sizeof(double)) &&
               acc_is_present(d_out, n * sizeof(double));
    }
};

// ---------------------------------------------------------------------------
// 方式 2：成员函数里的 pragma（被测对象 —— CLAUDE.md 只说了 lambda 会挂）
// ---------------------------------------------------------------------------
struct ViaMemberPragma
{
    std::vector<double> p, out;
    double* d_p = 0;
    double* d_out = 0;
    int n = 0;

    void init(int n_)
    {
        n = n_;
        p.assign(n, 0); out.assign(n, 0);
        for (int i = 0; i < n; i++) { p[i] = (double)(i + 1) * 0.5; }
        d_p = p.data(); d_out = out.data();
    }

    void map()
    {
        #pragma acc enter data copyin(d_p[0:n]) create(d_out[0:n])
    }
    void download()
    {
        #pragma acc update host(d_out[0:n])
    }
    void unmap()
    {
        #pragma acc exit data delete(d_p[0:n], d_out[0:n])
    }
};

// ---------------------------------------------------------------------------
// 方式 3：自由函数 + 指针参数上的 pragma（退路 B，Force.cpp:31 的 present 已证明可靠）
// ---------------------------------------------------------------------------
static void free_map_create(double* p, double* out, int n)
{
    #pragma acc enter data copyin(p[0:n]) create(out[0:n])
}
static void free_update_host(double* out, int n)
{
    #pragma acc update host(out[0:n])
}
static void free_unmap(double* p, double* out, int n)
{
    #pragma acc exit data delete(p[0:n], out[0:n])
}

// ---------------------------------------------------------------------------
int main()
{
    const int N = 16;

    printf("=== spike_state_map：验证 FpdState 的映射机制 ===\n\n");

    // ---- 方式 1：运行时 API ----
    printf("[方式 1] 运行时 API 建映射 + 另一 TU 的 present() 认账\n");
    {
        ViaRuntimeAPI s;
        s.init(N);
        s.map();
        check(s.is_present(), "Q1/Q4: acc_copyin/acc_create 后 acc_is_present == true");

        // 先不映射就调 kernel 会 crash（present 找不到），所以这里已经 map 了
        kernel_times2(s.d_p, s.d_out, N);
        s.download();
        check(verify_doubled(s.p, s.out, N), "Q1: kernel(present) 认账运行时 API 建的映射，out == 2*p 逐位");

        // Q4：重复 create 同地址应递增引用计数，delete 一次后应仍在 present 表
        acc_create(s.d_p, N * sizeof(double));           // 引用计数 +1
        acc_delete(s.d_p, N * sizeof(double));           // 引用计数 -1，回到 1
        check(acc_is_present(s.d_p, N * sizeof(double)),
              "Q4: create 两次 + delete 一次后仍在 present 表（引用计数对称）");
        s.unmap();
        check(!s.is_present(), "Q1: acc_delete 后 acc_is_present == false");
        printf("\n");
    }

    // ---- 方式 2：成员函数里的 pragma ----
    printf("[方式 2] 成员函数里的 #pragma acc enter data / update host\n");
    {
        ViaMemberPragma s;
        s.init(N);
        s.map();
        kernel_times2(s.d_p, s.d_out, N);
        s.download();
        check(verify_doubled(s.p, s.out, N), "Q2: 成员函数 pragma 建映射，present() 认账，out == 2*p 逐位");
        s.unmap();
        printf("\n");
    }

    // ---- 方式 3：自由函数 + 指针参数的 pragma ----
    printf("[方式 3] 自由函数 + 指针参数上的 pragma（退路 B）\n");
    {
        std::vector<double> p(N), out(N);
        for (int i = 0; i < N; i++) { p[i] = (double)(i + 1) * 0.5; }
        double* dp = p.data(); double* dout = out.data();
        free_map_create(dp, dout, N);
        kernel_times2(dp, dout, N);
        free_update_host(dout, N);
        check(verify_doubled(p, out, N), "Q3: 自由函数 pragma 建映射，present() 认账，out == 2*p 逐位");
        free_unmap(dp, dout, N);
        check(!acc_is_present(dp, N * sizeof(double)),
              "Q3: 自由函数 unmap 后 acc_is_present == false");
        printf("\n");
    }

    printf("===============================================\n");
    printf("%s\n", g_fail == 0 ? "全部通过" : "有失败项");
    return g_fail;
}
