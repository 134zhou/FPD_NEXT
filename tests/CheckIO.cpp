#include "IOBin.h"
#include <cstdio>
#include <cstring>
#include <vector>

// 与 Python 独立构造已知数据，检查所有保留字段的轴序、偏移和逐位读写。
// 目录内 python_N.fpd 由 baseline/ckpt_v4_20261025/interop.py 生成。
int run_io_check(const char* directory)
{
    int fails = 0;
    for (int N : {0, 2})
    {
        CkptHeader h;
        h.Nx = 8; h.Ny = 4; h.Nz = 2; h.N = N; h.step = 7;
        std::vector<double> expected[12], actual[12];
        CkptArrays a;
        double** slots[] = {&a.vx, &a.vy, &a.vz, &a.Rx, &a.Ry, &a.Rz,
                           &a.Rux, &a.Ruy, &a.Ruz, &a.Vx, &a.Vy, &a.Vz};
        for (int t = 0; t < 12; t++)
        {
            const size_t count = t < 3 ? h.size() : (size_t)N;
            expected[t].resize(count);
            actual[t].resize(count);
            *slots[t] = actual[t].data();
            for (size_t q = 0; q < count; q++)
            {
                expected[t][q] = t < 3 ? (q % 8) * 1e6 + ((q / 8) % 4) * 1e3 + q / 32 + t * 0.25
                                        : t * 100.0 + q * 0.5;
            }
        }
        const std::string suffix = std::to_string(N) + ".fpd";
        const std::string input = std::string(directory) + "/python_" + suffix;
        const std::string output = std::string(directory) + "/cpp_" + suffix;
        CkptHeader got;
        std::string err;
        bool ok = load_checkpoint(input.c_str(), got, a, err);
        ok = ok && got.version == 4 && got.Nx == 8 && got.Ny == 4 && got.Nz == 2
                && got.N == N && got.step == 7;
        for (int t = 0; t < 12 && ok; t++)
        {
            if (!expected[t].empty())
            {
                ok = memcmp(expected[t].data(), actual[t].data(), expected[t].size() * sizeof(double)) == 0;
            }
        }
        if (ok) { ok = save_checkpoint(output.c_str(), h, a, err); }
        printf("%s IO N=%d：Python 输入全部字段逐位正确，C++ 写出 %s\n",
               ok ? "PASS" : "FAIL", N, err.c_str());
        fails += !ok;

        if (N == 2)
        {
            for (const char* name : {"old", "corrupt", "truncated", "trailing"})
            {
                const std::string path = std::string(directory) + "/" + name + ".fpd";
                err.clear();
                const bool rejected = !load_checkpoint(path.c_str(), got, a, err);
                printf("%s IO 拒绝 %s：%s\n", rejected ? "PASS" : "FAIL", name, err.c_str());
                fails += !rejected;
            }
        }
    }
    return fails;
}
