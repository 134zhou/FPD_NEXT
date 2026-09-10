// spike：实测 Thomas 算法解 z 向三对角（壁面泊松的每 (kx,ky) 列），
// 对照稠密 Gauss 消元，并验证 (0,0) 奇异列的 d_0 定规。
//
// 三对角矩阵（负半定的离散 Neumann 拉普拉斯，λ_xy ≤ 0）：
//   次/超对角 a_k = c_k = 1
//   对角   d_k = λ_xy - 1   (k = 0 或 k = Nz-1)
//          d_k = λ_xy - 2   (其余)
// 奇异列 λ_xy = 0 时 A·1 = 0（行和为 0），定规是 d_0 -= 1（即 -1 → -2）。
//
// 判据（纯 CPU，无卡可跑）：
//   T1  Thomas vs 稠密 Gauss，Nz ∈ {2,6,32,64} × 若干 λ_xy<0：相对差 < 1e-12
//   T2  奇异列 λ_xy=0 加定规后：x_0 == 0（到 1e-16）且残差 |A0 x - b| < 1e-13
//   T3  相容性被破坏（人为 Σb≠0）时解会漂移 —— 证明「减均值」不是多余的
//
// 编译（任意 C++ 编译器即可，无 CUDA 依赖）：
//   g++ -O2 tests/spike/spike_tridiag.cpp -o /tmp/st && /tmp/st

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

// 前推系数 w_k（只依赖 λ_xy 与 k，与右端无关）。
// w_0 = 1/d_0,  w_k = 1/(d_k - w_{k-1})
static void build_w(int Nz, double lam, std::vector<double>& w)
{
    w.assign(Nz, 0.0);
    double d0 = lam - 1.0;           // 顶层
    w[0] = 1.0 / d0;
    for (int k = 1; k < Nz; k++)
    {
        const double dk = (k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
        w[k] = 1.0 / (dk - w[k-1]);
    }
}

// Thomas 就地解（实部；复数右端用同一组 w 各扫一遍）。返回 x 覆盖 b。
static void thomas(const std::vector<double>& w, std::vector<double>& b)
{
    const int Nz = (int)b.size();
    // 前代
    b[0] *= w[0];
    for (int k = 1; k < Nz; k++)
    {
        b[k] = (b[k] - b[k-1]) * w[k];
    }
    // 回代
    for (int k = Nz - 2; k >= 0; k--)
    {
        b[k] = b[k] - w[k] * b[k+1];
    }
}

// 稠密 Gauss 消元（带部分选主元），O(Nz^3)，只在小 Nz 上跑作独立参考。
static void dense_gauss(int Nz, double lam, const std::vector<double>& b, std::vector<double>& x)
{
    std::vector<double> A((size_t)Nz * Nz, 0.0);
    // 组装 A
    for (int k = 0; k < Nz; k++)
    {
        const double dk = (k == 0 || k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
        A[(size_t)k * Nz + k] = dk;
        if (k > 0)       { A[(size_t)k * Nz + (k-1)] = 1.0; }
        if (k < Nz - 1) { A[(size_t)k * Nz + (k+1)] = 1.0; }
    }
    std::vector<double> M = A;
    x = b;

    for (int col = 0; col < Nz; col++)
    {
        // 部分选主元
        int piv = col;
        double best = fabs(M[(size_t)col * Nz + col]);
        for (int r = col + 1; r < Nz; r++)
        {
            const double v = fabs(M[(size_t)r * Nz + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (piv != col)
        {
            for (int c = 0; c < Nz; c++) { std::swap(M[(size_t)col*Nz+c], M[(size_t)piv*Nz+c]); }
            std::swap(x[col], x[piv]);
        }
        const double pivv = M[(size_t)col * Nz + col];
        for (int r = col + 1; r < Nz; r++)
        {
            const double f = M[(size_t)r * Nz + col] / pivv;
            for (int c = col; c < Nz; c++) { M[(size_t)r * Nz + c] -= f * M[(size_t)col * Nz + c]; }
            x[r] -= f * x[col];
        }
    }
    for (int r = Nz - 1; r >= 0; r--)
    {
        for (int c = r + 1; c < Nz; c++) { x[r] -= M[(size_t)r * Nz + c] * x[c]; }
        x[r] /= M[(size_t)r * Nz + r];
    }
}

// 残差 |A x - b|（A 用同一个组装函数）
static double residual(int Nz, double lam, const std::vector<double>& x, const std::vector<double>& b)
{
    double mx = 0.0;
    for (int k = 0; k < Nz; k++)
    {
        const double dk = (k == 0 || k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
        double lhs = dk * x[k];
        if (k > 0)       { lhs += x[k-1]; }
        if (k < Nz - 1) { lhs += x[k+1]; }
        const double r = fabs(lhs - b[k]);
        if (r > mx) { mx = r; }
    }
    return mx;
}

// 确定性伪随机实向量
static void fill_b(std::vector<double>& b, unsigned seed)
{
    unsigned s = seed;
    for (size_t t = 0; t < b.size(); t++)
    {
        s = s * 1664525u + 1013904223u;
        b[t] = ((double)(s & 0xFFFFFF) / 16777216.0) * 2.0 - 1.0;
    }
}

int main()
{
    int fails = 0;
    printf("=== spike_tridiag：Thomas vs 稠密 Gauss + 奇异列定规 ===\n\n");

    // ---------- T1：Thomas vs 稠密 Gauss ----------
    {
        const int Nzs[] = {2, 6, 32, 64};
        const double lams[] = {-0.1, -1.0, -3.0, -8.0};   // 只测非奇异列（λ_xy<0）
        for (int a = 0; a < 4; a++)
        {
            const int Nz = Nzs[a];
            for (int b = 0; b < 4; b++)
            {
                const double lam = lams[b];
                std::vector<double> w, bt, bg;
                build_w(Nz, lam, w);
                bt.assign(Nz, 0.0);
                fill_b(bt, 9000u + (unsigned)Nz * 7u + (unsigned)b);
                bg = bt;
                thomas(w, bt);
                std::vector<double> xg;
                dense_gauss(Nz, lam, bg, xg);

                double maxrel = 0.0, maxmag = 0.0;
                for (int k = 0; k < Nz; k++)
                {
                    const double m = fabs(xg[k]);
                    if (m > maxmag) { maxmag = m; }
                    const double d = fabs(bt[k] - xg[k]);
                    if (d > maxrel) { maxrel = d; }
                }
                maxrel = maxmag > 0.0 ? maxrel / maxmag : maxrel;
                const bool ok = maxrel < 1e-12;
                if (!ok) { fails++; }
                printf("  T1  Nz=%2d  lam=%+.2f   Thomas vs Gauss  max 相对差 = %.3e   %s\n",
                       Nz, lam, maxrel, ok ? "PASS" : "FAIL");
            }
        }
    }

    // ---------- T2：奇异列 λ_xy=0 加定规 ----------
    {
        const int Nz = 32;
        // 先构造一个相容的 b（Σb = 0）：取一个随机 p，作用 A0 得 b，再令 p0=0
        std::vector<double> pref(Nz), b(Nz);
        fill_b(pref, 4242u);
        const double lam = 0.0;
        // 计算 b = A0 pref（用 residual 的组装，但这里手动算）
        for (int k = 0; k < Nz; k++)
        {
            const double dk = (k == 0 || k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
            double lhs = dk * pref[k];
            if (k > 0)       { lhs += pref[k-1]; }
            if (k < Nz - 1) { lhs += pref[k+1]; }
            b[k] = lhs;
        }
        // 定规：d_0 -= 1（-1 → -2）
        std::vector<double> w(Nz);
        w[0] = 1.0 / (lam - 2.0);          // d_0 已定规
        for (int k = 1; k < Nz; k++)
        {
            const double dk = (k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
            w[k] = 1.0 / (dk - w[k-1]);
        }
        std::vector<double> x = b;
        thomas(w, x);

        // x 应精确满足原方程 A0 x = b，且 x_0 = 0（浮点下 ~机器精度，非 O(1)）
        const double res = residual(Nz, lam, x, b);
        const double x0 = fabs(x[0]);
        const bool ok_res = res < 1e-13;
        const bool ok_x0 = x0 < 1e-13;
        printf("  T2  奇异列定规  Nz=%d   |A0 x - b| = %.3e (%s)   |x_0| = %.3e (%s)\n",
               Nz, res, ok_res ? "PASS" : "FAIL", x0, ok_x0 ? "PASS" : "FAIL");
        if (!ok_res || !ok_x0) { fails++; }

        // 与无定规的稠密 Gauss（最小二乘/伪逆）对比：去掉 x_0 后应逐点一致
        // （这里只验残差与 x_0，证明定规不是罚函数近似）
    }

    // ---------- T3：相容性被破坏时的漂移 ----------
    {
        const int Nz = 32;
        const double lam = 0.0;
        std::vector<double> w(Nz);
        w[0] = 1.0 / (lam - 2.0);
        for (int k = 1; k < Nz; k++)
        {
            const double dk = (k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
            w[k] = 1.0 / (dk - w[k-1]);
        }
        // 不相容右端：全 1（Σb = Nz ≠ 0）
        std::vector<double> b(Nz, 1.0);
        std::vector<double> x = b;
        thomas(w, x);

        // 不相容时解会有线性漂移：相邻差分近似恒定且非零。
        // 减均值后 Σx = 0；不减均值 Σx = Nz·(大数)。这里直接验「解不是常数」（常数
        // 是齐次解，会被不相容分量放大）。用 |x[Nz-1] - x[0]| 衡量漂移。
        const double drift = fabs(x[Nz-1] - x[0]);
        // 减掉均值后（生产路径会做），残差应小；这里展示不减均值的后果
        printf("  T3  不相容右端 Σb=%d≠0  解漂移 |x[Nz-1]-x[0]| = %.3e   %s\n",
               Nz, drift, drift > 1e-3 ? "PASS（证明减均值非多余）" : "FAIL（漂移没测出来）");
        if (!(drift > 1e-3)) { fails++; }
    }

    printf("\n%s\n", fails == 0 ? "全部通过" : "有失败项");
    return fails;
}
