// fpd_tool —— 纯 CPU 的 .fpd 检查工具，不需要 GPU
//
//   fpd_tool --dump-ckpt <file> [--at i j k]
//   fpd_tool --diff-ckpt <a> <b>
//   fpd_tool --verify-forces <ckpt> <config.used>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include "IOBin.h"
#include "Config.h"
#include "Potential.h"

static void usage()
{
    printf("用法:\n");
    printf("  fpd_tool --dump-ckpt <file> [--at i j k]   打印头部；--at 打印指定格点的 v\n");
    printf("  fpd_tool --diff-ckpt <a> <b>              逐位比较两个检查点\n");
    printf("  fpd_tool --verify-forces <ckpt> <config>  用配置重算粒子力，与文件里的 F 对照\n");
}

// 按 header 分配好数组的持有者
struct Holder
{
    std::vector<double> vx, vy, vz, p;
    std::vector<double> Rx, Ry, Rz, Rux, Ruy, Ruz, Vx, Vy, Vz, Fx, Fy, Fz;
    CkptArrays a;

    void alloc(const CkptHeader& h)
    {
        const size_t size = h.size();
        const size_t N = (size_t)h.N;
        vx.assign(size, 0); vy.assign(size, 0); vz.assign(size, 0);
        if (h.has_pressure()) { p.assign(size, 0); }
        Rx.assign(N,0); Ry.assign(N,0); Rz.assign(N,0);
        Rux.assign(N,0); Ruy.assign(N,0); Ruz.assign(N,0);
        Vx.assign(N,0); Vy.assign(N,0); Vz.assign(N,0);
        Fx.assign(N,0); Fy.assign(N,0); Fz.assign(N,0);

        a.vx = vx.data(); a.vy = vy.data(); a.vz = vz.data();
        a.p  = h.has_pressure() ? p.data() : 0;
        a.Rx = Rx.data(); a.Ry = Ry.data(); a.Rz = Rz.data();
        a.Rux = Rux.data(); a.Ruy = Ruy.data(); a.Ruz = Ruz.data();
        a.Vx = Vx.data(); a.Vy = Vy.data(); a.Vz = Vz.data();
        a.Fx = Fx.data(); a.Fy = Fy.data(); a.Fz = Fz.data();
    }
};

static int cmd_dump(int argc, char** argv)
{
    const char* path = argv[2];
    CkptHeader h;
    std::string err;

    if (!read_ckpt_header(path, h, err)) { fprintf(stderr, "错误: %s\n", err.c_str()); return 1; }

    printf("=== %s ===\n", path);
    printf("  version      %u\n", h.version);
    printf("  网格         %d x %d x %d   (size = %zu)\n", h.Nx, h.Ny, h.Nz, h.size());
    printf("  粒子数 N     %d\n", h.N);
    printf("  step         %lld\n", (long long)h.step);
    printf("  dt           %.17g\n", h.dt);
    printf("  kT           %.17g\n", h.kT);
    printf("  radius       %.17g\n", h.radius);
    printf("  xi           %.17g\n", h.xi);
    printf("  ratio_eta    %.17g\n", h.ratio_eta);
    printf("  noise_on     %d\n", h.noise_on);
    printf("  seed         %llu\n", (unsigned long long)h.seed);
    printf("  rng_draws    %llu\n", (unsigned long long)h.rng_draws);
    printf("  has_pressure %s\n", h.has_pressure() ? "是" : "否");

    // 派生量：由 Common.h 的工厂算出，不从文件读
    NS_Config ns = make_ns_config(h.Nx, h.Ny, h.Nz, h.dt, h.kT, h.noise_on != 0);
    PhiParams pp = make_phi_params(h.radius, h.xi, h.ratio_eta);
    printf("  --- 派生量（由 make_ns_config/make_phi_params 重算）---\n");
    printf("  W            %.17g\n", ns.W);
    printf("  inv_dt       %.17g\n", ns.inv_dt);
    printf("  n_range      %d\n", pp.n_range);

    Holder H; H.alloc(h);
    CkptHeader h2;
    if (!load_checkpoint(path, h2, H.a, err))
    { fprintf(stderr, "错误: %s\n", err.c_str()); return 1; }
    printf("  校验和       通过\n");

    // --at i j k：判据 3 用，验证 C++/Python 的轴序约定一致
    for (int i = 3; i + 3 < argc; i++)
    {
        if (strcmp(argv[i], "--at") == 0)
        {
            const int ii = atoi(argv[i+1]), jj = atoi(argv[i+2]), kk = atoi(argv[i+3]);
            if (ii < 0 || ii >= h.Nx || jj < 0 || jj >= h.Ny || kk < 0 || kk >= h.Nz)
            { fprintf(stderr, "错误: --at 下标越界\n"); return 1; }
            const int Nx = h.Nx, Ny = h.Ny;
            const size_t ijk = IDX(ii, jj, kk);
            printf("  --- v at (i=%d, j=%d, k=%d)  IDX=%zu ---\n", ii, jj, kk, ijk);
            printf("  vx = %.17g\n", H.vx[ijk]);
            printf("  vy = %.17g\n", H.vy[ijk]);
            printf("  vz = %.17g\n", H.vz[ijk]);
        }
    }

    if (h.N > 0)
    {
        printf("  --- 粒子 0 ---\n");
        printf("  R  = (%.17g, %.17g, %.17g)\n", H.Rx[0], H.Ry[0], H.Rz[0]);
        printf("  Ru = (%.17g, %.17g, %.17g)\n", H.Rux[0], H.Ruy[0], H.Ruz[0]);
        printf("  V  = (%.17g, %.17g, %.17g)\n", H.Vx[0], H.Vy[0], H.Vz[0]);
    }
    return 0;
}

// 逐位比较；返回 0 表示完全相同
static int cmd_diff(const char* pa, const char* pb)
{
    CkptHeader ha, hb;
    std::string err;
    if (!read_ckpt_header(pa, ha, err)) { fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }
    if (!read_ckpt_header(pb, hb, err)) { fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }

    if (ha.Nx != hb.Nx || ha.Ny != hb.Ny || ha.Nz != hb.Nz || ha.N != hb.N)
    { printf("维度不同\n"); return 1; }

    Holder A, B;
    A.alloc(ha); B.alloc(hb);
    if (!load_checkpoint(pa, ha, A.a, err)) { fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }
    if (!load_checkpoint(pb, hb, B.a, err)) { fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }

    if (ha.step != hb.step)
    { printf("[注意] step 不同: %lld vs %lld\n", (long long)ha.step, (long long)hb.step); }

    struct Item { const char* name; const double* x; const double* y; size_t n; };
    const size_t size = ha.size(), N = (size_t)ha.N;
    const Item items[] = {
        {"vx", A.vx.data(), B.vx.data(), size},
        {"vy", A.vy.data(), B.vy.data(), size},
        {"vz", A.vz.data(), B.vz.data(), size},
        {"Rx", A.Rx.data(), B.Rx.data(), N},
        {"Ry", A.Ry.data(), B.Ry.data(), N},
        {"Rz", A.Rz.data(), B.Rz.data(), N},
        {"Rux", A.Rux.data(), B.Rux.data(), N},
        {"Ruy", A.Ruy.data(), B.Ruy.data(), N},
        {"Ruz", A.Ruz.data(), B.Ruz.data(), N},
        {"Vx", A.Vx.data(), B.Vx.data(), N},
        {"Vy", A.Vy.data(), B.Vy.data(), N},
        {"Vz", A.Vz.data(), B.Vz.data(), N},
        // Fx/Fy/Fz 曾经漏掉。加粒子间力/壁面势之后，力是逐位重启对照里最该看的量
        // —— 它是 R 的【原因】，R 逐位相同但 F 不同说明「同一构型下力算得不一样」，
        // 那比轨迹漂移更早、更清楚地暴露问题。
        {"Fx", A.Fx.data(), B.Fx.data(), N},
        {"Fy", A.Fy.data(), B.Fy.data(), N},
        {"Fz", A.Fz.data(), B.Fz.data(), N},
    };

    int nbad = 0;
    printf("=== diff %s  vs  %s ===\n", pa, pb);
    for (size_t t = 0; t < sizeof(items)/sizeof(items[0]); t++)
    {
        const Item& it = items[t];
        size_t ndiff = 0; double maxabs = 0.0, maxrel = 0.0;
        for (size_t i = 0; i < it.n; i++)
        {
            // 逐位比较：用 memcmp 而非 ==，这样 NaN 和 ±0 也能区分
            if (memcmp(&it.x[i], &it.y[i], sizeof(double)) != 0)
            {
                ndiff++;
                const double d = fabs(it.x[i] - it.y[i]);
                const double s = fabs(it.x[i]) > fabs(it.y[i]) ? fabs(it.x[i]) : fabs(it.y[i]);
                if (d > maxabs) { maxabs = d; }
                if (s > 0 && d / s > maxrel) { maxrel = d / s; }
            }
        }
        if (ndiff == 0)
        {
            printf("  %-4s 逐位相同 (%zu 个)\n", it.name, it.n);
        }
        else
        {
            nbad++;
            printf("  %-4s 有 %zu/%zu 个不同   最大绝对差 %.3e   最大相对差 %.3e\n",
                   it.name, ndiff, it.n, maxabs, maxrel);
        }
    }
    printf("%s\n", nbad == 0 ? "全部逐位相同" : "存在差异");
    return nbad == 0 ? 0 : 1;
}

// J7：检查点自洽性。读 .fpd 里的 R，用配置里的势参数重算 F，与文件里的 F 对照。
// 势参数不进 .fpd header（见 PROGRESS.md 决策记录），这条判据把「不自描述」的代价
// 从静默变成可检测 —— 配错势参数会在这里 FAIL。
static int cmd_verify_forces(const char* ckpt_path, const char* cfg_path)
{
    std::string err;

    CkptHeader h;
    if (!read_ckpt_header(ckpt_path, h, err))
    { fprintf(stderr, "错误: %s\n", err.c_str()); return 1; }

    Holder holder;
    holder.alloc(h);
    if (!load_checkpoint(ckpt_path, h, holder.a, err))
    { fprintf(stderr, "错误: %s\n", err.c_str()); return 1; }

    FpdConfig fc;
    if (!load_config(cfg_path, fc, err))
    { fprintf(stderr, "错误: %s\n", err.c_str()); return 1; }

    NS_Config cfg = make_ns_config(fc);
    PotentialParams pot = make_potential_params(fc);
    WallParams wall = make_wall_params(fc);
    ExternalField ext = make_external_field(fc);

    std::vector<double> Fx(h.N), Fy(h.N), Fz(h.N);
    compute_particle_forces_cpu(cfg, pot, ext, wall, h.N,
                                holder.Rx.data(), holder.Ry.data(), holder.Rz.data(),
                                Fx.data(), Fy.data(), Fz.data());

    // 对比（容差：文件里的 F 是 GPU 全矩阵算的，这里是 CPU 串行 i<j，浮点顺序不同，
    // 差 ~1e-15 相对；配错参数会差 O(1)）
    double maxabs = 0.0, fmag = 0.0;
    for (int n = 0; n < h.N; n++)
    {
        double d1 = fabs(Fx[n] - holder.Fx[n]);
        double d2 = fabs(Fy[n] - holder.Fy[n]);
        double d3 = fabs(Fz[n] - holder.Fz[n]);
        maxabs = fmax(maxabs, fmax(d1, fmax(d2, d3)));
        fmag = fmax(fmag, fmax(fabs(holder.Fx[n]), fmax(fabs(holder.Fy[n]), fabs(holder.Fz[n]))));
    }
    const double tol = 1e-10 * (fmag > 0.0 ? fmag : 1.0);
    const bool ok = maxabs < tol;

    printf("=== verify-forces %s ===\n", ckpt_path);
    printf("  势 = %s   壁面势 = %s   外场 = (%g, %g, %g)   N = %d\n",
           fc.potential.c_str(), fc.wall_pot.c_str(), ext.gx, ext.gy, ext.gz, h.N);
    printf("  重算 F vs 文件 F: max |ΔF| = %.3e   阈值 %.3e   %s\n",
           maxabs, tol, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

int main(int argc, char** argv)
{
    if (argc < 3) { usage(); return 2; }

    if (strcmp(argv[1], "--dump-ckpt") == 0) { return cmd_dump(argc, argv); }
    if (strcmp(argv[1], "--diff-ckpt") == 0)
    {
        if (argc < 4) { usage(); return 2; }
        return cmd_diff(argv[2], argv[3]);
    }
    if (strcmp(argv[1], "--verify-forces") == 0)
    {
        if (argc < 4) { usage(); return 2; }
        return cmd_verify_forces(argv[2], argv[3]);
    }
    usage();
    return 2;
}
