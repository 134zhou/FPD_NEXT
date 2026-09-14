// fpd_tool —— 纯 CPU 的 .fpd 检查工具，不需要 GPU
//
//   fpd_tool --dump-ckpt <file> [--at i j k]
//   fpd_tool --diff-ckpt <a> <b>
//   fpd_tool --verify-forces <ckpt> <config.used>
//   fpd_tool --sed-stats <f1> <f2> ... [--window A B] [--bins K]
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include "IOBin.h"
#include "Config.h"
#include "Potential.h"
#include "Analysis.h"     // 分块平均与平台判据：沉降统计与 fpd_check 【共用同一份】

static void usage()
{
    printf("用法:\n");
    printf("  fpd_tool --dump-ckpt <file> [--at i j k]   打印头部；--at 打印指定格点的 v\n");
    printf("  fpd_tool --diff-ckpt <a> <b>              逐位比较两个检查点\n");
    printf("  fpd_tool --verify-forces <ckpt> <config>  用配置重算粒子力，与文件里的 F 对照\n");
    printf("  fpd_tool --sed-stats <f1> <f2> ...        沉降统计：<Vz>(t)、质心、最小壁面间隙、φ(z)\n");
    printf("      [--window A B]  只统计 step ∈ [A,B] 的检查点（缺省用后一半）\n");
    printf("      [--bins K]      φ(z) 的分箱数（缺省 16）\n");
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

// ---------------------------------------------------------------------------
// --sed-stats：从一串 .fpd 检查点重建沉降统计（纯 CPU，不碰 GPU）
//
// 只用 .fpd 里的【粒子块】（R/Ru/V），流体场一个字节都不看 —— 沉降的平均速度、
// 质心高度、壁面间隙、浓度剖面全都只需要粒子坐标。
//
// ⚠️ 分块平均【复用 tests/Analysis.cpp】。本工具不另写一份 —— 把相关样本当独立
//    样本会让误差棒低估约 7.5 倍，而「两处独立定义同一个统计量」是个没有判据能
//    发现的漂移点。找不到平台时按设计拒绝背书（pick_plateau 返回 -1）。
//
// ⚠️ φ(z) 是【按粒子中心的计数估计】，不是相场在格子上的求和。两者在有界面宽度的
//    体系里不等（中心法与 diffusion 界面的差别见 make_init.py 的 24% 提示）。
//    要精确的局域体积分数得把 φ 在网格上求和，那是另一件事（需要 config.used）。
// ---------------------------------------------------------------------------
struct SedRow
{
    long   step;
    double t;
    double Vz_mean, Vz_std, Vz_min, Vz_max;
    double z_cm, z_min, z_max, gap_min;
};

// ⚠️ zs 里累积的是【窗口内全部检查点】的粒子 z，所以除以 bin 体积【还必须再除以
//    帧数】。漏掉这一步会得到 φ > 1 这种一眼假但没人会去核的数（帧越多越离谱）。
static void print_sed_profile(const std::vector<double>& zs, int n_frames, int K,
                              double zlo, double zhi, double radius, double xi, int NxNy)
{
    if (n_frames <= 0) { return; }
    const double dz = (zhi - zlo) / (double)K;
    const double sphere = 4.0 / 3.0 * M_PI * radius * radius * radius;
    const double diffuse = sphere + 8.0 * M_PI * radius * xi * xi * M_PI * M_PI / 24.0;
    const double bin_vol = (double)NxNy * dz * (double)n_frames;

    std::vector<double> cnt(K, 0.0);
    for (size_t n = 0; n < zs.size(); n++)
    {
        int b = (int)std::floor((zs[n] - zlo) / dz);
        if (b < 0) { b = 0; }
        if (b >= K) { b = K - 1; }
        cnt[b] += 1.0;
    }
    printf("\n  φ(z) 剖面（按粒子中心计数，%d 帧平均；两个约定都给）\n", n_frames);
    printf("    %-18s %8s %10s %10s\n", "z 区间", "每帧计数", "φ(解析球)", "φ(扩散)");
    for (int b = 0; b < K; b++)
    {
        if (cnt[b] == 0.0) { continue; }
        const double z1 = zlo + b * dz, z2 = z1 + dz;
        printf("    [%7.2f, %7.2f) %8.2f %10.4f %10.4f\n",
               z1, z2, cnt[b] / (double)n_frames,
               cnt[b] * sphere / bin_vol, cnt[b] * diffuse / bin_vol);
    }
}

static int cmd_sed_stats(int argc, char** argv)
{
    std::vector<const char*> files;
    long win_lo = -1, win_hi = -1;
    int  K = 16;
    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--window") == 0 && i + 2 < argc)
        { win_lo = atol(argv[++i]); win_hi = atol(argv[++i]); }
        else if (strcmp(argv[i], "--bins") == 0 && i + 1 < argc) { K = atoi(argv[++i]); }
        else if (argv[i][0] == '-') { fprintf(stderr, "错误: 未知选项 %s\n", argv[i]); return 2; }
        else { files.push_back(argv[i]); }
    }
    if (files.empty()) { fprintf(stderr, "错误: --sed-stats 需要至少一个 .fpd\n"); return 2; }
    if (K < 1) { K = 1; }

    std::vector<SedRow> rows;
    std::vector<double> vz_series, step_series;
    std::vector<double> z_all;          // 统计窗口内累积的粒子 z
    int N_used = 0, NxNy = 0;
    double radius = 0.0, xi = 0.0, dt = 0.0, Nz = 0.0;

    for (size_t f = 0; f < files.size(); f++)
    {
        CkptHeader h;
        std::string err;
        if (!read_ckpt_header(files[f], h, err))
        { fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }

        Holder H;
        H.alloc(h);
        if (!load_checkpoint(files[f], h, H.a, err))
        { fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }

        const int N = h.N;
        if (N <= 0) { continue; }
        radius = h.radius; xi = h.xi; dt = h.dt; Nz = (double)h.Nz;
        N_used = N; NxNy = h.Nx * h.Ny;

        SedRow r;
        r.step = (long)h.step;
        r.t = (double)h.step * h.dt;
        double sz = 0.0, sv = 0.0, sv2 = 0.0;
        r.z_min = 1e300; r.z_max = -1e300;
        r.Vz_min = 1e300; r.Vz_max = -1e300;
        r.gap_min = 1e300;
        for (int n = 0; n < N; n++)
        {
            sz += H.Rz[n]; sv += H.Vz[n]; sv2 += H.Vz[n] * H.Vz[n];
            r.z_min = std::fmin(r.z_min, H.Rz[n]);
            r.z_max = std::fmax(r.z_max, H.Rz[n]);
            r.Vz_min = std::fmin(r.Vz_min, H.Vz[n]);
            r.Vz_max = std::fmax(r.Vz_max, H.Vz[n]);
            // 壁面间隙 h = (z ± 1/2) - a —— 与 Potential.h 的 wall_gap_* 同一约定。
            // 这是「壁面势有没有托住粒子」的【唯一在线体检】。
            const double hb = (H.Rz[n] + 0.5) - h.radius;
            const double ht = ((double)h.Nz - 0.5 - H.Rz[n]) - h.radius;
            r.gap_min = std::fmin(r.gap_min, std::fmin(hb, ht));
        }
        r.z_cm = sz / N;
        r.Vz_mean = sv / N;
        r.Vz_std = N > 1 ? std::sqrt(std::fmax((sv2 - N * r.Vz_mean * r.Vz_mean) / (N - 1), 0.0)) : 0.0;
        rows.push_back(r);
    }

    if (rows.empty()) { fprintf(stderr, "错误: 没有可用的检查点（N=0？）\n"); return 2; }

    // ⚠️ 按 step 排序，【不信任参数顺序】。shell glob 会把 init 文件（step=0）也带
    //    进来，且排在最后；不排序的话「后一半窗口」会变成「init 那一个样本」，
    //    而时间序列还会出现 step 回跳 —— 两种都是静默错统计。
    for (size_t i = 1; i < rows.size(); i++)
    {
        SedRow key = rows[i];
        size_t j = i;
        while (j > 0 && rows[j - 1].step > key.step) { rows[j] = rows[j - 1]; j--; }
        rows[j] = key;
    }
    for (size_t i = 1; i < rows.size(); i++)
    {
        if (rows[i].step == rows[i - 1].step)
        {
            fprintf(stderr,
                    "错误: 有两个文件的 step 都是 %ld。最可能的原因是 glob 把初始构型\n"
                    "      （如 out/sed_init.fpd，step=0）和同名检查点（out/sed_0000000.fpd）\n"
                    "      一起带了进来。请收紧通配，例如 out/sed_[0-9]*.fpd\n", rows[i].step);
            return 2;
        }
    }

    // 缺省窗口：后一半（【排序之后】才定，且必须落在 rows 的 step 上，
    // 否则会出现 win_hi < win_lo 这种「窗口里一个样本都没有」的静默空统计）
    if (win_lo < 0)
    {
        win_lo = rows[rows.size() / 2].step;
        win_hi = rows.back().step;
    }

    // 时间序列在【排序 + 定窗口之后】统一构建，保证顺序与窗口一致
    for (size_t i = 0; i < rows.size(); i++)
    {
        if (rows[i].step >= win_lo && rows[i].step <= win_hi)
        { vz_series.push_back(rows[i].Vz_mean); step_series.push_back(rows[i].t); }
    }

    printf("=== 沉降统计（%zu 个检查点，N=%d）===\n", rows.size(), N_used);
    printf("  %9s %8s %11s %11s %11s %9s %9s %9s\n",
           "step", "t", "<Vz>", "Vz_std", "z_cm", "z_min", "z_max", "min_gap");
    for (size_t i = 0; i < rows.size(); i++)
    {
        const SedRow& r = rows[i];
        printf("  %9ld %8.1f %11.5f %11.5f %11.4f %9.4f %9.4f %9.5f%s\n",
               r.step, r.t, r.Vz_mean, r.Vz_std, r.z_cm, r.z_min, r.z_max, r.gap_min,
               (r.step >= win_lo && r.step <= win_hi) ? "  *" : "");
    }
    printf("  （* = 参与统计窗口 step ∈ [%ld, %ld]）\n", win_lo, win_hi);

    // 采样间隔必须按【物理时间】记账，否则 dt 越小样本越相关（CLAUDE.md 数值定则）。
    // ⚠️ 从【全部检查点】的间距取，不从窗口序列取：窗口里只剩 1 个样本时会退化成
    //    返回 dt，而真实间隔是 interval_ckpt*dt —— 那会把 tau_int 报小几个量级。
    double samp_dt = 0.0;
    if (rows.size() >= 2)
    {
        samp_dt = rows[1].t - rows[0].t;
        for (size_t i = 2; i < rows.size(); i++)
        {
            const double d = rows[i].t - rows[i - 1].t;
            if (std::fabs(d - samp_dt) > 1e-9 * std::fmax(samp_dt, 1.0))
            { printf("  [注意] 检查点间隔不均匀（%.4f vs %.4f）——分块平均仍可用，"
                     "但 tau_int 的解释要小心\n", samp_dt, d); break; }
        }
    }
    printf("\n  采样间隔 = %.4f 物理时间单位（%zu 个样本，步长 %g）\n",
           samp_dt, vz_series.size(), dt);

    if (vz_series.size() < 20)
    {
        printf("  *** 样本数 %zu < 20，拒绝背书（n_b 不足时块方差本身相对误差 >16%%）***\n",
               vz_series.size());
    }
    else
    {
        BlockStat bs[64];
        const int ns = blocking_analysis(vz_series.data(), (int)vz_series.size(), bs, 64);
        const BlockStat pk = pick_plateau(bs, ns, 20);
        print_blocking(bs, ns, pk);
        if (pk.block_len < 0)
        {
            printf("  *** 未找到分块平台 → 误差棒不可信，必须延长模拟（按设计拒绝背书）***\n");
        }
        else
        {
            printf("  <Vz>_plateau = %.6f ± %.6f   (b=%d, n_b=%d, tau_int=%.1f 采样间隔)\n",
                   pk.mean, pk.stderr_mean, pk.block_len, pk.n_block, pk.tau_int_samp);
            printf("  即 tau_int = %.2f 物理时间单位；b >= 10*tau_int 的硬约束 %s\n",
                   pk.tau_int_samp * samp_dt,
                   (pk.block_len >= 10.0 * pk.tau_int_samp) ? "满足" : "*** 不满足 ***");
        }
    }

    // 全窗口的粒子 z 收集（用于 φ(z)）：重读一遍太浪费，这里只对窗口内的文件累加
    // ⚠️ 位置分布是【稳态量】，必须只统计平台段；混入瞬态会让剖面失真。
    int z_frames = 0;
    for (size_t f = 0; f < files.size(); f++)
    {
        CkptHeader h;
        std::string err;
        if (!read_ckpt_header(files[f], h, err)) { continue; }
        if ((long)h.step < win_lo || (long)h.step > win_hi) { continue; }
        Holder H;
        H.alloc(h);
        if (!load_checkpoint(files[f], h, H.a, err)) { continue; }
        for (int n = 0; n < h.N; n++) { z_all.push_back(H.Rz[n]); }
        z_frames++;
    }
    if (!z_all.empty())
    {
        print_sed_profile(z_all, z_frames, K, -0.5, Nz - 0.5, radius, xi, NxNy);
    }

    // 间隙体检：壁面势托没托住
    double gmin = 1e300;
    for (size_t i = 0; i < rows.size(); i++) { gmin = std::fmin(gmin, rows[i].gap_min); }
    printf("\n  全程最小壁面间隙 = %.5f", gmin);
    if (gmin <= 0.0)      { printf("   *** 粒子已侵入壁面（壁面势失效）***\n"); }
    else if (gmin < 0.5)  { printf("   [注意] 几乎贴壁，壁面势处于强排斥区\n"); }
    else                  { printf("   壁面势托住了\n"); }

    return 0;
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
    if (strcmp(argv[1], "--sed-stats") == 0) { return cmd_sed_stats(argc, argv); }
    usage();
    return 2;
}
