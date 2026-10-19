// fpd_calib —— 变密度压力求解的容差标定（独立程序，跑 fpd 之前手动运行）
//
//   fpd_calib <配置文件> [--set key=value ...] [--state ckpt.fpd]
//             [--warm N] [--gap G] [--meas M] [--target 1e-3] [--margin 10]
//
// 读与 fpd 完全相同的配置与 .fpd 状态，用生产路径推进轨迹，在若干测量点上量
//   ε_V(τ) = max_n |残差导致的粒子速度误差| / rms_n |同一步热噪声导致的粒子速度增量|
// （以及 ρ 加权能量比 ε_E），给出满足 ε <= target 的最松容差 τ*，推荐 τ*/margin。
// 结果只打印，不改配置、不写场数据。核心算法见 PcgCalib.h。
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>

#include "Config.h"
#include "IOBin.h"
#include "State.h"
#include "Step.h"
#include "PcgCalib.h"

static void usage()
{
    std::printf("用法: fpd_calib <配置文件> [--set key=value ...] [--state ckpt.fpd]\n"
                "                [--warm N] [--gap G] [--meas M] [--target 1e-3] [--margin 10]\n"
                "  --state   起始状态（缺省用配置的 init_file）；用运行中的检查点最能代表生产\n"
                "  --warm    首个测量点前推进的步数（缺省 200；从已平衡的检查点起可设 0）\n"
                "  --gap     测量点间隔步数（缺省 50）\n"
                "  --meas    测量点个数（缺省 4）\n"
                "  --target  ε 门槛（缺省 1e-3）\n"
                "  --margin  推荐值 = 最松合格档 / margin（缺省 10）\n");
}

int main(int argc, char** argv)
{
    if (argc < 2) { usage(); return 2; }

    FpdConfig fc;
    std::string err;
    if (!load_config(argv[1], fc, err)) { std::fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }

    CalibOpts o;
    std::string state;
    for (int i = 2; i < argc; i++)
    {
        const bool more = (i + 1 < argc);
        if (std::strcmp(argv[i], "--set") == 0 && more)
        {
            if (!apply_override(fc, argv[++i], err)) { std::fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }
        }
        else if (std::strcmp(argv[i], "--state") == 0 && more)  { state = argv[++i]; }
        else if (std::strcmp(argv[i], "--warm") == 0 && more)   { o.warm = std::atol(argv[++i]); }
        else if (std::strcmp(argv[i], "--gap") == 0 && more)    { o.gap = std::atol(argv[++i]); }
        else if (std::strcmp(argv[i], "--meas") == 0 && more)   { o.meas = std::atoi(argv[++i]); }
        else if (std::strcmp(argv[i], "--target") == 0 && more) { o.target = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--margin") == 0 && more) { o.margin = std::atof(argv[++i]); }
        else { std::fprintf(stderr, "错误: 无法识别的参数 '%s'\n", argv[i]); usage(); return 2; }
    }
    if (!state.empty()) { fc.init_file = state; }
    if (!validate_config(fc, err)) { std::fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }
    if (o.warm < 0 || o.gap < 1 || o.meas < 1 || o.target <= 0.0 || o.margin < 1.0)
    { std::fprintf(stderr, "错误: --warm>=0、--gap>=1、--meas>=1、--target>0、--margin>=1\n"); return 2; }

    if (fc.ratio_rho == 1.0)
    {
        std::printf("ratio_rho = 1：压力走常系数直解（机器精度），不使用 pcg_tol，无需标定。\n");
        return 0;
    }
    if (fc.noise_on == 0 || fc.kT == 0.0)
    {
        std::printf("无热噪声（noise_on=0 或 kT=0）：「误差比热噪声小」的判据不适用，\n"
                    "确定性算例请按所测量的精度保留紧容差（默认 1e-12）。\n");
        return 0;
    }

    CkptHeader h;
    CkptReader reader;
    if (!open_checkpoint(fc.init_file.c_str(), h, reader, err))
    { std::fprintf(stderr, "错误: %s\n", err.c_str()); return 2; }
    if (h.N < 1)
    {
        close_checkpoint(reader);
        std::fprintf(stderr, "错误: %s 没有粒子，无法定义粒子速度判据\n", fc.init_file.c_str());
        return 2;
    }

    FpdConfig effective = fc;
    effective.Nx = h.Nx; effective.Ny = h.Ny; effective.Nz = h.Nz;
    const StepParams sp = make_step_params(effective);

    FpdState st;
    st.init(sp.cfg, h.N, ST_FULL, fc.seed);
    CkptArrays arr = st.ckpt_arrays();
    if (!read_checkpoint(reader, h, arr, err))
    { std::fprintf(stderr, "错误: %s\n", err.c_str()); st.finish(); return 2; }
    st.upload(ST_PARTICLE | ST_VELOCITY);

    std::printf("=== fpd_calib：PCG 容差标定 ===\n");
    std::printf("  状态 %s（step=%lld）  %dx%dx%d  N=%d\n", fc.init_file.c_str(), (long long)h.step,
                h.Nx, h.Ny, h.Nz, h.N);
    std::printf("  ratio_rho=%g  ratio_eta=%g  kT=%g  dt=%g  g=(%g, %g, %g)  轨迹 pcg_tol=%.0e\n",
                fc.ratio_rho, fc.ratio_eta, fc.kT, fc.dt, sp.ext.gx, sp.ext.gy, sp.ext.gz, sp.cfg.pcg_tol);
    std::printf("  预热 %ld 步，之后每 %ld 步测一次，共 %d 次；门槛 %.0e，余量 %.0f 倍\n",
                o.warm, o.gap, o.meas, o.target, o.margin);
    std::printf("  分母：同一步热噪声导致的粒子速度增量；分子：PCG 残差导致的粒子速度误差\n\n");
    std::fflush(stdout);

    std::vector<CalibPoint> pts;
    const int rc = calib_run(sp, st, (long)h.step, o, pts);
    st.finish();
    if (rc != 0) { return rc; }

    double tau_star = 0.0;
    const double rec = calib_report(pts, o, sp.cfg.pcg_tol, tau_star);
    if (rec <= 0.0) { return 1; }
    std::printf("\n建议写入配置：\npcg_tol = %.0e\n", rec);
    return 0;
}
