#include "Config.h"
#include "Tests.h"

#include <cstdio>
#include <fstream>
#include <iterator>
#include <unistd.h>

// 直接调用配置层，无需初始化 GPU；临时文件只承载解析输入与归档往返。
int run_config_check()
{
    char path[] = "/tmp/fpd_config_XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) { std::perror("mkstemp"); return 1; }
    close(fd);
    const std::string base =
        "init_file=test.fpd\ndt=0.002\nn_steps=200\n"
        "kT=0.25\nradius=3.2\nxi=1\nratio_eta=50\n";
    std::string err;
    int fails = 0;
    auto check = [&](bool ok, const char* name)
    {
        std::printf("  %s  %s\n", name, ok ? "PASS" : "FAIL");
        if (!ok) { fails++; }
    };
    auto load = [&](const std::string& text, FpdConfig& c)
    {
        std::ofstream out(path);
        out << text;
        out.close();
        c = FpdConfig();
        err.clear();
        return load_config(path, c, err);
    };
    auto contents = [&]()
    {
        std::ifstream in(path);
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };

    FpdConfig c;
    check(load(base.substr(base.find('\n') + 1), c) && !validate_config(c, err) &&
          err.find("init_file") != std::string::npos,
          "init_file 缺失报错");
    check(load(base + "dt=0.003\n", c) && c.dt == 0.003 &&
          apply_override(c, "dt=0.004", err) && c.dt == 0.004,
          "重复键最后生效，--set 最后覆盖");
    check(!load(base + "radius_typo=3\n", c) && err.find("未知") != std::string::npos &&
          err.find("你是不是") == std::string::npos, "未知键报错且无拼写建议");
    check(!apply_override(c, "radius_typo=3", err) && err.find("你是不是") == std::string::npos,
          "--set 未知键无拼写建议");
    check(!load(base + "dt=0.002abc\n", c), "非法数值字符串报错");
    check(!load(base + "dt 0.002\n", c) && !apply_override(c, "dt", err),
          "配置和 --set 缺少等号报错");
    check(load(base + "potential=unknown\n", c) && !validate_config(c, err), "非法势名称报错");
    check(load(base + "wall_shift=unknown\n", c) && !validate_config(c, err), "非法移位名称报错");
    check(load(base + "pot_eps=10\nwall_alpha=1\n", c) && validate_config(c, err),
          "未使用的已知势参数允许存在");
    // 只测试配置层接受，不把不适合模拟的参数送入求解器。
    check(load(base + "dt=-1\nratio_eta=0\n", c) && validate_config(c, err),
          "数值范围由使用者负责");

    for (const char* key : {"Nx", "Ny", "Nz"})
    {
        const std::string kv = std::string(key) + "=32";
        check(load(base + kv + "\n", c) && !validate_config(c, err),
              (kv + " 在配置文件中被拒绝").c_str());
        check(load(base, c) && apply_override(c, kv.c_str(), err) &&
              !validate_config(c, err), (kv + " 在 --set 中被拒绝").c_str());
    }

    const char* models[] = {"none", "wca", "morse", "lj126"};
    const char* params[][4] = {
        {0, 0, 0, 0}, {"eps=1", "sigma=2", 0, 0},
        {"De=1", "alpha=1", "r_eq=2", "rcut=5"}, {"eps=1", "sigma=2", "rcut=5", 0}
    };
    for (int family = 0; family < 2; family++)
    {
        const std::string key = family ? "wallpotential" : "potential";
        const std::string prefix = family ? "wall_" : "pot_";
        for (int model = 0; model < 4; model++)
        {
            const std::string label = key + "=" + models[model];
            bool ok = load(base + label + "\n", c);
            // 每补一个参数都再次验证缺项，防止只检查势参数中的某一项。
            for (int j = 0; j < 4 && params[model][j]; j++)
            {
                ok = ok && !validate_config(c, err);
                ok = ok && apply_override(c, (prefix + params[model][j]).c_str(), err);
            }
            ok = ok && validate_config(c, err);
            check(ok, (label + " 必填参数可由 --set 补齐").c_str());

            // 归档全部字段后再次读入，两个归档必须完全一致，含实际派生参数。
            ok = ok && dump_config(c, path, err);
            const std::string first = contents();
            FpdConfig restored;
            ok = ok && load_config(path, restored, err) && validate_config(restored, err);
            ok = ok && dump_config(restored, path, err) && first == contents();
            ok = ok && first.find("wall_alpha = ") != std::string::npos &&
                 first.find("pot_De = ") != std::string::npos &&
                 first.find("h_rest") == std::string::npos && first.find("z_rest") == std::string::npos;
            check(ok, (label + " 全字段归档往返一致且无平衡位置预测").c_str());
        }
    }
    std::remove(path);
    std::printf("配置回归：%s\n", fails ? "失败" : "全部通过");
    return fails;
}
