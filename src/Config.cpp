#include "./include/Config.h"

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <fstream>
#include <sstream>

// 三者共用同一张表：解析、--set 覆盖、dump。
std::vector<FieldDesc> config_fields(FpdConfig& c)
{
    return {
        {"Nx",             F_INT,    &c.Nx,             true,  "网格 x 方向格子数"},
        {"Ny",             F_INT,    &c.Ny,             true,  "网格 y 方向格子数"},
        {"Nz",             F_INT,    &c.Nz,             true,  "网格 z 方向格子数"},
        {"dt",             F_DOUBLE, &c.dt,             true,  "时间步长"},
        {"n_steps",        F_LONG,   &c.n_steps,        true,  "绝对目标步数（不是增量）"},
        {"kT",             F_DOUBLE, &c.kT,             true,  "温度 k_B T"},
        {"radius",         F_DOUBLE, &c.radius,         true,  "粒子半径 a"},
        {"xi",             F_DOUBLE, &c.xi,             true,  "界面宽度 xi"},
        {"ratio_eta",      F_DOUBLE, &c.ratio_eta,      true,  "eta_c / eta_l"},

        {"gravity_z",      F_DOUBLE, &c.gravity_z,      false, "z 向外场力（每粒子）"},
        {"noise_on",       F_INT,    &c.noise_on,       false, "0/1，关掉则 W=0"},

        {"init_file",      F_STRING, &c.init_file,      false, ".fpd 初始构型或续跑文件；留空用内置默认"},
        {"out_dir",        F_STRING, &c.out_dir,        false, "输出目录"},
        {"run_name",       F_STRING, &c.run_name,       false, "输出文件名前缀"},
        {"interval_ckpt",  F_LONG,   &c.interval_ckpt,  false, "每多少步写一个 .fpd"},
        {"interval_log",   F_LONG,   &c.interval_log,   false, "每多少步打印一行进度"},
        {"save_pressure",  F_INT,    &c.save_pressure,  false, "0/1，把 p 存进 .fpd 供可视化"},
        {"seed",           F_ULL,    &c.seed,           false, "随机数种子"},
    };
}

static std::string trim(const std::string& s)
{
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) { return ""; }
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// 严格转换：整个字符串必须被消费掉，否则 "1.5abc" 这类会被静默截断
static bool set_field(const FieldDesc& f, const std::string& val, std::string& err)
{
    const char* s = val.c_str();
    char* end = 0;
    errno = 0;

    switch (f.type)
    {
        case F_STRING:
            *(std::string*)f.ptr = val;
            return true;

        case F_INT:
        {
            long v = std::strtol(s, &end, 10);
            if (end == s || *end != '\0' || errno == ERANGE)
            { err = std::string("键 '") + f.key + "' 的值不是合法整数: '" + val + "'"; return false; }
            *(int*)f.ptr = (int)v;
            return true;
        }
        case F_LONG:
        {
            long v = std::strtol(s, &end, 10);
            if (end == s || *end != '\0' || errno == ERANGE)
            { err = std::string("键 '") + f.key + "' 的值不是合法整数: '" + val + "'"; return false; }
            *(long*)f.ptr = v;
            return true;
        }
        case F_ULL:
        {
            unsigned long long v = std::strtoull(s, &end, 10);
            if (end == s || *end != '\0' || errno == ERANGE)
            { err = std::string("键 '") + f.key + "' 的值不是合法无符号整数: '" + val + "'"; return false; }
            *(unsigned long long*)f.ptr = v;
            return true;
        }
        case F_DOUBLE:
        {
            double v = std::strtod(s, &end);
            if (end == s || *end != '\0' || errno == ERANGE)
            { err = std::string("键 '") + f.key + "' 的值不是合法实数: '" + val + "'"; return false; }
            *(double*)f.ptr = v;
            return true;
        }
    }
    err = "内部错误：未知字段类型";
    return false;
}

// 未知 key 时给出「你是不是想写 X」的提示（前缀匹配 + 长度接近）
static std::string suggest(const std::vector<FieldDesc>& tab, const std::string& key)
{
    std::string best;
    for (size_t i = 0; i < tab.size(); i++)
    {
        const std::string k = tab[i].key;
        size_t n = k.size() < key.size() ? k.size() : key.size();
        size_t common = 0;
        while (common < n && k[common] == key[common]) { common++; }
        if (common >= 3 && common > best.size()) { best = k; }
    }
    return best.empty() ? "" : ("  你是不是想写 '" + best + "'？");
}

bool load_config(const char* path, FpdConfig& c, std::string& err)
{
    std::ifstream ifs(path);
    if (!ifs.is_open())
    {
        err = std::string("打不开配置文件: ") + path;
        return false;
    }

    std::vector<FieldDesc> tab = config_fields(c);
    std::vector<bool> seen(tab.size(), false);

    std::string line;
    int lineno = 0;
    while (std::getline(ifs, line))
    {
        lineno++;
        // 去掉 # 之后的内容
        size_t h = line.find('#');
        if (h != std::string::npos) { line = line.substr(0, h); }
        line = trim(line);
        if (line.empty()) { continue; }

        size_t eq = line.find('=');
        if (eq == std::string::npos)
        {
            std::ostringstream o;
            o << path << ":" << lineno << " 不是 'key = value' 形式: '" << line << "'";
            err = o.str();
            return false;
        }

        const std::string key = trim(line.substr(0, eq));
        const std::string val = trim(line.substr(eq + 1));

        size_t hit = tab.size();
        for (size_t i = 0; i < tab.size(); i++)
        {
            if (key == tab[i].key) { hit = i; break; }
        }
        // 未知 key 必须报错退出 —— 静默忽略拼错的键是最常见的浪费一天的事故
        if (hit == tab.size())
        {
            std::ostringstream o;
            o << path << ":" << lineno << " 未知的配置项 '" << key << "'" << suggest(tab, key);
            err = o.str();
            return false;
        }
        if (seen[hit])
        {
            std::ostringstream o;
            o << path << ":" << lineno << " 配置项 '" << key << "' 重复出现";
            err = o.str();
            return false;
        }
        if (!set_field(tab[hit], val, err))
        {
            std::ostringstream o;
            o << path << ":" << lineno << " " << err;
            err = o.str();
            return false;
        }
        seen[hit] = true;
    }

    // 必填项缺失也必须报错 —— 不给默认值兜底
    std::string missing;
    for (size_t i = 0; i < tab.size(); i++)
    {
        if (tab[i].required && !seen[i])
        {
            missing += (missing.empty() ? "" : ", ");
            missing += tab[i].key;
        }
    }
    if (!missing.empty())
    {
        err = std::string("配置文件缺少必填项: ") + missing;
        return false;
    }

    return true;
}

bool apply_override(FpdConfig& c, const char* kv, std::string& err)
{
    const std::string s = kv;
    size_t eq = s.find('=');
    if (eq == std::string::npos)
    {
        err = std::string("--set 需要 key=value 形式，收到: '") + kv + "'";
        return false;
    }
    const std::string key = trim(s.substr(0, eq));
    const std::string val = trim(s.substr(eq + 1));

    std::vector<FieldDesc> tab = config_fields(c);
    for (size_t i = 0; i < tab.size(); i++)
    {
        if (key == tab[i].key) { return set_field(tab[i], val, err); }
    }
    err = std::string("--set 用了未知的配置项 '") + key + "'" + suggest(tab, key);
    return false;
}

bool validate_config(const FpdConfig& c, std::string& err)
{
    if (c.Nx <= 0 || c.Ny <= 0 || c.Nz <= 0)
    { err = "Nx/Ny/Nz 必须为正"; return false; }
    if (c.dt <= 0.0)   { err = "dt 必须为正"; return false; }
    if (c.kT < 0.0)    { err = "kT 不能为负"; return false; }
    if (c.radius <= 0.0) { err = "radius 必须为正"; return false; }
    if (c.xi <= 0.0)   { err = "xi 必须为正"; return false; }
    if (c.ratio_eta < 1.0) { err = "ratio_eta 必须 >= 1（粒子不能比溶剂稀）"; return false; }
    if (c.n_steps < 0) { err = "n_steps 不能为负"; return false; }
    if (c.interval_ckpt <= 0) { err = "interval_ckpt 必须为正"; return false; }
    if (c.interval_log  <= 0) { err = "interval_log 必须为正"; return false; }

    // 相场支撑域不能大到自己绕一圈碰到自己
    PhiParams pp = make_phi_params(c);
    const int nmin = c.Nx < c.Ny ? (c.Nx < c.Nz ? c.Nx : c.Nz) : (c.Ny < c.Nz ? c.Ny : c.Nz);
    if (pp.n_range > nmin)
    {
        std::ostringstream o;
        o << "盒子太小：相场模板盒边长 " << pp.n_range
          << " 超过了最短的网格边 " << nmin << "（粒子会通过周期边界与自己作用）";
        err = o.str();
        return false;
    }

    // 显式粘性项的稳定性上界，只警告不拒绝（用户可能在做 dt 扫描）
    const double dt_max = 1.0 / (2.0 * 3.0 * c.ratio_eta);
    if (c.dt > dt_max)
    {
        fprintf(stderr,
                "[警告] dt = %g 超过显式粘性稳定上界 rho*dx^2/(2*d*eta_c) = %g，很可能发散\n",
                c.dt, dt_max);
    }
    return true;
}

bool dump_config(const FpdConfig& c, const char* path, std::string& err)
{
    std::ofstream ofs(path);
    if (!ofs.is_open())
    {
        err = std::string("无法写入 ") + path;
        return false;
    }

    // 注意：这里必须传一个非 const 引用给 config_fields，用一份拷贝即可
    FpdConfig tmp = c;
    std::vector<FieldDesc> tab = config_fields(tmp);

    ofs << "# 本次运行实际使用的配置（自动生成，可直接作为输入重跑）\n";
    for (size_t i = 0; i < tab.size(); i++)
    {
        ofs << tab[i].key << " = ";
        switch (tab[i].type)
        {
            case F_INT:    ofs << *(int*)tab[i].ptr; break;
            case F_LONG:   ofs << *(long*)tab[i].ptr; break;
            case F_ULL:    ofs << *(unsigned long long*)tab[i].ptr; break;
            case F_DOUBLE: ofs.precision(17); ofs << *(double*)tab[i].ptr; break;
            case F_STRING: ofs << *(std::string*)tab[i].ptr; break;
        }
        ofs << "\n";
    }

    // 派生量只作记录，带 # 前缀，重新读入时会被当注释跳过
    NS_Config ns = make_ns_config(c);
    PhiParams pp = make_phi_params(c);
    ofs << "\n# ---- 以下为派生量，仅供记录，不是输入项 ----\n";
    ofs.precision(17);
    ofs << "# derived inv_dt   = " << ns.inv_dt << "\n";
    ofs << "# derived W        = " << ns.W      << "   (= sqrt(2*kT/dt))\n";
    ofs << "# derived inv_xi   = " << pp.inv_xi << "\n";
    ofs << "# derived range    = " << pp.range  << "\n";
    ofs << "# derived n_range  = " << pp.n_range << "\n";
    ofs << "# derived range2   = " << pp.range2 << "\n";
    return true;
}

void print_config_help()
{
    FpdConfig c;
    std::vector<FieldDesc> tab = config_fields(c);
    printf("配置项（必填项标 *）：\n");
    for (size_t i = 0; i < tab.size(); i++)
    {
        printf("  %-16s %s %s\n", tab[i].key, tab[i].required ? "*" : " ", tab[i].help);
    }
}

NS_Config make_ns_config(const FpdConfig& c)
{
    return make_ns_config(c.Nx, c.Ny, c.Nz, c.dt, c.kT, c.noise_on != 0);
}

PhiParams make_phi_params(const FpdConfig& c)
{
    return make_phi_params(c.radius, c.xi, c.ratio_eta);
}
