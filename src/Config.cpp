#include "Config.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <fstream>
#include <sstream>

// 三者共用同一张表：解析、--set 覆盖、dump。
std::vector<FieldDesc> config_fields(FpdConfig& c)
{
    return {
        {"Nx",             F_INT,    &c.Nx,             false, "网格 x 方向格子数（由 init_file 提供，输入时拒绝）"},
        {"Ny",             F_INT,    &c.Ny,             false, "网格 y 方向格子数（由 init_file 提供，输入时拒绝）"},
        {"Nz",             F_INT,    &c.Nz,             false, "网格 z 方向格子数（由 init_file 提供，输入时拒绝）"},
        {"dt",             F_DOUBLE, &c.dt,             true,  "时间步长"},
        {"n_steps",        F_LONG,   &c.n_steps,        true,  "绝对目标步数（不是增量）"},
        {"kT",             F_DOUBLE, &c.kT,             true,  "温度 k_B T"},
        {"radius",         F_DOUBLE, &c.radius,         true,  "粒子半径 a"},
        {"xi",             F_DOUBLE, &c.xi,             true,  "界面宽度 xi"},
        {"ratio_eta",      F_DOUBLE, &c.ratio_eta,      true,  "eta_c / eta_l"},

        {"gravity_z",      F_DOUBLE, &c.gravity_z,      false, "z 向外场力（每粒子）"},
        {"gravity_x",      F_DOUBLE, &c.gravity_x,      false, "x 向外场力（每粒子）"},
        {"gravity_y",      F_DOUBLE, &c.gravity_y,      false, "y 向外场力（每粒子）"},
        {"noise_on",       F_INT,    &c.noise_on,       false, "0/1，关掉则 W=0"},

        {"potential",      F_STRING, &c.potential,      false, "粒子间势: none|wca|morse|lj126"},
        {"pot_eps",        F_DOUBLE, &c.pot_eps,        false, "WCA/LJ 的 ε"},
        {"pot_sigma",      F_DOUBLE, &c.pot_sigma,      false, "WCA/LJ 的 σ"},
        {"pot_De",         F_DOUBLE, &c.pot_De,         false, "Morse 阱深 De"},
        {"pot_alpha",      F_DOUBLE, &c.pot_alpha,      false, "Morse 宽度参数 alpha"},
        {"pot_r_eq",       F_DOUBLE, &c.pot_r_eq,       false, "Morse 平衡距离 r_eq"},
        {"pot_rcut",       F_DOUBLE, &c.pot_rcut,       false, "LJ/Morse 截断半径（WCA 派生，不应给）"},
        {"pot_shift",      F_STRING, &c.pot_shift,      false, "截断移位: none|energy|force"},

        {"wallpotential",   F_STRING, &c.wallpotential,   false, "粒子-壁面排斥势: none|wca|morse|lj126"},
        {"wall_eps",       F_DOUBLE, &c.wall_eps,       false, "壁面势 WCA/LJ 的 ε"},
        {"wall_sigma",     F_DOUBLE, &c.wall_sigma,     false, "壁面势 WCA/LJ 的 σ（自变量是表面间隙 h）"},
        {"wall_De",        F_DOUBLE, &c.wall_De,        false, "壁面势 Morse 阱深 De"},
        {"wall_alpha",     F_DOUBLE, &c.wall_alpha,     false, "壁面势 Morse 宽度参数 alpha"},
        {"wall_r_eq",      F_DOUBLE, &c.wall_r_eq,      false, "壁面势 Morse 平衡间隙 r_eq"},
        {"wall_rcut",      F_DOUBLE, &c.wall_rcut,      false, "壁面势 LJ/Morse 截断间隙（WCA 派生，不应给）"},
        {"wall_shift",     F_STRING, &c.wall_shift,     false, "壁面势截断移位: none|energy|force"},

        {"init_file",      F_STRING, &c.init_file,      false, ".fpd 初始构型或续跑文件；留空用内置默认"},
        {"out_dir",        F_STRING, &c.out_dir,        false, "输出目录"},
        {"run_name",       F_STRING, &c.run_name,       false, "输出文件名前缀"},
        {"interval_ckpt",  F_LONG,   &c.interval_ckpt,  false, "每多少步写一个 .fpd"},
        {"interval_log",   F_LONG,   &c.interval_log,   false, "每多少步打印一行进度"},
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

// 配置文件与 --set 共用字段查找和赋值，返回字段下标供必填项记录使用。
static int assign_field(const std::vector<FieldDesc>& tab, const std::string& key,
                        const std::string& val, std::string& err)
{
    for (size_t i = 0; i < tab.size(); i++)
    {
        if (key == tab[i].key)
        {
            return set_field(tab[i], val, err) ? (int)i : -1;
        }
    }
    err = "未知的配置项 '" + key + "'";
    return -1;
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

        const int hit = assign_field(tab, key, val, err);
        if (hit < 0)
        {
            err = std::string(path) + ":" + std::to_string(lineno) + " " + err;
            return false;
        }
        c.keys_given.push_back(key);
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

    if (assign_field(config_fields(c), key, val, err) < 0) { return false; }
    c.keys_given.push_back(key);
    return true;
}

// 用户是否真的明确提供过这个参数。
static bool was_given(const std::vector<std::string>& keys, const std::string& k)
{
    for (size_t i = 0; i < keys.size(); i++) { if (keys[i] == k) { return true; } }
    return false;
}

// 只检查势名称和必填参数；数值范围、截断与稳定性由使用者保证。
bool validate_config(const FpdConfig& c, std::string& err)
{
    if (c.init_file.empty()) { err = "必须提供 init_file"; return false; }
    for (const char* key : {"Nx", "Ny", "Nz"})
    {
        if (was_given(c.keys_given, key))
        { err = std::string("网格尺寸由 init_file 提供，不允许设置 ") + key; return false; }
    }
    FpdConfig tmp = c;
    for (const FieldDesc& f : config_fields(tmp))
    {
        if (f.required && !was_given(c.keys_given, f.key))
        { err = std::string("缺少必填项 ") + f.key; return false; }
    }
    struct PotFamily
    {
        const char* prefix;        // "pot_" / "wall_"
        const char* potential_key; // "potential" / "wallpotential"
        const char* shift_key;   // "pot_shift" / "wall_shift"
        const std::string* pot;  // 当前取值
        const std::string* shift;
    };
    const PotFamily fams[2] = {
        {"pot_",  "potential", "pot_shift",  &c.potential,  &c.pot_shift},
        {"wall_", "wallpotential", "wall_shift", &c.wallpotential, &c.wall_shift},
    };
    for (int f = 0; f < 2; f++)
    {
        const PotFamily& F = fams[f];
        const std::string& pn = *F.pot;
        if (pn != "none" && pn != "wca" && pn != "morse" && pn != "lj126")
        {
            err = std::string(F.potential_key) + " 的合法取值是 none | wca | morse | lj126，收到 '" + pn + "'";
            return false;
        }
        if (*F.shift != "none" && *F.shift != "energy" && *F.shift != "force")
        {
            err = std::string(F.shift_key) + " 的合法取值是 none | energy | force，收到 '" + *F.shift + "'";
            return false;
        }

        // 需要的势参数必须给（不给默认值兜底，与必填项一致）
        {
            auto p = [&F](const char* suffix) { return std::string(F.prefix) + suffix; };
            std::string needed[4];
            int nn = 0;
            if (pn == "wca")   { needed[nn++] = p("eps"); needed[nn++] = p("sigma"); }
            if (pn == "morse") { needed[nn++] = p("De");  needed[nn++] = p("alpha");
                                 needed[nn++] = p("r_eq"); needed[nn++] = p("rcut"); }
            if (pn == "lj126") { needed[nn++] = p("eps"); needed[nn++] = p("sigma"); needed[nn++] = p("rcut"); }
            std::string missing;
            for (int j = 0; j < nn; j++)
            {
                if (!was_given(c.keys_given, needed[j]))
                { missing += (missing.empty() ? "" : ", "); missing += needed[j]; }
            }
            if (!missing.empty())
            {
                err = std::string(F.potential_key) + " = " + pn + " 需要 " + missing +
                      "（势参数不给默认值兜底）";
                return false;
            }
        }
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
        if (std::strcmp(tab[i].key, "Nx") == 0 ||
            std::strcmp(tab[i].key, "Ny") == 0 ||
            std::strcmp(tab[i].key, "Nz") == 0) { continue; }
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
    ofs << "# z 向边界        = noslip (x/y 周期，z 上下无滑移硬壁)\n";
    ofs << "# 背景力密度补偿  = 0 (无滑移壁面本身就是真实的动量汇)\n";
    ofs << "# derived inv_dt   = " << ns.inv_dt << "\n";
    ofs << "# derived W        = " << ns.W      << "   (= sqrt(2*kT/dt))\n";
    ofs << "# derived inv_xi   = " << pp.inv_xi << "\n";
    ofs << "# derived range    = " << pp.range  << "\n";
    ofs << "# derived n_range  = " << pp.n_range << "\n";
    ofs << "# derived range2   = " << pp.range2 << "\n";
    if (c.potential != "none")
    {
        PotentialParams potp = make_potential_params(c);
        ofs << "# potential       = " << c.potential << "   shift = " << c.pot_shift
            << "   rcut = " << potp.rcut << "\n";
        ofs << "#   U(rcut) = " << potp.u_at_rc << "   U'(rcut) = " << potp.dudr_at_rc << "\n";
    }
    if (c.wallpotential != "none")
    {
        WallParams wp = make_wall_params(c);
        ofs << "# wallpotential   = " << c.wallpotential << "   shift = " << c.wall_shift
            << "   rcut = " << wp.pot.rcut << "   a = " << wp.radius << "\n";
        ofs << "#   U(rcut) = " << wp.pot.u_at_rc << "   U'(rcut) = " << wp.pot.dudr_at_rc << "\n";
    }
    ExternalField ext = make_external_field(c);
    if (ext.gx != 0.0 || ext.gy != 0.0 || ext.gz != 0.0)
    {
        ofs << "# external_field  = (" << ext.gx << ", " << ext.gy << ", " << ext.gz << ")\n";
    }
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

PotentialParams make_potential_params(const FpdConfig& c)
{
    int type;
    if      (c.potential == "wca")   { type = POT_WCA;   }
    else if (c.potential == "morse") { type = POT_MORSE; }
    else if (c.potential == "lj126") { type = POT_LJ126; }
    else                             { type = POT_NONE;  }

    int shift;
    if      (c.pot_shift == "none")   { shift = SHIFT_NONE;   }
    else if (c.pot_shift == "force")  { shift = SHIFT_FORCE;  }
    else                              { shift = SHIFT_ENERGY; }

    return make_potential_params(type, shift, c.pot_eps, c.pot_sigma,
                                 c.pot_De, c.pot_alpha, c.pot_r_eq, c.pot_rcut);
}

// 与 make_potential_params 逐字同构，只是换成 wall_* 字段。
// ⚠️ 半径 a 取自 radius（同一个粒子半径）：间隙 h 的定义依赖它，必须同源。
WallParams make_wall_params(const FpdConfig& c)
{
    int type;
    if      (c.wallpotential == "wca")   { type = POT_WCA;   }
    else if (c.wallpotential == "morse") { type = POT_MORSE; }
    else if (c.wallpotential == "lj126") { type = POT_LJ126; }
    else                            { type = POT_NONE;  }

    int shift;
    if      (c.wall_shift == "none")   { shift = SHIFT_NONE;   }
    else if (c.wall_shift == "force")  { shift = SHIFT_FORCE;  }
    else                               { shift = SHIFT_ENERGY; }

    WallParams wp;
    wp.pot    = make_potential_params(type, shift, c.wall_eps, c.wall_sigma,
                                      c.wall_De, c.wall_alpha, c.wall_r_eq, c.wall_rcut);
    wp.radius = c.radius;
    return wp;
}

ExternalField make_external_field(const FpdConfig& c)
{
    ExternalField ext;
    ext.gx = c.gravity_x;
    ext.gy = c.gravity_y;
    ext.gz = c.gravity_z;
    return ext;
}
