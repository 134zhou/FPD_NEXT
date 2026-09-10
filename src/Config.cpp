#include "Config.h"

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
        {"gravity_x",      F_DOUBLE, &c.gravity_x,      false, "x 向外场力（每粒子）"},
        {"gravity_y",      F_DOUBLE, &c.gravity_y,      false, "y 向外场力（每粒子）"},
        {"gravity_compensate", F_INT, &c.gravity_compensate, false,
                                                   "-1=auto（periodic→1, noslip→0）| 0/1 显式"},
        {"boundary_z",     F_STRING, &c.boundary_z,     false, "z 向边界: periodic|noslip（壁面）"},
        {"noise_on",       F_INT,    &c.noise_on,       false, "0/1，关掉则 W=0"},

        {"potential",      F_STRING, &c.potential,      false, "粒子间势: none|wca|morse|lj126"},
        {"pot_eps",        F_DOUBLE, &c.pot_eps,        false, "WCA/LJ 的 ε"},
        {"pot_sigma",      F_DOUBLE, &c.pot_sigma,      false, "WCA/LJ 的 σ"},
        {"pot_De",         F_DOUBLE, &c.pot_De,         false, "Morse 阱深 De"},
        {"pot_alpha",      F_DOUBLE, &c.pot_alpha,      false, "Morse 宽度参数 alpha"},
        {"pot_r_eq",       F_DOUBLE, &c.pot_r_eq,       false, "Morse 平衡距离 r_eq"},
        {"pot_rcut",       F_DOUBLE, &c.pot_rcut,       false, "LJ/Morse 截断半径（WCA 派生，不应给）"},
        {"pot_shift",      F_STRING, &c.pot_shift,      false, "截断移位: none|energy|force"},

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
        c.keys_given.push_back(key);
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
        if (key == tab[i].key)
        {
            if (!set_field(tab[i], val, err)) { return false; }
            c.keys_given.push_back(key);
            return true;
        }
    }
    err = std::string("--set 用了未知的配置项 '") + key + "'" + suggest(tab, key);
    return false;
}

// 势参数查表的 helper：防「用不到的参数被静默忽略」。
static bool is_pot_param(const std::string& key)
{
    return key.size() > 4 && key.compare(0, 4, "pot_") == 0;
}
static bool potential_uses(const std::string& pot, const std::string& key)
{
    if (pot == "wca")   { return key == "pot_eps" || key == "pot_sigma"; }
    if (pot == "morse") { return key == "pot_De" || key == "pot_alpha" || key == "pot_r_eq" || key == "pot_rcut"; }
    if (pot == "lj126") { return key == "pot_eps" || key == "pot_sigma" || key == "pot_rcut"; }
    return false;
}
static bool was_given(const std::vector<std::string>& keys, const std::string& k)
{
    for (size_t i = 0; i < keys.size(); i++) { if (keys[i] == k) { return true; } }
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

    // --- z 向边界 ---
    if (c.boundary_z != "periodic" && c.boundary_z != "noslip")
    {
        err = "boundary_z 的合法取值是 periodic | noslip，收到 '" + c.boundary_z + "'";
        return false;
    }
    const bool wall = (c.boundary_z == "noslip");
    if (wall && c.Nz < 2)
    {
        err = "boundary_z=noslip 要求 Nz >= 2（两面壁之间至少要有 1 个自由 z 面）";
        return false;
    }
    // 壁面提供真实动量汇，再叠加背景力密度等于双重扣除。显式写 1 直接拒绝
    // （而不是静默改成 0）—— 「不做静默行为变更」是本项目的一贯规矩。
    if (wall && c.gravity_compensate == 1)
    {
        err = "boundary_z=noslip 与 gravity_compensate=1 冲突：壁面已经提供真实的动量汇，"
              "再叠加背景力密度 bg=-ΣF/size 等于双重扣除。"
              "请删掉 gravity_compensate 这一行（默认 auto，壁面下取 0）或显式设为 0";
        return false;
    }

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

    // --- 势函数检查：合法性 ---
    if (c.potential != "none" && c.potential != "wca" && c.potential != "morse" && c.potential != "lj126")
    {
        err = "potential 的合法取值是 none | wca | morse | lj126，收到 '" + c.potential + "'";
        return false;
    }
    if (c.pot_shift != "none" && c.pot_shift != "energy" && c.pot_shift != "force")
    {
        err = "pot_shift 的合法取值是 none | energy | force，收到 '" + c.pot_shift + "'";
        return false;
    }

    // --- 参数查表：用不到的势参数报错（静默忽略是浪费一天的事故）---
    for (size_t i = 0; i < c.keys_given.size(); i++)
    {
        const std::string& k = c.keys_given[i];
        if (is_pot_param(k) && k != "pot_shift" && !potential_uses(c.potential, k))
        {
            err = "potential = " + c.potential + " 用不到 " + k + "，请删掉它或改 potential";
            return false;
        }
    }
    // 需要的势参数必须给（不给默认值兜底，与必填项一致）
    {
        const char* needed[4] = {0, 0, 0, 0};
        if (c.potential == "wca")   { needed[0]="pot_eps";  needed[1]="pot_sigma"; }
        if (c.potential == "morse") { needed[0]="pot_De";   needed[1]="pot_alpha"; needed[2]="pot_r_eq"; needed[3]="pot_rcut"; }
        if (c.potential == "lj126") { needed[0]="pot_eps";  needed[1]="pot_sigma"; needed[2]="pot_rcut"; }
        std::string missing;
        for (int j = 0; j < 4; j++)
        {
            if (needed[j] && !was_given(c.keys_given, needed[j]))
            { missing += (missing.empty() ? "" : ", "); missing += needed[j]; }
        }
        if (!missing.empty())
        {
            err = "potential = " + c.potential + " 需要 " + missing + "（势参数不给默认值兜底）";
            return false;
        }
    }

    // --- 最小镜像硬约束：rcut < min(N)/2（严格小于，等号会双重计数镜像）---
    // ⚠️ 壁面模式下 z 不再周期，最小镜像【不作用于 z】⇒ Nz 不该参与这个上界。
    //    用 min(Nx,Ny) 而不是 min(Nx,Ny,Nz)：否则 z 向薄盒子会被无理由拒绝。
    if (c.potential != "none")
    {
        PotentialParams potp = make_potential_params(c);
        const int nmin_eff = wall ? (c.Nx < c.Ny ? c.Nx : c.Ny) : nmin;
        const double half = 0.5 * (double)nmin_eff;
        if (potp.rcut >= half)
        {
            std::ostringstream o;
            o << "势截断 rcut = " << potp.rcut << " 不小于最短（周期）网格边的一半 " << half
              << "（周期方向取 min(" << (wall ? "Nx,Ny" : "Nx,Ny,Nz") << ") = " << nmin_eff
              << "），最小镜像约定不成立"
              << "（旧代码 cutoff=22.2 配 Nz=32 正是这个错误）";
            err = o.str();
            return false;
        }

        // 残余力比：|U'(rcut)| / max|U'|，> 1e-6 报警（shift=force 时力在 rcut 恒 0，跳过）。
        // 峰值在【物理可达区】[2a, rcut] 上取 —— 下界是粒子直径 2a，更近就重叠了，
        // 那里 U' 可大几个量级但物理不可达，会稀释残余力比。
        if (c.pot_shift != "force")
        {
            const double r_lo = 2.0 * c.radius;
            double max_du = 0.0;
            for (int k = 0; k <= 500; k++)
            {
                const double r  = r_lo + (potp.rcut - r_lo) * (double)k / 500.0;
                const double du = std::fabs(pair_force_over_r_bare(potp, r) * r);
                if (du > max_du) { max_du = du; }
            }
            const double du_rc   = std::fabs(pair_force_over_r_bare(potp, potp.rcut) * potp.rcut);
            const double rho_res = max_du > 0.0 ? du_rc / max_du : 0.0;
            if (rho_res > 1e-6)
            {
                fprintf(stderr,
                        "[警告] 势在 rcut=%g 处残余力比 |U'(rcut)|/max|U'| = %.3e > 1e-6，"
                        "能量移位后力仍跳变。可改 pot_shift=force（改势形状）、调大盒子或调 pot_alpha\n",
                        potp.rcut, rho_res);
            }
        }
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
        // 跳过当前 potential 用不到的势参数：dump 出「用不到」会让 config.used
        // 重跑时 validate_config 报错（完整 dump vs 用不到报错的矛盾）。
        const std::string k = tab[i].key;
        if (is_pot_param(k) && k != "potential")
        {
            if (k == "pot_shift") { if (c.potential == "none") { continue; } }
            else if (!potential_uses(c.potential, k)) { continue; }
        }
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
    ofs << "# derived wall_z            = " << ns.wall_z
        << "   (boundary_z=" << c.boundary_z << ")\n";
    ofs << "# derived gravity_compensate = " << gravity_compensate_of(c)
        << "   (配置值 " << c.gravity_compensate << "，-1 表示 auto)\n";
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
    ExternalField ext = make_external_field(c);
    if (ext.gx != 0.0 || ext.gy != 0.0 || ext.gz != 0.0)
    {
        ofs << "# external_field  = (" << ext.gx << ", " << ext.gy << ", " << ext.gz
            << ")   compensate = " << c.gravity_compensate << "\n";
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

// boundary_z 字符串 -> wall_z 整数。这是【唯一】的转换点。
int wall_z_of(const FpdConfig& c)
{
    return (c.boundary_z == "noslip") ? 1 : 0;
}

// gravity_compensate 的 auto 解析（-1 -> 由 boundary_z 决定）。唯一转换点。
// 调用方【必须】把它打印出来 —— 静默是本项目反复点名的事故类别。
int gravity_compensate_of(const FpdConfig& c)
{
    if (c.gravity_compensate >= 0) { return c.gravity_compensate; }
    return wall_z_of(c) ? 0 : 1;
}

NS_Config make_ns_config(const FpdConfig& c)
{
    return make_ns_config(c.Nx, c.Ny, c.Nz, c.dt, c.kT, c.noise_on != 0, wall_z_of(c));
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

ExternalField make_external_field(const FpdConfig& c)
{
    ExternalField ext;
    ext.gx = c.gravity_x;
    ext.gy = c.gravity_y;
    ext.gz = c.gravity_z;
    return ext;
}
