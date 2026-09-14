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
        {"boundary_z",     F_STRING, &c.boundary_z,     false, "z 向边界: noslip（唯一合法值）"},
        {"noise_on",       F_INT,    &c.noise_on,       false, "0/1，关掉则 W=0"},

        {"potential",      F_STRING, &c.potential,      false, "粒子间势: none|wca|morse|lj126"},
        {"pot_eps",        F_DOUBLE, &c.pot_eps,        false, "WCA/LJ 的 ε"},
        {"pot_sigma",      F_DOUBLE, &c.pot_sigma,      false, "WCA/LJ 的 σ"},
        {"pot_De",         F_DOUBLE, &c.pot_De,         false, "Morse 阱深 De"},
        {"pot_alpha",      F_DOUBLE, &c.pot_alpha,      false, "Morse 宽度参数 alpha"},
        {"pot_r_eq",       F_DOUBLE, &c.pot_r_eq,       false, "Morse 平衡距离 r_eq"},
        {"pot_rcut",       F_DOUBLE, &c.pot_rcut,       false, "LJ/Morse 截断半径（WCA 派生，不应给）"},
        {"pot_shift",      F_STRING, &c.pot_shift,      false, "截断移位: none|energy|force"},

        {"wall_pot",       F_STRING, &c.wall_pot,       false, "粒子-壁面排斥势: none|wca|morse|lj126"},
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
//
// ⚠️ 粒子间势（pot_* / potential）与壁面势（wall_* / wall_pot）【共用这三个函数】，
//    靠 prefix + selector 参数区分。曾经想过复制一份 wall_ 版，否掉了 ——
//    复制粘贴后各自漂移正是 C1/C2 的成因。
//
// selector 必须显式传：'potential' 不以 'pot_' 开头（天然被前缀判据排除），
// 但 'wall_pot' 恰好以 'wall_' 开头，不排除就会被当成「用不到的势参数」报错。
static bool is_pot_param(const std::string& key, const char* prefix, const char* selector)
{
    if (key == selector) { return false; }
    const size_t n = std::strlen(prefix);
    return key.size() > n && key.compare(0, n, prefix) == 0;
}
static bool potential_uses(const std::string& pot, const std::string& key,
                           const char* prefix, const char* selector)
{
    if (key == selector) { return false; }
    auto p = [prefix](const char* suffix) { return std::string(prefix) + suffix; };
    if (pot == "wca")   { return key == p("eps") || key == p("sigma"); }
    if (pot == "morse") { return key == p("De") || key == p("alpha") || key == p("r_eq") || key == p("rcut"); }
    if (pot == "lj126") { return key == p("eps") || key == p("sigma") || key == p("rcut"); }
    return false;
}
static bool was_given(const std::vector<std::string>& keys, const std::string& k)
{
    for (size_t i = 0; i < keys.size(); i++) { if (keys[i] == k) { return true; } }
    return false;
}

// ---------------------------------------------------------------------------
// 壁面势的平衡间隙 h_rest：二分反解 |F_wall(h)| = |g_z|。
//
// 两个用处：
//   1. dump_config 的派生量行 —— 让人一眼看出「粒子会停在哪」（配错了立刻可见）
//   2. 残余力比的采样下界（见 validate_config 里的注释）
//
// 用【单壁近似】：平衡点离壁必 < rcut <= Nz/2，所以对面那面壁的贡献恒为 0。
// gz = 0 时没有平衡点（势把粒子推回体相），返回 false，调用方退回约定下界。
// 势太软、在最深可达间隙仍压不住重力时也返回 false —— 那是【配置错误】，不是数值问题。
// ---------------------------------------------------------------------------
static bool wall_rest_gap(WallParams wp, double gz, double& h_rest)
{
    if (wp.pot.type == POT_NONE) { return false; }
    const double target = std::fabs(gz);
    if (target <= 0.0) { return false; }
    if (wp.pot.rcut <= 0.0) { return false; }

    // F(h) = pair_force_over_r(h²)·h 在 (0, rcut) 上单调递减：+∞ → 0
    auto F = [&wp](double h) { return pair_force_over_r(wp.pot, h * h) * h; };
    double lo = 1e-3 * wp.pot.rcut, hi = wp.pot.rcut;
    if (F(lo) <= target) { return false; }   // 势太软：重力能压穿

    for (int it = 0; it < 200; it++)
    {
        const double mid = 0.5 * (lo + hi);
        if (F(mid) > target) { lo = mid; } else { hi = mid; }
    }
    h_rest = 0.5 * (lo + hi);
    return true;
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

    // --- z 向边界：只有 noslip ---
    // periodic 给出【专门的】「它被删了」而不是「取值不合法」—— 静默或含糊都是
    // 本项目反复点名的事故类别。
    if (c.boundary_z == "periodic")
    {
        err = "boundary_z=periodic 已于 Phase 8-A 删除（三维 FFT 全周期路径随 z 周期"
              "边界一起移除）。请改用 boundary_z = noslip（z 上下无滑移硬壁）或删掉这一行";
        return false;
    }
    if (c.boundary_z != "noslip")
    {
        err = "boundary_z 的合法取值只有 noslip，收到 '" + c.boundary_z + "'";
        return false;
    }
    if (c.Nz < 2)
    {
        err = "boundary_z=noslip 要求 Nz >= 2（两面壁之间至少要有 1 个自由 z 面）";
        return false;
    }

    // 相场支撑域：x/y 是周期的，粒子仍会通过周期边界与自己作用；z 向不周期，
    // 但同样需要 Nz 装得下「离两壁都 >= range」的粒子。
    PhiParams pp = make_phi_params(c);
    const int nmin = c.Nx < c.Ny ? (c.Nx < c.Nz ? c.Nx : c.Nz) : (c.Ny < c.Nz ? c.Ny : c.Nz);
    if (pp.n_range > nmin)
    {
        std::ostringstream o;
        o << "盒子太小：相场模板盒边长 " << pp.n_range
          << " 超过了最短的网格边 " << nmin
          << "（x/y 周期方向粒子会与自己作用；z 向也放不下离两壁各 "
          << (pp.n_range / 2) << " 的粒子）";
        err = o.str();
        return false;
    }

    // --- 势函数检查：合法性 ---
    // 粒子间势与壁面势共用同一套 enum 与同一段逻辑，只换 prefix/selector。
    // 「加新势族要改三处、漏一处就静默」正是这个 helper 要消灭的东西。
    struct PotFamily
    {
        const char* prefix;      // "pot_" / "wall_"
        const char* selector;    // "potential" / "wall_pot"
        const char* shift_key;   // "pot_shift" / "wall_shift"
        const std::string* pot;  // 当前取值
        const std::string* shift;
    };
    const PotFamily fams[2] = {
        {"pot_",  "potential", "pot_shift",  &c.potential,  &c.pot_shift},
        {"wall_", "wall_pot",  "wall_shift", &c.wall_pot,   &c.wall_shift},
    };
    for (int f = 0; f < 2; f++)
    {
        const PotFamily& F = fams[f];
        const std::string& pn = *F.pot;
        if (pn != "none" && pn != "wca" && pn != "morse" && pn != "lj126")
        {
            err = std::string(F.selector) + " 的合法取值是 none | wca | morse | lj126，收到 '" + pn + "'";
            return false;
        }
        if (*F.shift != "none" && *F.shift != "energy" && *F.shift != "force")
        {
            err = std::string(F.shift_key) + " 的合法取值是 none | energy | force，收到 '" + *F.shift + "'";
            return false;
        }

        // 用不到的势参数报错（静默忽略是浪费一天的事故）
        for (size_t i = 0; i < c.keys_given.size(); i++)
        {
            const std::string& k = c.keys_given[i];
            if (is_pot_param(k, F.prefix, F.selector) && k != F.shift_key &&
                !potential_uses(pn, k, F.prefix, F.selector))
            {
                err = std::string(F.selector) + " = " + pn + " 用不到 " + k +
                      "，请删掉它或改 " + F.selector;
                return false;
            }
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
                err = std::string(F.selector) + " = " + pn + " 需要 " + missing +
                      "（势参数不给默认值兜底）";
                return false;
            }
        }
    }

    // --- 最小镜像硬约束：rcut < min(N)/2（严格小于，等号会双重计数镜像）---
    // ⚠️ z 不再周期，最小镜像【不作用于 z】⇒ Nz 不该参与这个上界。
    //    用 min(Nx,Ny) 而不是 min(Nx,Ny,Nz)：否则 z 向薄盒子会被无理由拒绝。
    if (c.potential != "none")
    {
        PotentialParams potp = make_potential_params(c);
        const int nmin_eff = c.Nx < c.Ny ? c.Nx : c.Ny;
        const double half = 0.5 * (double)nmin_eff;
        if (potp.rcut >= half)
        {
            std::ostringstream o;
            o << "势截断 rcut = " << potp.rcut << " 不小于最短周期网格边的一半 " << half
              << "（周期方向只有 x/y，取 min(Nx,Ny) = " << nmin_eff
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

    // --- 壁面排斥势 ---
    // 壁面是【非周期】的，所以上面那条最小镜像硬约束不适用（没有镜像可双重计数）。
    if (c.wall_pot != "none")
    {
        WallParams wp = make_wall_params(c);

        // 体相区长度 = Nz - 2*wall_rcut。小于 0 意味着每个粒子都同时受两面壁作用。
        // 强受限体系可能【有意】如此，所以只警告不拒绝（与 dt 稳定性上界同一原则）。
        if (2.0 * wp.pot.rcut > (double)c.Nz)
        {
            fprintf(stderr,
                    "[警告] 2*wall_rcut = %g > Nz = %d：盒子里没有「同时离两壁都超过截断」"
                    "的体相区，每个粒子都同时受两面壁作用。若非有意，请调小 wall_rcut "
                    "或加大 Nz\n", 2.0 * wp.pot.rcut, c.Nz);
        }

        // 残余力比：|U'(rcut)| / max|U'|，> 1e-6 报警。
        // ⚠️ 采样下界用【平衡间隙 h_rest】而不是 0。壁面势的物理可达区是
        //    [h_rest, rcut]（粒子停在 h_rest 附近）。若照抄粒子间势用 [0, rcut]，
        //    h→0 处 U' 大几个量级会把比值稀释到永远低于报警线，判据变成哑的。
        if (c.wall_shift != "force")
        {
            double h_rest = 0.0;
            const bool have_rest = wall_rest_gap(wp, c.gravity_z, h_rest);
            if (c.gravity_z != 0.0 && !have_rest)
            {
                fprintf(stderr,
                        "[警告] 壁面势在最深可达间隙处的排斥力仍小于 |gravity_z| = %g，"
                        "粒子会在重力下压穿壁面势。请调大 wall_eps 或减小 gravity_z\n",
                        std::fabs(c.gravity_z));
            }
            const double r_lo = have_rest ? h_rest : 0.1 * wp.pot.rcut;
            double max_du = 0.0;
            for (int k = 0; k <= 500; k++)
            {
                const double r  = r_lo + (wp.pot.rcut - r_lo) * (double)k / 500.0;
                const double du = std::fabs(pair_force_over_r_bare(wp.pot, r) * r);
                if (du > max_du) { max_du = du; }
            }
            const double du_rc   = std::fabs(pair_force_over_r_bare(wp.pot, wp.pot.rcut) * wp.pot.rcut);
            const double rho_res = max_du > 0.0 ? du_rc / max_du : 0.0;
            if (rho_res > 1e-6)
            {
                fprintf(stderr,
                        "[警告] 壁面势在 wall_rcut=%g 处残余力比 |U'(rcut)|/max|U'| = %.3e "
                        "> 1e-6（采样区间 [%g, %g]）。可改 wall_shift=force、调大 wall_rcut "
                        "或调 wall_alpha\n", wp.pot.rcut, rho_res, r_lo, wp.pot.rcut);
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
        // 粒子间势与壁面势各跑一遍，同一段逻辑（见 is_pot_param 的注释）。
        const std::string k = tab[i].key;
        {
            struct { const char* prefix; const char* selector;
                     const char* shift_key; const std::string* pot; } fams[2] = {
                {"pot_",  "potential", "pot_shift",  &c.potential },
                {"wall_", "wall_pot",  "wall_shift", &c.wall_pot  },
            };
            bool skip = false;
            for (int f = 0; f < 2 && !skip; f++)
            {
                if (!is_pot_param(k, fams[f].prefix, fams[f].selector)) { continue; }
                if (k == fams[f].shift_key) { skip = (*fams[f].pot == "none"); }
                else { skip = !potential_uses(*fams[f].pot, k, fams[f].prefix, fams[f].selector); }
            }
            if (skip) { continue; }
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
    if (c.wall_pot != "none")
    {
        WallParams wp = make_wall_params(c);
        double h_rest = 0.0;
        const bool have_rest = wall_rest_gap(wp, c.gravity_z, h_rest);
        ofs << "# wall_pot        = " << c.wall_pot << "   shift = " << c.wall_shift
            << "   rcut = " << wp.pot.rcut << "   a = " << wp.radius << "\n";
        ofs << "#   U(rcut) = " << wp.pot.u_at_rc << "   U'(rcut) = " << wp.pot.dudr_at_rc << "\n";
        // 平衡间隙是运维时唯一能一眼看出「壁面势配错了」的数（配软了压在壁上、
        // 配硬了粒子悬在半空）。配不出平衡点时要明说，不能只留个空白。
        ofs << "#   h_rest  = ";
        if (have_rest)
        {
            ofs << h_rest << "   (解 |F_wall(h)| = |gravity_z| = "
                << std::fabs(c.gravity_z) << ")\n";
            ofs << "#   粒子中心静止高度 z_rest = " << (-0.5 + wp.radius + h_rest) << "\n";
        }
        else
        {
            ofs << "无（gravity_z = 0，或势太软压不住重力）\n";
        }
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
    if      (c.wall_pot == "wca")   { type = POT_WCA;   }
    else if (c.wall_pot == "morse") { type = POT_MORSE; }
    else if (c.wall_pot == "lj126") { type = POT_LJ126; }
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
