#!/usr/bin/env python3
"""
生成 .fpd 初始构型（step=0、速度场为零）。

**纯 stdlib**，在裸 conda base 里直接能跑，无需 numpy。

  python3 tools/make_init.py --grid 128 64 32 --single              -o init/single.fpd
  python3 tools/make_init.py --grid 128 128 64 --lattice 363        -o init/n363.fpd
  python3 tools/make_init.py --grid 64 64 64 --random 50 --min-sep 8 -o init/n50.fpd
  python3 tools/make_init.py --grid 8 4 2 --pattern index           -o /tmp/t.fpd

  python3 tools/make_init.py --grid 128 64 32 --phi 0.05              -o init/phi05.fpd

--pattern index 把 vx/vy/vz 填成 i*1e6 + j*1e3 + k，用于验证 C++ 与 Python
两侧的轴序约定一致（判据 3）。轴序搞反不报错，只会得到一个转置的场。

z 向边界只有一种：z 上下无滑移硬壁，x/y 周期。粒子中心自动被约束在离两壁
>= --wall-margin 的区间内，越界的构型会直接报错退出。

体积分数 φ 有【两个】常用约定，本脚本两个都打印，不替你选：
  - 标称值   用解析球体积 (4/3)πa³          —— 文献里最常见的约定
  - 扩散界面 用 ∫φ dV = (4/3)πa³ + 8πaξ²π²/24 —— 相场真值，比解析球大 24%
误用哪一个会让 M_i / λ^T 差 24%（README 记录的陷阱），所以两个数都摆在台面上。
"""

import argparse
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fpd_format as F


def place_single(Nx, Ny, Nz):
    # 放在盒心的一般亚格点位置，沿用 --check 的约定（避开整数与半格这两个特殊点）
    return [(Nx / 2.0 + 0.3, Ny / 2.0 + 0.7, Nz / 2.0 + 0.5)]


def place_lattice(Nx, Ny, Nz, n, min_sep):
    """规则简单立方堆积，尽量接近 n 个且间距 >= min_sep。"""
    best = None
    # 找一个能容纳 n 个粒子的最小分格数
    for nc in range(1, max(Nx, Ny, Nz) + 1):
        cx, cy, cz = Nx / nc, Ny / nc, Nz / nc
        if min(cx, cy, cz) < min_sep:
            break
        if nc ** 3 >= n:
            best = nc
            break
    if best is None:
        raise SystemExit("盒子放不下 %d 个间距 >= %g 的粒子" % (n, min_sep))

    pts = []
    for k in range(best):
        for j in range(best):
            for i in range(best):
                if len(pts) >= n:
                    break
                pts.append(((i + 0.5) * Nx / best,
                            (j + 0.5) * Ny / best,
                            (k + 0.5) * Nz / best))
    return pts


def minsq_with(p, q, Nx, Ny, Nz):
    """p 与 q 的距离平方，带最小镜像。x/y 周期，z 向【不】折叠（无滑移壁面）。"""
    dx = p[0] - q[0]; dy = p[1] - q[1]; dz = p[2] - q[2]
    dx -= Nx * round(dx / Nx)
    dy -= Ny * round(dy / Ny)
    return dx * dx + dy * dy + dz * dz, dz


def place_random(Nx, Ny, Nz, n, min_sep, seed, z_margin=0.0, max_try=200000):
    """拒绝采样，保证间距 >= min_sep。z 向不折叠，且离两壁 >= z_margin。"""
    rng = random.Random(seed)
    pts = []
    tries = 0
    s2 = min_sep * min_sep
    zlo, zhi = -0.5 + z_margin, Nz - 0.5 - z_margin
    if zlo >= zhi:
        raise SystemExit(
            "盒子太薄：Nz=%d 放不下离两壁各 %g 的粒子" % (Nz, z_margin))
    while len(pts) < n:
        tries += 1
        if tries > max_try:
            raise SystemExit(
                "放了 %d/%d 个粒子后 %d 次尝试仍失败：盒子太小或 min_sep 太大"
                % (len(pts), n, max_try))
        p = (rng.uniform(0, Nx), rng.uniform(0, Ny), rng.uniform(zlo, zhi))
        ok = True
        for q in pts:
            dsq, _ = minsq_with(p, q, Nx, Ny, Nz)
            if dsq < s2:
                ok = False
                break
        if ok:
            pts.append(p)
    return pts


def phi_to_n(phi, Nx, Ny, Nz, radius, diffuse=False, xi=1.0):
    """由目标体积分数反推粒子数。

    diffuse=False 用解析球 (4/3)πa³（文献常用约定）；
    diffuse=True  用扩散界面真值 ∫φ = (4/3)πa³ + 8πaξ²π²/24 + O(ξ³)。
    两者差 24%（a=3.2, ξ=1 时）—— 这正是 README 点名的「用解析球算 M_i 会毁掉验证」
    那件事，所以【调用方必须知道自己在用哪一个】，本函数不给默认值兜底之外的隐含选择。
    """
    sphere = 4.0 / 3.0 * math.pi * radius ** 3
    vol = sphere + (8.0 * math.pi * radius * xi * xi * math.pi ** 2 / 24.0) if diffuse else sphere
    return int(round(phi * Nx * Ny * Nz / vol))


def main():
    ap = argparse.ArgumentParser(description="生成 .fpd 初始构型")
    ap.add_argument("--grid", nargs=3, type=int, required=True, metavar=("NX", "NY", "NZ"))
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--single", action="store_true", help="单粒子放盒心")
    g.add_argument("--lattice", type=int, metavar="N", help="规则堆积 N 个粒子")
    g.add_argument("--random", type=int, metavar="N", help="随机不重叠 N 个粒子")
    g.add_argument("--phi", type=float, metavar="F",
                   help="随机不重叠，粒子数由目标体积分数 F 反推（标称解析球约定）")
    g.add_argument("--empty", action="store_true", help="无粒子（纯流体测试用）")

    ap.add_argument("--min-sep", type=float, default=None,
                    help="最小间距，默认 2*radius+xi（刚好不重叠）")
    ap.add_argument("--radius", type=float, default=3.2)
    ap.add_argument("--xi", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--pattern", choices=["zero", "index"], default="zero",
                    help="速度场填法；index 用于轴序互操作判据")
    ap.add_argument("--wall-margin", type=float, default=None,
                    help="粒子中心离两壁的最小距离，默认 ceil(radius+3*xi)+1"
                         "（≈相场支撑域半径 range）。有壁面势时可显式调小")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    Nx, Ny, Nz = a.grid
    min_sep = a.min_sep if a.min_sep is not None else (2.0 * a.radius + a.xi)

    # 粒子离两壁的边距。默认 ceil(radius+3*xi)+1 —— 这是【相场支撑域】的经验尺寸
    # （a+3ξ 之内 φ 还有 O(1e-3) 的量级），不是 C++ 的 range。
    #   C++: pp.range = (int)(2*(radius+xi))            （Common.h make_phi_params）
    #   此处: ceil(radius+3*xi)+1
    # 默认参数下两者都给 8，所以差异看不出来；改参数时会分叉。
    # ⚠️ 别再引入第三个「边距该多大」的定义 —— 要么用这个，要么用 --wall-margin 显式给。
    z_margin = a.wall_margin if a.wall_margin is not None \
               else float(math.ceil(a.radius + 3.0 * a.xi) + 1)

    # 有壁面势时的硬下界：C++ 的 check_wall_bounds 要求 h = (z+1/2)-a > 0
    #（wallpotential=none 时只要求中心在盒内）。低于它构型会被直接中止。
    h_floor = -0.5 + a.radius
    if z_margin < h_floor:
        raise SystemExit(
            "--wall-margin %g 会让粒子表面越过壁面（需要 >= radius = %g，"
            "因为壁面势的自变量是表面间隙 h = (z+1/2)-a）" % (z_margin, h_floor))

    if a.single:
        pts = place_single(Nx, Ny, Nz)
        if not (-0.5 + z_margin <= pts[0][2] <= Nz - 0.5 - z_margin):
            raise SystemExit("盒心放不下离两壁各 %g 的粒子（Nz=%d）" % (z_margin, Nz))
    elif a.empty:
        pts = []
    elif a.lattice is not None:
        if a.lattice <= 0:
            raise SystemExit("--lattice 必须为正")
        pts = place_lattice(Nx, Ny, Nz, a.lattice, min_sep)
        bad = [p for p in pts if not (-0.5 + z_margin <= p[2] <= Nz - 0.5 - z_margin)]
        if bad:
            raise SystemExit("规则堆积有 %d 个粒子离壁 < %g，请改用 --random"
                             % (len(bad), z_margin))
    else:
        n = a.random
        if a.phi is not None:
            if not (0.0 < a.phi < 1.0):
                raise SystemExit("--phi 必须在 (0, 1) 之间")
            n = phi_to_n(a.phi, Nx, Ny, Nz, a.radius)
            if n <= 0:
                raise SystemExit("--phi %g 在 %dx%dx%d 的盒子里不足 1 个粒子"
                                 % (a.phi, Nx, Ny, Nz))
        pts = place_random(Nx, Ny, Nz, n, min_sep, a.seed, z_margin=z_margin)

    N = len(pts)
    size = Nx * Ny * Nz

    header = F.default_header(Nx, Ny, Nz, N)

    if a.pattern == "index":
        # v[IDX(i,j,k)] = i*1e6 + j*1e3 + k，IDX = i + j*Nx + k*Nx*Ny
        vals = [0.0] * size
        for k in range(Nz):
            for j in range(Ny):
                for i in range(Nx):
                    vals[i + j * Nx + k * Nx * Ny] = i * 1e6 + j * 1e3 + k
        fields = {"vx": vals, "vy": vals, "vz": vals}
    else:
        fields = {"vx": None, "vy": None, "vz": None}

    particles = {
        "Rx":  [p[0] for p in pts],
        "Ry":  [p[1] for p in pts],
        "Rz":  [p[2] for p in pts],
        # 不折叠位置的初值就是折叠位置（MSD 从这里起算）
        "Rux": [p[0] for p in pts],
        "Ruy": [p[1] for p in pts],
        "Ruz": [p[2] for p in pts],
    }

    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)

    h = F.write_fpd(a.out, header, fields, particles)

    print("已写出 %s" % a.out)
    print("  网格 %dx%dx%d   粒子 %d   校验和 0x%016x" % (Nx, Ny, Nz, N, h))
    print("  体积 %.2f MB" % (os.path.getsize(a.out) / 1024.0 / 1024.0))
    if N > 1:
        dmin = min(minsq_with(pts[i], pts[j], Nx, Ny, Nz)[0] ** 0.5
                   for i in range(N) for j in range(i + 1, N))
        print("  最小粒子间距 %.3f  (2a+xi = %.3f)" % (dmin, 2 * a.radius + a.xi))
        if dmin < 2 * a.radius:
            print("  [警告] 有粒子重叠，重叠区 eta 会叠加超过 eta_c")
    print("  z 向无滑移壁面；粒子离壁 >= %.2f（壁面在 z = -0.5 与 %g）"
          % (z_margin, Nz - 0.5))
    if any(p[2] < 0.5 or p[2] > Nz - 1.5 for p in pts):
        print("  [注意] 有粒子的中心落在最外两个格胞内（0.5 或 Nz-1.5 以内）")

    # 体积分数：两个约定都给。文献里的 φ 用哪个不统一，只给一个才是陷阱。
    if N > 0:
        sphere = 4.0 / 3.0 * math.pi * a.radius ** 3
        diffuse = sphere + 8.0 * math.pi * a.radius * a.xi ** 2 * math.pi ** 2 / 24.0
        print("  体积分数 φ（标称解析球 4/3·πa³ = %.2f） = %.4f"
              % (sphere, N * sphere / size))
        print("  体积分数 φ（扩散界面 ∫φ   = %.2f） = %.4f   ← 比标称大 %.1f%%"
              % (diffuse, N * diffuse / size,
                 100.0 * (diffuse / sphere - 1.0)))
        print("  ⚠️ 报 φ 时必须写明用的哪个约定（差 %.0f%% 会让 M_i/λ^T 直接对不上）"
              % (100.0 * (diffuse / sphere - 1.0)))
        if a.phi is not None:
            print("  （--phi %g 按【标称解析球】约定反推，得到 %d 个粒子）"
                  % (a.phi, N))


if __name__ == "__main__":
    main()
