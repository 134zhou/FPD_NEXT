#!/usr/bin/env python3
"""生成 .fpd 初态。修改下方数据后直接运行本文件。"""

import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fpd_format as F


# 初态数据：当前默认生成 config/sed.cfg 使用的沉降构型。
NX, NY, NZ = 128, 64, 32
RADIUS, XI = 3.2, 1.0
PHI = 0.05                    # 标称解析球体积分数
SEED = 20260914
MIN_SEP = 2.0 * RADIUS + XI
WALL_MARGIN = 8.0
VELOCITY_PATTERN = "zero"    # 互操作判据改为 "index"
OUT = "out/sed_init.fpd"
POINTS = None                 # None: 生成沉降构型；也可直接写 [] 或坐标列表


def make_particles():
    if POINTS is not None:
        return POINTS

    sphere = 4.0 / 3.0 * math.pi * RADIUS ** 3
    count = round(PHI * NX * NY * NZ / sphere)
    rng = random.Random(SEED)
    points = []
    zlo, zhi = -0.5 + WALL_MARGIN, NZ - 0.5 - WALL_MARGIN

    for _ in range(200000):
        if len(points) == count:
            return points
        p = (rng.uniform(0, NX), rng.uniform(0, NY), rng.uniform(zlo, zhi))
        for q in points:
            dx = p[0] - q[0]
            dy = p[1] - q[1]
            dz = p[2] - q[2]
            dx -= NX * round(dx / NX)
            dy -= NY * round(dy / NY)
            if dx * dx + dy * dy + dz * dz < MIN_SEP ** 2:
                break
        else:
            points.append(p)

    raise SystemExit("无法生成 %d 个互不重叠的粒子" % count)


def make_fields():
    if VELOCITY_PATTERN == "zero":
        return {"vx": None, "vy": None, "vz": None}

    values = [i * 1e6 + j * 1e3 + k
              for k in range(NZ) for j in range(NY) for i in range(NX)]
    return {"vx": values, "vy": values, "vz": values}


def main():
    points = make_particles()
    particles = {}
    if points:
        rx, ry, rz = (list(x) for x in zip(*points))
        particles = {"Rx": rx, "Ry": ry, "Rz": rz,
                     "Rux": rx, "Ruy": ry, "Ruz": rz}

    directory = os.path.dirname(os.path.abspath(OUT))
    os.makedirs(directory, exist_ok=True)
    checksum = F.write_fpd(
        OUT, F.default_header(NX, NY, NZ, len(points)), make_fields(), particles)

    sphere = 4.0 / 3.0 * math.pi * RADIUS ** 3
    diffuse = sphere + math.pi ** 3 * RADIUS * XI ** 2 / 3.0
    size = NX * NY * NZ
    print("已写出 %s：%dx%dx%d，N=%d，校验和 0x%016x"
          % (OUT, NX, NY, NZ, len(points), checksum))
    print("体积分数：标称 %.4f，扩散界面 %.4f"
          % (len(points) * sphere / size, len(points) * diffuse / size))


if __name__ == "__main__":
    main()
