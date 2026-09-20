#!/usr/bin/env python3
"""
往返校验：用 ParaView 读回转换产物，与源 .fpd 逐点比对。

  $PV/bin/pvpython tools/verify_vtk.py vis/A.pvd out/A_0000005.fpd --config config/smoke.cfg

这比「检查 XML 是否合法」强得多 —— 它验证的是 ParaView 实际看到的数值
与模拟真正写出的数值一致，且插值方向、轴序、时间轴都对。
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fpd_format as F
import fpd2vtk
from fpd_config import read_values


def main():
    ap = argparse.ArgumentParser(description="验证 .pvd 与对应 .fpd 的数组和时间轴")
    ap.add_argument("pvd")
    ap.add_argument("fpd")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    pvd, fpd = args.pvd, args.fpd

    from paraview.simple import OpenDataFile, UpdatePipeline
    from paraview import servermanager as sm

    header = F.read_header(fpd)
    Nx, Ny, Nz = header["Nx"], header["Ny"], header["Nz"]
    step = header["step"]
    dt = read_values(args.config, ("dt",))["dt"]
    t_want = step * dt

    src = OpenDataFile(pvd)
    src.UpdatePipeline()
    times = list(src.TimestepValues) if hasattr(src, "TimestepValues") else []
    print("时间轴: %d 个时间步  %s" % (len(times), times[:6]))

    ok = True
    if not times:
        print("  [FAIL] .pvd 没有时间轴")
        ok = False
    elif not any(abs(t - t_want) < 1e-12 * max(1.0, abs(t_want)) for t in times):
        print("  [FAIL] 时间轴里没有 %.17g（step=%d, dt=%g）" % (t_want, step, dt))
        ok = False
    else:
        print("  [PASS] 时间轴含 t=%.17g" % t_want)

    UpdatePipeline(time=t_want, proxy=src)
    data = sm.Fetch(src)

    nblk = data.GetNumberOfBlocks() if hasattr(data, "GetNumberOfBlocks") else 0
    print("块数: %d  %s" % (nblk, "[PASS]" if nblk == 2 else "[FAIL] 期望 2"))
    if nblk != 2:
        return 1

    fluid, parts = data.GetBlock(0), data.GetBlock(1)

    # --- 源数据 ---
    want = ["vx", "vy", "vz", "Rx", "Ry", "Rz", "Vx", "Vy", "Vz", "Rux", "Ruy", "Ruz"]
    arrs = fpd2vtk.read_arrays(fpd, header, want)

    from vtk.util import numpy_support as ns

    # --- 流体：速度（应等于源数据插值到格心）---
    vel = ns.vtk_to_numpy(fluid.GetCellData().GetArray("velocity"))
    cx, cy, cz = fpd2vtk.stagger_to_center(arrs["vx"], arrs["vy"], arrs["vz"], Nx, Ny, Nz)
    ref = np.stack([cx.ravel(), cy.ravel(), cz.ravel()], axis=1).astype(np.float32)
    d = np.max(np.abs(vel - ref))
    scale = max(1e-30, float(np.max(np.abs(ref))))
    rel = d / scale
    print("velocity: %d 个格子  最大相对差 %.3e  %s"
          % (vel.shape[0], rel, "[PASS]" if rel < 1e-6 else "[FAIL]"))
    ok &= (rel < 1e-6)

    # --- 粒子 ---
    N = header["N"]
    pts = ns.vtk_to_numpy(parts.GetPoints().GetData())
    R = np.stack([arrs["Rx"], arrs["Ry"], arrs["Rz"]], axis=1)
    d = np.max(np.abs(pts - R))
    print("粒子位置: %d 个  最大绝对差 %.3e  %s"
          % (N, d, "[PASS]" if d < 1e-9 else "[FAIL]"))
    ok &= (d < 1e-9)

    for name, keys in (("velocity", ("Vx", "Vy", "Vz")),
                       ("unwrapped", ("Rux", "Ruy", "Ruz"))):
        a = ns.vtk_to_numpy(parts.GetPointData().GetArray(name))
        r = np.stack([arrs[keys[0]], arrs[keys[1]], arrs[keys[2]]], axis=1).astype(np.float32)
        sc = max(1e-30, float(np.max(np.abs(r))))
        d = np.max(np.abs(a - r)) / sc
        print("粒子 %-10s 最大相对差 %.3e  %s"
              % (name + ":", d, "[PASS]" if d < 1e-6 else "[FAIL]"))
        ok &= (d < 1e-6)

    print("\n%s" % ("全部通过" if ok else "有失败项"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
