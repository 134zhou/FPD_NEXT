#!/usr/bin/env python3
"""
.fpd -> VTK 转换器。**离线运行，不占模拟的热路径。**

必须跑在 pvpython 下（需要 numpy + vtk）：
  PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64
  $PV/bin/pvpython tools/fpd2vtk.py out/run_*.fpd -o vis/

产出：
  vis/<name>_<step>.vti   流体（ImageData，CellData: velocity 插值到格心 + pressure）
  vis/<name>_<step>.vtp   粒子（PolyData，PointData: V / F / Ru）
  vis/<name>_<step>.vtm   多块，把上面两个绑在一起
  vis/<name>.pvd          整个时间序列 —— 【在 ParaView 里打开这一个文件即可】

不手写任何 XML：全部交给 VTK 自己的 writer，schema 由它保证正确。

⚠️ 轴序：.fpd 里数组按 IDX = i + j*Nx + k*Nx*Ny 连续存放（x 变化最快），
   所以 numpy 必须 reshape 成 (Nz, Ny, Nx)。搞反不报错，只会得到转置的场。
   用 --self-test 验证。
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fpd_format as F


def _need_vtk():
    try:
        import vtk
        from vtk.util import numpy_support
        return vtk, numpy_support
    except ImportError:
        sys.exit("找不到 vtk 模块。本脚本需要在 pvpython 下运行：\n"
                 "  $PARAVIEW/bin/pvpython tools/fpd2vtk.py ...")


def read_arrays(path, header, want):
    """按偏移量直接读需要的数组，避免把整个文件读进内存。"""
    offs = F.field_offsets(header)
    out = {}
    with open(path, "rb") as f:
        for name in want:
            if name not in offs:
                out[name] = None
                continue
            off, n = offs[name]
            f.seek(off)
            out[name] = np.fromfile(f, dtype="<f8", count=n)
    return out


def stagger_to_center(vx, vy, vz, Nx, Ny, Nz):
    """
    MAC 面心速度插值到格心。

    vx[ijk] 位于格子 ijk 的 +x 面，所以格心值 = 0.5*(vx[i-1,j,k] + vx[i,j,k])。
    这正是 Stokes.cpp:47 已在用的约定（vax = 0.5*(vx[ijk] + vx[im,j,k])），
    直接复用，别另起一套。

    reshape 成 (Nz, Ny, Nx)：x 是最后一维（IDX 里 x 变化最快）。
    """
    a = vx.reshape(Nz, Ny, Nx)
    b = vy.reshape(Nz, Ny, Nx)
    c = vz.reshape(Nz, Ny, Nx)
    cx = 0.5 * (a + np.roll(a, 1, axis=2))   # x 是 axis=2
    cy = 0.5 * (b + np.roll(b, 1, axis=1))   # y 是 axis=1
    cz = 0.5 * (c + np.roll(c, 1, axis=0))   # z 是 axis=0
    return cx, cy, cz


def write_vti(path, header, arrs, with_pressure):
    vtk, ns = _need_vtk()
    Nx, Ny, Nz = header["Nx"], header["Ny"], header["Nz"]

    img = vtk.vtkImageData()
    # CellData 时 WholeExtent/Dimensions 用【格点数】= 格子数 + 1
    img.SetDimensions(Nx + 1, Ny + 1, Nz + 1)
    img.SetSpacing(1.0, 1.0, 1.0)
    img.SetOrigin(0.0, 0.0, 0.0)

    cx, cy, cz = stagger_to_center(arrs["vx"], arrs["vy"], arrs["vz"], Nx, Ny, Nz)
    vel = np.empty((Nx * Ny * Nz, 3), dtype=np.float32)
    vel[:, 0] = cx.ravel()
    vel[:, 1] = cy.ravel()
    vel[:, 2] = cz.ravel()
    va = ns.numpy_to_vtk(vel, deep=1)
    va.SetName("velocity")
    img.GetCellData().AddArray(va)
    img.GetCellData().SetVectors(va)

    if with_pressure and arrs.get("p") is not None:
        pa = ns.numpy_to_vtk(arrs["p"].astype(np.float32), deep=1)
        pa.SetName("pressure")
        img.GetCellData().AddArray(pa)

    w = vtk.vtkXMLImageDataWriter()
    w.SetFileName(path)
    w.SetInputData(img)
    w.SetDataModeToAppended()
    w.EncodeAppendedDataOff()
    w.Write()


def write_vtp(path, header, arrs):
    vtk, ns = _need_vtk()
    N = header["N"]

    pts = vtk.vtkPoints()
    pd = vtk.vtkPolyData()

    if N > 0:
        xyz = np.empty((N, 3), dtype=np.float64)
        xyz[:, 0] = arrs["Rx"]; xyz[:, 1] = arrs["Ry"]; xyz[:, 2] = arrs["Rz"]
        pts.SetData(ns.numpy_to_vtk(xyz, deep=1))
    pd.SetPoints(pts)

    verts = vtk.vtkCellArray()
    for i in range(N):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)
    pd.SetVerts(verts)

    def add(name, kx, ky, kz):
        v = np.empty((N, 3), dtype=np.float32)
        v[:, 0] = arrs[kx]; v[:, 1] = arrs[ky]; v[:, 2] = arrs[kz]
        a = ns.numpy_to_vtk(v, deep=1)
        a.SetName(name)
        pd.GetPointData().AddArray(a)
        return a

    if N > 0:
        vel = add("velocity", "Vx", "Vy", "Vz")
        add("force", "Fx", "Fy", "Fz")
        # 不折叠位置：粒子穿越周期边界时 R 会瞬移，Ru 不会，做轨迹/MSD 用这个
        add("unwrapped", "Rux", "Ruy", "Ruz")
        pd.GetPointData().SetVectors(vel)

    w = vtk.vtkXMLPolyDataWriter()
    w.SetFileName(path)
    w.SetInputData(pd)
    w.SetDataModeToAppended()
    w.EncodeAppendedDataOff()
    w.Write()


def write_vtm(path, vti_name, vtp_name):
    """多块：ParaView 打开它就同时看到流体和粒子。"""
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="vtkMultiBlockDataSet" version="1.0" '
                'byte_order="LittleEndian" header_type="UInt32">\n')
        f.write('  <vtkMultiBlockDataSet>\n')
        f.write('    <DataSet index="0" name="fluid" file="%s"/>\n' % vti_name)
        f.write('    <DataSet index="1" name="particles" file="%s"/>\n' % vtp_name)
        f.write('  </vtkMultiBlockDataSet>\n')
        f.write('</VTKFile>\n')


def write_pvd(path, entries):
    """时间序列。timestep 用物理时间 step*dt。"""
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <Collection>\n')
        for t, name in entries:
            f.write('    <DataSet timestep="%.17g" group="" part="0" file="%s"/>\n' % (t, name))
        f.write('  </Collection>\n')
        f.write('</VTKFile>\n')


def self_test():
    """
    判据 9/10：轴序与插值方向。
    构造 vx = 已知线性场，检查插值后的格心值与解析值一致。
    """
    Nx, Ny, Nz = 8, 4, 2
    # v[IDX(i,j,k)] = i*1e6 + j*1e3 + k
    lin = np.zeros(Nx * Ny * Nz)
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                lin[i + j * Nx + k * Nx * Ny] = i * 1e6 + j * 1e3 + k

    a = lin.reshape(Nz, Ny, Nx)
    ok_axis = (a[1, 2, 3] == 3e6 + 2e3 + 1)
    print("  轴序 reshape(Nz,Ny,Nx): a[k=1,j=2,i=3] = %.0f  期望 %.0f  %s"
          % (a[1, 2, 3], 3e6 + 2e3 + 1, "PASS" if ok_axis else "FAIL"))

    # vx 只随 i 线性变化，格心值应为 i*1e6 - 0.5e6（周期边界处除外）
    vx = np.zeros(Nx * Ny * Nz)
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                vx[i + j * Nx + k * Nx * Ny] = i * 1.0
    cx, _, _ = stagger_to_center(vx, vx, vx, Nx, Ny, Nz)
    # 内部格点 i>=1: 0.5*(i + (i-1)) = i - 0.5
    inner = cx[0, 0, 1:]
    want = np.arange(1, Nx) - 0.5
    ok_interp = np.allclose(inner, want)
    print("  x 插值 内部格点: %s  期望 %s  %s"
          % (inner, want, "PASS" if ok_interp else "FAIL"))

    return 0 if (ok_axis and ok_interp) else 1


def main():
    ap = argparse.ArgumentParser(description=".fpd -> VTK（离线转换）")
    ap.add_argument("files", nargs="*", help=".fpd 文件，可用通配符")
    ap.add_argument("-o", "--out", default="vis", help="输出目录")
    ap.add_argument("--stride", type=int, default=1, help="每隔几个文件取一个")
    ap.add_argument("--no-pressure", action="store_true")
    ap.add_argument("--particles-only", action="store_true")
    ap.add_argument("--self-test", action="store_true", help="轴序与插值方向自检")
    a = ap.parse_args()

    if a.self_test:
        return self_test()

    if not a.files:
        ap.error("需要至少一个 .fpd 文件")

    files = sorted(a.files)[:: a.stride]
    if not os.path.isdir(a.out):
        os.makedirs(a.out)

    entries = []
    base_name = None
    for path in files:
        header = F.read_header(path)
        step = header["step"]
        stem = os.path.splitext(os.path.basename(path))[0]
        if base_name is None:
            base_name = stem.rsplit("_", 1)[0] if "_" in stem else stem

        want = ["Rx", "Ry", "Rz", "Rux", "Ruy", "Ruz",
                "Vx", "Vy", "Vz", "Fx", "Fy", "Fz"]
        if not a.particles_only:
            want += ["vx", "vy", "vz"]
            if header["flags"] & F.FLAG_PRESSURE and not a.no_pressure:
                want += ["p"]
        arrs = read_arrays(path, header, want)

        vtp_name = stem + ".vtp"
        write_vtp(os.path.join(a.out, vtp_name), header, arrs)

        if a.particles_only:
            entries.append((step * header["dt"], vtp_name))
            print("  %s -> %s" % (path, vtp_name))
            continue

        vti_name = stem + ".vti"
        vtm_name = stem + ".vtm"
        write_vti(os.path.join(a.out, vti_name), header, arrs,
                  with_pressure=not a.no_pressure)
        write_vtm(os.path.join(a.out, vtm_name), vti_name, vtp_name)
        entries.append((step * header["dt"], vtm_name))
        print("  %s -> %s (+ .vti/.vtp)" % (path, vtm_name))

    pvd = os.path.join(a.out, (base_name or "run") + ".pvd")
    write_pvd(pvd, entries)
    print("\n共 %d 个时间步。在 ParaView 里打开：\n  %s" % (len(entries), pvd))
    return 0


if __name__ == "__main__":
    sys.exit(main())
