"""用独立构造的 v3 文件验证历史转换入口；在 pvpython 下运行。"""

from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import fpd_format
import fpd_format_v3 as F
import fpd2vtk_v3


def fixture(path, N, step):
    # 不调用格式模块的写入函数，固定历史字段顺序及已知值。
    fields = [np.arange(64, dtype=float) + t * 100 for t in range(3)]
    particles = [np.arange(N, dtype=float) * 0.25 + t for t in range(12)]
    body = struct.pack("<8sIIiiiiq", b"FPDCKPT\0", 0x01020304, 3, 8, 4, 2, N, step)
    for a in fields + particles:
        body += struct.pack("<%dd" % len(a), *a)
    checksum = 14695981039346656037
    for b in body:
        checksum = ((checksum ^ b) * 1099511628211) & ((1 << 64) - 1)
    path.write_bytes(body + struct.pack("<Q", checksum))
    return fields, particles


with tempfile.TemporaryDirectory() as tmp:
    directory = Path(tmp)
    config = directory / "run.config.used"
    config.write_text("dt = 0.25\n")
    paths = []
    for step, N in ((0, 0), (1, 2), (2, 2)):
        path = directory / ("run_%07d.fpd" % step)
        fields, particles = fixture(path, N, step)
        paths.append(path)
        header = F.read_header(path)
        offsets = F.field_offsets(header)
        assert path.stat().st_size == 48 + 24 * 64 + 96 * N
        assert offsets["_checksum"][0] == path.stat().st_size - 8
        arrays = fpd2vtk_v3.converter.read_arrays(path, header, ["vx", "vy", "vz"] + F.PARTICLE_ARRAYS)
        for name, expected in zip(["vx", "vy", "vz"] + F.PARTICLE_ARRAYS, fields + particles):
            np.testing.assert_array_equal(arrays[name], expected)
        try:
            fpd_format.read_header(path)
        except ValueError:
            pass
        else:
            raise AssertionError("当前 v4 模块意外接受 v3")
    v4 = directory / "v4.fpd"
    v4.write_bytes(struct.pack("<8sIIiiiiq", b"FPDCKPT\0", 0x01020304, 4, 8, 4, 2, 0, 0))
    assert fpd_format.read_header(v4)["version"] == 4
    try:
        F.read_header(v4)
    except ValueError:
        pass
    else:
        raise AssertionError("v3 模块意外接受 v4")
    print("PASS v3 N=0/2 全字段、偏移、长度及 v3/v4 模块隔离", flush=True)

    command = [sys.executable, str(ROOT / "tools/fpd2vtk_v3.py"),
               *map(str, paths), "--config", str(config), "-o", str(directory / "vis"), "-j", "2"]
    subprocess.run(command, check=True)
    for step, N in ((0, 0), (1, 2), (2, 2)):
        stem = "run_%07d" % step
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(directory / "vis" / (stem + ".vtp")))
        reader.Update()
        data = reader.GetOutput()
        assert data.GetNumberOfPoints() == N
        assert data.GetPointData().GetArray("force") is None
        if N:
            for name, start in (("velocity", 6), ("unwrapped", 3)):
                want = np.column_stack([np.arange(N) * 0.25 + t for t in range(start, start + 3)])
                np.testing.assert_array_equal(vtk_to_numpy(data.GetPointData().GetArray(name)), want)
            want = np.column_stack([np.arange(N) * 0.25 + t for t in range(3)])
            np.testing.assert_array_equal(vtk_to_numpy(data.GetPoints().GetData()), want)
        fluid = vtk.vtkXMLImageDataReader()
        fluid.SetFileName(str(directory / "vis" / (stem + ".vti")))
        fluid.Update()
        field = np.arange(64, dtype=float).reshape(2, 4, 8)
        want = np.column_stack([((field + t * 100 + np.roll(field + t * 100, 1, axis=axis)) / 2).ravel()
                                for t, axis in enumerate((2, 1, 0))]).astype(np.float32)
        np.testing.assert_array_equal(vtk_to_numpy(fluid.GetOutput().GetCellData().GetArray("velocity")), want)
        root = ET.parse(directory / "vis" / (stem + ".vtm")).getroot()
        assert [d.get("file") for d in root.iter("DataSet")] == [stem + ".vti", stem + ".vtp"]
    pvd = directory / "vis/run.pvd"
    entries = [(d.get("timestep"), d.get("file")) for d in ET.parse(pvd).getroot().iter("DataSet")]
    assert entries == [("0", "run_0000000.vtm"), ("0.25", "run_0000001.vtm"), ("0.5", "run_0000002.vtm")]
    print("PASS 双进程完整转换，流体/粒子数值、空粒子、无 force、VTM/PVD 时间轴", flush=True)
    subprocess.run(command, check=True)
    assert [(d.get("timestep"), d.get("file")) for d in ET.parse(pvd).getroot().iter("DataSet")] == entries
    print("PASS 重跑增量跳过及索引合并", flush=True)
