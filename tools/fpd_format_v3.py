"""历史 .fpd v3 的只读布局，固定自提交 42c279a。

头部 40 字节，随后为 3 个流体数组、12 个粒子数组和 8 字节校验和。
总长度为 48 + 24*size + 96*N；末尾 Fx/Fy/Fz 仅用于定位，不导出到 VTK。
本模块不用于当前生产程序和初态生成。
"""

import struct

VERSION = 3
PARTICLE_ARRAYS = ["Rx", "Ry", "Rz", "Rux", "Ruy", "Ruz",
                   "Vx", "Vy", "Vz", "Fx", "Fy", "Fz"]


def default_header(Nx, Ny, Nz, N, step=0):
    return {"version": VERSION, "Nx": Nx, "Ny": Ny, "Nz": Nz,
            "N": N, "step": step}


def read_header(path):
    with open(path, "rb") as f:
        magic, tag, version, Nx, Ny, Nz, N, step = struct.unpack(
            "<8sIIiiiiq", f.read(40))
    if magic != b"FPDCKPT\0":
        raise ValueError("%s 不是 .fpd 文件（magic 不匹配）" % path)
    if tag != 0x01020304:
        raise ValueError("%s 字节序标记不匹配: 0x%08x" % (path, tag))
    if version != VERSION:
        raise ValueError("%s 版本 %d，本工具只认 %d" % (path, version, VERSION))
    header = default_header(Nx, Ny, Nz, N, step)
    header["_data_offset"] = 40
    return header


def size_of(header):
    return header["Nx"] * header["Ny"] * header["Nz"]


def field_offsets(header):
    """v3 的粒子数组顺序包含末尾受力，避免误用 v4 的文件长度。"""
    off = header["_data_offset"]
    out = {}
    for name in ("vx", "vy", "vz"):
        out[name] = (off, size_of(header))
        off += 8 * size_of(header)
    for name in PARTICLE_ARRAYS:
        out[name] = (off, header["N"])
        off += 8 * header["N"]
    out["_checksum"] = (off, 1)
    return out
