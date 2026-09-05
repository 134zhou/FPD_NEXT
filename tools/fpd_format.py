"""
.fpd 格式的 Python 侧读写。**纯 stdlib**，在裸 conda base 里直接能跑。

布局定义见 src/include/IOBin.h —— 那里是规范，这里是它的镜像实现。

⚠️ 两处独立实现同一个格式和同一个 FNV-1a 哈希，是典型的漂移点。
   改动任一侧后必须重跑互操作判据：
       python3 tools/make_init.py --grid 8 4 2 --pattern index -o /tmp/t.fpd
       ./build/fpd_tool --dump-ckpt /tmp/t.fpd --at 3 1 1

⚠️ 数组按线性 IDX = i + j*Nx + k*Nx*Ny 顺序连续存放，即【x 变化最快】。
   numpy 侧 reshape 必须是 (Nz, Ny, Nx)，x 是最后一维。搞反不报错，
   只会得到一个转置的场。
"""

import struct

MAGIC       = b"FPDCKPT\0"
VERSION     = 1
ENDIAN_TAG  = 0x01020304
FLAG_PRESSURE = 0x1

FNV_OFFSET = 14695981039346656037
FNV_PRIME  = 1099511628211
MASK64     = (1 << 64) - 1

# header 字段：(名字, struct 格式)。顺序必须与 IOBin.h 的 FPD_HEADER_FIELDS 完全一致。
HEADER_FIELDS = [
    ("version",   "<I"),
    ("flags",     "<I"),
    ("Nx",        "<i"),
    ("Ny",        "<i"),
    ("Nz",        "<i"),
    ("N",         "<i"),
    ("step",      "<q"),
    ("dt",        "<d"),
    ("kT",        "<d"),
    ("radius",    "<d"),
    ("xi",        "<d"),
    ("ratio_eta", "<d"),
    ("noise_on",  "<i"),
    ("seed",      "<Q"),
    ("rng_draws", "<Q"),
]

# 粒子数组的顺序，同样必须与 C++ 侧一致
PARTICLE_ARRAYS = ["Rx", "Ry", "Rz", "Rux", "Ruy", "Ruz",
                   "Vx", "Vy", "Vz", "Fx", "Fy", "Fz"]


class Fnv:
    """FNV-1a 64。与 IOBin.cpp 的 Hasher 逐位等价。"""

    def __init__(self):
        self.h = FNV_OFFSET

    def feed(self, buf):
        h = self.h
        for b in buf:
            h = ((h ^ b) * FNV_PRIME) & MASK64
        self.h = h


def _fnv_bytes(h, buf):
    """比逐字节循环快得多的等价实现（大数组时用这个）。"""
    for b in buf:
        h = ((h ^ b) * FNV_PRIME) & MASK64
    return h


def default_header(Nx, Ny, Nz, N, dt=0.002, kT=0.25, radius=3.2, xi=1.0,
                   ratio_eta=50.0, noise_on=1, seed=1234, step=0,
                   rng_draws=0, has_pressure=False):
    return {
        "version": VERSION,
        "flags": FLAG_PRESSURE if has_pressure else 0,
        "Nx": Nx, "Ny": Ny, "Nz": Nz, "N": N,
        "step": step,
        "dt": dt, "kT": kT, "radius": radius, "xi": xi, "ratio_eta": ratio_eta,
        "noise_on": noise_on, "seed": seed, "rng_draws": rng_draws,
    }


def write_fpd(path, header, fields, particles, chunk_doubles=1 << 20):
    """
    fields    : dict，键 vx/vy/vz（必需）与 p（当 has_pressure 时必需）。
                值可以是 None 表示全零（大网格时避免在内存里造 3 个巨型列表），
                也可以是长度为 size 的可迭代 float 序列。
    particles : dict，键见 PARTICLE_ARRAYS，值为长度 N 的序列；缺的按全零处理。
    """
    size = header["Nx"] * header["Ny"] * header["Nz"]
    N = header["N"]
    has_p = bool(header["flags"] & FLAG_PRESSURE)

    h = FNV_OFFSET

    with open(path, "wb") as f:

        def emit(buf):
            nonlocal h
            f.write(buf)
            h = _fnv_bytes(h, buf)

        emit(MAGIC)
        emit(struct.pack("<I", ENDIAN_TAG))
        for name, fmt in HEADER_FIELDS:
            emit(struct.pack(fmt, header[name]))

        zero_chunk = bytes(8 * chunk_doubles)

        def emit_field(vals):
            if vals is None:
                # 全零：分块写，512^3 也不会爆内存
                left = size
                while left > 0:
                    n = min(left, chunk_doubles)
                    emit(zero_chunk[: 8 * n] if n != chunk_doubles else zero_chunk)
                    left -= n
            else:
                if len(vals) != size:
                    raise ValueError("场长度 %d != size %d" % (len(vals), size))
                for i in range(0, size, chunk_doubles):
                    part = vals[i: i + chunk_doubles]
                    emit(struct.pack("<%dd" % len(part), *part))

        emit_field(fields.get("vx"))
        emit_field(fields.get("vy"))
        emit_field(fields.get("vz"))
        if has_p:
            emit_field(fields.get("p"))

        for name in PARTICLE_ARRAYS:
            vals = particles.get(name)
            if vals is None:
                vals = [0.0] * N
            if len(vals) != N:
                raise ValueError("粒子数组 %s 长度 %d != N %d" % (name, len(vals), N))
            emit(struct.pack("<%dd" % N, *vals))

        # 校验和本身不参与哈希
        f.write(struct.pack("<Q", h))

    return h


def read_header(path):
    with open(path, "rb") as f:
        magic = f.read(8)
        if magic[:7] != MAGIC[:7]:
            raise ValueError("%s 不是 .fpd 文件（magic 不匹配）" % path)
        (tag,) = struct.unpack("<I", f.read(4))
        if tag != ENDIAN_TAG:
            raise ValueError("%s 字节序标记不匹配: 0x%08x" % (path, tag))
        header = {}
        for name, fmt in HEADER_FIELDS:
            n = struct.calcsize(fmt)
            (header[name],) = struct.unpack(fmt, f.read(n))
        if header["version"] != VERSION:
            raise ValueError("%s 版本 %d，本工具只认 %d"
                             % (path, header["version"], VERSION))
        header["_data_offset"] = f.tell()
    return header


def size_of(header):
    return header["Nx"] * header["Ny"] * header["Nz"]


def field_offsets(header):
    """返回 {名字: (字节偏移, 元素个数)}，供 numpy 的 fromfile/memmap 直接定位。"""
    off = header["_data_offset"]
    size = size_of(header)
    N = header["N"]
    out = {}
    for name in ("vx", "vy", "vz"):
        out[name] = (off, size)
        off += 8 * size
    if header["flags"] & FLAG_PRESSURE:
        out["p"] = (off, size)
        off += 8 * size
    for name in PARTICLE_ARRAYS:
        out[name] = (off, N)
        off += 8 * N
    out["_checksum"] = (off, 1)
    return out


def verify_checksum(path):
    """整文件重算 FNV，与尾部 8 字节比对。返回 (是否一致, 期望值, 实际值)。"""
    with open(path, "rb") as f:
        data = f.read()
    if len(data) < 8:
        raise ValueError("%s 太短" % path)
    body, tail = data[:-8], data[-8:]
    (want,) = struct.unpack("<Q", tail)
    got = _fnv_bytes(FNV_OFFSET, body)
    return (want == got, want, got)
