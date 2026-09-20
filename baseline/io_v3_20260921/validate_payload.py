#!/usr/bin/env python3
"""比较 v2/v3 共用的数值载荷，并检查 v3 逐位重启。"""
from pathlib import Path
import struct

base = Path(__file__).resolve().parent
old_header = 104  # v2: magic/endian/version/flags/dims/step/物理量/RNG
new_header = 40   # v3: magic/endian/version/dims/step
size = 32 * 32 * 32 * 8
for name in ("bulk", "wall"):
    for step in (0, 50):
        old = (base / "before" / f"{name}_{step:07d}.fpd").read_bytes()
        new = (base / "after" / f"{name}_{step:07d}.fpd").read_bytes()
        assert struct.unpack_from("<I", old, 12)[0] == 2
        assert struct.unpack_from("<I", new, 12)[0] == 3
        assert old[old_header:old_header + 3 * size] == new[new_header:new_header + 3 * size]
        # v2 检查点含压力，v3 删除该数组；其后是同序粒子数组。
        assert old[old_header + 4 * size:-8] == new[new_header + 3 * size:-8]
        print(name, step, "共用载荷逐位相同")
    continuous = (base / "after" / f"{name}_0000050.fpd").read_bytes()
    restart = (base / "after" / f"{name}_restart_0000050.fpd").read_bytes()
    assert continuous == restart
    print(name, "50 步与 25+25 步重启整文件逐位相同")
