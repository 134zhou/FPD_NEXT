#!/usr/bin/env python3
"""将 v4 checkpoint 复制为 step=0 的初态，保留全部场与粒子数据。"""

import argparse
import os
from pathlib import Path
import struct

from fpd_format import Fnv, HEADER_FIELDS, field_offsets, read_header


def checkpoint_to_init(source, destination):
    header = read_header(source)
    checksum_offset = field_offsets(header)["_checksum"][0]
    # 偏移从格式字段表推导，避免另存一份头部布局。
    step_offset = 12
    for name, fmt in HEADER_FIELDS:
        if name == "step":
            break
        step_offset += struct.calcsize(fmt)

    with open(source, "rb") as src:
        if os.fstat(src.fileno()).st_size != checksum_offset + 8:
            raise ValueError("文件大小与 v4 头部不符")
        # 独占创建，禁止覆盖源文件或已有输出。
        with open(destination, "xb") as dst:
            try:
                old_hash, new_hash = Fnv(), Fnv()
                original = src.read(header["_data_offset"])
                modified = bytearray(original)
                struct.pack_into("<q", modified, step_offset, 0)
                old_hash.feed(original)
                new_hash.feed(modified)
                dst.write(modified)
                remaining = checksum_offset - len(original)
                while remaining:
                    block = src.read(min(1 << 20, remaining))
                    if not block:
                        raise ValueError("文件数据被截断")
                    old_hash.feed(block)
                    new_hash.feed(block)
                    dst.write(block)
                    remaining -= len(block)
                tail = src.read(8)
                if len(tail) != 8 or src.read(1):
                    raise ValueError("文件尾部长度不符")
                if struct.unpack("<Q", tail)[0] != old_hash.h:
                    raise ValueError("原文件 FNV 校验失败")
                dst.write(struct.pack("<Q", new_hash.h))
                dst.flush()
            except BaseException:
                # 不留下看似可用的半成品。
                Path(destination).unlink()
                raise
    return header["step"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", help="输入 v4 checkpoint (.fpd)")
    parser.add_argument("destination", help="输出初态 (.fpd)，须为新文件")
    args = parser.parse_args()
    print("正在复制并校验文件（大文件的 FNV 校验需要时间）…", flush=True)
    try:
        step = checkpoint_to_init(args.source, args.destination)
    except (OSError, ValueError, struct.error) as exc:
        parser.exit(1, "错误：%s\n" % exc)
    print("%s -> %s：step %d -> 0，场与粒子数据保持不变"
          % (args.source, args.destination, step))


if __name__ == "__main__":
    main()
