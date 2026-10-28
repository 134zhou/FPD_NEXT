"""从仓库根目录运行：python3 baseline/checkpoint_to_init_20261028/check.py。"""
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from checkpoint_to_init import checkpoint_to_init
from fpd_format import default_header, write_fpd, read_header, verify_checksum, PARTICLE_ARRAYS


with tempfile.TemporaryDirectory() as directory:
    work = Path(directory)
    for n, step in [(0, 12345), (2, 9876543210), (2, 0)]:
        src, dst = work / "source.fpd", work / "init.fpd"
        fields = {name: [float(i + shift) for i in range(64)]
                  for name, shift in [("vx", 0), ("vy", 100), ("vz", 200)]}
        particles = {name: [float(i + j) for i in range(n)]
                     for j, name in enumerate(PARTICLE_ARRAYS)}
        write_fpd(src, default_header(8, 4, 2, n, step), fields, particles)
        original = src.read_bytes()
        subprocess.run([sys.executable, str(ROOT / "tools/checkpoint_to_init.py"),
                        str(src), str(dst)], check=True)
        result = dst.read_bytes()
        assert src.read_bytes() == original
        assert result[:32] == original[:32] and result[40:-8] == original[40:-8]
        assert len(result) == len(original)
        assert read_header(dst)["step"] == 0 and verify_checksum(dst)[0]
        dump = subprocess.run([str(ROOT / "build/fpd_tool"), "--dump-ckpt", str(dst),
                               "--at", "3", "1", "1"], check=True,
                              capture_output=True, text=True).stdout
        print(dump, end="")
        for target in (src, dst):
            before = target.read_bytes()
            try:
                checkpoint_to_init(src, target)
            except FileExistsError:
                pass
            else:
                raise AssertionError("覆盖已有文件")
            assert target.read_bytes() == before
        dst.unlink()
        print("PASS N=%d step=%d 字节一致、校验和、C++ 读取、覆盖保护" % (n, step))

    for label, damaged in [
        ("corrupt", original[:40] + bytes([original[40] ^ 1]) + original[41:]),
        ("truncated", original[:-1]),
        ("trailing", original + b"x"),
        ("v3", original[:12] + (3).to_bytes(4, "little") + original[16:]),
    ]:
        src.write_bytes(damaged)
        try:
            checkpoint_to_init(src, dst)
        except ValueError:
            pass
        else:
            raise AssertionError(label)
        assert not dst.exists()
        print("PASS 拒绝 " + label)

    # C++ 写出的真实中间 checkpoint 也必须能转换并回读。
    src = ROOT / "baseline/ckpt_v4_20261025/after/pair_disjoint_0000010.fpd"
    checkpoint_to_init(src, dst)
    original, result = src.read_bytes(), dst.read_bytes()
    assert original[:32] == result[:32] and original[40:-8] == result[40:-8]
    assert read_header(dst)["step"] == 0 and verify_checksum(dst)[0]
    subprocess.run([str(ROOT / "build/fpd_tool"), "--dump-ckpt", str(dst)], check=True)
    print("PASS C++ checkpoint -> Python 转换 -> C++ 读取")
