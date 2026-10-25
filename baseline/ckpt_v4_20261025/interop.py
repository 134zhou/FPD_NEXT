"""从仓库根目录运行：python3 baseline/ckpt_v4_20261025/interop.py。"""
from pathlib import Path
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
import fpd_format as F

DIR = Path(__file__).resolve().parent / 'interop'
DIR.mkdir(exist_ok=True)
names = ['vx', 'vy', 'vz', 'Rx', 'Ry', 'Rz', 'Rux', 'Ruy', 'Ruz', 'Vx', 'Vy', 'Vz']
for n in (0, 2):
    arrays = {}
    for t, name in enumerate(names):
        arrays[name] = ([i*1e6+j*1e3+k+t*.25 for k in range(2) for j in range(4) for i in range(8)]
                        if t < 3 else [t*100.0+q*.5 for q in range(n)])
    F.write_fpd(str(DIR / ('python_%d.fpd' % n)), F.default_header(8, 4, 2, n, 7), arrays, arrays)
raw = (DIR / 'python_2.fpd').read_bytes()
old = bytearray(raw[:-8] + bytes(48))
struct.pack_into('<I', old, 12, 3)
old += struct.pack('<Q', F._fnv_bytes(F.FNV_OFFSET, old))
corrupt = bytearray(raw)
corrupt[45] ^= 1
for name, data in [('old', old), ('corrupt', corrupt), ('truncated', raw[:-17]), ('trailing', raw+b'x')]:
    (DIR / (name+'.fpd')).write_bytes(data)
subprocess.run([str(ROOT/'build/fpd_check'), '--check-io', str(DIR)], check=True)
for n in (0, 2):
    cpp = DIR / ('cpp_%d.fpd' % n)
    header = F.read_header(str(cpp))
    assert header['N'] == n and header['step'] == 7
    assert cpp.stat().st_size == 48+24*64+72*n
    assert F.verify_checksum(str(cpp))[0]
    assert cpp.read_bytes() == (DIR / ('python_%d.fpd' % n)).read_bytes()
    assert set(F.field_offsets(header)) == set(names+['_checksum'])
    print('PASS Python 读取 C++ N=%d：文件长度、全部字节、校验和及字段集' % n)
try:
    F.read_header(str(DIR/'old.fpd'))
except ValueError:
    print('PASS Python 拒绝 v3')
else:
    raise AssertionError('接受了 v3')
for name in ('corrupt', 'truncated', 'trailing'):
    assert not F.verify_checksum(str(DIR/(name+'.fpd')))[0]
    print('PASS Python 校验拒绝 '+name)
subprocess.run([str(ROOT/'build/fpd_tool'), '--dump-ckpt', str(DIR/'cpp_0.fpd'), '--at', '3', '1', '1'], check=True)
