"""生成小型回归初态；仅测试使用的 v3 样本，不是用户数据转换工具。"""
from pathlib import Path
import struct
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import fpd_format as F
p = Path(__file__).resolve().parent
for case, points in [('wall', [(16.3,16.7,4.7)]), ('pair', [(12.,16.,16.),(18.,16.,16.)]), ('pair_disjoint', [(12.,16.,16.),(32.,16.,16.)])]:
    parts = dict(zip(('Rx','Ry','Rz'), map(list, zip(*points))))
    parts.update({u:parts[r] for u,r in zip(('Rux','Ruy','Ruz'),('Rx','Ry','Rz'))})
    path = p / (case+'_v4.fpd')
    F.write_fpd(str(path), F.default_header(64 if case == "pair_disjoint" else 32,32,32,len(points)), {}, parts)
    body = bytearray(path.read_bytes()[:-8] + bytes(24*len(points)))
    struct.pack_into('<I', body, 12, 3)
    (p/(case+'_v3.fpd')).write_bytes(body+struct.pack('<Q', F._fnv_bytes(F.FNV_OFFSET, body)))
