"""比较 v3/v4 的共同状态；旧受力只用于确认基准确实覆盖非零受力。"""
from pathlib import Path
import struct
p = Path(__file__).resolve().parent
for case, n in [('wall', 1), ('pair_disjoint', 2)]:
    for step in (0, 5, 10):
        name = '%s_%07d.fpd' % (case, step)
        old, new = (p/'before'/name).read_bytes(), (p/'after'/name).read_bytes()
        assert struct.unpack_from('<I', old, 12)[0] == 3
        assert struct.unpack_from('<I', new, 12)[0] == 4
        assert old[:12] == new[:12] and old[16:40] == new[16:40]
        assert len(old)-len(new) == 24*n
        off = 40
        for field in ['vx','vy','vz','Rx','Ry','Rz','Rux','Ruy','Ruz','Vx','Vy','Vz']:
            size = (64 if case == 'pair_disjoint' else 32)*32*32 if field[0].islower() else n
            a, b = old[off:off+8*size], new[off:off+8*size]
            assert a == b, (case, step, field, sum(x!=y for x,y in zip(a,b)))
            off += 8*size
        forces = struct.unpack_from('<%dd' % (3*n), old, off)
        assert max(map(abs, forces)) > 0
        print('PASS %s step=%d：12 个共同数组逐位一致，节省 %d 字节，旧 max|F|=%.9g' %
              (case, step, 24*n, max(map(abs, forces))))
a = (p/'after/wall_0000010.fpd').read_bytes()
b = (p/'restart/wall_0000010.fpd').read_bytes()
assert a == b
print('PASS 单粒子连续 10 步 vs 5+5 重启：整个文件逐位一致')

assert (p/'zero_v3.log').read_bytes() == (p/'zero_v4.log').read_bytes()
assert not list((p/'zero_v3').glob('*.fpd')) and not list((p/'zero_v4').glob('*.fpd'))
print('PASS 零步运行：均提示已达到目标，不生成检查点')
for case, step in [('failed', 1)]:
    name = 'wall_%07d.fpd' % step
    old = (p/(case+'_v3')/name).read_bytes()
    new = (p/(case+'_v4')/name).read_bytes()
    assert old[16:40] == new[16:40] and old[40:-32] == new[40:-8]
    assert len(old)-len(new) == 24
    print('PASS %s 落盘：共同状态逐位一致，v4 无受力' % case)
