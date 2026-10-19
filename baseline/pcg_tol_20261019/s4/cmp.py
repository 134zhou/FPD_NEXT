"""比较 .fpd 粒子位置/速度：python3 cmp.py ref.fpd a.fpd [b.fpd ...]"""
import struct, sys, math
sys.path.insert(0, '/home/doll/Documents/C++/FPD_NEXT/tools')
import fpd_format as F

def read(path, names):
    h = F.read_header(path)
    offs = F.field_offsets(h)
    out = {}
    with open(path, 'rb') as f:
        for n in names:
            off, cnt = offs[n]
            f.seek(off)
            out[n] = struct.unpack('<%dd' % cnt, f.read(8 * cnt))
    return h, out

names = ['Rux', 'Ruy', 'Ruz', 'Vx', 'Vy', 'Vz']
_, ref = read(sys.argv[1], names)
N = len(ref['Rux'])
for p in sys.argv[2:]:
    _, d = read(p, names)
    dR = [math.sqrt(sum((d[c][i] - ref[c][i]) ** 2 for c in ('Rux', 'Ruy', 'Ruz'))) for i in range(N)]
    dV = [math.sqrt(sum((d[c][i] - ref[c][i]) ** 2 for c in ('Vx', 'Vy', 'Vz'))) for i in range(N)]
    vz = sum(d['Vz']) / N
    print('%-40s max|dR| %.3e  rms|dR| %.3e  max|dV| %.3e  <Vz> %.6e'
          % (p.split('/')[-1], max(dR), math.sqrt(sum(x * x for x in dR) / N), max(dV), vz))
