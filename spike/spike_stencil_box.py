#!/usr/bin/env python3
"""
模板盒尺寸的判决性检验（Phase 7-B）。

复现并验证修复：stencil_point 的模板盒原本是 [kn-range+1, kn+range]（宽 2*range），
但半径 range 的球心落在格胞【任意】位置时，能触及的整数层最多有 2*range+1 个 ——
负方向会静默漏掉一层。

本脚本做两件事：
  1. 对每个 Loc × 每个位置，把「模板盒遍历」与「无盒暴力枚举整个空间」对照，
     列出漏点；
  2. 用 n_range = 2*range+1、盒 [kn-range, kn+range] 重跑同一对照，确认零漏点。

纯 stdlib，任何 python3 可跑，秒级。

对照结论（修复前）：
  FACE_X  Rz=16.5  盒 2140 点 vs 暴力 2152 点  差 -9.64e-04
  FACE_Z  Rz=27.0  盒 1787 点 vs 暴力 1811 点  差 -2.79e-03
  EDGE_ZX Rz=16.5  盒 2132 点 vs 暴力 2142 点  差 -8.33e-04
  CELL 各位置均无漏（体心层的偏移是 0，只在整数位置、距离恰为 range 时才会漏）

这解释了 Phase 0 spike 里 ∫φ = 170.3051225223 与连续球坐标积分 170.306 的 5e-4
偏差 —— 那不是「格点离散化误差」，是漏点。修复后 ∫φ = 170.3060863408，差 4.7e-7。
"""
import math
radius, xi, Nx, Ny, Nz = 3.2, 1.0, 32, 32, 32
RANGE = int(2*(radius+xi)); RANGE_M1 = RANGE-1; NR = 2*RANGE; R2 = RANGE*RANGE
OFF = {"CELL":(0.0,0.0,0.0), "FACE_X":(0.5,0,0), "FACE_Y":(0,0.5,0),
       "FACE_Z":(0,0,0.5), "EDGE_XY":(0.5,0.5,0), "EDGE_YZ":(0,0.5,0.5), "EDGE_ZX":(0.5,0,0.5)}
# 壁面模式下的有效 kr 范围
VALID = {"CELL":(0,Nz-1), "FACE_X":(0,Nz-1), "FACE_Y":(0,Nz-1), "EDGE_XY":(0,Nz-1),
         "FACE_Z":(0,Nz-2), "EDGE_YZ":(-1,Nz-1), "EDGE_ZX":(-1,Nz-1)}

def order(dx,dy,dz):
    return 0.5*(math.tanh((radius - math.sqrt(dx*dx+dy*dy+dz*dz))*1.0)+1.0)

def sum_box(loc, Rx,Ry,Rz):
    ox,oy,oz = OFF[loc]; lo,hi = VALID[loc]
    inx,iny,inz = math.floor(Rx), math.floor(Ry), math.floor(Rz)
    s = 0.0; n = 0
    for li in range(NR):
        for lj in range(NR):
            for lk in range(NR):
                ir = li + inx - RANGE_M1; jr = lj + iny - RANGE_M1; kr = lk + inz - RANGE_M1
                dx = ir-Rx+ox; dy = jr-Ry+oy; dz = kr-Rz+oz
                if dx*dx+dy*dy+dz*dz > R2: continue
                if not (lo <= kr <= hi): continue
                s += order(dx,dy,dz); n += 1
    return s, n

def sum_brute(loc, Rx,Ry,Rz):
    ox,oy,oz = OFF[loc]; lo,hi = VALID[loc]
    W = RANGE + 3
    c = math.floor(Rz)
    s = 0.0; n = 0
    for ir in range(math.floor(Rx)-W, math.floor(Rx)+W+1):
        for jr in range(math.floor(Ry)-W, math.floor(Ry)+W+1):
            for kr in range(c-W, c+W+1):
                dx = ir-Rx+ox; dy = jr-Ry+oy; dz = kr-Rz+oz
                if dx*dx+dy*dy+dz*dz > R2: continue
                if not (lo <= kr <= hi): continue
                s += order(dx,dy,dz); n += 1
    return s, n

print("Rx=16.3 Ry=16.7（与实际自检一致）；比较『模板盒』与『无盒暴力枚举』")
for loc in ("CELL","FACE_X","FACE_Z","EDGE_ZX"):
    for Rz in (4.0, 27.0, 16.5):
        a,na = sum_box(loc, 16.3, 16.7, Rz)
        b,nb = sum_brute(loc, 16.3, 16.7, Rz)
        flag = "" if abs(a-b) < 1e-15 else "  <== 盒漏点"
        print("  %-8s Rz=%5.1f  盒=%12.7f (%3d 点)  暴力=%12.7f (%3d 点)  差=%9.2e%s"
              % (loc, Rz, a, na, b, nb, a-b, flag))

print()
print("=== 验证修法：n_range = 2*range+1，盒 = [in-range, in+range] ===")
NR2 = 2*RANGE+1
def sum_box2(loc, Rx,Ry,Rz):
    ox,oy,oz = OFF[loc]; lo,hi = VALID[loc]
    inx,iny,inz = math.floor(Rx), math.floor(Ry), math.floor(Rz)
    s=0.0; n=0
    for li in range(NR2):
        for lj in range(NR2):
            for lk in range(NR2):
                ir = li + inx - RANGE; jr = lj + iny - RANGE; kr = lk + inz - RANGE
                dx = ir-Rx+ox; dy = jr-Ry+oy; dz = kr-Rz+oz
                if dx*dx+dy*dy+dz*dz > R2: continue
                if not (lo <= kr <= hi): continue
                s += order(dx,dy,dz); n += 1
    return s,n
bad = 0
for loc in OFF:
    for Rz in (4.0, 27.0, 16.5, 16.0, 0.5, 31.0, 7.3):
        for Rx in (16.3, 16.0, 4.5):
            a,_ = sum_box2(loc, Rx, 16.7, Rz)
            b,_ = sum_brute(loc, Rx, 16.7, Rz)
            if abs(a-b) > 1e-14:
                bad += 1
                print("  仍漏: %s Rz=%.1f Rx=%.1f  差=%.2e" % (loc,Rz,Rx,a-b))
print("  全部覆盖检查：%s" % ("通过（无漏点）" if bad==0 else "%d 处仍有漏点" % bad))
