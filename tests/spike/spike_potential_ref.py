#!/usr/bin/env python3
# spike：独立实现三种势的 U(r) 与 g(r) = -(1/r)dU/dr，产出黄金表。
#
# 这份代码【不看 C++ 源码】，只按势函数的解析公式独立实现，
# 用作 CheckPotential.cpp 里 J2 判据的黄金值。这正是项目的一贯做法：
# spike_template_routine.cpp 曾用独立 Python 复核 ∫φ = 170.3051225223 到全位一致。
#
# 纯 stdlib，任何 python3 都能跑：
#   python3 tests/spike/spike_potential_ref.py > tests/spike/potential_golden.txt
#
# 约定：g(r) 是标量因子，C++ 侧 F_i = g * (R_i - R_j)。g > 0 表示排斥。
#       能量用 shift=energy（默认），即 U(r) = U_bare(r) - U_bare(rcut)，U(rcut)=0。

import math

def lj_bare_U(r, eps, sigma):
    s6 = (sigma / r) ** 6
    return 4.0 * eps * (s6 * s6 - s6)

def lj_bare_g(r, eps, sigma):
    # g = -(1/r)dU/dr = 24 eps / r^2 * s^6 (2 s^6 - 1)
    s6 = (sigma / r) ** 6
    return 24.0 * eps * s6 * (2.0 * s6 - 1.0) / (r * r)

def morse_bare_U(r, De, alpha, r_eq):
    E = math.exp(-alpha * (r - r_eq))
    return De * (E * E - 2.0 * E)

def morse_bare_g(r, De, alpha, r_eq):
    # g = -(1/r)dU/dr = 2 alpha De E (E-1) / r
    E = math.exp(-alpha * (r - r_eq))
    return 2.0 * alpha * De * E * (E - 1.0) / r

def wca_rcut(sigma):
    return 2.0 ** (1.0 / 6.0) * sigma

def emit(name, r, U, g):
    # %.17g 保证与 IEEE double 往返无损
    print("%-6s %.17g %.17g %.17g" % (name, r, U, g))

def table_lj(name, eps, sigma, rcut, rs):
    U_rc = lj_bare_U(rcut, eps, sigma) if rcut > 0 else 0.0
    for r in rs:
        if r >= rcut:
            emit(name, r, 0.0, 0.0)
        else:
            emit(name, r, lj_bare_U(r, eps, sigma) - U_rc, lj_bare_g(r, eps, sigma))

def table_morse(De, alpha, r_eq, rcut, rs):
    U_rc = morse_bare_U(rcut, De, alpha, r_eq)
    for r in rs:
        if r >= rcut:
            emit("morse", r, 0.0, 0.0)
        else:
            emit("morse", r, morse_bare_U(r, De, alpha, r_eq) - U_rc,
                 morse_bare_g(r, De, alpha, r_eq))

def main():
    # 采样点：覆盖 [0.6*sigma, 1.4*rcut] 两侧，含特征点 r_eq 与 2^{1/6}*sigma
    rs = [2.0, 3.0, 4.0, 5.0, 6.0, 6.5, 7.0, 7.4, 8.0, 8.30621915748936,
          9.0, 10.0, 12.0, 14.0, 15.0, 16.0, 18.0, 20.0]

    print("# spike_potential_ref 黄金表   （r, U, g）  g = -(1/r)dU/dr，shift=energy")
    print("# 参数见各段注释；%.17g 与 IEEE double 往返无损")

    # WCA：eps=1, sigma=7.4（简单参数，供 N=3 解析对照 J4）
    print("# wca  eps=1 sigma=7.4  rcut=2^(1/6)*sigma=8.30621915748936")
    rcut_wca = wca_rcut(7.4)
    table_lj("wca", 1.0, 7.4, rcut_wca, rs)

    # LJ 12-6：旧代码参数 eps=57.1428571429, sigma=7.4, rcut=15
    print("# lj126  eps=57.1428571429 sigma=7.4  rcut=15")
    table_lj("lj126", 57.1428571429, 7.4, 15.0, rs)

    # Morse：De=50, alpha=1, r_eq=7.4, rcut=15
    print("# morse  De=50 alpha=1 r_eq=7.4  rcut=15")
    table_morse(50.0, 1.0, 7.4, 15.0, rs)

    # 手算锚点（供肉眼核对，也写进 README）
    print()
    print("# 锚点：lj126(eps=57.1428571429,sigma=7.4,r=8) U=%.6f |F|=g*r=%.6f" %
          (lj_bare_U(8.0, 57.1428571429, 7.4) - lj_bare_U(15.0, 57.1428571429, 7.4),
           lj_bare_g(8.0, 57.1428571429, 7.4) * 8.0))
    print("# 锚点：morse(De=50,alpha=1,r_eq=7.4,r=8) g=%.6f (吸引, g<0)" %
          morse_bare_g(8.0, 50.0, 1.0, 7.4))

if __name__ == "__main__":
    main()
