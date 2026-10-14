import numpy as np
from collections import defaultdict
from fractions import Fraction as F
exec(open('f4gamma5.py').read().split("def run(label")[0])
def moments(label, Tf, ours_ex, cds_ex):
    T = val(Tf, w0); P = (2*Tf[0], 2*F(Tf[1]))          # period 2T = P[0] w + P[1]
    def canon(c):                                         # canonical rep: w-coefficient in [-Tf0, Tf0)
        a, b = F(c[0]), F(c[1])
        while a >= Tf[0]: a, b = a - P[0], b - P[1]
        while a < -Tf[0]: a, b = a + P[0], b + P[1]
        return (a, b)
    m = defaultdict(int)                                  # needed exponents of s(c) in W_needed = CDS / ours ; block (y) = s(y)/s(-y)
    for y in cds_ex[0]: m[canon(y)] += 1; m[canon((-y[0], -y[1]))] -= 1
    for y in ours_ex[0]: m[canon(y)] -= 1; m[canon((-y[0], -y[1]))] += 1
    lines = defaultdict(list)
    for (a, b), v in m.items():
        if v: lines[(a, b - (b // 1))].append((b, v))   # line: fixed w-coefficient a and b mod 1 (step-1 shifts move along b)
    print(f"{label}: T = {Tf[0]}w + {Tf[1]};  needed correction has {sum(1 for v in m.values() if v)} nonzero sinh-factor exponents on {len(lines)} lines")
    kinds = defaultdict(int)
    for key, pts in sorted(lines.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        m0 = sum(v for _, v in pts); m1 = sum(b*v for b, v in pts)
        kind = "finite (CDD)" if m0 == 0 and m1 == 0 else ("Gamma-type tower" if m0 == 0 else "Barnes double-Gamma (linear growth)")
        kinds[kind] += 1
        print(f"   line a = {str(key[0]):>3s}, b = {str(key[1]):>4s} mod 1:  sum m = {m0:+d},  sum b m = {str(m1):>5s}   ->  {kind}")
    print("   summary:", dict(kinds))
Ta = (6, F(9, 2)); s = (F(3), (F(9, 2) - 1)/2)
oursA = normalize([(y[0] + s[0], y[1] + s[1]) for y in SBa] + [(y[0] - s[0], y[1] - s[1]) for y in SBa], 1, Ta)
moments("(a)", Ta, oursA, cds_exact(Ta, 1))
Tb = (6, F(-3, 2))
moments("(b)", Tb, normalize(BBb, int(sgn_b), Tb), cds_exact(Tb, 0))
