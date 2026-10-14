import numpy as np
from collections import Counter
from fractions import Fraction as F
w0 = 2.37
def run(name, Tf, c378, ours_blocks, ours_sign):
    T = float(Tf[0])*w0 + float(Tf[1])
    R = lambda v: round(v, 7)
    def red(c):
        s = 1
        while c > T + 1e-9: c -= 2*T; s = -s
        while c <= -T + 1e-9: c += 2*T; s = -s
        return R(c), s
    class Amp:
        def __init__(self, num=(), den=(), sign=1):
            n_, d_, sg = Counter(), Counter(), sign
            for c in num: c2, s = red(c); n_[c2] += 1; sg *= s
            for c in den: c2, s = red(c); d_[c2] += 1; sg *= s
            common = n_ & d_; n_ -= common; d_ -= common
            self.num, self.den, self.sign = n_, d_, sg
        @staticmethod
        def blocks(ys, sign=1): return Amp(list(ys), [-y for y in ys], sign)
        def shift(self, s): return Amp([c + s for c in self.num.elements()], [c + s for c in self.den.elements()], self.sign)
        def __mul__(self, o): return Amp(list(self.num.elements()) + list(o.num.elements()), list(self.den.elements()) + list(o.den.elements()), self.sign*o.sign)
        def __eq__(self, o): return self.num == o.num and self.den == o.den and self.sign == o.sign
    H = 2*T/(w0 + c378); Bv = (H - 12)/3
    def brk(x): return [(x - 1, 1), (x + 1, 1), (x - 1 + Bv, -1), (x + 1 - Bv, -1)]
    ys = []
    for x in (1, H/3 + 1, H - 1, H - H/3 - 1):
        for z, s in brk(x): ys.append(z*T/H if s > 0 else -z*T/H)
    S11 = Amp.blocks(ys)
    yb = -Bv*T/H                                           # CDS pole 1+1 -> 1' (the coupling-dependent pole)
    S1p1 = S11.shift(yb/2)*S11.shift(-yb/2)                # S_{1',1}
    S1p1p = S1p1.shift(yb/2)*S1p1.shift(-yb/2)             # S_{1',1'}
    ours = Amp.blocks([float(a)*w0 + float(b) for a, b in ours_blocks], ours_sign)
    print(f"{name}: T = {T:.3f}, H = {H:.5f}, CDS pole 1+1->1' at y = {yb:.4f}")
    print(f"   our breather amplitude == CDS bootstrap S_(1',1'):  {ours == S1p1p}")
    print(f"   our breather amplitude == CDS S_11:                 {ours == S11}")
    return ours, S1p1p
SBa = [(0, F(3, 4)), (1, F(5, 4)), (3, F(7, 4)), (3, F(11, 4)), (5, F(13, 4)), (6, F(15, 4)), (-2, F(-3, 4)), (-4, F(-15, 4))]
Ta = (6, F(9, 2)); s = (F(3), (F(9, 2) - 1)/2)
BBa = [(y[0] + s[0], y[1] + s[1]) for y in SBa] + [(y[0] - s[0], y[1] - s[1]) for y in SBa]
run("(a)", Ta, 1, BBa, 1)
BBn, sgn_b, Tb, Ab = np.load("f4BB_0,-0.5,-1,-1.5.npy", allow_pickle=True)
run("(b)", (6, F(-3, 2)), 0, [(0, y) for y in BBn], int(sgn_b))
BBo, sgo, To, Ao = np.load("f4BB_1,1.5,1,1.5.npy", allow_pickle=True)
run("(original)", (6, F(3, 2)), 1, [(0, y) for y in BBo], int(sgo))
