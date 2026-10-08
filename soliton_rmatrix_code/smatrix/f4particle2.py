import numpy as np
from collections import Counter
BBs, BBsign, T, A0 = np.load("f4BB_1,1.5,1,1.5.npy", allow_pickle=True); BBs = list(BBs); T = float(T)
omega = 2.37; H = 2*T/(omega + 1)
R = lambda v: round(v, 7)
def red(c):
    s = 1
    while c > T + 1e-9: c -= 2*T; s = -s
    while c <= -T + 1e-9: c += 2*T; s = -s
    return R(c), s
class Amp:
    """product of sinh(theta/2 + i pi c/(2T)) factors: num / den, times sign"""
    def __init__(self, num=(), den=(), sign=1):
        n_, d_, sg = Counter(), Counter(), sign
        for c in num: c2, s = red(c); n_[c2] += 1; sg *= s
        for c in den: c2, s = red(c); d_[c2] += 1; sg *= s
        common = n_ & d_; n_ -= common; d_ -= common
        self.num, self.den, self.sign = n_, d_, sg
    @staticmethod
    def blocks(ys, sign=1): return Amp([y for y in ys], [-y for y in ys], sign)
    def shift(self, s): return Amp([c + s for c in self.num.elements()], [c + s for c in self.den.elements()], self.sign)
    def __mul__(self, o): return Amp(list(self.num.elements()) + list(o.num.elements()), list(self.den.elements()) + list(o.den.elements()), self.sign*o.sign)
    def __eq__(self, o): return self.num == o.num and self.den == o.den and self.sign == o.sign
    def val(self, th):
        f = lambda c: np.sinh(th/2 + 1j*np.pi*c/(2*T))
        return self.sign*np.prod([f(c) for c in self.num.elements()])/np.prod([f(c) for c in self.den.elements()])
def cds_blocks(a, b, H):
    Bv = (H - 12)/3
    def brk(x): return [(x - 1, 1), (x + 1, 1), (x - 1 + Bv, -1), (x + 1 - Bv, -1)]
    def sq(x): return brk(x) + brk(H - x)
    spec = {(1, 1): [1, H/3+1], (1, 2): [H/6+2, H/2+2], (1, 3): [2, H/3, H/3+2], (2, 2): [1, H/3-3, H/3+1, H/3+3],
            (2, 3): [H/6+1, H/6+3, H/2+1, H/2+3], (3, 3): [1, 3, H/3-1, H/3+1, H/3+1, H/3+3]}[(a, b)]
    ys = []
    for x in spec:
        for (z, s) in sq(x): ys.append(z*T/H if s > 0 else -z*T/H)
    return Amp.blocks(ys)
Y = lambda U: U*T/H
def sym(S, s): return S.shift(s)*S.shift(-s)
def closure(S11, lab, S12_ref=None, S13_ref=None):
    b111 = sym(S11, T/3) == S11
    S13 = sym(S11, Y(2)/2); S12a = sym(S11, Y(H/3 + 2)/2)
    ub_a, ub_b = T - Y(H/2 - 1), T - Y(5*H/6)          # 2 = 1 + 3 : S_{2,1}(th) = S_{1,1}(th + i ubar_{12}^3) S_{3,1}(th - i ubar_{32}^1)
    S12b = S11.shift(ub_a)*S13.shift(-ub_b)
    th = 0.37+0.21j
    out = f"{lab:24s}: 111 bootstrap {'holds' if b111 else 'FAILS'};  S_12 (1+1->2) == S_12 (1+3->2): {S12a == S12b}"
    if S12_ref is not None: out += f";  S_12 == CDS S_12: {S12a == S12_ref};  S_13 == CDS S_13: {S13 == S13_ref}"
    print(out)
    return S12a, S13
ours = Amp.blocks(BBs, BBsign); c11 = cds_blocks(1, 1, H)
print(f"T = {T:.3f}, H = {H:.5f}")
closure(c11, "CDS S_11 (control)", cds_blocks(1, 2, H), cds_blocks(1, 3, H))
o12, o13 = closure(ours, "soliton-derived S_11", cds_blocks(1, 2, H), cds_blocks(1, 3, H))
# how do the two differ?  ratio ours / CDS for S_11
ratio = ours*Amp(list(c11.den.elements()), list(c11.num.elements()), c11.sign)
print("S_11(ours)/S_11(CDS): numerator sinh-args", sorted(ratio.num.elements()), " denominator", sorted(ratio.den.elements()), " sign", ratio.sign)
