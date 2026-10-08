import numpy as np
from fractions import Fraction as F
exec(open('f4cand.py').read().split("def test(name")[0])
w0 = 2.37
def red_ab(y, Tf):
    a, b = F(y[0]), F(y[1]); TT = (F(Tf[0]), F(Tf[1])); s = 1
    while val((a, b), w0) > val(TT, w0) + 1e-9: a, b = a - 2*TT[0], b - 2*TT[1]; s = -s
    while val((a, b), w0) <= -val(TT, w0) + 1e-9: a, b = a + 2*TT[0], b + 2*TT[1]; s = -s
    return (a, b), s
def normalize(ys, sgn, Tf):
    out = []; s_ = sgn
    for y in ys:
        y2, s = red_ab(y, Tf); s_ *= s
        if y2 == (0, 0): continue
        if y2 == (F(Tf[0]), F(Tf[1])): s_ *= -1; continue
        out.append(y2)
    res = []
    for y in out:
        m = [z for z in res if z == (-y[0], -y[1])]
        if m: res.remove(m[0])
        else: res.append(y)
    return sorted(res, key=lambda y: val(y, w0)), s_
def cds_exact(Tf, c378):
    """CDS S_11 at H = 2T/(w + c378), in exact y = A w + B form (y = x T/H)"""
    # x T/H = x (w + c378)/2 ;  B = (H - 12)/3 -> in y units B*T/H = (2T - 12(w + c378))/6
    def y_of(cx, dB):   # x = cx + dB*B (cx may contain H: handled by caller)
        return None
    # blocks of S_11 = {1}{H-1}{H/3+1}{2H/3-1} with {x} = (x-1)(x+1)/((x-1+B)(x+1-B))
    Hy = (F(Tf[0]), F(Tf[1]))                                  # H in y units is T
    one = (F(1, 2), F(c378, 2))                                # x = 1  -> y = (w + c378)/2
    By = ((F(2)*Hy[0] - 12*one[0]*2)/6, (F(2)*Hy[1] - 12*one[1]*2)/6)   # B in y units = (2T - 12*(w+c378)/... )/6
    # recompute B_y properly: B = (H - 12)/3 ; y(B) = B * T/H = (T - 12 T/H)/3 = (T - 6 (w + c378))/3
    By = ((Hy[0] - 6)/3, (Hy[1] - 6*F(c378))/3)
    def add(*terms):
        a = sum(t[0]*k for t, k in terms); b = sum(t[1]*k for t, k in terms); return (a, b)
    ys = []
    for (base, sgn_) in ((add((one, 1)), 1), (add((Hy, 1), (one, -1)), 1), (add((Hy, F(1, 3)), (one, 1)), 1), (add((Hy, F(2, 3)), (one, -1)), 1)):
        ys += [add((base, 1), (one, -1)), add((base, 1), (one, 1))]
        ys += [(-v[0], -v[1]) for v in (add((base, 1), (one, -1), (By, 1)), add((base, 1), (one, 1), (By, -1)))]
    return normalize(ys, 1, Tf)
SBa = [(0, F(3, 4)), (1, F(5, 4)), (3, F(7, 4)), (3, F(11, 4)), (5, F(13, 4)), (6, F(15, 4)), (-2, F(-3, 4)), (-4, F(-15, 4))]
def exact(y, w=2.37):
    for a in range(-24, 25):
        b = y - a*w
        if abs(b*4 - round(b*4)) < 1e-6 and abs(b) < 30: return (a, F(round(b*4), 4))
    raise ValueError(y)
BBn, sgn_b, Tb, Ab = np.load("f4BB_0,-0.5,-1,-1.5.npy", allow_pickle=True)
BBb = [exact(y) for y in BBn]
show = lambda y: f"{y[0]}w{'+' if y[1] >= 0 else ''}{y[1]}"
for name, SB, sgn, Tf, c in (("(a)", SBa, 1, (6, F(9, 2)), 1), ("(b)", None, None, (6, F(-3, 2)), 0)):
    if SB is not None:
        s = (F(Tf[0], 2), (F(Tf[1]) - 1)/2)
        ours = normalize([(y[0] + s[0], y[1] + s[1]) for y in SB] + [(y[0] - s[0], y[1] - s[1]) for y in SB], 1, Tf)
    else:
        ours = normalize(BBb, int(sgn_b), Tf)
    cds = cds_exact(Tf, c)
    ratio = normalize(list(cds[0]) + [(-y[0], -y[1]) for y in ours[0]], cds[1]*ours[1], Tf)
    print(f"{name}  T = {Tf[0]}w + {Tf[1]},  H = 2T/(w + {c})")
    print(f"   ours: {[show(y) for y in ours[0]]}")
    print(f"   CDS : {[show(y) for y in cds[0]]}")
    print(f"   CDS/ours = sign {ratio[1]:+d} x blocks {[show(y) for y in ratio[0]]}")
