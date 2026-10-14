import numpy as np, time
from collections import defaultdict
from fractions import Fraction as F
exec(open('f4gamma.py').read().split("Ta = (6, F(9, 2)); s =")[0])
def Wa_factors(a, Tf):
    """exact: W_a = [s(-1-a)/s(T-1-a)] [s(a)/s(T+a)] [s(T-a)/s(-a)] [s(T+1+a)/s(1+a)],  s(c) = sinh(theta/2 + i pi c/2T)"""
    T_ = (F(Tf[0]), F(Tf[1])); a = (F(a[0]), F(a[1]))
    add = lambda *xs: (sum(x[0]*k for x, k in xs), sum(x[1]*k for x, k in xs))
    one = (F(0), F(1))
    num = [add((one, -1), (a, -1)), a, add((T_, 1), (a, -1)), add((T_, 1), (one, 1), (a, 1))]
    den = [add((T_, 1), (one, -1), (a, -1)), add((T_, 1), (a, 1)), add((a, -1),), add((one, 1), (a, 1))]
    return num, den
def s_val(c, T, th): return np.sinh(th/2 + 1j*np.pi*val(c, w0)/(2*T))
def run(label, Tf, ours_ex, cds_ex):
    t0 = time.time(); T = val(Tf, w0); bet = np.pi*(T - 1)/T
    # 1. verify the exact formula against the infinite-product evaluation
    devs = []
    for e in ((1, F(1, 2)), (2, F(3, 2)), (0, F(5, 2))):
        Tm = (F(Tf[0]) - e[0], F(Tf[1]) - e[1])
        Fp = Fmaker([val(e, w0), val(Tm, w0)], T, J=4000)
        for th in (0.31+0.17j, -0.4+0.5j):
            W_num = Fp(th + 1j*bet)*Fp(th)**2*Fp(th - 1j*bet)
            W_ex = 1.0
            for a in (e, Tm):
                n_, d_ = Wa_factors(a, Tf)
                W_ex *= np.prod([s_val(c, T, th) for c in n_])/np.prod([s_val(c, T, th) for c in d_])
            devs.append(abs(W_num/W_ex - 1))
    print(f"{label}: exact W_e formula vs infinite product (J = 4000): max deviation {max(devs):.1e}")
    # 2. exact lattice problem in sinh-factor exponents (positions canonical mod 2T; the sign is checked separately)
    def canon(c):
        a, b = F(c[0]), F(c[1])
        while val((a, b), w0) > T + 1e-9: a, b = a - 2*Tf[0], b - 2*F(Tf[1])
        while val((a, b), w0) <= -T + 1e-9: a, b = a + 2*Tf[0], b + 2*F(Tf[1])
        return (a, b)
    need = defaultdict(int)
    for y in cds_ex[0]: need[canon(y)] += 1; need[canon((-y[0], -y[1]))] -= 1
    for y in ours_ex[0]: need[canon(y)] -= 1; need[canon((-y[0], -y[1]))] += 1
    grid = [(al, F(b4, 4)) for al in range(0, 13) for b4 in range(-24, 25) if 0 < val((al, F(b4, 4)), w0) < T/2 + 1e-9]
    cols = []
    for e in grid:
        Tm = (F(Tf[0]) - e[0], F(Tf[1]) - e[1]); col = defaultdict(int)
        for a in (e, Tm):
            n_, d_ = Wa_factors(a, Tf)
            for c in n_: col[canon(c)] += 1
            for c in d_: col[canon(c)] -= 1
        cols.append({k: v for k, v in col.items() if v})
    keys = sorted(set(k for k, v in need.items() if v) | set(k for c in cols for k in c), key=lambda y: val(y, w0))
    M = np.array([[c.get(k, 0) for c in cols] for k in keys], float)
    rhs = np.array([need.get(k, 0) for k in keys], float)
    sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    res = np.abs(M @ sol - rhs).max()
    print(f"   {len(grid)} candidate pairs (quarter-step grid), {len(keys)} sinh-factor positions, rank {np.linalg.matrix_rank(M)}: best residual {res:.3f}  [{time.time()-t0:.0f}s]")
    return res
Ta = (6, F(9, 2)); s = (F(3), (F(9, 2) - 1)/2)
oursA = normalize([(y[0] + s[0], y[1] + s[1]) for y in SBa] + [(y[0] - s[0], y[1] - s[1]) for y in SBa], 1, Ta)
run("(a)", Ta, oursA, cds_exact(Ta, 1))
Tb = (6, F(-3, 2))
run("(b)", Tb, normalize(BBb, int(sgn_b), Tb), cds_exact(Tb, 0))
