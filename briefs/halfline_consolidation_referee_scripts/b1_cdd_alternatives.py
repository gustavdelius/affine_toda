"""Do the CDD alternatives to Xi that reproduce Kim-Cho's one loop (b1_dn_scalar.py, part B) also satisfy
unitarity, crossing-unitarity and the bootstrap closure, and do they differ from K^N_1 = K_1(2-B) beyond one loop?"""
import mpmath as mp
from b1_dn_scalar import blkI, rhs_plusXi, Xi, cont, axioms, K1
mp.mp.dps = 30
half = mp.mpf(1)/2
def P(xf):
    return lambda m, T, o: blkI(T/2 - xf(T, o), m, T)*blkI(T/2 + xf(T, o), m, T)
x0, xw, xT = (lambda T, o: 0), (lambda T, o: o/2), (lambda T, o: (T - o)/2)
xh, xwh = (lambda T, o: half), (lambda T, o: half + o/2)
alts = {'P1/2 Pw/2 P(T-w)/2': (xh, xw, xT), 'P0 P1/2+w/2 P(T-w)/2': (x0, xwh, xT), 'P1/2 P1/2+w/2 P(T-w)/2': (xh, xwh, xT)}
pts = [mp.mpc('0.41', '0.17'), mp.mpc('-1.1', '0.33')]
ok = True
for name, xs in alts.items():
    def f(m, T, o, xs=xs):
        ma, mb = m - (T/2 - half), m + (T/2 - half)
        c = lambda mm: P(xs[0])(mm, T, o)*P(xs[1])(mm, T, o)*P(xs[2])(mm, T, o)
        return rhs_plusXi(m, T, o)/(Xi(ma, T, o)*Xi(mb, T, o))*c(ma)*c(mb)
    for n in (4, 5, 6):
        h = 2*n - 2
        B = mp.mpf('0.6')
        ax = axioms(lambda t: cont(f, t, n, B), n, B, pts)
        dual = max(abs(cont(f, t, n, B)/K1(t, h, 2 - B) - 1) for t in pts)
        print(f"{name:24s} n={n}: axioms residual {mp.nstr(ax, 2)}, |K/K_1(2-B) - 1| at B=0.6: {mp.nstr(dual, 3)}")
        ok &= ax < 1e-20 and dual > 1e-3
print("PASS" if ok else "FAIL", "alternative CDD factors satisfy all axioms, match one loop (d_4), and are not the dual of K_1")
