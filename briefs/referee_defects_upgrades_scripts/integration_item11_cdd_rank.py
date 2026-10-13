"""Item 11 (ch. 43): (a) the simpler condition f(D) f(D + i pi) = 1 versus the string condition eq-dpt-soliton-cdd,
prod_{k<h} f(D - i pi + 2 pi i k/h) = 1: sufficient for even h (pairs k, k + h/2 differ by i pi), not necessary
(T^cl_1), not sufficient for odd h (sinh example); (b) the residue of T_Vbar at X^h = q is rank one on each charge
sector, M_Vbar = q u w^T exactly (sympy, symbolic q-powers), so on the infinite direct sum of sectors it has infinite
rank.  Written from the chapter's definitions (eq-dpt-particle, eq-dt-tvbar); no book or referee code."""
import numpy as np
import sympy as sp

def string(f, D, h):
    return np.prod([f(D - 1j*np.pi + 2j*np.pi*k/h) for k in range(h)], axis=0)

def Tcl1(D, h):          # eq-dpt-particle, a = 1: (e^{-D} + e^{i pi/h})/(e^{-D} + e^{-i pi/h})
    return (np.exp(-D) + np.exp(1j*np.pi/h))/(np.exp(-D) + np.exp(-1j*np.pi/h))

sinhf = lambda D: (np.sinh(D) + 0.3j)/(np.sinh(D) - 0.3j)
Ds = np.array([0.41, -0.7 + 0.2j, 1.3 - 0.1j, 2.2 + 0.05j])
ok = True
print('(a) |string - 1| and |f(D) f(D + i pi) - 1|:')
for h in range(2, 9):
    st_s = abs(string(sinhf, Ds, h) - 1).max()
    st_T = abs(string(lambda D: Tcl1(D, h), Ds, h) - 1).max()
    sim_T = abs(Tcl1(Ds, h)*Tcl1(Ds + 1j*np.pi, h) - 1).max()
    print(f'   h = {h}: sinh-type string {st_s:.2e}; T^cl_1 string {st_T:.2e}, T^cl_1 simple {sim_T:.2e}')
    ok &= st_T < 1e-12
    if h % 2 == 0:
        ok &= st_s < 1e-12
    if h > 2:
        ok &= sim_T > 0.1
    if h in (3, 5):
        ok &= st_s > 0.1
sim_s = abs(sinhf(Ds)*sinhf(Ds + 1j*np.pi) - 1).max()
print(f'   sinh-type simple condition: {sim_s:.1e}; its string value at h = 3, D = 0.41: {string(sinhf, np.array([0.41]), 3)[0]:.10f}')
ok &= sim_s < 1e-14
# even h, by hand: pairing k <-> k + h/2 (arguments differ by i pi) -- check the pairing identity symbolically
for h in (4, 6):
    args = [sp.Symbol('D') - sp.I*sp.pi + 2*sp.I*sp.pi*sp.Rational(k, h) for k in range(h)]
    ok &= all(sp.simplify(args[k + h//2] - args[k] - sp.I*sp.pi) == 0 for k in range(h//2))
# (b) rank one at X^h = q on a sector of Vbar (eq-dt-tvbar): entry (i+j <- i) = X^j q prod_{l=i}^{i+j} z_l, z_l = q^{(Qb+s)_l},
# prod_l z_l = q^{-1} (sum Qb = -1, sum s = 0)
for h in (3, 4, 5):
    X = sp.Symbol('X', nonzero=True)
    zs = sp.symbols(f'z0:{h-1}', nonzero=True)
    q = X**h                                   # the species-r pole X^h = q
    z = list(zs) + [1/(q*sp.prod(zs))]
    M = sp.zeros(h, h)
    for i in range(h):
        for j in range(h):
            M[(i + j) % h, i] = X**j*q*sp.prod([z[(i + l) % h] for l in range(j + 1)])
    minors_zero = all(sp.simplify(M[a, c]*M[b, d] - M[a, d]*M[b, c]) == 0
                      for a in range(h) for b in range(h) for c in range(h) for d in range(h))
    nonzero = all(sp.simplify(M[a, c]) != 0 for a in range(h) for c in range(h))
    print(f'(b) h = {h}: all 2x2 minors of the Vbar sector block vanish at X^h = q: {minors_zero}; all entries nonzero: {nonzero}')
    ok &= minors_zero and nonzero
print('PASS' if ok else 'FAIL', ': item 11 - even h: simpler condition sufficient (not necessary); odd h: not sufficient; rank 1 per sector')
