"""eq-hirota-bilinear for non-simply-laced / twisted algebras, by folding explicit parent solutions.
Book's equation (10-solitons.qmd:48):
  (1/2)(D_x^2 - D_t^2) tau_j.tau_j = m^2 n_j (alpha_j^2/2) (tau_j^2 - prod_{l!=j} tau_l^{-a_lj}),  a_lj = 2 alpha_l.alpha_j/alpha_l^2,
with the ansatz phi = (i/beta) sum_j (2 alpha_j/alpha_j^2) ln tau_j.  Folding (sec-affine-folding): alpha_O = mean of orbit,
label n_O = |O| n_j.  Matching the ansatz: tau_O = tau_j^{|O| alpha_O^2/2}  (=tau_j for orthogonal orbits)."""
import numpy as np, mpmath as mp
from lie import Algebra
from hirota import soliton, MpCtx
mp.mp.dps = 40
def report(n, ok): print(('PASS ' if ok else 'FAIL ') + n)

def fold_data(alpha, n, perm):
    """alpha: (r+1) x r parent affine roots, n labels, perm: node permutation. returns orbits, alpha_O, n_O, a_PO"""
    seen = set(); orbits = []
    for j in range(len(n)):
        if j in seen: continue
        O = [j]; k = perm[j]
        while k != j: O.append(k); k = perm[k]
        seen |= set(O); orbits.append(O)
    aO = np.array([alpha[O].mean(0) for O in orbits])
    nO = np.array([len(O)*n[O[0]] for O in orbits])
    Cf = np.array([[2*aO[P] @ aO[Q]/(aO[P] @ aO[P]) for Q in range(len(orbits))] for P in range(len(orbits))])
    assert np.allclose(nO @ aO, 0)
    from fractions import Fraction
    expo = [mp.mpf(Fraction(float(len(O)*(aO[i] @ aO[i])/2)).limit_denominator(100).numerator)/Fraction(float(len(O)*(aO[i] @ aO[i])/2)).limit_denominator(100).denominator for i, O in enumerate(orbits)]
    return orbits, aO, nO, np.round(Cf, 10), expo

def residual(tauO, aO, nO, Cf, pts, mm=1):
    """max relative residual of the book's bilinear equation; tauO(x,t) -> list"""
    worst = 0
    for x0, t0 in pts:
        T = tauO(x0, t0)
        for J in range(len(nO)):
            f = lambda x, t: tauO(x, t)[J]
            tx = mp.diff(lambda x: f(x, t0), x0); txx = mp.diff(lambda x: f(x, t0), x0, 2)
            tt = mp.diff(lambda t: f(x0, t), t0); ttt = mp.diff(lambda t: f(x0, t), t0, 2)
            tj = T[J]
            lhs = (tj*txx - tx**2) - (tj*ttt - tt**2)            # (1/2)(D_x^2 - D_t^2) tau.tau
            prod = mp.mpf(1)
            for L in range(len(nO)):
                if L != J: prod *= T[L]**(-int(round(Cf[L, J])))
            from fractions import Fraction
            cf = Fraction(float(nO[J]*(aO[J] @ aO[J])/2)).limit_denominator(100)   # exact n_j alpha_j^2/2
            rhs = mm**2*mp.mpf(cf.numerator)/cf.denominator*(tj**2 - prod)
            worst = max(worst, abs(lhs - rhs)/(abs(tj)**2 + abs(prod)))
    return worst

pts = [(mp.mpf('0.3'), mp.mpf('-0.2')), (mp.mpf('-0.7'), mp.mpf('0.5'))]

# ---------- c_2^(1) from a_3^(1): reflection j -> -j (fixes 0, 2)
alpha = np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [-1, 0, 0, 1]], float)  # alpha_j = e_j - e_{j+1}, j=0..3 (in R^4)
n = np.array([1, 1, 1, 1]); perm = [0, 3, 2, 1]
orbits, aO, nO, Cf, expo = fold_data(alpha, n, perm)
print('c2 folding: orbits', orbits, 'alpha_O^2', [round(float(a @ a), 6) for a in aO], 'labels', nO.tolist())
print('   folded Cartan a_PO =\n', Cf)
h = 4; om = mp.e**(2j*mp.pi/h)
def an_tau(spec, ths, xis, h):
    om = mp.e**(2j*mp.pi/h)
    def taus(x, t):
        Es = [mp.e**(2*mp.sin(mp.pi*a/h)*(x*mp.cosh(th) - t*mp.sinh(th)) + xi) for a, th, xi in zip(spec, ths, xis)]
        out = []
        for j in range(h):
            s = 0
            for k in range(len(spec)+1):
                import itertools
                for S in itertools.combinations(range(len(spec)), k):
                    term = 1
                    for i in S: term *= om**(j*spec[i])*Es[i]
                    for i1, i2 in itertools.combinations(S, 2):
                        A_, B_ = mp.pi*spec[i1]/h, mp.pi*spec[i2]/h; c = mp.cosh(ths[i1]-ths[i2])
                        term *= (c - mp.cos(A_-B_))/(c - mp.cos(A_+B_))
                    s += term
            out.append(s)
        return out
    return taus
th = mp.mpf('0.6'); xi = mp.mpc('0.2', '0.4')
for name, spec, ths, xis in [('c2 soliton = coincident pair (1,3) of a3', [1, 3], [th, th], [xi, xi]),
                             ('c2 soliton = species 2 of a3', [2], [th], [xi])]:
    par = an_tau(spec, ths, xis, h)
    tauO = lambda x, t, par=par: [par(x, t)[O[0]]**expo[i] for i, O in enumerate(orbits)]
    chk = max(abs(par(*p)[1] - par(*p)[3]) for p in pts)
    res = residual(tauO, aO, nO, Cf, pts)
    report(f'{name}: tau_1=tau_3 ({float(chk):.0e}); book bilinear eq for c2 (labels 1,2,1; alpha^2 2,1,2), moving, res {float(res):.1e}', res < 1e-25)

# ---------- a_2^(2) from a_2^(1) by j -> -j: non-orthogonal orbit {1,2}, tau_O = tau_1^(1/2)
alpha = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]], float); n = np.array([1, 1, 1]); perm = [0, 2, 1]
orbits, aO, nO, Cf, expo = fold_data(alpha, n, perm)
print('a2^(2) folding: orbits', orbits, 'alpha_O^2', [round(float(a @ a), 6) for a in aO], 'labels', nO.tolist(), 'tau exponents', expo)
print('   folded Cartan a_PO =\n', Cf)
par = an_tau([1, 2], [th, th], [xi, xi], 3)
tauO = lambda x, t: [par(x, t)[O[0]]**expo[i] for i, O in enumerate(orbits)]
res = residual(tauO, aO, nO, Cf, pts)
report(f'a2^(2) (Kac labels 1 on long, 2 on short): coincident (1,2) pair of a2, tau_O = sqrt(tau_1) = 1 - omega-stuff; res {float(res):.1e}', res < 1e-25)
E = mp.mpf('0.37'); t1 = 1 - E + E**2/4
print('   check: a2 pair tau_1 = (1 - E/2)^2 :', mp.nstr(an_tau([1, 2], [0, 0], [mp.log(E), mp.log(E)], 3)(0, 0)[1] - (1 - E/2)**2, 5))

# ---------- d_4^(1) foldings via static Hirota recursion (foundations_code/hirota/hirota.py)
Al = Algebra('d', 4); ctx = MpCtx(40)
def static_tau(delta, mu):
    t, tail = soliton(Al, delta, mu, ctx=ctx, Kmax=12, tol=mp.mpf(10)**-30)
    return t
def run(name, perm, delta, mu):
    orbits, aO, nO, Cf, expo = fold_data(Al.alpha, Al.n.astype(int), perm)
    t = static_tau(delta, mu)
    degs = [len(c)-1 for c in t]
    xi0 = mp.mpc('0.1', '0.3')
    def par(x, tt):
        E = mp.e**(mu*x + xi0)
        return [sum(c*E**p for p, c in enumerate(cs)) for cs in t]
    inv = max(abs(par(mp.mpf('0.4'), 0)[j] - par(mp.mpf('0.4'), 0)[perm[j]]) for j in range(5))
    tauO = lambda x, tt: [par(x, tt)[O[0]]**expo[i] for i, O in enumerate(orbits)]
    res = residual(tauO, aO, nO, Cf, [(mp.mpf('0.4'), mp.mpf(0)), (mp.mpf('-0.3'), mp.mpf(0))])
    print(f'   {name}: orbits {orbits}, alpha_O^2 {[round(float(a @ a), 4) for a in aO]}, labels {nO.tolist()}, parent deg tau_j {degs}')
    report(f'{name}: parent solution sigma-invariant ({float(inv):.0e}); book bilinear eq of folded algebra, res {float(res):.1e}', res < 1e-25 and inv < 1e-25)
mu1 = mp.sqrt(2); mu2 = mp.sqrt(6)
d = {a: [mp.mpf(int(round(complex(x).real))) for x in Al.delta[a]] for a in Al.delta}   # exact integer first-order data for d4
assert all(np.allclose(Al.delta[a], [float(v) for v in d[a]]) for a in d)
run('d3^(2) <- d4 (centre element (01)(34)), species 1', [1, 0, 2, 4, 3], d[1], mu1)
run('d3^(2) <- d4 (centre element (01)(34)), species 2', [1, 0, 2, 4, 3], d[2], mu2)
run('g2^(1) <- d4 (triality fixing 0,2), coincident orbit {1,3,4}', [0, 3, 2, 4, 1], [d[1][j]+d[3][j]+d[4][j] for j in range(5)], mu1)
run('g2^(1) <- d4 (triality), species 2', [0, 3, 2, 4, 1], d[2], mu2)
