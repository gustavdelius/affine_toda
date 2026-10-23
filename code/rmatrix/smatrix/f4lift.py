import numpy as np, itertools
from scipy.special import loggamma
omega = 2.37
BBs, BBsign, T = np.load("f4BB_m1.npy", allow_pickle=True); BBs = list(BBs)
A0 = [omega + 1, 3*omega + 0.5 + 1, 4*omega + 1, 6*omega + 1.5]
def Fmaker(A, J=3000):
    js = np.arange(1, J+1); Aj = 2*T*(js-1)
    def lf1(z):
        z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A)
    def logf(t):
        d = lambda z, z0: lf1(z) - lf1(z0); r = 0.371
        return np.sum(d(t + Aj, Aj + r) - d(t + Aj + T, Aj + T + r) + d(-t + Aj + T, Aj + T + r) - d(-t + Aj + 2*T, Aj + 2*T + r))
    c_of = lambda t: np.prod([np.sin(np.pi*(t - a)) for a in A])
    K = c_of(0.123)*c_of(-0.123)*np.exp(logf(0.123) + logf(-0.123)); nrm = np.sqrt(K + 0j)
    Fr = lambda t: c_of(t)*np.exp(logf(t))/nrm; sg = np.sign(Fr(1e-6).real)
    return lambda theta: sg*Fr(2*T*theta/(2j*np.pi))
blk = lambda y, th: np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T))
F0 = Fmaker(A0); ths = [0.37+0.21j, -0.8+1.3j, 1.1-0.4j]
X = {}
for j in range(3):                                      # shift a_j by +1 (j = 0, 1, 2; the last zero is T itself)
    A1 = list(A0); A1[j] += 1; F1 = Fmaker(A1)
    ratio = np.array([F1(th)/F0(th) for th in ths]); a = A0[j]
    for cand in ([a + 1, T - a - 1], [-(a + 1), -(T - a - 1)], [a, T - a], [-a, -(T - a)], [a + 1, -(T - a)], [a + 1, T - a], [-(a + 1), T - a - 1]):
        for sg in (1, -1):
            mdl = np.array([sg*np.prod([blk(y, th) for y in cand]) for th in ths])
            if np.allclose(ratio, mdl, rtol=2e-3): X[j] = (cand, sg)
    print(f"shift a_{j} -> a_{j}+1: F changes by {X.get(j)}  (ratio at test points {np.round(ratio, 4)})")
def red(y):
    y = (y + T) % (2*T) - T; return y if abs(y + T) > 1e-9 else T
def norm_blocks(yl, sgn):
    out = []; s_ = sgn
    for y in yl:
        y = red(y)
        if abs(y) < 1e-9: continue
        if abs(abs(y) - T) < 1e-9: s_ *= -1; continue
        out.append(y)
    res = []
    for y in out:
        m = [z for z in res if abs(z + y) < 1e-9]
        if m: res.remove(m[0])
        else: res.append(y)
    return sorted(res), s_
exec(open('e6sol.py').read().split("def cds(a, b, H):")[1].split("Hc = 2*T")[0].join(["def cds(a, b, H):", ""]))
tB = T - 1
def BB_for(shifts):
    yl = list(BBs); sg = BBsign
    for j, k in enumerate(shifts):
        if k == 0 or j not in X: continue
        cand, s = X[j]
        for _ in range(abs(k)):
            for off, mult in ((tB, 1), (0, 2), (-tB, 1)):
                for y in cand: yl += [(y + off) if k > 0 else -(y + off)]*mult
            sg *= s**4
    return norm_blocks(yl, sg)
from scipy.optimize import minimize_scalar
def best_match(BB):
    out = []
    for a in (1, 2, 3, 4):
        def mm(H):
            D_ = cds(a, a, H)
            if len(D_[0]) != len(BB[0]) or D_[1] != BB[1]: return np.inf
            return np.abs(np.array(D_[0]) - np.array(BB[0])).max()
        Hs = np.linspace(4, 30, 13001); v = np.array([mm(h) for h in Hs]); i = np.argmin(v)
        if np.isfinite(v[i]):
            rr = minimize_scalar(mm, bounds=(Hs[max(i-2, 0)], Hs[min(i+2, 13000)]), method='bounded', options={'xatol': 1e-12})
            out.append((rr.fun, a, rr.x))
    return min(out) if out else (np.inf, None, None)
print("searching integer lift shifts (d_a2, d_a6, d_a8) in {-2..2}^3:")
hits = []
for sh in itertools.product(range(-2, 3), repeat=3):
    mis, a, H = best_match(BB_for(sh))
    if mis < 1e-7: hits.append((sh, a, H)); print(f"   shifts {sh}: IDENTICAL to CDS S_{a}{a} at H = {H:.8f}  (2T/H = {2*T/H:.6f}, B = {(H-12)/3:.6f})")
if not hits: print("   no exact match in this range")
