import numpy as np, time
from allS import *
import allS as AS
rng = np.random.default_rng(7)
tB = {3: 3*omega, 1: 3*omega + 0.5, 2: 3*omega + 0.5}          # lowest-breather poles in S_aa
def make_breather(a):
    t0 = tB[a]; D = dims[a]**2; X = rng.normal(size=(D, 6)) + 1j*rng.normal(size=(D, 6)); d = 1e-7
    Y, _ = apply_S(a, a, 1j*np.pi*(t0 + d)/T, X); U_, sv, _ = np.linalg.svd(d*Y)
    assert sv[1] < 1e-4*sv[0]
    sig = U_[:, 0]                                                 # singlet in W_a (x) W_a (fused basis)
    beta = np.pi*t0/T
    off = [o - beta/2 for o in part[a]['off']] + [o + beta/2 for o in part[a]['off']]
    vec = np.kron(part[a]['B'], part[a]['B']) @ sig
    part[f'B{a}'] = dict(off=off, B=vec[:, None], hw=np.array([1.0+0j])); dims[f'B{a}'] = 1
    return beta
beta = {a: make_breather(a) for a in (1, 2, 3)}
mB = {a: 2*({1: 2*np.cos(np.pi/4 + np.pi/(2*(6+2/omega))), 2: 2*np.cos(np.pi/(6+2/omega)), 3: 1.0}[a])*np.cos(beta[a]/2) for a in (1, 2, 3)}
print("lowest breather masses (M3=1):", {a: round(v, 5) for a, v in mB.items()})
# superset of real pole positions of S_33 (t')
cand33 = set()
for j in (1, 2, 3):
    Aj = 2*T*(j-1)
    for a_ in (omega, 2*omega + 0.5, 3*omega):
        for k in range(0, int(6*T)):
            cand33 |= {a_ - k - Aj, -1 - a_ - k - Aj, Aj + T - a_ + k, 1 + Aj + T + a_ + k, -(a_ - k - Aj), -(Aj + T - a_ + k)}
def scalar(A_, b, theta):
    K, res = K_top(A_, b, theta); return K
def blocks(A_, b):
    offs = [o for _, o in swaps(A_, b)]
    cands = sorted({round(tp - T*o/np.pi, 9) for tp in cand33 for o in offs})
    cands = [c for c in cands if -T + 1e-6 < c <= T + 1e-6]
    poles = []
    for t0 in cands:
        n = [abs(scalar(A_, b, 1j*np.pi*(t0 + d)/T)) for d in (1e-4, 1e-6)]
        od = np.log10(n[1]/n[0])/2
        if od > 0.5: poles += [t0]*int(round(od))
    return poles
def model(poles, theta):
    out = 1.0
    for y in poles: out *= np.sinh(theta/2 + 1j*np.pi*y/(2*T))/np.sinh(theta/2 - 1j*np.pi*y/(2*T))
    return out
results = {}
for a in (3, 1, 2):
    for b in (3, 1, 2):
        t1 = time.time(); A_ = f'B{a}'
        pls = blocks(A_, b)
        r = [scalar(A_, b, th)/model(pls, th) for th in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)]
        sign = np.mean(r); spread = np.abs(np.array(r) - sign).max()
        cross = scalar(A_, b, 1j*np.pi - (0.3+0.2j))/scalar(A_, b, 0.3+0.2j)
        results[(a, b)] = (pls, sign)
        print(f"S[B{a}, soliton {b}] = {sign.real:+.3f} x prod (y)_T over y = {[round(y, 4) for y in pls]}"
              f"   (match {spread:.0e}; S(i pi - th)/S(th) = {cross.real:+.4f}{cross.imag:+.4f}i) [{time.time()-t1:.0f}s]")
np.save('breather_soliton.npy', np.array([results, beta], dtype=object), allow_pickle=True)
