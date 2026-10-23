import numpy as np, time
from fractions import Fraction as Fr
from spinorn import *
exec(open('bootstrap_exact.py').read().split("# numerically identified S[B^a, soliton 3]")[0].replace("T = (3, Fr(1)); twoT = (6, Fr(2))", "")
     .replace("w = 2.37", f"w = {omega}").replace("3*w + 1", f"{N}*w + {Fr(N-1, 2)}").replace("(3, 1)", f"({N}, {Fr(N-1, 2)})").replace("(-3, -1)", f"(-{N}, -{Fr(N-1, 2)})")
     .replace("a - 6, b - 2", f"a - {2*N}, b - {N-1}").replace("a + 6, b + 2", f"a + {2*N}, b + {N-1}"))
w_ = lambda y: float(y[0]*omega + y[1])
tB = T - 1                                           # lowest breather of soliton N: t = T - 1 = N w + (N-3)/2
def singlet(t0):
    d = 1e-7; U_, sv, _ = np.linalg.svd(d*S(1j*np.pi*(t0 + d)/T)); assert sv[1] < 1e-5*sv[0]; return U_[:, 0]
sig = singlet(tB); beta = np.pi*tB/T; hw = np.eye(D)[0]
def SBn(theta):                                       # breather B^N scattering soliton N : scalar
    v = np.kron(sig, hw).reshape(D, D, D)
    def app(v, M, pos):
        w = np.moveaxis(v, [pos, pos+1], [0, 1]); sh = w.shape
        return np.moveaxis((M @ w.reshape(D*D, -1)).reshape(sh), [0, 1], [pos, pos+1])
    v = app(v, S(theta + 1j*beta/2), 1); v = app(v, S(theta - 1j*beta/2), 0)
    tgt = np.kron(hw, sig); return np.vdot(tgt, v.ravel())/np.vdot(tgt, tgt)
# candidates: shifted Gamma-pole superset of F and reflections
cand = set()
for j in (1, 2):
    Aj_ = 2*T*(j-1)
    for a in A:
        for k in range(0, int(4*T)):
            cand |= {a - k - Aj_, -1 - a - k - Aj_, Aj_ + T - a + k, 1 + Aj_ + T + a + k}
cands = sorted({round(s*(tp + o), 9) for tp in cand for o in (tB/2, -tB/2) for s in (1, -1)})
cands = [c for c in cands if 1e-6 < c < T - 1e-6]
t1 = time.time(); pz = []
for t0 in cands:
    n_ = [abs(SBn(1j*np.pi*(t0 + d)/T)) for d in (1e-3, 1e-5)]
    od = np.log10(n_[1]/n_[0])/2
    if abs(od) > 0.5: pz.append((t0, int(round(od))))
def exact(y):
    for A2 in range(-40, 41):
        for C4 in range(-60, 61):
            if abs(A2/2*omega + C4/4 - y) < 1e-6: return (Fr(A2, 2), Fr(C4, 4))
    return None
blocks = []
for y, od in pz:
    e = exact(y if od > 0 else -y)
    blocks += [e]*abs(od)
def model(bl, th, sg=1):
    out = sg
    for (a_, c_) in bl:
        y = float(a_)*omega + float(c_); out *= np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T))
    return out
r = np.array([SBn(th)/model(blocks, th) for th in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)]); sg = int(np.sign(r.mean().real))
print(f"S[B{N},{N}] = {'+' if sg > 0 else '-'} " + " ".join(f"({show(y)})" for y in blocks) + f"   (match {np.abs(r/sg - 1).max():.1e}; {len(cands)} candidates, {time.time()-t1:.0f}s)")
bl, sgn = normalize(blocks, sg)
BB = shift(bl, sgn, (Fr(N, 2), Fr(N - 3, 4)))            # S[B^N, B^N]: fuse soliton N + N -> B^N (half angle tB/2 = N w/2 + 1/4)
print(f"S[B{N},B{N}] (bootstrap) = {'+' if BB[1] > 0 else '-'} " + " ".join(f"({show(y)})" for y in BB[0]))
# DGZ c_N^(1) particle S-matrix, H = 2N + B, mapped with y = x(omega+1/2)/2, B = -1/(omega+1/2)
def to_y(c, d): c, d = Fr(c), Fr(d); return (c/2, c/4 - d/2)
def dgz(a, b):
    if a < b: a, b = b, a
    num, den = [], []
    for p in range(1, b + 1):
        for (c, d) in [(a - b - 1 + 2*p, 0), (2*N - a + b + 1 - 2*p, 1)]:
            num += [(c - 1, d), (c + 1, d)]; den += [(c - 1, d + 1), (c + 1, d - 1)]
    return normalize([to_y(c, d) for c, d in num] + [tuple(-v for v in to_y(c, d)) for c, d in den], 1)
for a in range(1, N+1):
    if dgz(a, a) == BB: print(f"   -> IDENTICAL to the DGZ c_{N}^(1) particle amplitude S_{a}{a}")
for a_ in range(1, N+1): print(f"   DGZ S_{a_}{a_} = {'+' if dgz(a_, a_)[1] > 0 else '-'} " + " ".join(f"({show(y)})" for y in dgz(a_, a_)[0]))
if False: print(f"   DGZ S_{N}{N} = {'+' if dgz(N, N)[1] > 0 else '-'} " + " ".join(f"({show(y)})" for y in dgz(N, N)[0]))
np.save(f'SBn_n{N}.npy', np.array([blocks, sg], dtype=object), allow_pickle=True)
