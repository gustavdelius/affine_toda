import numpy as np, os, itertools, time
from fractions import Fraction as Fr
from scipy.special import loggamma
from bspin import *
Pk = list(np.load(f'bspin_P_n{N}_w{omega}.npy'))
def br(x, l): return (x - q**l)/(1 - x*q**l)
rho = lambda x: [1.0, br(x, 2), br(x, 6), br(x, 2)*br(x, 10)]          # n = 3 (from bspin_struct.py)
Rmat = lambda x: sum(r*P_ for r, P_ in zip(rho(x), Pk))
x1 = 0.61-0.77j; Rn = intertwiner(rep(x1), rep(1.0), range(N+1), Wn, Wn)[0]; Rn = Rn/Rn[0, 0]
closed_err = np.abs(Rmat(x1) - Rn).max()
lifts = [int(v) for v in os.environ.get('LIFTS', '0,0,0').split(',')]
A = [(2*j - 1)*omega + lifts[j-1] for j in range(1, N+1)]; T = A[-1]; mu = 2*T; J = int(os.environ.get('JTRUNC', '2000'))
def c_of(t): return np.prod([np.sin(np.pi*(t - a)) for a in A], axis=0)
def lf1(z):
    z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A) - len(A)*np.log(np.pi)
js = np.arange(1, J+1); Aj = 2*T*(js-1)
def logf_rel(t):
    d = lambda z, z0: lf1(z) - lf1(z0); r = 0.371
    return np.sum(d(t + Aj, Aj + r) - d(t + Aj + T, Aj + T + r) + d(-t + Aj + T, Aj + T + r) - d(-t + Aj + 2*T, Aj + 2*T + r))
_K = c_of(0.123)*c_of(-0.123)*np.exp(logf_rel(0.123) + logf_rel(-0.123)); _norm = np.sqrt(_K + 0j)
_F = lambda t: c_of(t)*np.exp(logf_rel(t))/_norm; _sg = np.sign(_F(1e-6).real)
def F(theta): return _sg*_F(mu*theta/(2j*np.pi))
def S(theta): return F(theta)*Rmat(np.exp(mu*theta))
# crossing-sign test (same convention as for c_n)
C = np.load(f'bspin_C_n{N}_w{omega}.npy'); tc = np.load(f'bspin_tc_n{N}_w{omega}.npy')[0]; Ci = np.linalg.inv(C)
P = np.zeros((D*D, D*D))
for i in range(D):
    for j in range(D): P[j*D+i, i*D+j] = 1
def pt1(M_): return M_.reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(D*D, D*D)
th = 0.31+0.17j
lhs = P @ S(1j*np.pi - th); rhs = np.kron(C, np.eye(D)) @ pt1(P @ (P @ S(th)) @ P) @ np.kron(Ci, np.eye(D))
cross = np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel())
print(f"lifts {lifts}: T = {T:.3f}; closed-form R err {closed_err:.0e}; unitarity {np.abs(S(th) @ S(-th) - np.eye(D*D)).max():.0e};"
      f" crossing test factor {cross.real:+.4f}{cross.imag:+.4f}i (resid {np.abs(lhs - cross*rhs).max()/np.abs(lhs).max():.0e}); x(i pi)/x_c = {np.exp(2j*np.pi*T)/tc:.4f}")
# fusion poles
rng = np.random.default_rng(0); X = rng.normal(size=(D*D, 40)) + 1j*rng.normal(size=(D*D, 40))
def probe(t0):
    n_ = [np.linalg.norm(S(1j*np.pi*(t0 + d)/T) @ X) for d in (1e-5, 1e-7)]
    sv = np.linalg.svd(1e-7*S(1j*np.pi*(t0 + 1e-7)/T), compute_uv=False); return np.log10(n_[1]/n_[0])/2, int(np.sum(sv > 1e-6*sv[0]))
for t0, lab in [(A[0], "3+3 -> 2 (22)"), (A[1], "3+3 -> 1 (7)"), (T - 1, "lowest breather (1)")]:
    od, rk = probe(t0); print(f"   t = {t0:7.3f} {lab:22s} order {od:.2f} rank {rk}")
# lowest breather and S[B,3]
tB = T - 1
d_ = 1e-7; U_, sv, _ = np.linalg.svd(d_*S(1j*np.pi*(tB + d_)/T)); sig = U_[:, 0]; beta = np.pi*tB/T; hw = np.eye(D)[0]
def SB(theta):
    v = np.kron(sig, hw).reshape(D, D, D)
    def app(v, M, pos):
        w = np.moveaxis(v, [pos, pos+1], [0, 1]); sh = w.shape
        return np.moveaxis((M @ w.reshape(D*D, -1)).reshape(sh), [0, 1], [pos, pos+1])
    v = app(v, S(theta + 1j*beta/2), 1); v = app(v, S(theta - 1j*beta/2), 0)
    tgt = np.kron(hw, sig); return np.vdot(tgt, v.ravel())/np.vdot(tgt, tgt)
cand = set()
for j in (1, 2):
    Ajj = 2*T*(j-1)
    for a in A:
        for k in range(int(4*T)): cand |= {a - k - Ajj, -1 - a - k - Ajj, Ajj + T - a + k, 1 + Ajj + T + a + k}
cands = sorted({round(s*(tp + o), 9) for tp in cand for o in (tB/2, -tB/2) for s in (1, -1)})
cands = [c_ for c_ in cands if 1e-6 < c_ < T - 1e-6]
pz = []
for t0 in cands:
    n_ = [abs(SB(1j*np.pi*(t0 + d)/T)) for d in (1e-3, 1e-5)]; od = np.log10(n_[1]/n_[0])/2
    if abs(od) > 0.5: pz.append((t0, int(round(od))))
ys = [];  [ys.extend([y if o > 0 else -y]*abs(o)) for y, o in pz]
def model(yl, th_, sg=1):
    out = sg
    for y in yl: out *= np.sinh(th_/2 + 1j*np.pi*y/(2*T))/np.sinh(th_/2 - 1j*np.pi*y/(2*T))
    return out
r = np.array([SB(th_)/model(ys, th_) for th_ in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)]); sg = int(np.sign(r.mean().real))
print(f"   S[B,3]: {len(ys)} blocks y = {[round(y, 4) for y in ys]}, sign {sg:+d}, match {np.abs(r/sg - 1).max():.0e}")
# S[B,B] by bootstrap (numeric block arithmetic, period 2T)
def red(y):
    y = (y + T) % (2*T) - T
    return y if abs(y + T) > 1e-9 else T
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
BB = norm_blocks([y + tB/2 for y in ys] + [y - tB/2 for y in ys], sg*sg)
print(f"   S[B,B] = {BB[1]:+d} x blocks {[round(y, 4) for y in BB[0]]}")
# DGZ a_{2n-1}^(2) S_nn = prod_{a=0}^{n-1} (2a)(H-2a)/((2a+B)(H-2a-B)),  H = h + B/2,  h = 2n-1 ; try to find s = T/H
h = 2*N - 1
def dgz_blocks(Hv):
    Bv = 2*(Hv - h); s_ = T/Hv; yl = []
    for a in range(N):
        yl += [2*a*s_, (Hv - 2*a)*s_, -(2*a + Bv)*s_, -(Hv - 2*a - Bv)*s_]
    return norm_blocks(yl, 1)
best = None
poles = [y for y in BB[0] if y > 0]
for y1 in poles:
    Hv = 2*T/y1                                  # assume the block (2) of DGZ is our pole y1
    db = dgz_blocks(Hv)
    if len(db[0]) == len(BB[0]) and np.allclose(db[0], BB[0], atol=1e-7) and db[1] == BB[1]:
        best = Hv; break
if best: print(f"   -> IDENTICAL to DGZ a_{2*N-1}^(2) particle S_{N}{N} with H = {best:.6f}, B = {2*(best - h):.6f};  (omega+?) check: T/H = {T/best:.6f}")
else: print("   -> no DGZ S_nn identification found for this lift choice")
