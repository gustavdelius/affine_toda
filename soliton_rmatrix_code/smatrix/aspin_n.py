import numpy as np, os, math, time
from scipy.special import loggamma
from bspin import *
t0_ = time.time()
def br(x, l): return (x - q**l)/(1 - x*q**l)
def rho(x):                       # conjecture: Lambda^{n-k} eigenvalue = prod_{i=1}^{ceil(k/2)} <4k-8i+6>
    out = []
    for k in range(N+1):
        v = 1.0
        for i in range(1, (k+1)//2 + 1): v = v*br(x, 4*k - 8*i + 6)
        out.append(v)
    return out
cache = f"bspinN_P_n{N}_w{omega}.npy"
x0, x1 = 1.37+0.41j, 0.61-0.77j
if os.path.exists(cache): Pk = list(np.load(cache)); R1 = None
else:
    R0 = intertwiner(rep(x0), rep(1.0), range(N+1), Wn, Wn)[0]; R0 = R0/R0[0, 0]; lam = rho(x0); Pk = []
    for kk in range(N+1):
        M = np.eye(D*D, dtype=complex)
        for jj in range(N+1):
            if jj != kk: M = M @ (R0 - lam[jj]*np.eye(D*D))/(lam[kk] - lam[jj])
        Pk.append(M)
    np.save(cache, np.array(Pk))
    R1 = intertwiner(rep(x1), rep(1.0), range(N+1), Wn, Wn)[0]; R1 = R1/R1[0, 0]
Rmat = lambda x: sum(r*P_ for r, P_ in zip(rho(x), Pk))
ranks = [int(round(np.trace(P_).real)) for P_ in Pk]
print(f"b_{N}^(1) spinor: channel ranks {ranks} (Lambda^k: {[math.comb(2*N+1, N-k) for k in range(N+1)]});"
      + (f" conjectured eigenvalues reproduce an independent solve at x1 to {np.abs(Rmat(x1) - R1).max():.0e};" if R1 is not None else "")
      + f" projector idempotence {max(np.abs(P_ @ P_ - P_).max() for P_ in Pk):.0e}  [{time.time()-t0_:.0f}s]")
A = [(2*j - 1)*omega + (j - 1) for j in range(1, N+1)]; T = A[-1]; mu = 2*T; J = int(os.environ.get('JTRUNC', '2000'))
def c_of(t): return np.prod([np.sin(np.pi*(t - a)) for a in A], axis=0)
def lf1(z):
    z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A) - len(A)*np.log(np.pi)
js = np.arange(1, J+1); Aj = 2*T*(js-1)
def logf_rel(t):
    d = lambda z, z0: lf1(z) - lf1(z0); r = 0.371
    return np.sum(d(t + Aj, Aj + r) - d(t + Aj + T, Aj + T + r) + d(-t + Aj + T, Aj + T + r) - d(-t + Aj + 2*T, Aj + 2*T + r))
_K = c_of(0.123)*c_of(-0.123)*np.exp(logf_rel(0.123) + logf_rel(-0.123)); _norm = np.sqrt(_K + 0j)
_F = lambda t: c_of(t)*np.exp(logf_rel(t))/_norm; _sg = np.sign(_F(1e-6).real)
def S(theta): return _sg*_F(mu*theta/(2j*np.pi))*Rmat(np.exp(mu*theta))
th = 0.31+0.17j; print(f"T = {T:.3f}; unitarity {np.abs(S(th) @ S(-th) - np.eye(D*D)).max():.0e}")
rng = np.random.default_rng(0); X = rng.normal(size=(D*D, 60)) + 1j*rng.normal(size=(D*D, 60))
def probe(t0):
    n_ = [np.linalg.norm(S(1j*np.pi*(t0 + d)/T) @ X) for d in (1e-5, 1e-7)]
    sv = np.linalg.svd(1e-7*S(1j*np.pi*(t0 + 1e-7)/T), compute_uv=False); return np.log10(n_[1]/n_[0])/2, int(np.sum(sv > 1e-6*sv[0]))
for j in range(1, N):
    od, rk = probe(A[j-1]); print(f"   t = a_{j} = {A[j-1]:7.3f}: {N}+{N} -> soliton {N-j}   order {od:.2f} rank {rk}")
od, rk = probe(T - 1); print(f"   t = T-1 = {T-1:7.3f}: lowest breather   order {od:.2f} rank {rk}")
H = T/(omega + 0.5)
Ms = {N - j: 2*np.cos(np.pi*A[j-1]/(2*T)) for j in range(1, N)}
print("   soliton masses M_a/M_n:", {a: round(v, 6) for a, v in sorted(Ms.items())}, " vs 2 sin(a pi/H):", {a: round(2*np.sin(a*np.pi/H), 6) for a in sorted(Ms)})
# lowest breather amplitudes
tB = T - 1; d_ = 1e-7; U_, sv, _ = np.linalg.svd(d_*S(1j*np.pi*(tB + d_)/T)); sig = U_[:, 0]; beta = np.pi*tB/T; hw = np.eye(D)[0]
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
cands = [c_ for c_ in cands if 1e-6 < c_ < T - 1e-6]; pz = []
for t0 in cands:
    n_ = [abs(SB(1j*np.pi*(t0 + d)/T)) for d in (1e-3, 1e-5)]; od = np.log10(n_[1]/n_[0])/2
    if abs(od) > 0.5: pz.append((t0, int(round(od))))
ys = []; [ys.extend([y if o > 0 else -y]*abs(o)) for y, o in pz]
def model(yl, th_, sg=1):
    out = sg
    for y in yl: out *= np.sinh(th_/2 + 1j*np.pi*y/(2*T))/np.sinh(th_/2 - 1j*np.pi*y/(2*T))
    return out
r = np.array([SB(th_)/model(ys, th_) for th_ in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)]); sg = int(np.sign(r.mean().real))
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
BB = norm_blocks([y + tB/2 for y in ys] + [y - tB/2 for y in ys], sg*sg)
h = 2*N - 1; Bv = 2*(H - h); s_ = T/H; yl = []
for a in range(N): yl += [2*a*s_, (H - 2*a)*s_, -(2*a + Bv)*s_, -(H - 2*a - Bv)*s_]
DG = norm_blocks(yl, 1)
same = len(DG[0]) == len(BB[0]) and np.allclose(DG[0], BB[0], atol=1e-7) and DG[1] == BB[1]
print(f"   S[B,{N}] fit match {np.abs(r/sg - 1).max():.0e};  S[B,B] vs DGZ a_{2*N-1}^(2) S_{N}{N} (H = T/(w+1/2) = {H:.5f}, B = -1/(w+1/2) = {Bv:.5f}): {'IDENTICAL' if same else 'DIFFERENT'}  [{time.time()-t0_:.0f}s]")
if not same: print("     ours", np.round(BB[0], 4), BB[1], "\n     DGZ ", np.round(DG[0], 4), DG[1])
