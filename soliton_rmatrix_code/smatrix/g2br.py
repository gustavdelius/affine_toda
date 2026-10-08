import numpy as np, os, time
from scipy.special import loggamma
exec(open('g2sol.py').read().split("def label(z):")[0])
offs = [float(v) for v in os.environ.get('AOFF', '0,2,3').split(',')]
A = [omega + offs[0], 4*omega + offs[1], 6*omega + offs[2]]; T = A[-1]; mu = 2*T; J = 3000
def c_of(t): return np.prod([np.sin(np.pi*(t - a)) for a in A], axis=0)
def lf1(z):
    z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A)
js = np.arange(1, J+1); Aj = 2*T*(js-1)
def logf(t):
    d = lambda z, z0: lf1(z) - lf1(z0); r = 0.371
    return np.sum(d(t + Aj, Aj + r) - d(t + Aj + T, Aj + T + r) + d(-t + Aj + T, Aj + T + r) - d(-t + Aj + 2*T, Aj + 2*T + r))
_K = c_of(0.123)*c_of(-0.123)*np.exp(logf(0.123) + logf(-0.123)); _nm = np.sqrt(_K + 0j)
_F = lambda t: c_of(t)*np.exp(logf(t))/_nm; _sg = np.sign(_F(1e-6).real)
def S(theta): return _sg*_F(mu*theta/(2j*np.pi))*Rcoef(np.exp(mu*theta))
th = 0.31+0.17j; print(f"{ALG}, zeros {np.round(A, 3)}, T = {T:.3f}: unitarity {np.abs(S(th) @ S(-th) - np.eye(N2)).max():.0e}")
rng = np.random.default_rng(0); X = rng.normal(size=(N2, 20)) + 1j*rng.normal(size=(N2, 20))
for t0_, lab in [(A[0], "fusion 1"), (A[1], "fusion 2 (self?)"), (T - 1, "lowest breather")]:
    nn = [np.linalg.norm(S(1j*np.pi*(t0_ + d)/T) @ X) for d in (1e-5, 1e-7)]
    sv = np.linalg.svd(1e-7*S(1j*np.pi*(t0_ + 1e-7)/T), compute_uv=False)
    print(f"   t = {t0_:7.3f} {lab:18s} order {np.log10(nn[1]/nn[0])/2:.2f} rank {int(np.sum(sv > 1e-6*sv[0]))}")
tB = T - 1; d_ = 1e-7; U_, sv, _ = np.linalg.svd(d_*S(1j*np.pi*(tB + d_)/T)); sig = U_[:, 0]; beta = np.pi*tB/T
hw = np.zeros(n); hw[np.argmax(W @ np.array([1.0, 0.3]))] = 1
def SB(theta):
    v = np.kron(sig, hw).reshape(n, n, n)
    def app(v, M, pos):
        w = np.moveaxis(v, [pos, pos+1], [0, 1]); sh = w.shape
        return np.moveaxis((M @ w.reshape(N2, -1)).reshape(sh), [0, 1], [pos, pos+1])
    v = app(v, S(theta + 1j*beta/2), 1); v = app(v, S(theta - 1j*beta/2), 0)
    tgt = np.kron(hw, sig); return np.vdot(tgt, v.ravel())/np.vdot(tgt, tgt)
cand = set()
for j in (1, 2):
    Ajj = 2*T*(j-1)
    for a in A:
        for k in range(int(3*T)): cand |= {a - k - Ajj, -1 - a - k - Ajj, Ajj + T - a + k, 1 + Ajj + T + a + k}
cands = sorted({round(s*(tp + o), 9) for tp in cand for o in (tB/2, -tB/2) for s in (1, -1)})
cands = [c_ for c_ in cands if 1e-6 < c_ < T - 1e-6]; pz = []
for t0_ in cands:
    nn = [abs(SB(1j*np.pi*(t0_ + d)/T)) for d in (1e-3, 1e-5)]; od = np.log10(nn[1]/nn[0])/2
    if abs(od) > 0.5: pz.append((t0_, int(round(od))))
ys = []; [ys.extend([y if o > 0 else -y]*abs(o)) for y, o in pz]
blkf = lambda y, th_: np.sinh(th_/2 + 1j*np.pi*y/(2*T))/np.sinh(th_/2 - 1j*np.pi*y/(2*T))
r = np.array([SB(th_)/np.prod([blkf(y, th_) for y in ys]) for th_ in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)]); sg = int(np.sign(r.mean().real))
print(f"   S[B,sol]: {len(ys)} blocks, sign {sg:+d}, fit {np.abs(r/sg - 1).max():.0e}  [{time.time()-t0:.0f}s]")
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
def cds_g2(a, H):                                   # CDS eq. (3.10), H = 6 + 3B
    Bv = (H - 6)/3
    def brace(x, nu): return [(x - nu*Bv - 1, 1), (x + nu*Bv + 1, 1), (x + nu*Bv + Bv - 1, -1), (x - nu*Bv - Bv + 1, -1)]
    spec = {11: [(1, 0), (H/2, 0.5), (H - 1, 0)], 12: [(H/3, 1), (2*H/3, 1)], 22: [(H/3 - 1, 1), (H/3 + 1, 1), (2*H/3 - 1, 1), (2*H/3 + 1, 1)]}[a]
    yl = []
    for x, nu in spec:
        for z, s in brace(x, nu): yl.append(z*T/H if s > 0 else -z*T/H)
    return norm_blocks(yl, 1)
Hpred = 12*T/(T - 3) if ALG == 'g2' else 6*T/(T + 3)
poles_ = [y for y in BB[0] if y > 0]
Hc = sorted({c*T/y for y in poles_ for c in (1, 2, 3, 4, 5, 6, 0.5, 1.5)} | {Hpred})
for a in (11, 22):
    best = (np.inf, None)
    for H in Hc:
        D_ = cds_g2(a, H)
        if len(D_[0]) == len(BB[0]) and D_[1] == BB[1]:
            mis = np.abs(np.array(D_[0]) - np.array(BB[0])).max()
            if mis < best[0]: best = (mis, H)
    tag = "IDENTICAL" if best[0] < 1e-7 else "different"
    print(f"   S[B,B] vs CDS (g2,d4^(3)) S_{a}: {tag}" + (f" at H = {best[1]:.6f} (predicted {Hpred:.6f})" if best[1] else ""))
print("   ours:", np.round(BB[0], 3), BB[1])
