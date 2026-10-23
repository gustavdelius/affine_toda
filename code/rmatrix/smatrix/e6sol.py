import numpy as np, time, os
from scipy.special import loggamma
src = open('f4R2.py').read().split("def label(z):")[0]
exec(src)
t0_ = time.time()
Ball = np.hstack([np.hstack([cs[k] for k in sorted(cs)]) for cs in CS]); sizes = dimsC
Binv = np.linalg.inv(Ball); Pk = []; off = 0
for j, cs in enumerate(CS):
    Bj = np.hstack([cs[k] for k in sorted(cs)]); d = Bj.shape[1]; Pk.append(Bj @ Binv[off:off+d]); off += d
print("projectors:", [int(round(np.trace(P_).real)) for P_ in Pk], f"idempotence {max(np.abs(P_@P_-P_).max() for P_ in Pk):.0e}")
def br(x, l): return (x - q**l)/(1 - x*q**l)
forms = {324: lambda x: 1.0, 273: lambda x: br(x, 2), 52: lambda x: br(x, 8), 26: lambda x: br(x, 2)*br(x, 12), 1: lambda x: br(x, 8)*br(x, 18)}
Rmat = lambda x: sum(forms[d](x)*P_ for d, P_ in zip(dimsC, Pk))
lifts = [int(v) for v in os.environ.get('LIFTS', '0,1,2,3').split(',')]
A = [omega + lifts[0], 4*omega + lifts[1], 6*omega + lifts[2], 9*omega + lifts[3]]; T = A[-1]; mu = 2*T; J = 1500
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
th = 0.31+0.17j; print(f"lifts {lifts}, T = {T:.3f}: unitarity {np.abs(S(th) @ S(-th) - np.eye(n*n)).max():.0e}")
rng = np.random.default_rng(0); X = rng.normal(size=(n*n, 60)) + 1j*rng.normal(size=(n*n, 60))
def probe(t0):
    nn_ = [np.linalg.norm(S(1j*np.pi*(t0 + d)/T) @ X) for d in (1e-5, 1e-7)]
    sv = np.linalg.svd(1e-7*S(1j*np.pi*(t0 + 1e-7)/T), compute_uv=False); return np.log10(nn_[1]/nn_[0])/2, int(np.sum(sv > 1e-6*sv[0]))
for t0, lab in [(A[0], "26+26 -> 299"), (A[1], "26+26 -> 53"), (A[2], "26+26 -> 26"), (T - 1, "lowest breather")]:
    od, rk = probe(t0); print(f"   t = {t0:7.3f} {lab:16s} order {od:.2f} rank {rk}")
tB = T - 1; d_ = 1e-7; U_, sv, _ = np.linalg.svd(d_*S(1j*np.pi*(tB + d_)/T)); sig = U_[:, 0]; beta = np.pi*tB/T
hw = np.zeros(n); hw[np.argmax(W @ np.array([8, 4, 2, 1.]))] = 1                  # a highest-weight vector of the 26
def SB(theta):
    v = np.kron(sig, hw).reshape(n, n, n)
    def app(v, M, pos):
        w = np.moveaxis(v, [pos, pos+1], [0, 1]); sh = w.shape
        return np.moveaxis((M @ w.reshape(n*n, -1)).reshape(sh), [0, 1], [pos, pos+1])
    v = app(v, S(theta + 1j*beta/2), 1); v = app(v, S(theta - 1j*beta/2), 0)
    tgt = np.kron(hw, sig); return np.vdot(tgt, v.ravel())/np.vdot(tgt, tgt)
cand = set()
for j in (1, 2):
    Ajj = 2*T*(j-1)
    for a in A:
        for k in range(int(3*T)): cand |= {a - k - Ajj, -1 - a - k - Ajj, Ajj + T - a + k, 1 + Ajj + T + a + k}
cands = sorted({round(s*(tp + o), 9) for tp in cand for o in (tB/2, -tB/2) for s in (1, -1)})
cands = [c_ for c_ in cands if 1e-6 < c_ < T - 1e-6]; pz = []
for t0 in cands:
    nn_ = [abs(SB(1j*np.pi*(t0 + d)/T)) for d in (1e-3, 1e-5)]; od = np.log10(nn_[1]/nn_[0])/2
    if abs(od) > 0.5: pz.append((t0, int(round(od))))
ys = []; [ys.extend([y if o > 0 else -y]*abs(o)) for y, o in pz]
def model(yl, th_, sg=1):
    out = sg
    for y in yl: out *= np.sinh(th_/2 + 1j*np.pi*y/(2*T))/np.sinh(th_/2 - 1j*np.pi*y/(2*T))
    return out
r = np.array([SB(th_)/model(ys, th_) for th_ in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)]); sg = int(np.sign(r.mean().real))
print(f"   S[B,26]: {len(ys)} blocks, sign {sg:+d}, fit match {np.abs(r/sg - 1).max():.0e}  ({len(cands)} candidates) [{time.time()-t0_:.0f}s]")
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
# CDS (f4^(1), e6^(2)) particle S-matrices, eq. (4.6) first form, B = (H-12)/3, [x]_0 = {x}_0 {H-x}_0
def cds(a, b, H):
    Bv = (H - 12)/3
    def brk(x): return [(x - 1, 1), (x + 1, 1), (x - 1 + Bv, -1), (x + 1 - Bv, -1)]
    def sq(x): return brk(x) + brk(H - x)
    spec = {(1, 1): [1, H/3+1], (1, 2): [H/6+2, H/2+2], (1, 3): [2, H/3, H/3+2], (1, 4): [3, H/3-1, H/3+1, H/3+3],
            (2, 2): [1, H/3-3, H/3+1, H/3+3], (2, 3): [H/6+1, H/6+3, H/2+1, H/2+3], (2, 4): [H/6, H/6+2, H/6+4, H/2, H/2+2, H/2+4],
            (3, 3): [1, 3, H/3-1, H/3+1, H/3+1, H/3+3], (3, 4): [2, 4, H/3-2, H/3, H/3, H/3+2, H/3+2, H/3+4],
            (4, 4): [1, 3, 5, H/3-3, H/3-1, H/3-1, H/3+1, H/3+1, H/3+1, H/3+3, H/3+3, H/3+5]}[(a, b)]
    yl = []
    for x in spec:
        for (z, s) in sq(x): yl.append(z*T/H if s > 0 else -z*T/H)
    return norm_blocks(yl, 1)
Hc = 2*T/(omega + lifts[0])
for a in (1, 2, 3, 4):
    D_ = cds(a, a, Hc); same = len(D_[0]) == len(BB[0]) and np.allclose(D_[0], BB[0], atol=1e-7) and D_[1] == BB[1]
    print(f"   S[B,B] vs CDS S_{a}{a} at H = 2T/(w+{lifts[0]}) = {Hc:.5f}: {'IDENTICAL' if same else 'different'}")
print("   ours:", np.round(BB[0], 3), BB[1])
