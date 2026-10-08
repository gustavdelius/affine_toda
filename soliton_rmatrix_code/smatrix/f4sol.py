import numpy as np, time, os
from scipy.special import loggamma
exec(open('e6R.py').read().split("u, res = solve_R(0.83+0.29j)")[0])
t0_ = time.time()
def Rapply(x, V):
    u, _ = solve_R(x); C_ = Binv @ V; out = np.zeros_like(V, dtype=complex)
    for ti, t in enumerate(types):
        for a in range(t['mult']):
            s, e_ = offs[ti][a]
            for b in range(t['mult']): out += u[uidx[(ti, b, a)]]*(t['B'][b] @ C_[s:e_])
    return out
m1 = int(os.environ.get('M1', '1'))
_o = [float(v) for v in os.environ.get('AOFF', f'{m1},{0.5+m1},1,1.5').split(',')]
A = [omega + _o[0], 3*omega + _o[1], 4*omega + _o[2], 6*omega + _o[3]]; T = A[-1]; mu = 2*T; J = 1500
def c_of(t): return np.prod([np.sin(np.pi*(t - a)) for a in A], axis=0)
def lf1(z):
    z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A) - len(A)*np.log(np.pi)
js = np.arange(1, J+1); Aj = 2*T*(js-1)
def logf_rel(t):
    d = lambda z, z0: lf1(z) - lf1(z0); r = 0.371
    return np.sum(d(t + Aj, Aj + r) - d(t + Aj + T, Aj + T + r) + d(-t + Aj + T, Aj + T + r) - d(-t + Aj + 2*T, Aj + 2*T + r))
_K = c_of(0.123)*c_of(-0.123)*np.exp(logf_rel(0.123) + logf_rel(-0.123)); _norm = np.sqrt(_K + 0j)
_F = lambda t: c_of(t)*np.exp(logf_rel(t))/_norm; _sg = np.sign(_F(1e-6).real)
def Fsc(theta): return _sg*_F(mu*theta/(2j*np.pi))
def Sapply(theta, V): return Fsc(theta)*Rapply(np.exp(mu*theta), V)
# crossing point and crossing-sign test
def mats(x):
    return [((x if i == 0 else 1)*E[i], F[i]/(x if i == 0 else 1), Kd[i]) for i in range(5)]
xa = 0.83+0.29j; tc = None
for eps in (1, -1):
    for l in range(-14, 15):
        tt = eps*q**l; rows = []
        for (Ei, Fi, Kk), (Ej, Fj, Kj) in zip(mats(xa), mats(xa*tt)):
            Kii = np.linalg.inv(Kk)
            for A_, B_ in (((-Kii@Ei).T, Ej), ((-Fi@Kk).T, Fj), (Kii.T, Kj)):
                rows.append(np.kron(np.eye(n), A_.T) - np.kron(B_, np.eye(n)))
        G = sum(r_.conj().T @ r_ for r_ in rows); del rows
        ev_, V_ = np.linalg.eigh(G)
        if ev_[0] < 1e-16*ev_[-1]: tc = tt; Cc = V_[:, 0].reshape(n, n); print(f"crossing point: V*(x) = V(x * {'+' if eps > 0 else '-'}q^{l});  x(i pi)/x_c = {np.exp(2j*np.pi*T)/tt:.4f}")
Pm = np.zeros((N2, N2))
for i in range(n):
    for j in range(n): Pm[j*n+i, i*n+j] = 1
def pt1(M_): return M_.reshape(n, n, n, n).transpose(2, 1, 0, 3).reshape(N2, N2)
th = 0.31+0.17j; Ci = np.linalg.inv(Cc)
Sfull = lambda z: Fsc(z)*Rmatrix(solve_R(np.exp(mu*z))[0])
lhs = Pm @ Sfull(1j*np.pi - th); rhs = np.kron(Cc, np.eye(n)) @ pt1(Pm @ (Pm @ Sfull(th)) @ Pm) @ np.kron(Ci, np.eye(n))
cr = np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel())
Uu = Sfull(th) @ Sfull(-th)
print(f"m1 = {m1}, T = {T:.3f}: unitarity {np.abs(Uu - np.eye(N2)).max():.0e}; crossing test factor {cr.real:+.4f}{cr.imag:+.4f}i (resid {np.abs(lhs - cr*rhs).max()/np.abs(lhs).max():.0e})  [{time.time()-t0_:.0f}s]")
rng = np.random.default_rng(0); X = rng.normal(size=(N2, 30)) + 1j*rng.normal(size=(N2, 30))
def probe(t0):
    nn_ = [np.linalg.norm(Sapply(1j*np.pi*(t0 + d)/T, X)) for d in (1e-5, 1e-7)]
    sv = np.linalg.svd(1e-7*Sfull(1j*np.pi*(t0 + 1e-7)/T), compute_uv=False); return np.log10(nn_[1]/nn_[0])/2, int(np.sum(sv > 1e-6*sv[0]))
for t0, lab in [(A[0], "27+27 -> 378"), (A[1], "27+27 -> 79"), (A[2], "27+27 -> 27"), (T - 1, "lowest breather")]:
    od, rk = probe(t0); print(f"   t = {t0:7.3f} {lab:16s} order {od:.2f} rank {rk}")
np.save("f4sol_state.npy", np.array([m1, T]), allow_pickle=True)
print(f"[{time.time()-t0_:.0f}s]")
