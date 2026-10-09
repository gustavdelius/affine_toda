"""Shared tools: fast projectors via highest-weight vectors, spinor reps, particle-labelled eigenvalue scans."""
import numpy as np, itertools
from scipy.linalg import null_space, orth
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def site(op, i, N):
    mats = [I2]*N; mats[i] = op; out = mats[0]
    for m in mats[1:]: out = np.kron(out, m)
    return out
def spin_weights(N): return np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=N)])
def spinor_rep(N, q, x, alg='c'):
    """alg='c': U_q(d_{N+1}^(2)) spinor (c_N^(1) solitons); alg='a': U_q(b_N^(1)) spinor (a_{2N-1}^(2) solitons)"""
    Wn = spin_weights(N); e, f, k, qi = {}, {}, {}, {}
    for i in range(1, N):
        e[i] = site(sp, i-1, N) @ site(sm, i, N); f[i] = e[i].T.copy(); k[i] = np.diag(q**(2*(Wn[:, i-1] - Wn[:, i]))).astype(complex); qi[i] = q**2
    e[N] = site(sp, N-1, N); f[N] = e[N].T.copy(); k[N] = np.diag(q**(2*Wn[:, N-1])).astype(complex); qi[N] = q
    if alg == 'c':
        e[0] = x*site(sm, 0, N); f[0] = site(sp, 0, N)/x; k[0] = np.diag(q**(-2*Wn[:, 0])).astype(complex); qi[0] = q
    else:
        e[0] = x*site(sm, 0, N) @ site(sm, 1, N); f[0] = site(sp, 0, N) @ site(sp, 1, N)/x
        k[0] = np.diag(q**(-2*(Wn[:, 0] + Wn[:, 1]))).astype(complex); qi[0] = q**2
    return e, f, k, qi
def flip(Da, Db):
    """P: V_a (x) V_b -> V_b (x) V_a"""
    P = np.zeros((Da*Db, Da*Db))
    for i in range(Da):
        for j in range(Db): P[j*Da + i, i*Db + j] = 1
    return P
def components(E, F, K, W, nodes, simple):
    """multiplicity-free decomposition of V (x) V under the finite algebra; returns list of (dim, projector, hw weight)"""
    n = W.shape[0]; I = np.eye(n)
    dE = {i: np.kron(E[i], I) + np.kron(K[i], E[i]) for i in nodes}
    dF = {i: np.kron(F[i], np.linalg.inv(K[i])) + np.kron(I, F[i]) for i in nodes}
    WW = (W[:, None, :] + W[None, :, :]).reshape(n*n, -1)
    key = lambda w: tuple(np.round(np.asarray(w)*4).astype(int))
    wsp = {}
    for idx, w in enumerate(WW): wsp.setdefault(key(w), []).append(idx)
    comps = []
    for kk, idx in wsp.items():
        M = np.vstack([dE[i][:, idx] for i in nodes]); ns = null_space(M, rcond=1e-10)
        for c in range(ns.shape[1]):
            v = np.zeros(n*n, complex); v[idx] = ns[:, c]; comps.append((kk, v))
    sk = {i: key(simple[i]) for i in nodes}
    spaces_all = []
    for k0, v in comps:
        spaces = {k0: v[:, None]/np.linalg.norm(v)}; frontier = [k0]
        while frontier:
            new = []
            for kk in frontier:
                B = spaces[kk]
                for i in nodes:
                    k2 = tuple(a - b for a, b in zip(kk, sk[i]))
                    if k2 not in wsp: continue
                    Y = dF[i] @ B
                    if np.linalg.norm(Y) < 1e-12: continue
                    cur = spaces.get(k2); Z = Y if cur is None else np.hstack([cur, Y]); O = orth(Z, rcond=1e-10)
                    if cur is None or O.shape[1] > cur.shape[1]: spaces[k2] = O; new.append(k2)
            frontier = list(set(new))
        spaces_all.append(spaces)
    Bs = [np.hstack([s[kk] for kk in sorted(s)]) for s in spaces_all]
    Ball = np.hstack(Bs); Binv = np.linalg.inv(Ball); out = []; off = 0
    for (k0, v), B in zip(comps, Bs):
        d = B.shape[1]; out.append((d, B @ Binv[off:off+d], k0)); off += d
    return out
def br(x, l, s, q): return (x - s*q**l)/(1 - s*x*q**l)
def maxdev_scan(Mfun, lxs, signs=(1,)):
    """max_x ||s|-1| of eig(Mfun(x)) over x = sign*exp(lx)"""
    best = (0, None, None)
    for sg in signs:
        for lx in lxs:
            e = np.linalg.eigvals(Mfun(sg*np.exp(lx))); d = np.abs(np.abs(e) - 1).max()
            if d > best[0]: best = (d, sg, lx)
    return best
def krein_err(e):
    return max(np.min(np.abs(e - 1/np.conj(s))) for s in e)
def transfer3(Rpl, D, l12, l13):
    """T_1 = R_13(x13) R_12(x12) on V^{x3}; Rpl(x) particle-labelled D^2 x D^2; l = log x differences"""
    I = np.eye(D); P23 = np.kron(I, flip(D, D))
    R12 = np.kron(Rpl(np.exp(l12)), I)
    R13 = P23 @ np.kron(Rpl(np.exp(l13)), I) @ P23
    return R13 @ R12
def scan3(Rpl, D, L=6.0, n=400, seed=0):
    rng = np.random.default_rng(seed); worst = (0, None); kre = 0
    for _ in range(n):
        l12, l13 = rng.uniform(-L, L, 2)
        e = np.linalg.eigvals(transfer3(Rpl, D, l12, l13)); d = np.abs(np.abs(e) - 1).max()
        if d > worst[0]: worst = (d, (l12, l13)); kre = krein_err(e)
    return worst, kre
def jimbo(n, q, x):
    """Jimbo U_q(sl_n^(1)) R on C^n x C^n homogeneous gradation, normalised; particle-labelled (R, not P R)"""
    E = lambda i, j: (lambda M: (M.__setitem__((i, j), 1), M)[1])(np.zeros((n, n)))
    R = np.zeros((n*n, n*n), complex)
    for i in range(n): R += (x*q - 1/(x*q))*np.kron(E(i, i), E(i, i))
    for i in range(n):
        for j in range(n):
            if i != j:
                R += (x - 1/x)*np.kron(E(i, i), E(j, j)) + (q - 1/q)*x**np.sign(j - i)*np.kron(E(i, j), E(j, i))
    return R/(x*q - 1/(x*q))
def apply_pair(M, v, i, j, D, N):
    """apply D^2 x D^2 matrix M to tensor factors (i, j) of v (shape D^N x m)"""
    m = v.shape[1]; t = v.reshape((D,)*N + (m,))
    t = np.moveaxis(t, [i, j], [0, 1]); sh = t.shape
    t = (M @ t.reshape(D*D, -1)).reshape(sh)
    return np.moveaxis(t, [0, 1], [i, j]).reshape(D**N, m)
def transferN(Rpl, D, ls):
    """T_1 = R_{1N}(x_{1N}) ... R_{12}(x_{12}); ls = log x_{1j}, j=2..N"""
    N = len(ls) + 1; T = np.eye(D**N, dtype=complex)
    for j, l in enumerate(ls, start=1): T = apply_pair(Rpl(np.exp(l)), T, 0, j, D, N)
    return T
def scanN(Rpl, D, N, L=6.0, n=200, seed=0):
    rng = np.random.default_rng(seed); worst = (0, None); kre = 0
    for _ in range(n):
        ls = rng.uniform(-L, L, N - 1)
        e = np.linalg.eigvals(transferN(Rpl, D, ls)); d = np.abs(np.abs(e) - 1).max()
        if d > worst[0]: worst = (d, ls); kre = krein_err(e)
    return worst, kre
