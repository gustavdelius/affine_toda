"""U_q(d_4^(2)) on the 8 = 7 + 1 (c_3^(1) soliton 1), built from the relations; R-matrix poles and crossing point"""
import numpy as np, os, time
from scipy.optimize import least_squares
omega = float(os.environ.get('OMEGA', '2.37')); q = float(__import__("os").environ.get("QSIGN", "1"))*np.exp(-1j*np.pi*omega)
e3 = np.eye(3)
simple = {1: e3[0]-e3[1], 2: e3[1]-e3[2], 3: e3[2], 0: -e3[0]}
W = np.array([s*e3[i] for i in range(3) for s in (1, -1)] + [np.zeros(3), np.zeros(3)]); n = 8
nrm = {i: np.dot(simple[i], simple[i]) for i in simple}; qi = {i: q**nrm[i] for i in simple}
def K(i): return np.diag(q**(2*W @ simple[i])).astype(complex)
def a_ij(i, j): return int(round(2*np.dot(simple[i], simple[j])/nrm[i]))
P = {i: [(t, s) for s in range(n) for t in range(n) if np.allclose(W[t], W[s] + simple[i])] for i in simple}
def qb(m, Q): return (Q**m - Q**(-m))/(Q - 1/Q)
def serre(Ei, Ej, m, Q):
    def qf(r): return np.prod([qb(s, Q) for s in range(1, r+1)]) if r > 0 else 1.0
    return sum((-1)**r*qf(m)/(qf(r)*qf(m-r))*np.linalg.matrix_power(Ei, m-r) @ Ej @ np.linalg.matrix_power(Ei, r) for r in range(m+1))
def build(v, nodes):
    E, F = {}, {}; k = 0
    for i in nodes:
        E[i] = np.zeros((n, n), complex); F[i] = np.zeros((n, n), complex)
        for (t, s) in P[i]: E[i][t, s] = v[k] + 1j*v[k+1]; k += 2
        for (t, s) in P[i]: F[i][s, t] = v[k] + 1j*v[k+1]; k += 2
    return E, F
def resid(EE, FF, nodes, active):
    res = []
    for i in nodes:
        for j in nodes:
            if i not in active and j not in active: continue
            c = EE[i] @ FF[j] - FF[j] @ EE[i]
            if i == j: c = c - (K(i) - np.linalg.inv(K(i)))/(qi[i] - 1/qi[i])
            res.append(c.ravel())
            if i != j:
                m = 1 - a_ij(i, j); res.append(serre(EE[i], EE[j], m, qi[i]).ravel()); res.append(serre(FF[i], FF[j], m, qi[i]).ravel())
    r = np.concatenate(res); return np.concatenate([r.real, r.imag])
def solve(active, fixed=({}, {})):
    rng = np.random.default_rng(0); m = sum(4*len(P[i]) for i in active); best = None
    for t in range(12):
        f = lambda v: resid({**fixed[0], **build(v, active)[0]}, {**fixed[1], **build(v, active)[1]}, sorted(set(active) | set(fixed[0])), active)
        sol = least_squares(f, rng.normal(size=m), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        r = np.abs(sol.fun).max()
        if best is None or r < best[0]: best = (r, sol.x)
        if r < 1e-12: break
    return best
cache = f"d42vec_w{omega}.npy"
if os.path.exists(cache): E, F = np.load(cache, allow_pickle=True)
else:
    r, v = solve([1, 2, 3]); E, F = build(v, [1, 2, 3])
    r0, v0 = solve([0], (E, F)); E0, F0 = build(v0, [0]); E.update(E0); F.update(F0)
    print(f"U_q(d_4^(2)) vector: finite residual {r:.1e}, node-0 residual {r0:.1e}; Cartan row 0 {[a_ij(0, j) for j in range(4)]}; singlet coupled: {np.abs(E[0][7]).max() + np.abs(E[0][:, 7]).max() + np.abs(E[0][6]).max() > 1e-8}")
    np.save(cache, np.array([E, F], dtype=object), allow_pickle=True)
def rep(x):
    return ({i: (x if i == 0 else 1)*E[i] for i in range(4)}, {i: F[i]/(x if i == 0 else 1) for i in range(4)}, {i: K(i) for i in range(4)}, dict(qi))
if __name__ == "__main__":
    from rsolve import intertwiner
    t0 = time.time()
    def Rn(x):
        R_ = intertwiner(rep(x), rep(1.0), range(4), W, W)[0]
        hwi = int(np.argmax(W @ np.array([4, 2, 1.]))); top = hwi*n + hwi
        return R_/R_[top, top]
    found = []
    for kph in range(4):
        eps = np.exp(1j*np.pi*kph/2)
        for l in range(0, 13):
            xs = eps*q**(-l)
            try:
                n1 = np.linalg.norm(Rn(xs*(1 + 1e-5))); n2 = np.linalg.norm(Rn(xs*(1 + 1e-7)))
            except Exception: continue
            if n2/n1 > 30:
                sv = np.linalg.svd(1e-7*Rn(xs*(1 + 1e-7)), compute_uv=False); found.append((['+', '+i', '-', '-i'][kph] + f"q^-{l}", int(np.sum(sv > 1e-6*sv[0]))))
    print("R_11 (top-normalised) poles and residue ranks:", found, f"[{time.time()-t0:.0f}s]")
    def mats(x): return [((x if i == 0 else 1)*E[i], F[i]/(x if i == 0 else 1), K(i)) for i in range(4)]
    xa = 0.83+0.29j
    for kph in range(4):
        eps = np.exp(1j*np.pi*kph/2)
        for l in range(-12, 13):
            tt = eps*q**l; G = 0
            for (Ei, Fi, Kk), (Ej, Fj, Kj) in zip(mats(xa), mats(xa*tt)):
                Kii = np.linalg.inv(Kk)
                for A_, B_ in (((-Kii@Ei).T, Ej), ((-Fi@Kk).T, Fj), (Kii.T, Kj)):
                    r_ = np.kron(np.eye(n), A_.T) - np.kron(B_, np.eye(n)); G = G + r_.conj().T @ r_
            ev_ = np.linalg.eigvalsh(G)
            if ev_[0] < 1e-13*ev_[-1]: print(f"crossing point: V*(x) = V(x * {['+', '+i', '-', '-i'][kph]}q^{l})")
