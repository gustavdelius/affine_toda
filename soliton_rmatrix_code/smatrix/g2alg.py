"""U_q(g_2^(1)) on 7 (ALG=g2: d_4^(3) solitons) and U_q(d_4^(3)) on 8 = 7+1 (ALG=d43: g_2^(1) solitons)"""
import numpy as np, os, time
from scipy.optimize import least_squares
ALG = os.environ.get('ALG', 'g2'); omega = float(os.environ.get('OMEGA', '2.37')); q = np.exp(-1j*np.pi*omega)
a1 = np.array([1.0, 0.0]); a2 = np.array([-1.5, np.sqrt(3)/2])                     # alpha1 short (norm 1), alpha2 long (norm 3)
shorts = [s*v for v in (a1, a1 + a2, 2*a1 + a2) for s in (1, -1)]
theta, theta_s = 3*a1 + 2*a2, 2*a1 + a2
simple = {1: a1, 2: a2, 0: -theta if ALG == 'g2' else -theta_s}
W = np.array(shorts + [np.zeros(2)] + ([np.zeros(2)] if ALG == 'd43' else [])); n = len(W)
nrm = {i: np.dot(simple[i], simple[i]) for i in simple}; qi = {i: q**nrm[i] for i in simple}
def K(i): return np.diag(q**(2*W @ simple[i])).astype(complex)
def a_ij(i, j): return int(round(2*np.dot(simple[i], simple[j])/nrm[i]))
def pairs(i): return [(t, s) for s in range(n) for t in range(n) if np.allclose(W[t], W[s] + simple[i])]
P = {i: pairs(i) for i in simple}
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
def solve(active, fixed=({}, {}), seed=0):
    rng = np.random.default_rng(seed); m = sum(4*len(P[i]) for i in active); best = None
    for t in range(12):
        f = lambda v: resid({**fixed[0], **build(v, active)[0]}, {**fixed[1], **build(v, active)[1]}, sorted(set(active) | set(fixed[0])), active)
        sol = least_squares(f, rng.normal(size=m), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        r = np.abs(sol.fun).max()
        if best is None or r < best[0]: best = (r, sol.x)
        if r < 1e-12: break
    return best
cache = f"{ALG}rep_w{omega}.npy"
if os.path.exists(cache): E, F = np.load(cache, allow_pickle=True)
else:
    r, v = solve([1, 2]); E, F = build(v, [1, 2])
    r0, v0 = solve([0], (E, F)); E0, F0 = build(v0, [0]); E.update(E0); F.update(F0)
    print(f"{ALG}: {n}-dim rep; finite residual {r:.1e}, affine node residual {r0:.1e}; Cartan row 0 {[a_ij(0, j) for j in range(3)]}")
    np.save(cache, np.array([E, F], dtype=object), allow_pickle=True)
def rep(x):
    return ({i: (x if i == 0 else 1)*E[i] for i in range(3)}, {i: F[i]/(x if i == 0 else 1) for i in range(3)}, {i: K(i) for i in range(3)}, dict(qi))
if __name__ == "__main__":
    from rsolve import intertwiner
    t0 = time.time(); xs = 1.3*np.exp(2j*np.pi*(np.arange(24) + 0.17)/24)
    Rs = []
    for x in xs:
        R_ = intertwiner(rep(x), rep(1.0), range(3), W, W)[0]; Rs.append(R_/R_.ravel()[np.argmax(np.abs(R_.ravel()))])
    U_, sv, _ = np.linalg.svd(np.array([R_.ravel() for R_ in Rs]).T, full_matrices=False)
    dimC = int(np.sum(sv > 1e-9*sv[0])); Bc = [U_[:, k].reshape(n*n, n*n) for k in range(dimC)]
    print(f"commutant dimension {dimC}  [{time.time()-t0:.0f}s]")
    np.save(f"{ALG}_comm_w{omega}.npy", np.array(Bc))
