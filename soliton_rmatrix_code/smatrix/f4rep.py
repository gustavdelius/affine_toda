"""26-dim representation of U_q(f_4^(1)) (symmetry of imaginary e_6^(2) Toda solitons), solved numerically from the relations"""
import numpy as np, itertools, os, time
from scipy.optimize import least_squares
omega = float(os.environ.get('OMEGA', '2.37')); q = float(__import__("os").environ.get("QSIGN", "1"))*np.exp(-1j*np.pi*omega)
e = np.eye(4)
simple = {1: e[1]-e[2], 2: e[2]-e[3], 3: e[3], 4: 0.5*(e[0]-e[1]-e[2]-e[3])}     # Bourbaki F4: 1-2=>3-4, alpha1,2 long
theta = e[0] + e[1]; simple[0] = -theta
short = [s*e[i] for i in range(4) for s in (1, -1)] + [0.5*np.array(sg) for sg in itertools.product((1, -1), repeat=4)]
W = np.array(short + [np.zeros(4), np.zeros(4)]); n = len(W)                     # 24 short roots + 2 zero weights
zero = [24, 25]
nrm = {i: np.dot(simple[i], simple[i]) for i in simple}
qi = {i: q**nrm[i] for i in simple}
def K(i): return np.diag(q**(2*W @ simple[i])).astype(complex)
def a_ij(i, j): return int(round(2*np.dot(simple[i], simple[j])/nrm[i]))
def pairs(i):          # (target, source): weight(target) = weight(source) + alpha_i
    return [(t, s) for s in range(n) for t in range(n) if np.allclose(W[t], W[s] + simple[i])]
P = {i: pairs(i) for i in simple}
print("allowed entries per generator:", {i: len(P[i]) for i in simple})
def build(vec, nodes):
    E, F = {}, {}; k = 0
    for i in nodes:
        E[i] = np.zeros((n, n), complex); F[i] = np.zeros((n, n), complex)
        for (t, s) in P[i]: E[i][t, s] = vec[k] + 1j*vec[k+1]; k += 2
        for (t, s) in P[i]: F[i][s, t] = vec[k] + 1j*vec[k+1]; k += 2
    return E, F
def qb(m, Q): return (Q**m - Q**(-m))/(Q - 1/Q)
def serre(Ei, Ej, m, Q):
    from math import comb
    def qf(r): return np.prod([qb(s, Q) for s in range(1, r+1)]) if r > 0 else 1.0
    return sum((-1)**r*qf(m)/(qf(r)*qf(m-r))*np.linalg.matrix_power(Ei, m-r) @ Ej @ np.linalg.matrix_power(Ei, r) for r in range(m+1))
def residual(E, F, nodes, fixedE=None, fixedF=None):
    EE = dict(E); FF = dict(F)
    if fixedE: EE.update(fixedE); FF.update(fixedF)
    res = []
    allnodes = sorted(EE)
    for i in allnodes:
        for j in allnodes:
            if i not in nodes and j not in nodes: continue
            c = EE[i] @ FF[j] - FF[j] @ EE[i]
            if i == j: c = c - (K(i) - np.linalg.inv(K(i)))/(qi[i] - 1/qi[i])
            res.append(c.ravel())
            if i != j:
                m = 1 - a_ij(i, j)
                res.append(serre(EE[i], EE[j], m, qi[i]).ravel()); res.append(serre(FF[i], FF[j], m, qi[i]).ravel())
    r = np.concatenate(res); return np.concatenate([r.real, r.imag])
def solve(nodes, fixedE=None, fixedF=None, tries=8, seed=0):
    rng = np.random.default_rng(seed); m = sum(4*len(P[i]) for i in nodes); best = None
    for t in range(tries):
        x0 = rng.normal(size=m)
        sol = least_squares(lambda v: residual(*build(v, nodes), nodes, fixedE, fixedF), x0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        r = np.abs(sol.fun).max()
        if best is None or r < best[0]: best = (r, sol.x)
        if r < 1e-12: break
    return best
if __name__ == "__main__":
    t0 = time.time()
    r, v = solve([1, 2, 3, 4])
    E, F = build(v, [1, 2, 3, 4])
    print(f"U_q(f_4) on 26: max relation residual {r:.1e}  [{time.time()-t0:.0f}s]")
    r0, v0 = solve([0], E, F)
    E0, F0 = build(v0, [0])
    print(f"affine node 0 (alpha_0 = -theta, Cartan row 0 = {[a_ij(0, j) for j in range(5)]}): max residual {r0:.1e}  [{time.time()-t0:.0f}s]")
    E.update(E0); F.update(F0)
    np.save(f"f4rep_w{omega}.npy", np.array([E, F], dtype=object), allow_pickle=True)
