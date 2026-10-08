"""27-dim representation of U_q(e_6^(2)) (symmetry of imaginary f_4^(1) Toda solitons): f_4 part on 26 + singlet; node 0 = -theta_short"""
import numpy as np, itertools, os, time
from scipy.optimize import least_squares
omega = float(os.environ.get('OMEGA', '2.37')); q = np.exp(-1j*np.pi*omega)
E26, F26 = np.load(f"f4rep_w{omega}.npy", allow_pickle=True)
e = np.eye(4)
simple = {1: e[1]-e[2], 2: e[2]-e[3], 3: e[3], 4: 0.5*(e[0]-e[1]-e[2]-e[3]), 0: -e[0]}
short = [s*e[i] for i in range(4) for s in (1, -1)] + [0.5*np.array(sg) for sg in itertools.product((1, -1), repeat=4)]
W = np.array(short + [np.zeros(4)]*3); n = 27
nrm = {i: np.dot(simple[i], simple[i]) for i in simple}; qi = {i: q**nrm[i] for i in simple}
def K(i): return np.diag(q**(2*W @ simple[i])).astype(complex)
def a_ij(i, j): return int(round(2*np.dot(simple[i], simple[j])/nrm[i]))
E, F = {}, {}
for i in range(1, 5):
    E[i] = np.zeros((n, n), complex); F[i] = np.zeros((n, n), complex); E[i][:26, :26] = E26[i]; F[i][:26, :26] = F26[i]
P0 = [(t, s) for s in range(n) for t in range(n) if np.allclose(W[t], W[s] + simple[0])]
def build0(v):
    E0 = np.zeros((n, n), complex); F0 = np.zeros((n, n), complex); k = 0
    for (t, s) in P0: E0[t, s] = v[k] + 1j*v[k+1]; k += 2
    for (t, s) in P0: F0[s, t] = v[k] + 1j*v[k+1]; k += 2
    return E0, F0
def qb(m, Q): return (Q**m - Q**(-m))/(Q - 1/Q)
def serre(Ei, Ej, m, Q):
    def qf(r): return np.prod([qb(s, Q) for s in range(1, r+1)]) if r > 0 else 1.0
    return sum((-1)**r*qf(m)/(qf(r)*qf(m-r))*np.linalg.matrix_power(Ei, m-r) @ Ej @ np.linalg.matrix_power(Ei, r) for r in range(m+1))
def residual_full(EE, FF, nodes=range(5)):
    res = []
    for i in nodes:
        for j in nodes:
            c = EE[i] @ FF[j] - FF[j] @ EE[i]
            if i == j: c = c - (K(i) - np.linalg.inv(K(i)))/(qi[i] - 1/qi[i])
            res.append(c.ravel())
            if i != j:
                m = 1 - a_ij(i, j); res.append(serre(EE[i], EE[j], m, qi[i]).ravel()); res.append(serre(FF[i], FF[j], m, qi[i]).ravel())
            res.append((K(i) @ EE[j] @ np.linalg.inv(K(i)) - qi[i]**a_ij(i, j)*EE[j]).ravel())
    r = np.concatenate(res); return np.concatenate([r.real, r.imag])
def res0(v):
    E0, F0 = build0(v); EE = dict(E); FF = dict(F); EE[0] = E0; FF[0] = F0
    rr = []
    for j in range(5):
        c = EE[0] @ FF[j] - FF[j] @ EE[0]
        if j == 0: c = c - (K(0) - np.linalg.inv(K(0)))/(qi[0] - 1/qi[0])
        rr.append(c.ravel())
        if j != 0:
            rr.append((EE[j] @ FF[0] - FF[0] @ EE[j]).ravel())
            m = 1 - a_ij(0, j); rr.append(serre(EE[0], EE[j], m, qi[0]).ravel()); rr.append(serre(FF[0], FF[j], m, qi[0]).ravel())
            m2 = 1 - a_ij(j, 0); rr.append(serre(EE[j], EE[0], m2, qi[j]).ravel()); rr.append(serre(FF[j], FF[0], m2, qi[j]).ravel())
    r = np.concatenate(rr); return np.concatenate([r.real, r.imag])
if __name__ == "__main__":
    t0 = time.time(); rng = np.random.default_rng(1); best = None
    print("Cartan row 0:", [a_ij(0, j) for j in range(5)], " row 4:", [a_ij(4, j) for j in range(5)], f"; {len(P0)} entries for e_0")
    for tr in range(10):
        sol = least_squares(res0, rng.normal(size=4*len(P0)), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        r = np.abs(sol.fun).max()
        if best is None or r < best[0]: best = (r, sol.x)
        if r < 1e-12: break
    E[0], F[0] = build0(best[1])
    full = np.abs(residual_full(E, F)).max()
    sing = np.abs(E[0][26]).max() + np.abs(E[0][:, 26]).max()
    print(f"U_q(e_6^(2)) on 27: node-0 residual {best[0]:.1e}; all relations {full:.1e}; singlet coupled by e_0: {sing > 1e-8}  [{time.time()-t0:.0f}s]")
    np.save(f"e6rep_w{omega}.npy", np.array([E, F], dtype=object), allow_pickle=True)
