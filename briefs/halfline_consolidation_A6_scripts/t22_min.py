"""Direct minimization of the (convex) C = +1 energy, discretized, in the full field space (no symmetry assumed) for d_4, d_5,
and in the invariant subspace for d_6; gives the vacuum energy and boundary values, and initial data for shooting."""
import sys
import numpy as np
from scipy.optimize import minimize
from cl_lin import droots

n = int(sys.argv[1])
full = len(sys.argv) > 2 and sys.argv[2] == 'full'
SAVEZ = True
A, nj = droots(n)
if full:
    P = np.eye(n)
else:
    N1 = n + 1
    orb = {}
    for j in range(N1):
        jj = min(j, n - j)
        orb.setdefault('end' if jj <= 1 else jj, []).append(j)
    Rm = np.array([A[o[0]] - A[j] for o in orb.values() for j in o[1:]])
    _, sv, Vt = np.linalg.svd(Rm)
    P = Vt[int(np.sum(sv > 1e-10)):].T
d = P.shape[1]
AP = A @ P
L, M = 16.0, 3200
# graded grid, finer near the wall
s = np.linspace(0, 1, M + 1)
xs = -L*(1 - s)**1.5
dx = np.diff(xs)
w = np.zeros(M + 1); w[:-1] += dx/2; w[1:] += dx/2
sq = np.sqrt(nj)


def energy(zf):
    z = zf.reshape(M + 1, d)
    dz = np.diff(z, axis=0)
    kin = 0.5*np.sum(dz**2/dx[:, None])
    e = np.exp(z @ AP.T)                 # (M+1, n+1)
    pot = np.sum(w*((e - 1) @ nj))
    eb = np.exp(AP @ z[-1]/2)
    bnd = 2*np.sum(sq*(eb - 1))
    g = np.zeros_like(z)
    g[:-1] -= dz/dx[:, None]; g[1:] += dz/dx[:, None]
    g += (w[:, None]*((e*nj) @ AP))
    g[-1] += (sq*eb) @ AP
    return kin + pot + bnd, g.ravel()


rng = np.random.default_rng(1)
res = []
for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 6):
    z0 = rng.normal(scale=0.5, size=(M + 1, d))*np.exp(xs/3)[:, None]
    r = minimize(energy, z0.ravel(), jac=True, method='L-BFGS-B', options=dict(maxiter=20000, ftol=1e-15, gtol=1e-10))
    z = r.x.reshape(M + 1, d)
    res.append((r.fun, P @ z[-1]))
    np.save(f'minprof_d{n}_{trial}.npy', np.c_[xs, z @ P.T])
    print(f'd_{n} ({"full" if full else "invariant"} space) start {trial}: E = {r.fun:.10f}, psi(0) = {np.round(P @ z[-1], 6)}, '
          f'|grad| = {np.abs(r.jac).max():.1e}', flush=True)
np.save(f'min_d{n}_{"full" if full else "inv"}.npy', np.array([res[int(np.argmin([t[0] for t in res]))][1]]))
