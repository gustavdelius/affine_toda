"""Lattice energy minimization for real-coupling d_n (or any simply-laced root data) Toda on x in [-X, 0], phi(-X) = 0.
E = sum 1/2 |dphi|^2/dx + sum w V(phi) + B(phi(0)), V = sum_j n_j (e^{a_j.phi} - 1), B = 2 sum_j C_j sqrt(n_j) (e^{a_j.phi/2} - 1)."""
import sys
import numpy as np
from scipy.optimize import minimize
from dn_classical import droots


def make_grid(X=14.0, N=1400, p=2.0):
    s = np.linspace(0, 1, N + 1)
    x = -X*(1 - s)**p          # dense near x = 0 ... actually dense near 0 for p > 1? (1-s)^p small near s=1
    return x


class Lattice:
    def __init__(self, al, nm, C, x):
        self.al, self.nm, self.C, self.x = al, nm, np.asarray(C, float), x
        self.dx = np.diff(x)
        w = np.zeros(len(x)); w[:-1] += self.dx/2; w[1:] += self.dx/2
        self.w = w
        self.r = al.shape[1]

    def unpack(self, z):
        phi = np.zeros((len(self.x), self.r)); phi[1:] = z.reshape(-1, self.r); return phi

    def energy(self, z):
        phi = self.unpack(z)
        d = np.diff(phi, axis=0)
        ap = phi @ self.al.T                       # (N+1, nodes)
        ex = np.exp(np.clip(ap, -700, 700))
        V = (ex - 1) @ self.nm
        eb = np.exp(ap[-1]/2)
        Bv = 2*np.sum(self.C*np.sqrt(self.nm)*(eb - 1))
        E = 0.5*np.sum(np.sum(d**2, axis=1)/self.dx) + np.sum(self.w*V) + Bv
        g = np.zeros_like(phi)
        gd = d/self.dx[:, None]
        g[:-1] -= gd; g[1:] += gd
        g += (self.w[:, None]*(ex*self.nm)) @ self.al
        g[-1] += (self.C*np.sqrt(self.nm)*eb) @ self.al
        return E, g[1:].ravel()


def run(n, C, starts=20, seed=0, X=14.0, N=1400, verbose=True, al=None, nm=None):
    if al is None:
        al, nm = droots(n)
    x = make_grid(X, N)
    lat = Lattice(al, nm, C, x)
    rng = np.random.default_rng(seed)
    res = []
    for k in range(starts):
        amp = rng.normal(size=al.shape[1])*rng.uniform(0.2, 3)
        prof = np.exp(1.5*x[1:])[:, None]*amp[None, :]
        r = minimize(lat.energy, prof.ravel(), jac=True, method='L-BFGS-B', options={'maxiter': 20000, 'gtol': 1e-10, 'maxfun': 40000})
        phi = lat.unpack(r.x)
        res.append((r.fun, phi[-1] @ al.T, phi))
        if verbose:
            print(f'  start {k}: E = {r.fun:.8f}  alpha.phi(0) = {np.round(phi[-1] @ al.T, 5)}', flush=True)
    res.sort(key=lambda t: t[0])
    return res, lat


if __name__ == '__main__':
    n = int(sys.argv[1]); pat = sys.argv[2]
    C = np.array([1.0 if ch == '+' else -1.0 if ch == '-' else 0.0 for ch in pat])
    res, lat = run(n, C, starts=int(sys.argv[3]) if len(sys.argv) > 3 else 12)
    print('best', res[0][0], res[0][1])
