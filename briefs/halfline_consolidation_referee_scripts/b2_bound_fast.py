"""Bound states of the linearized problem about the d_n C = +1 vacuum: zeros of
D(Omega^2) = det[eps'(0) + R eps(0)] over the solutions that decay at -infinity (all channels closed),
for -4 < Omega^2 < m_min^2.  Starts normalized to 1 at x = -L, moderate tolerances (faster version of
bound_det in b2_dn_classical.py).  Also: does the d_5 particle-3 pole at Omega = 1 come with a zero of D?"""
import numpy as np
from scipy.integrate import solve_ivp
import b2_dn_classical as b2

def D(n, psi, p0, Om2, L=14.0):
    al, mk = b2.dn_roots(n)
    H0 = (al.T*mk) @ al
    w, U = np.linalg.eigh(H0)
    R = 0.5*(al.T*(np.sqrt(mk)*np.exp(al @ p0/2))) @ al
    def rhs(x, y):
        M2 = (al.T*(mk*np.exp(al @ psi(x)[0]))) @ al
        return np.concatenate([y[n:], (M2 - Om2*np.eye(n)) @ y[:n]])
    cols = []
    for b in range(n):
        kap = np.sqrt(w[b] - Om2)
        sol = solve_ivp(rhs, (-L, 0.0), np.concatenate([U[:, b], kap*U[:, b]]), method='DOP853', rtol=1e-10, atol=1e-12)
        yf = sol.y[:, -1]
        v = yf[n:] + R @ yf[:n]
        cols.append(v/np.linalg.norm(yf))
    return np.linalg.det(np.array(cols).T)

if __name__ == '__main__':
  for n in (4, 5):
      psi, p0 = b2.vacuum(n)
      al, mk = b2.dn_roots(n)
      m2min = np.min(np.linalg.eigvalsh((al.T*mk) @ al))
      g = np.concatenate([np.linspace(-4.0, m2min - 1e-3, 45), np.linspace(0.95, 1.05, 11) if n == 5 else []])
      g = np.sort(g[g < m2min])
      d = np.array([D(n, psi, p0, x) for x in g])
      ch = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1])))
      print(f"d_{n}: {ch} sign changes of D on -4 < Omega^2 < m_min^2 = {m2min:.4f}; min |D| = {np.abs(d).min():.3e}"
            + (f"; D near Omega^2 = 1: {np.round(d[np.abs(g-1) < 0.06], 4)}" if n == 5 else ''), flush=True)
