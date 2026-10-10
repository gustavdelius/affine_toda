# N=5 kinks, Q=0 and Q=2 at alpha=pi (graded rule): follow the broken region from psi=0.42pi downwards by local maximisation
import sys, os, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from ik_common import *
Nn = 5; Qs = charges(Nn); m1 = m_first(Nn)
idx0 = np.where(Qs == 0)[0]; idx2 = np.where(Qs == 2)[0]
Om0 = np.exp(1j*np.pi*m1[idx0]); Om2 = np.exp(1j*np.pi*m1[idx2])
def dev(Rp, ls):
    T = transferN(Rp, 3, ls)
    d0 = np.abs(np.abs(np.linalg.eigvals(Om0[:, None]*T[np.ix_(idx0, idx0)])) - 1).max()
    d2 = np.abs(np.abs(np.linalg.eigvals(Om2[:, None]*T[np.ix_(idx2, idx2)])) - 1).max()
    return max(d0, d2)
best = np.array([-0.838, -1.012, -1.272, -0.871])
for psi_t in (0.42, 0.4175, 0.415, 0.4125, 0.41, 0.4075, 0.405, 0.4025, 0.401):
    xi_pi = 2/(3*psi_t) + 0.0013; psi = psi_of_xi(xi_pi); Rp = setup(q_of_xi(xi_pi))
    starts = [best] + [best + np.random.default_rng(k).normal(0, 0.15, 4) for k in range(4)]
    top = (dev(Rp, best), best)
    for s in starts:
        r = minimize(lambda v: -dev(Rp, v), s, method='Nelder-Mead', options={'maxfev': 400, 'xatol': 1e-4, 'fatol': 1e-10})
        if -r.fun > top[0]: top = (-r.fun, r.x)
    if top[0] > 1e-6: best = top[1]
    print(f"psi/pi={psi:+.4f}: max dev {top[0]:.2e} at rel. rapidities {np.round(top[1], 3)}", flush=True)
