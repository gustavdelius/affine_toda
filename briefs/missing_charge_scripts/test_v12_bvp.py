"""Inspect the lowest-Re(S) critical points in missing sector (1,2) for L=1, 2, 3 and test the continuum BVP in V_(1,2)."""
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../textbook_code/thimbles')))

import numpy as np
from scipy.integrate import solve_bvp
from lie import An
from model import Chain
from newton_batch import find_cps

g = An(3)
w_reps = g.fundamental_weights_rep(2)
lam12 = w_reps[(1, 2)]

# 1. Inspect lowest Re(S) critical points in sector (1,2) for L=1, 2, 3 at kappa=1.0
for L in (1, 2, 3):
    ch = Chain(g, L, 1.0, lam12)
    cps, cnt = find_cps(ch, nstart=8000, re_spread=3.5, im_spread=1.5)
    S = np.array([ch.action(u) for u in cps])
    o = np.argsort(S.real)
    print(f"\n=== Sector (1,2), L={L}, kappa=1.0: #cps={len(cps)} ===")
    for idx in o[:6]:
        U = cps[idx]
        ue = np.array([g.eps_coords(U[x]) / (2 * np.pi) for x in range(L)])
        print(f"  S = {S[idx].real:9.5f}{S[idx].imag:+9.5f}j")
        for x in range(L):
            print(f"    u[{x+1}]/(2pi) = {np.round(ue[x], 4)}")

# 2. Solve continuum BVP in V_(1,2) (shows mesh blowup because tau_1 = tau_3 in R must cross 0 on R):
print("\n=== Solving continuum BVP for (alpha, beta) in V_(1,2) ===")
ell = 16.0
x = np.linspace(-ell / 2, ell / 2, 401)
def ode(x, y):
    al, dal, be, dbe = y
    ddal = np.exp(2 * be) * np.sin(2 * al)
    ddbe = np.exp(2 * be) * np.cos(2 * al) - np.exp(-2 * be)
    return np.vstack([dal, ddal, dbe, ddbe])

def bc(ya, yb):
    return np.array([ya[0] - 0.0, yb[0] - np.pi, ya[2] - 0.0, yb[2] - 0.0])

for sign_be in (+1.0, -1.0, +0.3, -0.3):
    y_init = np.zeros((4, x.size))
    y_init[0] = 2.0 * np.arctan(np.exp(np.sqrt(2.0) * x))
    y_init[1] = 2.0 * np.sqrt(2.0) / (2.0 * np.cosh(np.sqrt(2.0) * x))
    y_init[2] = sign_be / np.cosh(x)
    res = solve_bvp(ode, bc, x, y_init, tol=1e-8, max_nodes=5000)
    print(f"  sign_be={sign_be:+.1f}: success={res.success} ({res.message})")
