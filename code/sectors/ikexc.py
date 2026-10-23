import numpy as np
from lib import *
exec(open('a22three.py').read().split('q = np.exp(-1j*0.77)')[0])
m = np.array([1, 0, -1]); g = np.abs(m)
Qg = np.diag([(-1.0)**(g[i]*g[j]) for i in range(3) for j in range(3)])
lx = np.linspace(0.02, 14, 300)
for xi_pi in (0.3, 0.45, 0.6, 0.9, 1.5, 3.0):
    xi = (xi_pi + 0.0013)*np.pi; q = np.conj(1j*np.exp(1j*np.pi**2/(3*xi))); Rp = setup(q)
    dp = max(np.abs(np.abs(np.linalg.eigvals(Rp(-np.exp(l)))) - 1).max() for l in lx)
    dg = max(np.abs(np.abs(np.linalg.eigvals(Qg @ Rp(-np.exp(l)))) - 1).max() for l in lx)
    e = np.linalg.eigvals(Qg @ Rp(-np.exp(1.0)))
    print(f"xi={xi_pi}pi: K0-K_odd ~ R(-x): plain {dp:.1e}, graded {dg:.1e}, Krein err (graded, logx=1) {krein_err(e):.0e}")
