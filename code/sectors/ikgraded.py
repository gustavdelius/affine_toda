import numpy as np
from lib import *
exec(open('a22three.py').read().split('q = np.exp(-1j*0.77)')[0])
m = np.array([1, 0, -1]); g = np.abs(m)
Qg = np.diag([(-1.0)**(g[i]*g[j]) for i in range(3) for j in range(3)])
Pg3 = Qg @ flip(3, 3)          # graded flip: charged kinks odd
def transfer_graded(Rpl, ls, D=3):
    Nn = len(ls) + 1; T = np.eye(D**Nn, dtype=complex)
    for j, l in enumerate(ls, start=1):
        v = T
        for k in range(j, 1, -1): v = apply_pair(Pg3, v, k-1, k, D, Nn)
        v = apply_pair(Rpl(np.exp(l)), v, 0, 1, D, Nn)
        for k in range(2, j+1): v = apply_pair(Pg3, v, k-1, k, D, Nn)
        T = v
    return T
lx = np.linspace(0.02, 14, 300)
for xi_pi in (0.30, 0.35, 0.45, 0.55, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3, 1.5, 2.0, 3.0, 6.0):
    xi = (xi_pi + 0.0013)*np.pi; q = np.conj(1j*np.exp(1j*np.pi**2/(3*xi))); Rp = setup(q)
    Rg = lambda x: Qg @ Rp(x)          # graded particle-labelled: Pg R-check = Qg P R-check
    psi = (2/(3*(xi_pi + 0.0013))) % 2; psi = min(psi, 2 - psi)
    d2p = max(np.abs(np.abs(np.linalg.eigvals(Rp(np.exp(l)))) - 1).max() for l in lx)
    d2g = max(np.abs(np.abs(np.linalg.eigvals(Rg(np.exp(l)))) - 1).max() for l in lx)
    rng = np.random.default_rng(2); d3g = 0
    for k in range(150):
        e = np.linalg.eigvals(transfer_graded(Rg, rng.uniform(-8, 8, 2))); d3g = max(d3g, np.abs(np.abs(e) - 1).max())
    print(f"xi={xi_pi:4.2f}pi |psi|/pi={psi:.3f}: 2-body plain {d2p:.1e}  graded {d2g:.1e};  3-body graded {d3g:.1e}", flush=True)
