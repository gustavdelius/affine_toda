"""Coordinator check of varsigma = (-1)^n s: zeros of the corner entry of K (book intertwiner code) and U_varsigma."""
import numpy as np
from boundary_common import dn_vector, q_phys, coideal_gen, intertwiners
def Kraw(n, om, th, s):
    q = q_phys(om); eps = np.ones(n + 1)/(2*np.cos(np.pi*om/2)); eps[0] *= s; y = np.exp(om*th)
    pairs = [(coideal_gen(dn_vector(n, q, y), j, 1/q, eps[j]), coideal_gen(dn_vector(n, q, 1/y), j, 1/q, eps[j])) for j in range(n + 1)]
    sols = intertwiners(pairs)[0]; assert len(sols) == 1
    return sols[0]
om = 0.37
for n in (4, 5, 6):
    h = 2*n - 2; T = (n - 1)*om
    for s in (1, -1):
        out = []
        for sg in (1, -1):   # theta with y^h = sg
            th = (0 if sg == 1 else 1j*np.pi)/(om*h) + 0j
            th = th if sg == -1 else 2j*np.pi/(om*h)
            K = Kraw(n, om, th + 1e-5, s)
            out.append(abs(K[-1, 0]/K[0, 0]) + abs(K[0, -1]/K[0, 0]))
        # unitarity with the book normalization
        th = 0.31 + 0.17j; mu = -1j*T*th/np.pi; vs = (-1)**n*s
        nrm = lambda t, m: (-np.sin(np.pi*T/2)/np.sin(np.pi*m) if vs == 1 else -1j*np.sin(np.pi*T/2)/np.cos(np.pi*m))
        Kp = Kraw(n, om, th, s); Km = Kraw(n, om, -th, s)
        P = (Km/Km[0, 0]*nrm(-th, -mu)) @ (Kp/Kp[0, 0]*nrm(th, mu))
        c2 = np.cos(2*np.pi*mu)
        U = vs*(c2 - vs*np.cos(np.pi*T))*(c2 - vs*np.cos(np.pi*om))/(2*np.sin(np.pi*om/2)**2*(c2 - vs))
        print(f'n={n} s={s:+d}: |corner/K11| at y^h=+1: {out[0]:.2e}, at y^h=-1: {out[1]:.2e};  '
              f'|K(-th)K(th) - U_vs 1| = {np.abs(P - U*np.eye(P.shape[0])).max():.1e}')
