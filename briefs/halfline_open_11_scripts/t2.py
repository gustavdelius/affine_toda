import numpy as np
from sp_core import *
from boundary_common import q_phys
def Kbar_num(n, om, l, eta, th1, th2):
    q = q_phys(om); N = n+1
    V = lambda th: Vy(n, q, th, om); W = lambda th: Vy(n, q, th, om, True)
    K1 = KV(n, om, th1, l, eta)
    A1 = R(W(-th2), V(-th1)) @ one(K1, N) @ R(V(th1), W(-th2))
    B1 = R(W(th2), V(-th1)); C1 = one(K1, N) @ R(V(th1), W(th2))
    M = []
    for k in range(N):
        E = np.zeros((N, N)); E[k, k] = 1
        M.append((A1 @ one(E, N) - one(E, N) @ B1 @ C1).ravel())
    u, s, vh = np.linalg.svd(np.array(M).T)
    return vh[-1].conj()
for n in (2, 3, 4):
    N = n+1
    for l in range(1, N):
        for om in (0.37, 0.81):
            for eta in (0.3+0.4j, -0.2+1.1j):
                D = Kbar_num(n, om, l, eta, 0.31+0.12j, -0.17+0.05j)
                P = D*np.diag(KV(n, om, -0.17+0.05j, l, eta))
                r = P[l]/P[0]
                x = 0.5j*np.pi*om*(N - 2*l)
                cands = {'s(e+x)/s(e-x)': np.sinh(eta+x)/np.sinh(eta-x), 's(e-x)/s(e+x)': np.sinh(eta-x)/np.sinh(eta+x)}
                print(n, l, om, eta, 'r=', np.round(r, 6), {k: np.round(v/r, 6) for k, v in cands.items()})
