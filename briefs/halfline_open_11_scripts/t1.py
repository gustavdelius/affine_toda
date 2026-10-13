# 1. check KV solves braid RE (VV); solve diagonal KVbar from mixed RE; check RE (VbarVbar)
import numpy as np
from sp_core import *
from boundary_common import q_phys
np.set_printoptions(precision=5, suppress=True, linewidth=200)
om = 0.37; q = q_phys(om)
for n in (2, 3):
    N = n + 1
    for l in range(0, N):
        eta = 0.3 + 0.4j
        th1, th2 = 0.31 + 0.12j, -0.17 + 0.05j
        V = lambda th: Vy(n, q, th, om); W = lambda th: Vy(n, q, th, om, True)
        K1, K2 = KV(n, om, th1, l, eta), KV(n, om, th2, l, eta)
        # braid RE: R_VV(th1-th2) K2(th1) R_VV(th1+th2) K2(th2) = K2(th2) R(th1+th2) K2(th1) R(th1-th2)
        lhs = R(V(-th2), V(-th1)) @ one(K1, N) @ R(V(th1), V(-th2)) @ one(K2, N)
        rhs = one(K2, N) @ R(V(th2), V(-th1)) @ one(K1, N) @ R(V(th1), V(th2))
        r, c = prop(lhs, rhs)
        # mixed: V(th1) (x) Vbar(th2): solve diagonal D = Kbar(th2)
        # Rbar_{V Vbar}... LHS: R_{Vbar V}(z1) (1(x)K_V(th1)) R_{V Vbar}(z2) (1(x)D) ; RHS: (1(x)D) R_{Vbar V}(z2)(1(x)K_V(th1)) R_{V Vbar}(z1)
        A1 = R(W(-th2), V(-th1)) @ one(K1, N) @ R(V(th1), W(-th2))
        B1 = R(W(th2), V(-th1)); C1 = one(K1, N) @ R(V(th1), W(th2))
        # equation: A1 (1(x)D) = (1(x)D) B1 C1 ; unknown diag D (N entries)
        M = []
        for k in range(N):
            E = np.zeros((N, N)); E[k, k] = 1
            M.append((A1 @ one(E, N) - one(E, N) @ B1 @ C1).ravel())
        M = np.array(M).T
        u, s, vh = np.linalg.svd(M)
        D = vh[-1].conj()
        print(n, l, 'RE_VV resid %.1e' % r, ' mixed sv', s[-1]/s[0], s[-2]/s[0], ' D ratio to KV:', D/np.diag(KV(n, om, th2, l, eta)), 'vs KV(-)', D*np.diag(KV(n,om,th2,l,eta)))
        # block constants: D*KV(th2) on S vs S^c, as function of eta
        P = D*np.diag(KV(n, om, th2, l, eta))
        if 0 < l < N:
            rat = P[l]/P[0]
            print('     ratio Sc/S', rat, ' candidates: -1', ' sinh-ish', np.sinh(eta)**0)
