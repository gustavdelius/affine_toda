import numpy as np
from t2 import Kbar_num
from sp_core import KV
from boundary_common import q_phys
om = 0.37; q = q_phys(om)
for n in (2, 3, 4):
    N = n+1; T = N*om/2
    for l in range(1, N):
        for eta in (0.3+0.4j, -0.5+0.9j):
            ths = [-0.17+0.05j, 0.25-0.1j, 0.5+0.2j, 0.1+0.33j]
            gs = []
            for th in ths:
                D = Kbar_num(n, om, l, eta, 0.31+0.12j, th)
                P = D*np.diag(KV(n, om, th, l, eta)); r = P[l]/P[0]
                gs.append(r*np.sinh(eta+T*th)/np.sinh(eta-T*th))
            best = None
            for c in range(-2*N, 2*N+1):
                h = [g*np.exp(-c*om*th) for g, th in zip(gs, ths)]
                X = [np.exp(2*T*th) for th in ths]
                # h = k (xb X - 1)/(xb - X): unknowns k, xb ; h(xb - X) = k xb X - k ; two eq from first two points: linear in (xb, k, k xb)?
                # write h xb - h X = k xb X - k -> unknowns u=xb, v=k, w=k*xb: h u - X w + v = h X
                M = np.array([[h[i], 1, -X[i]] for i in range(3)]); u, v, w = np.linalg.solve(M, [h[i]*X[i] for i in range(3)])
                res = abs(h[3]*u + v - X[3]*w - h[3]*X[3]) + abs(w - u*v)
                if best is None or res < best[0]: best = (res, c, u, v)
            res, c, xb, k = best
            etab = np.log(xb)/2
            print(n, l, eta, 'c=%d res %.1e' % (c, res), 'xibar', np.round(xb, 6), 'etabar-eta', np.round(etab - eta, 6), '/(i pi om)', np.round((etab-eta)/(1j*np.pi*om), 6), 'k', np.round(k, 6))
