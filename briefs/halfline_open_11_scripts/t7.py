import numpy as np
from sp_core import *
from boundary_common import q_phys
om = 0.37; q = q_phys(om)
for n in (1, 2, 3):
    N = n+1; T = N*om/2
    for l in range(0, N+1):
        eta = 0.3+0.4j
        V = lambda th: Vy(n, q, th, om); W = lambda th: Vy(n, q, th, om, True)
        # mixed RE check with KVbar formula
        th1, th2 = 0.31+0.12j, -0.17+0.05j
        K1 = KV(n, om, th1, l, eta); D = KVbar(n, om, th2, l, eta)
        lhs = R(W(-th2), V(-th1)) @ one(K1, N) @ R(V(th1), W(-th2)) @ one(D, N)
        rhs = one(D, N) @ R(W(th2), V(-th1)) @ one(K1, N) @ R(V(th1), W(th2))
        rm, _ = prop(lhs, rhs)
        out = []
        for th in (0.23+0.11j, -0.4+0.07j):
            a = 0.5j*np.pi
            kA = lambda t: sum(KVbar(n, om, a - t, l, eta)[i, i]*np.kron(np.eye(N)[i], np.eye(N)[i]) for i in range(N))  # V(-t)(x)Vbar(t)
            kB = lambda t: sum(KV(n, om, a - t, l, eta)[i, i]*np.kron(np.eye(N)[i], np.eye(N)[i]) for i in range(N))   # Vbar(-t)(x)V(t)
            v1 = R(V(th), W(-th)) @ kA(-th); r1, c1 = prop(kB(th), v1)
            v2 = R(W(th), V(-th)) @ kB(-th); r2, c2 = prop(kA(th), v2)
            out.append((r1, r2, c1, c2))
        print(n, l, 'mixedRE %.1e' % rm, ' crossing resid', ['%.1e %.1e' % (o[0], o[1]) for o in out], 'Phi', [np.round(o[2], 5) for o in out], [np.round(o[3], 5) for o in out])
