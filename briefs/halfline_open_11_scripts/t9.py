import numpy as np
from sp_core import *
for n in (1, 2, 3, 4):
    N = n+1
    for l in range(0, N+1):
        res = []
        for ch in (0, 1):
            ks = []
            for om in (0.37, 0.53):
                T = N*om/2; eta = 0.3+0.4j
                th = 0.21+0.13j
                c = Phis(n, om, l, eta, th)[ch]
                g = np.sinh(2*T*th - 1j*np.pi*T)/np.sinh(2*T*th + 1j*np.pi*T)
                for sgn in (1, -1):
                    cc = c/g**sgn
                    k = np.log(-cc+0j)/(-1j*np.pi*om)
                    # check other theta
                    th2 = -0.4+0.05j; c2 = Phis(n, om, l, eta, th2)[ch]
                    g2 = np.sinh(2*T*th2 - 1j*np.pi*T)/np.sinh(2*T*th2 + 1j*np.pi*T)
                    if abs(c2/g2**sgn/cc - 1) < 1e-9:
                        ks.append((sgn, k))
            # k defined mod 2/om; reconcile with two omegas
            res.append(ks)
        print(n, l, [[(s, np.round(k.real, 4)) for s, k in ks] for ks in res])
