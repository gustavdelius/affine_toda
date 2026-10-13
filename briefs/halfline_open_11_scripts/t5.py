import numpy as np
from t2 import Kbar_num
from sp_core import KV
from boundary_common import q_phys
for n in (2, 3):
    N = n+1; T = N*0.37/2
    for l in (1, 2):
        om = 0.37; eta = 0.3+0.4j
        for th2 in (-0.17+0.05j, 0.25-0.1j, 0.5+0.2j):
            D = Kbar_num(n, om, l, eta, 0.31+0.12j, th2)
            P = D*np.diag(KV(n, om, th2, l, eta)); r = P[l]/P[0]
            # candidate: r = sinh(eta+T th + x) sinh(eta - T th + x') ... print ratio to sinh(eta+T th)/sinh(eta-T th)
            g = np.sinh(eta + T*th2)/np.sinh(eta - T*th2)
            print(n, l, th2, np.round(r, 6), ' r/g', np.round(r/g, 6), ' r*g', np.round(r*g, 6))
