import numpy as np
from t2 import Kbar_num
from sp_core import KV
for n in (2, 3):
    N = n+1
    for l in range(0, N):
        om = 0.37; eta = 0.3+0.4j
        out = []
        for th1, th2 in [(0.31+0.12j, -0.17+0.05j), (0.6-0.1j, -0.17+0.05j), (-0.4+0.3j, -0.17+0.05j)]:
            D = Kbar_num(n, om, l, eta, th1, th2)
            P = D*np.diag(KV(n, om, th2, l, eta)); out.append(P/P[0])
        print(n, l, np.round(np.array(out), 5).tolist())
