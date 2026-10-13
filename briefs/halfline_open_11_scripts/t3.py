import numpy as np
from t2 import Kbar_num
from sp_core import KV
from boundary_common import q_phys
def rr(n, om, l, eta):
    D = Kbar_num(n, om, l, eta, 0.31+0.12j, -0.17+0.05j)
    P = D*np.diag(KV(n, om, -0.17+0.05j, l, eta)); return P[l]/P[0]
for n in (2, 3, 4):
    N = n+1
    for l in range(1, N):
        om = 0.37; q = q_phys(om)
        etas = [0.3+0.4j, -0.2+1.1j, 0.7-0.3j, 0.1+0.2j, 1.3+0.9j]
        xs = [np.exp(2*e) for e in etas]; rs = [rr(n, om, l, e) for e in etas]
        # Mobius r = (A x + B)/(C x + 1): r C x + r = A x + B -> A x + B - C r x = r
        M = np.array([[x, 1, -r*x] for x, r in zip(xs[:3], rs[:3])]); A, B, C = np.linalg.solve(M, rs[:3])
        res = max(abs((A*x+B)/(C*x+1) - r) for x, r in zip(xs[3:], rs[3:]))
        # express: fixed points / try A = -q^p etc
        w = -q  # e^{-i pi om}
        print(n, l, 'res %.1e' % res, 'A', np.round(A, 6), 'B', np.round(B, 6), 'C', np.round(C, 6),
              ' log(-A)/(-i pi om)=', np.round(np.log(-A+0j)/(-1j*np.pi*om), 4), ' log(-1/C)/(-i pi om)', np.round(np.log(-1/C+0j)/(-1j*np.pi*om), 4), 'B/A', np.round(B/A,6), 'AC', np.round(A*C,6))
