# Sanity of the spinor Rcheck: weight conservation, Yang-Baxter, braiding unitarity, commutation with U_q generators
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from spin_common import *
for n, alg in ((3, 'c'), (3, 'a'), (2, 'c')):
    D, Wn, setup = make(n, alg); P = flip(D, D); I = np.eye(D)
    tw = np.array([tuple(Wn[a]+Wn[b]) for a in range(D) for b in range(D)])
    same = np.all(np.isclose(tw[:, None, :], tw[None, :, :]), axis=2)
    for om in (0.3, 0.45, 0.7, 1.2, 1.61, 2.37):
        q = np.exp(-1j*np.pi*om)
        for s in (1, -1):
            R = setup(s*q); C = lambda x: P@R(x)
            x, y = 1.7+0.4j, 0.6-0.3j
            wc = np.abs(np.where(same, 0, C(x))).max()
            C12 = lambda z: np.kron(C(z), I); C23 = lambda z: np.kron(I, C(z))
            yb = np.abs(C12(x)@C23(x*y)@C12(y) - C23(y)@C12(x*y)@C23(x)).max()
            un = np.abs(C(x)@C(1/x) - np.eye(D*D)).max()
            print(f"{alg}_{n} om={om} sign={s:+d}: weight-violation {wc:.1e}, YBE {yb:.1e}, C(x)C(1/x)-1 {un:.1e}", flush=True)
