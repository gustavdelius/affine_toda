import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from spin_common import *
n, alg = 3, 'c'
D, Wn, setup = make(n, alg); P = flip(D, D)
Pi = np.diag([(-1.0)**sum(1 for w in Wn[i] if w < 0) for i in range(D)]); K = np.kron(Pi, np.eye(D))
for om in (0.3, 1.61):
    q = np.exp(-1j*np.pi*om); Rp, Rm = setup(q), setup(-q)
    for x in (0.37, 1.0, 1.7+0.4j, 1.9, 2.9, 7.3, 1.2+0.5j, 0.3-1.1j):
        Cp, Cm = P@Rp(x), P@Rm(x)
        print(f"om={om} x={x}: |Cm - K Cp K| = {np.abs(Cm - K@Cp@K).max():.1e}, |Cm|max {np.abs(Cm).max():.1e}, |Cp|max {np.abs(Cp).max():.1e}")
