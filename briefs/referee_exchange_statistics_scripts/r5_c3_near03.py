# c_3 spinors near omega = 0.3 (where the code's R-matrix is singular): R sanity and the Bloch-torus scan of t6 for Q=0
import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from spin_common import *
n, alg, steps = 3, sys.argv[1], int(sys.argv[2]); oms = [float(v) for v in sys.argv[3].split(',')]
D, Wn, setup = make(n, alg); P = flip(D, D); I = np.eye(D)
tw = np.array([tuple(Wn[a]+Wn[b]) for a in range(D) for b in range(D)])
same = np.all(np.isclose(tw[:, None, :], tw[None, :, :]), axis=2)
tot, first = weights_multi(Wn, 2)
lx = np.concatenate([np.linspace(0.01, 5, 120), np.linspace(5, 25, 30)])
grid = [np.array(a)*2*np.pi/steps for a in itertools.product(range(steps), repeat=n)]
for om in oms:
    q = np.exp(-1j*np.pi*om); Rp = setup(q); C = lambda x: P@Rp(x)
    x, y = 1.7+0.4j, 0.6-0.3j
    yb = np.abs(np.kron(C(x), I)@np.kron(I, C(x*y))@np.kron(C(y), I) - np.kron(I, C(y))@np.kron(C(x*y), I)@np.kron(I, C(x))).max()
    wc = np.abs(np.where(same, 0, C(x))).max()
    Rs = [Rp(np.exp(l)) for l in lx]
    idx = np.where(np.all(np.isclose(tot, 0), axis=1))[0]
    blocks = [R[np.ix_(idx, idx)] for R in Rs]; nunb = 0; best = 9e9
    for a in grid:
        Om = np.exp(1j*(first[idx] @ a)); d = 0
        for B in blocks:
            d = max(d, np.abs(np.abs(np.linalg.eigvals(Om[:, None]*B)) - 1).max())
            if d > 1e-6: break
        if d < 1e-6: nunb += 1
        best = min(best, d)
    print(f"{alg}_3 omega={om}: YBE {yb:.1e}, weight-violation {wc:.1e}; Q=0 unbroken grid points {nunb}/{len(grid)}", flush=True)
