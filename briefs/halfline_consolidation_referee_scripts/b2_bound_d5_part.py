import sys, numpy as np
import b2_dn_classical as b2
from b2_bound_fast import D
lo, hi, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
n = 5
psi, p0 = b2.vacuum(n)
g = np.linspace(lo, hi, N)
d = [D(n, psi, p0, x) for x in g]
for x, v in zip(g, d):
    print(f"{x:.5f} {v:.6e}", flush=True)
