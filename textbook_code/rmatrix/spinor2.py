import numpy as np
from spinor import spin_rep, spin_weights, Rgen, q, n
from dtw import rep, weights
Ws, Wv = spin_weights(), weights(n)
def scan(fun, D, phases):
    found = []
    for ph, lab in phases:
        for l in np.arange(-7, 7.25, 0.25):
            dims = []
            for d in (1e-5, 1e-7):
                R, _ = fun(ph*q**l*(1+d)); R = R/np.linalg.norm(R); s = np.linalg.svd(R, compute_uv=False)
                dims.append(int(np.sum(s > 1e-3*s[0])))
            if dims[0] == dims[1] < D: found.append((lab, l, dims[1]))
    return found
phases = [(1, "+"), (-1, "-"), (1j, "+i"), (-1j, "-i")]
for name, fun in [("R_33 spinor x spinor", lambda x: Rgen(spin_rep(x), spin_rep(1.0), Ws, Ws)),
                  ("R_13 vector x spinor", lambda x: Rgen(rep(n, q, x), spin_rep(1.0), Wv, Ws)),
                  ("R_11 vector x vector", lambda x: Rgen(rep(n, q, x), rep(n, q, 1.0), Wv, Wv))]:
    res = scan(fun, 64, phases)
    print(name + ":", "; ".join(f"x*={lab}q^{l:+.2f} -> {d}" for lab, l, d in res))
