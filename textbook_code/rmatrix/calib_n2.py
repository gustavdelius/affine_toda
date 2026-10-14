import numpy as np
from spinform import make
from dtw import rep, weights
from rsolve import intertwiner
q, n = 0.7, 2
spin_rep, Ws = make(n); Wv = weights(n)
def Rgen(A, B, wA, wB):
    Rc, null, _, _ = intertwiner(A, B, range(n+1), wA, wB); return Rc/np.abs(Rc).max(), null
funs = {"vector x spinor (heavy x light)": (lambda x: Rgen(rep(n, q, x), spin_rep(1.0), Wv, Ws), 24),
        "vector x vector (heavy x heavy)": (lambda x: Rgen(rep(n, q, x), rep(n, q, 1.0), Wv, Wv), 36),
        "spinor x spinor (light x light)": (lambda x: Rgen(spin_rep(x), spin_rep(1.0), Ws, Ws), 16)}
for name, (fun, D) in funs.items():
    print(name, "unique:", fun(1.3+0.4j)[1] == 1)
    hits = []
    for l in np.arange(-5, 5.5, 0.5):
        for ph, lab in ((1, '+'), (-1, '-'), (1j, '+i'), (-1j, '-i')):
            dims = []
            for d in (1e-5, 1e-7):
                R, _ = fun(ph*q**l*(1+d)); R = R/np.linalg.norm(R); s = np.linalg.svd(R, compute_uv=False); dims.append(int(np.sum(s > 1e-3*s[0])))
            if dims[0] == dims[1] < D: hits.append(f"{lab}q^{l:+.1f}->{dims[1]}")
    print("    special points:", "; ".join(hits))
