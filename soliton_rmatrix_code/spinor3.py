import numpy as np, sys
from spinor import spin_rep, spin_weights, Rgen, q, n
from dtw import rep, weights
Ws, Wv = spin_weights(), weights(n)
funs = {"R_33": lambda x: Rgen(spin_rep(x), spin_rep(1.0), Ws, Ws),
        "R_13": lambda x: Rgen(rep(n, q, x), spin_rep(1.0), Wv, Ws),
        "R_11": lambda x: Rgen(rep(n, q, x), rep(n, q, 1.0), Wv, Wv)}
name = sys.argv[1]; fun = funs[name]
out = []
for ph, lab in [(1j, "i"), (-1j, "-i")]:
    for l in np.arange(-6, 6.5, 0.5):
        R, _ = fun(ph*q**l*(1+1e-6)); R = R/np.linalg.norm(R); s = np.linalg.svd(R, compute_uv=False)
        d = int(np.sum(s > 1e-3*s[0]))
        if d < 64:
            R2, _ = fun(ph*q**l*(1+1e-8)); R2 = R2/np.linalg.norm(R2); s2 = np.linalg.svd(R2, compute_uv=False)
            d2 = int(np.sum(s2 > 1e-3*s2[0]))
            if d2 == d: out.append(f"x*={lab}q^{l:+.1f} -> {d}")
print(name + " (imaginary axis):", "; ".join(out) if out else "none", flush=True)
