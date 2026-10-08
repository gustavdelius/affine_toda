import numpy as np, time
from check22 import R22
t = time.time()
for l in (-6, -4, -2):
    dims = []
    for d in (1e-5, 1e-7):
        R, _ = R22(-(0.7**l)*(1+d)); R = R/np.linalg.norm(R)
        s = np.linalg.svd(R, compute_uv=False); dims.append(int(np.sum(s > 1e-3*s[0])))
    print(f"R_22 at x = -q^{l:+d}: image dimension {dims}  [{time.time()-t:.0f}s]", flush=True)
