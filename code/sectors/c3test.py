import numpy as np, time
from c3fuse import *
t0 = time.time()
th = 0.23 + 0.13j
for a, b in [(3, 1), (1, 1), (3, 2), (1, 2), (2, 2)]:
    Yab, r1 = apply(a, b, th, np.eye(dims(a)*dims(b), dtype=complex), S33)
    Yba, r2 = apply(b, a, -th, np.eye(dims(a)*dims(b), dtype=complex), S33)
    print(f"S_{a}{b}: closure {r1:.0e}; unitarity {np.abs(Yba @ Yab - np.eye(len(Yab))).max():.1e}  [{time.time()-t0:.0f}s]")
# scalar factor modulus on real axis for symmetric offsets
for a, b in [(3, 1), (2, 2)]:
    for thr in (0.07, 0.4):
        print(f"|prod F33| for S_{a}{b} at theta={thr}: {abs(np.exp(scalar_logF(a, b, thr))):.12f}")
