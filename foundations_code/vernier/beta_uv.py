"""First collapse threshold beta^2_UV (foundations Section 6.1) for a_{2n-1}^(2), longest roots |alpha|^2=2."""
import numpy as np, itertools
exec(open('oneloop_particles.py').read().split('print("Calibration')[0])
for n in [2, 3, 4]:
    roots, kac = a2n1_roots(n); R = np.array(roots); r = len(R)
    best = (np.inf, None)
    for size in range(2, 8):
        for S in itertools.combinations_with_replacement(range(r), size):
            v = R[list(S)]; den = (v**2).sum() - (v.sum(0)**2).sum()
            if den > 1e-12:
                b = 8*np.pi*(size - 1)/den
                if b < best[0] - 1e-12: best = (b, S)
    print(f"n={n}: beta^2_UV = {best[0]/np.pi:.6f} pi, attained by cluster {best[1]} (root indices 0..{r-1}, {r-1} = long root)")
