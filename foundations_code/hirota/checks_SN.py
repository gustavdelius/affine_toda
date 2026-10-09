"""(S) centre characters separate mass-degenerate species; (N) no resonances K m_a = m_c, 2 <= K <= max n_j."""
import numpy as np
from lie import Algebra
for nm in ['a3','a5','a7','d4','d5','d6','d7','d8','d9','d10','d13','e6','e7','e8']:
    A = Algebra(nm[0], int(nm[1:]))
    sp = sorted(A.mass)
    groups = {}
    for a in sp:
        groups.setdefault(round(A.mass[a], 8), []).append(a)
    sep = True
    for m, g in groups.items():
        if len(g) > 1:
            chars = [tuple(np.round([A.chi(z, a) for z in A.centre], 6)) for a in g]
            if len(set(chars)) < len(g):
                sep = False
    Kmax = int(max(A.n))
    gap = min(abs(K * A.mass[a] - A.mass[c]) for K in range(2, Kmax + 1) for a in sp for c in sp) if Kmax >= 2 else np.inf
    degs = [g for g in groups.values() if len(g) > 1]
    print(f"{nm}: degenerate sets {degs}; (S) centre separates: {sep}; (N) min_K>=2 |K m_a - m_c| = {gap:.4f}")
