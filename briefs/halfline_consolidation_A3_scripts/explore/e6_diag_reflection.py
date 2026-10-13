"""Any constant complex diagonal K_+ solves eq-classical-reflection (book RE_plus, r_toda, real lambda).
Log: ../logs/e6_diag_reflection.log (run inline on 2026-10-13 with identical code)."""
import numpy as np
from booklib import load_classical
C = load_classical(); RE_plus, r_toda = C['RE_plus'], C['r_toda']
rng = np.random.default_rng(4)
for N in (2, 3, 4):
    rf = lambda z, N=N: r_toda(N, z)
    w = 0
    for _ in range(5):
        D = np.diag(rng.normal(size=N) + 1j * rng.normal(size=N))
        for l1, l2 in [(0.7, 1.3), (1.9, 0.45), (2.2, -0.6)]:
            w = max(w, np.linalg.norm(RE_plus(lambda l: D, rf, l1, l2, N)) / np.linalg.norm(np.kron(D, D)))
    print(f'N={N}: random constant complex diagonal K_+: max relative residual of eq-classical-reflection {w:.1e}')
