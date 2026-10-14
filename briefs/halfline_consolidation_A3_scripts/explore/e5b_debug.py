"""a_2: are omega o tau-hat or tau-hat (diagram reflections, E_j -> E_tau(j) for all affine nodes) implementable by conjugation
in the defining representation?  Log: ../logs/e5b_debug.log (run inline on 2026-10-13 with identical code)."""
import numpy as np
from booklib import load_classical
C = load_classical(); algebra = C['algebra']
h = 3
roots, Es, Hs, nj, sub = algebra('a2(1)')
for c in range(h):
    M = np.zeros((h, h))
    for i in range(h):
        M[(c + 1 - i) % h, i] = -1
    tau = [next(k for k in range(h) if np.allclose(M @ roots[j], roots[k])) for j in range(h)]
    N = h; Mat = []; Mat2 = []
    for p in range(N):
        for q in range(N):
            K = np.zeros((N, N)); K[p, q] = 1
            Mat.append(np.concatenate([np.concatenate([(K @ Es[j] + Es[tau[j]].T @ K).ravel(), (K @ Es[j].T + Es[tau[j]] @ K).ravel()]) for j in range(h)]))
            Mat2.append(np.concatenate([np.concatenate([(K @ Es[j] - Es[tau[j]] @ K).ravel(), (K @ Es[j].T - Es[tau[j]].T @ K).ravel()]) for j in range(h)]))
    sv = np.linalg.svd(np.array(Mat).T, compute_uv=False); sv2 = np.linalg.svd(np.array(Mat2).T, compute_uv=False)
    print(f'c={c} tau={tau}: generators-only intertwiner for theta = omega o tau: nullity {(sv < 1e-12*sv[0]).sum()}, sv_min/max {sv[-1]/sv[0]:.2e}')
    print(f'      tau-hat alone (E_j -> E_tau(j)) implementable: nullity {(sv2 < 1e-12*sv2[0]).sum()}')
