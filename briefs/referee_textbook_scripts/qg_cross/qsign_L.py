"""Is R(x;-q) = L R(x;q) L^-1 with L a diagonal sign matrix on spinor(x)spinor (c_n^(1) solitons, Part VI conventions)?"""
import numpy as np
from rsolve import intertwiner
from qsign_cn import make
omega = 2.37
for n in (1, 2, 3):
    D = 2**n; Rs = {}
    for s in (1, -1):
        q = s*np.exp(-1j*np.pi*omega); rep, W = make(n, q)
        Rc = intertwiner(rep(1.37+0.41j), rep(1.0), range(n+1), W, W)[0]; Rs[s] = Rc/Rc[0, 0]
    A, B = Rs[1], Rs[-1]; nz = np.abs(A) > 1e-10; idx = np.argwhere(nz)
    l = {}
    for start in range(D*D):                      # propagate l_out/l_in = B/A over each connected block
        if start in l: continue
        l[start] = 1; changed = True
        while changed:
            changed = False
            for r, c in idx:
                ratio = np.round((B[r, c]/A[r, c]).real)
                if c in l and r not in l: l[r] = ratio*l[c]; changed = True
                elif r in l and c not in l: l[c] = l[r]/ratio; changed = True
    bad = sum(abs(l[r]/l[c] - B[r, c]/A[r, c]) > 1e-8 for r, c in idx)
    print(f"n={n}: same sparsity {np.array_equal(nz, np.abs(B) > 1e-10)}; entry ratios R(-q)/R(q) in {sorted(set(np.round((B[nz]/A[nz]).real).astype(int)))}; "
          f"diagonal sign L with R(-q) = L R(q) L^-1 exists: {bad == 0}")
