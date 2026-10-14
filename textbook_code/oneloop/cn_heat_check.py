"""Heat traces of the c_n^(1) solitons, n <= 4: prp-cn-threshold-blocks part 3 (Chapter 18).

Run from this directory:  python3 cn_heat_check.py     (about 15 min)

For every soliton a < n of c_2^(1), c_3^(1) and c_4^(1), realized as the sigma-invariant sector of a_{2n-1}^(1),
compare the Fourier-grid heat trace Tr(e^{-tA} - e^{-tA0}) with prp-net-count part 3 at t = 0.3, 0.6, 1, 1.5.
The double zero of X_n at lambda_0 = m_{n-a}^2 (a threshold) carries weight 2. With weight 1 instead, the
formula is off by e^{-t lambda_0}, which the grid rules out.
"""
import sys, os, warnings
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'hirota'))
from heat_trace_checks import levels_of, H_PT
from sectors import fold_sector

warnings.filterwarnings('ignore')
N, LP, TOL = 768, 48, 5e-5
ts = (0.3, 0.6, 1.0, 1.5)

for n in (2, 3, 4):
    for a in range(1, n):
        S = fold_sector(f'c{n}(1)<a{2*n-1}', a)
        lev = levels_of(S)
        double = [L['mu']**2*(1 - w*w) for L in lev.values() for w, o in L['ord'].items() if w > 0 and o == 2]
        lam0 = double[0]
        thresholds = sorted(L['mu']**2 for L in lev.values())
        at_threshold = any(abs(lam0 - m2) < 1e-9 for m2 in thresholds)
        Hm, H0, trV, kmax, x = S.grid(N, LP)
        ev = np.linalg.eigvals(Hm); ev0 = np.linalg.eigvalsh(H0.real)
        grid = [np.sum(np.exp(-t*ev)).real - np.sum(np.exp(-t*ev0)) for t in ts]
        err2 = max(abs(g - H_PT(lev, t)) for g, t in zip(grid, ts))
        err1 = min(abs(g - (H_PT(lev, t) - np.exp(-t*lam0))) for g, t in zip(grid, ts))
        ok = err2 < TOL and err1 > 100*TOL and len(double) == 1 and at_threshold
        print(f"{'PASS' if ok else 'FAIL'} c{n}(1) soliton {a}: double zero at lambda_0 = {lam0:.6f} "
              f"(a threshold: {at_threshold}); max |grid - formula| = {err2:.1e} with weight 2, "
              f"min = {err1:.1e} with weight 1", flush=True)
