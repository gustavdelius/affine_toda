"""Refine the coupling thresholds of the Theta-invariant sector (zeta = 1) for a_2: dense rapidity sampling
(including large rapidities and near-coincident ones), omega bracketing."""
import numpy as np, sys
from b_transfer import Tdr, Zshift, kron_all
rng = np.random.default_rng(21)
n = 2
def dev(cfg, types, om, nsamp=14):
    Np = len(types); s = -1 if cfg == 'MM' else 1
    Z = kron_all([Zshift(n, s)]*Np); worst = 0
    for i in range(nsamp):
        hi = [1.0, 3.0, 6.0][i % 3]
        th = np.sort(rng.uniform(0.02, hi, Np))[::-1]
        Ts = [Tdr(j, types, th, cfg[0], cfg[1], n, om) for j in range(Np)]
        G = sum(rng.normal()*T for T in Ts) + 0.37*Z
        w, P = np.linalg.eig(G); Pi = np.linalg.inv(P)
        sel = np.abs(np.diag(Pi @ Z @ P) - s**Np) < 1e-6
        for T in Ts:
            worst = max(worst, np.abs(np.abs(np.diag(Pi @ T @ P)[sel]) - 1).max())
    return worst
cases = {('NN', 'VV'): [0.62, 0.65, 0.66, 0.67, 0.68, 0.7, 0.75], ('UN', 'VV'): [0.62, 0.65, 0.66, 0.67, 0.68, 0.7, 0.75],
         ('NN', 'VVV'): [0.40, 0.43, 0.44, 0.45, 0.46, 0.48], ('NN', 'VW'): [0.5, 0.53, 0.55, 0.56, 0.57, 0.58]}
for (cfg, t), oms in cases.items():
    print(f'{cfg} {t}: ' + '  '.join(f'om={om}: {dev(cfg, list(t), om):.1e}' for om in oms)); sys.stdout.flush()
