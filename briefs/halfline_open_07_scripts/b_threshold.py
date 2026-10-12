"""Problem 7(b): coupling dependence of the Theta-invariant Z_3 sector (zeta = 1; zeta = -1 for s = -1, N odd).
Prints max_j max ||Lambda|-1| over random real rapidities in that sector, on a grid of omega."""
import numpy as np, sys
from b_transfer import Tdr, Zshift, kron_all
rng = np.random.default_rng(9)
n = 2
def dev_sector(cfg, types, om, nsamp=6):
    Np = len(types); s = -1 if cfg == 'MM' else 1
    target = s**Np
    Z = kron_all([Zshift(n, s)]*Np); worst = 0
    for _ in range(nsamp):
        th = np.sort(rng.uniform(0.05, 3.0, Np))[::-1]
        Ts = [Tdr(j, types, th, cfg[0], cfg[1], n, om) for j in range(Np)]
        G = sum(rng.normal()*T for T in Ts) + 0.37*Z
        w, P = np.linalg.eig(G); Pi = np.linalg.inv(P)
        zz = np.diag(Pi @ Z @ P); sel = np.abs(zz - target) < 1e-6
        for T in Ts:
            lam = np.diag(Pi @ T @ P)[sel]
            worst = max(worst, np.abs(np.abs(lam) - 1).max())
    return worst
grid = np.round(np.arange(0.1, 3.01, 0.1), 2)
cases = [('UU', 'VVV'), ('UU', 'VVW'), ('UN', 'VV'), ('NN', 'VV'), ('NN', 'VVV'), ('MM', 'VV'), ('UN', 'VW'), ('NN', 'VW')]
print('omega ' + ' '.join(f'{c}:{t:4s}' for c, t in cases))
for om in grid:
    vals = [dev_sector(c, list(t), om) for c, t in cases]
    print(f'{om:4.1f}  ' + ' '.join(f'{v:9.1e}' for v in vals)); sys.stdout.flush()
