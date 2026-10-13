"""Two solitons of a_2 on an interval: deviation max_j ||Lambda|-1| in the Theta-invariant sector (zeta = 1) on a
log-spaced grid of real rapidities theta_1 > theta_2 in [0.02, 8]; prints the max and where it occurs."""
import numpy as np, sys, time
from b_transfer import Tdr, Zshift, kron_all
rng = np.random.default_rng(2); n = 2
g = np.geomspace(0.02, 8, 13)
def dev_at(cfg, types, th, om):
    s = -1 if cfg == 'MM' else 1; Np = len(types)
    Z = kron_all([Zshift(n, s)]*Np)
    Ts = [Tdr(j, types, th, cfg[0], cfg[1], n, om) for j in range(Np)]
    G = sum(rng.normal()*T for T in Ts) + 0.37*Z
    w, P = np.linalg.eig(G); Pi = np.linalg.inv(P)
    sel = np.abs(np.diag(Pi @ Z @ P) - s**Np) < 1e-6
    return max(np.abs(np.abs(np.diag(Pi @ T @ P)[sel]) - 1).max() for T in Ts)
if __name__ == "__main__":
  for cfg, t in (('NN', 'VV'), ('UN', 'VV'), ('UU', 'VV'), ('UU', 'VW'), ('NN', 'VW')):
      for om in (0.25, 0.37, 0.55, 1.43):
          t0 = time.time(); best = (0, None); nbad = 0; tot = 0
          for a in range(len(g)):
              for b in range(a):
                  d = dev_at(cfg, list(t), [g[a], g[b]], om); tot += 1
                  nbad += d > 1e-8
                  if d > best[0]:
                      best = (d, (round(g[a], 3), round(g[b], 3)))
          print(f'{cfg} {t} om={om}: max dev {best[0]:.1e} at {best[1]}; broken at {nbad}/{tot} grid points '
                f'[{time.time()-t0:.0f}s]'); sys.stdout.flush()
