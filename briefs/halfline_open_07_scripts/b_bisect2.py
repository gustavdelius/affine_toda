"""NN two-soliton zeta = 1 threshold, finer bracket (fine rapidity grid incl. very small theta_2)."""
import numpy as np
from b_grid2 import dev_at
g = np.r_[np.geomspace(0.0005, 0.05, 8), np.geomspace(0.08, 12, 12)]
for om in (0.53, 0.54, 0.55):
    d = max((dev_at('NN', ['V', 'V'], [g[a], g[b]], om), g[a], g[b]) for a in range(len(g)) for b in range(a))
    print(f'NN VV om={om}: max dev {d[0]:.1e} at ({d[1]:.3g}, {d[2]:.3g})', flush=True)
