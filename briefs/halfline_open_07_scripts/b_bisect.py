"""Bracket the coupling thresholds of the zeta = 1 sector (a_2, two and three solitons) with a fine rapidity grid
that includes very small and large rapidities (where the breaking first appears)."""
import numpy as np, sys, itertools
from b_grid2 import dev_at
g = np.geomspace(0.001, 10, 17)
for cfg, t, oms in (('NN', 'VV', (0.40, 0.44, 0.48, 0.5, 0.52)), ('NN', 'VW', (0.48, 0.5, 0.51, 0.52)),
                    ('UN', 'VV', (0.68, 0.7, 0.72, 0.74)), ('UN', 'VW', (0.5, 0.6, 0.65, 0.68))):
    out = []
    for om in oms:
        d = max(dev_at(cfg, list(t), [g[a], g[b]], om) for a in range(len(g)) for b in range(a))
        out.append(f'om={om}: {d:.1e}')
    print(f'{cfg} {t}: ' + '  '.join(out)); sys.stdout.flush()
g3 = np.geomspace(0.003, 8, 8)
for cfg, oms in (('NN', (0.2, 0.3, 0.35)), ('UN', (0.2, 0.3, 0.37))):
    out = []
    for om in oms:
        d = max(dev_at(cfg, ['V']*3, [g3[a], g3[b], g3[c]], om) for a, b, c in itertools.combinations(range(len(g3))[::-1], 3))
        out.append(f'om={om}: {d:.1e}')
    print(f'{cfg} VVV: ' + '  '.join(out)); sys.stdout.flush()
