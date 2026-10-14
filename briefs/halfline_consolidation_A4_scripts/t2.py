import sys
from dr import *
g = np.geomspace(0.001, 10, 17)
def maxdev(cfg, types, om):
    best = (0, None)
    for a in range(len(g)):
        for b in range(a):
            dv = sector_dev(list(types), [g[a], g[b]], cfg[0], cfg[1], om)[0]
            if dv > best[0]:
                best = (dv, (round(g[a], 4), round(g[b], 4)))
    return best
for cfg, types, oms in (('UU', 'VV', (0.3, 0.7, 1.43)), ('UU', 'VW', (0.3, 0.7, 1.43)),
                        ('UN', 'VV', (0.6, 0.66, 0.68, 0.70, 0.72)), ('UN', 'VW', (0.6, 0.65, 0.68, 0.7)),
                        ('NN', 'VV', (0.5, 0.53, 0.54, 0.56, 0.67)), ('NN', 'VW', (0.48, 0.5, 0.51, 0.52))):
    for om in oms:
        dv, at = maxdev(cfg, types, om)
        print(cfg, types, om, f'{dv:.1e}', at); sys.stdout.flush()
