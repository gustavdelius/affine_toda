import sys
from dr import *
for cfg, types, oms in (('UN', 'VV', (0.66, 0.68, 0.69, 0.70)), ('NN', 'VV', (0.52, 0.53, 0.535, 0.54)),
                        ('NN', 'VV', (0.6667, 2/3))):
    for om in oms:
        row = []
        for t1 in (5, 10, 15, 20, 25):
            row.append(max(sector_dev(list(types), [t1, t2], cfg[0], cfg[1], om)[0] for t2 in (1e-4, 1e-3, 1e-2, 0.1)))
        print(cfg, types, round(om, 5), ' '.join(f'{x:.1e}' for x in row)); sys.stdout.flush()
