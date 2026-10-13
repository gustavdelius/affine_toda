from dr import *
for cfg in ('UU', 'UN', 'NN'):
    for types in ('VV', 'VW'):
        for om in (0.37, 0.83):
            dev, comm, k = sector_dev(list(types), [0.7, 0.25], cfg[0], cfg[1], om)
            print(cfg, types, om, f'dev={dev:.1e} comm={comm:.1e} dim={k}')
