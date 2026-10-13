from dr import *
g = np.geomspace(0.001, 12, 11)
for cfg in ('UU', 'UN', 'NN'):
    for types in ('VV', 'VW'):
        out = []
        for om in (0.1, 0.2, 0.37):
            out.append(max(sector_dev(list(types), [g[a], g[b]], cfg[0], cfg[1], om)[0] for a in range(len(g)) for b in range(a)))
        print(cfg, types, ' '.join(f'{x:.0e}' for x in out), flush=True)
