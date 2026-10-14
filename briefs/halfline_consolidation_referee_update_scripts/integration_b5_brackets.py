import numpy as np
from integration_b5_krein import deviation
t2s = np.geomspace(1e-4, 3, 8)
vw_grid = [(t1, t2) for t1 in (12.0, 20.0) for t2 in (1e-4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)]
for t, w, lo, hi in (('VV', 'UN', 0.695, 0.701), ('VW', 'UN', 0.664, 0.670), ('VV', 'NN', 0.533, 0.539), ('VW', 'NN', 0.497, 0.503)):
    grid = [(12.0, x) for x in t2s] if t == 'VV' else vw_grid
    d = [max(deviation(list(t), list(g), w, om) for g in grid) for om in (lo, hi)]
    print(t, w, lo, hi, [f'{x:.2e}' for x in d], flush=True)
