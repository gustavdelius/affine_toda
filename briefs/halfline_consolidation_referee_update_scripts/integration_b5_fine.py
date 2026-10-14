import numpy as np
from integration_b5_krein import deviation
t1s = (6.0, 8.0, 12.0, 16.0, 20.0, 30.0)
t2s = list(np.arange(0.1, 6.01, 0.1)) + [1e-3, 1e-2]
for walls, oms, types in (('UN', (0.6640, 0.6694), 'VW'), ('NN', (0.4970, 0.5030), 'VW'),
                          ('UN', (0.6950, 0.7010), 'VV'), ('NN', (0.5330, 0.5390), 'VV')):
    for om in oms:
        devs = {}
        for a in t1s:
            for b in t2s:
                devs[(a, b)] = max(deviation(list(types), [a, b], walls, om), deviation(list(types), [b, a], walls, om))
        k = max(devs, key=devs.get)
        print(f'{walls} {types} omega = {om}: max deviation over {len(devs)} rapidity pairs = {devs[k]:.2e} at {k}', flush=True)
