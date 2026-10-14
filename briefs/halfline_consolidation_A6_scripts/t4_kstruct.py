import numpy as np
from qk import K_raw, eps_star
np.set_printoptions(linewidth=220, precision=4, suppress=True)
n = 4; om = 0.37
for s in (1, -1):
    for y in (np.exp(om*(0.3+0.1j)), np.exp(om*(-0.7+0.2j))):
        K = K_raw(n, om, y, s)
        K = K/K[0, 0]
        print('s', s, 'y', y)
        print(np.abs(K))
