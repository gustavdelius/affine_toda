# d_4, d_5: symmetric reduction for C = +1, then full-space check and the linearized reflection
import numpy as np
from scipy.optimize import minimize
from dn_classical import *
np.set_printoptions(linewidth=180, precision=6, suppress=True)
for n in (4, 5):
    h = 2*n - 2
    C = np.ones(n + 1)
    st = Static(n, C)
    best = None
    rng = np.random.default_rng(1)
    for trial in range(40):
        c0 = rng.normal(size=n)*rng.choice([0.3, 1, 3])
        r = minimize(st.F, c0, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-13, 'maxiter': 20000, 'maxfev': 40000})
        if np.isfinite(r.fun) and (best is None or r.fun < best.fun - 1e-9):
            best = r
            print(n, trial, 'F =', r.fun, 'c =', r.x)
    s = st.solve(best.x, dense=True)
    Y = s.y[:, -1]
    print('  BC residual', np.abs(st.bc_res(Y[:n], Y[n:2*n])).max(), ' alpha.phi(0) =', st.al @ Y[:n])
    np.save(f'vac_d{n}_plus.npy', best.x)
