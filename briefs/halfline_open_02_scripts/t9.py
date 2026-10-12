# static solutions of d_n satisfying the C pattern BC: root finding from random starts; energies
import sys, time
import numpy as np
from scipy.optimize import fsolve
from dn_classical import *
np.set_printoptions(linewidth=180, precision=6, suppress=True)
n = int(sys.argv[1]); pat = sys.argv[2] if len(sys.argv) > 2 else '+'*(n + 1)
C = np.array([1.0 if ch == '+' else -1.0 for ch in pat])
st = Static(n, C, L=18.0)
def G(c):
    s = st.solve(c)
    if s is None:
        return np.ones(n)*1e3
    Y = s.y[:, -1]
    return st.bc_res(Y[:n], Y[n:2*n])
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
sols = []
t0 = time.time()
for trial in range(int(sys.argv[4]) if len(sys.argv) > 4 else 60):
    c0 = rng.normal(size=n)*10**rng.uniform(-1, 1.5)
    c, info, ier, msg = fsolve(G, c0, full_output=True, xtol=1e-12)
    if ier == 1 and np.abs(G(c)).max() < 1e-9:
        E = st.F(c)
        if not any(abs(E - e) < 1e-7 and np.allclose(c, cc, atol=1e-5) for e, cc in sols):
            sols.append((E, c))
            s = st.solve(c); phi0 = s.y[:n, -1]
            print(f'trial {trial}: E = {E:.10f}  c = {c}  alpha.phi(0) = {st.al @ phi0}', flush=True)
print('time', time.time() - t0)
sols.sort(key=lambda t: t[0])
np.save(f'sols_d{n}_{pat}.npy', np.array([np.r_[E, c] for E, c in sols]))
