import numpy as np, sys
from lib import *
exec(open('a22three.py').read().split('q = np.exp(-1j*0.77)')[0])
# consistency: transferN with N=3 equals transfer3
q = np.exp(-1j*0.77); Rpl = setup(q)
print("N=3 consistency", np.abs(transferN(Rpl, 3, [0.3, -0.8]) - transfer3(Rpl, 3, 0.3, -0.8)).max())
for xi_pi in (0.30, 0.35, 2.0, 3.0, 6.0, 12.0):
    xi = (xi_pi + 0.0013)*np.pi; qTW = 1j*np.exp(1j*np.pi**2/(3*xi)); Rpl = setup(np.conj(qTW))
    psi = (2*np.pi/(3*xi_pi + 0.0039)) % 2; psi = min(psi, 2 - psi)
    row = []
    for N in (3, 4, 5):
        (d, _), kr = scanN(Rpl, 3, N, L=6, n={3: 300, 4: 200, 5: 100}[N])
        row.append(f"N={N}: {d:.1e}")
    print(f"xi={xi_pi:5.2f}pi  |psi|/pi={psi:.3f} (1/N thresholds .5,.333,.25,.2):  " + "  ".join(row), flush=True)
