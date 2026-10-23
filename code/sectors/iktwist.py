import numpy as np
from lib import *
exec(open('a22three.py').read().split('q = np.exp(-1j*0.77)')[0])
m = np.array([1, 0, -1])
for xi_pi in (0.9, 1.1, 0.6):
    xi = xi_pi*np.pi; q = np.conj(1j*np.exp(1j*np.pi**2/(3*xi))); Rpl = setup(q)
    best = (9, None)
    for al in np.linspace(0, 2*np.pi, 73)[:-1]:
        Q = np.diag([np.exp(1j*al*m[i]*m[j]) for i in range(3) for j in range(3)])
        d = max(np.abs(np.abs(np.linalg.eigvals(Rpl(np.exp(l)) @ Q)) - 1).max() for l in np.linspace(0.05, 12, 120))
        if d < best[0]: best = (d, al)
    print(f"a_2^(2) xi={xi_pi}pi: best charge-bicharacter twist exp(i a m_i m_j): a={best[1]:.3f}, 2-body maxdev {best[0]:.2e}")
