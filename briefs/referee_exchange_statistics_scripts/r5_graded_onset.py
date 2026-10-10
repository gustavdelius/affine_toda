# Graded rule onset: N kinks, even Q at alpha=pi, odd Q at alpha=0; clustered line log x_{1j} = t incl. tiny |t|, plus random sets
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from ik_common import *
Nn = int(sys.argv[1]); Qs = charges(Nn); m1 = m_first(Nn); rng = np.random.default_rng(9)
secs = [(Q, 1.0 if Q % 2 == 0 else 0.0) for Q in range(Nn + 1)]
ts = np.concatenate([-np.logspace(-5, 0.4, 500), np.logspace(-5, 0.4, 200)])
rand = [rng.uniform(-3, 3, Nn-1) for _ in range(300 if Nn < 5 else 100)]
for psi_t in [float(v) for v in sys.argv[2].split(',')]:
    xi_pi = 2/(3*psi_t) + 0.0013; psi = psi_of_xi(xi_pi); Rp = setup(q_of_xi(xi_pi)); top = (0, None, None)
    for ls in [t*np.ones(Nn-1) for t in ts] + rand:
        T = transferN(Rp, 3, ls)
        for Q, a in secs:
            idx = np.where(Qs == Q)[0]; Om = np.exp(1j*np.pi*a*m1[idx])
            d = np.abs(np.abs(np.linalg.eigvals(Om[:, None]*T[np.ix_(idx, idx)])) - 1).max()
            if d > top[0]: top = (d, Q, ls)
    print(f"N={Nn} psi/pi={psi:+.4f} (2/N={2/Nn:.3f}): max dev {top[0]:.2e} in Q={top[1]} at {np.round(top[2], 5)}", flush=True)
