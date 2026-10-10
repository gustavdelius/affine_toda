# N=5 kinks, graded rule (even Q at alpha=pi, odd Q at alpha=0), denser rapidity sampling just above 2 pi/5
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from ik_common import *
Nn = 5; Qs = charges(Nn); m1 = m_first(Nn); rng = np.random.default_rng(11)
lss = [rng.uniform(-6, 6, 4) for _ in range(150)] + [rng.uniform(-1.5, 1.5, 4) for _ in range(150)] + [rng.uniform(-0.3, 0.3, 4) for _ in range(100)] \
      + [np.cumsum(rng.uniform(0, 2, 4)) for _ in range(100)]
for psi_t in [float(v) for v in sys.argv[1].split(',')]:
    xi_pi = 2/(3*psi_t) + 0.0013; psi = psi_of_xi(xi_pi); Rp = setup(q_of_xi(xi_pi))
    res = {}
    for ls in lss:
        T = transferN(Rp, 3, ls)
        for Q, a in ((0, 1.0), (2, 1.0), (4, 1.0), (1, 0.0), (3, 0.0)):
            idx = np.where(Qs == Q)[0]; Om = np.exp(1j*np.pi*a*m1[idx])
            d = np.abs(np.abs(np.linalg.eigvals(Om[:, None]*T[np.ix_(idx, idx)])) - 1).max()
            if d > res.get(Q, (0, None))[0]: res[Q] = (d, ls)
    print(f"psi/pi={psi:+.4f} (2pi/5 = 0.4): " + "; ".join(f"Q={Q}{'(a=pi)' if Q%2==0 else '(a=0)'}: {res[Q][0]:.1e}" for Q in sorted(res)), flush=True)
    print(f"     worst Q=0 at rapidities {np.round(res[0][1], 3)}", flush=True)
