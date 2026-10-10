# N=5, graded rule: targeted scan along clustered configurations ls = t + small spreads (soliton 1 against a cluster of 4)
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from ik_common import *
Nn = 5; Qs = charges(Nn); m1 = m_first(Nn); rng = np.random.default_rng(5)
cfgs = [t*np.ones(4) for t in np.linspace(-3, 3, 121)]
cfgs += [t*np.ones(4) + rng.normal(0, s, 4) for t in np.linspace(-2.5, 2.5, 51) for s in (0.05, 0.2) for _ in range(2)]
for psi_t in [float(v) for v in sys.argv[1].split(',')]:
    xi_pi = 2/(3*psi_t) + 0.0013; psi = psi_of_xi(xi_pi); Rp = setup(q_of_xi(xi_pi)); res = {}
    for ls in cfgs:
        T = transferN(Rp, 3, ls)
        for Q, a in ((0, 1.0), (2, 1.0), (1, 0.0)):
            idx = np.where(Qs == Q)[0]; Om = np.exp(1j*np.pi*a*m1[idx])
            dd = np.abs(np.abs(np.linalg.eigvals(Om[:, None]*T[np.ix_(idx, idx)])) - 1).max()
            if dd > res.get(Q, (0, None))[0]: res[Q] = (dd, ls)
    print(f"psi/pi={psi:+.4f}: " + "; ".join(f"Q={Q}: {res[Q][0]:.1e} at {np.round(res[Q][1], 2)}" for Q in sorted(res)), flush=True)
