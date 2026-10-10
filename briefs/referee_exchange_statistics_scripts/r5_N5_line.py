# N=5 kinks, Q=0 and Q=2 at alpha=pi: fine 1D scan along the clustered line log x_{1j} = t (j=2..5), plus small transverse offsets
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from ik_common import *
Nn = 5; Qs = charges(Nn); m1 = m_first(Nn)
idx0 = np.where(Qs == 0)[0]; idx2 = np.where(Qs == 2)[0]
Om0 = np.exp(1j*np.pi*m1[idx0]); Om2 = np.exp(1j*np.pi*m1[idx2])
ts = np.concatenate([-np.linspace(0.002, 2.5, 1250), np.linspace(0.002, 2.5, 300)])
offs = [np.zeros(4), np.array([0.05, -0.05, 0.02, -0.02]), np.array([0.15, 0.0, -0.15, 0.0])]
for psi_t in [float(v) for v in sys.argv[1].split(',')]:
    xi_pi = 2/(3*psi_t) + 0.0013; psi = psi_of_xi(xi_pi); Rp = setup(q_of_xi(xi_pi)); top = (0, None)
    for o in offs:
        for t in ts:
            T = transferN(Rp, 3, t + o)
            d = max(np.abs(np.abs(np.linalg.eigvals(Om0[:, None]*T[np.ix_(idx0, idx0)])) - 1).max(),
                    np.abs(np.abs(np.linalg.eigvals(Om2[:, None]*T[np.ix_(idx2, idx2)])) - 1).max())
            if d > top[0]: top = (d, t, o)
    print(f"psi/pi={psi:+.4f}: max dev {top[0]:.2e} at t={top[1]:.3f}, offset {top[2]}", flush=True)
