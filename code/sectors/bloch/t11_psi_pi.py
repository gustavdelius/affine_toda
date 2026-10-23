"""psi = pi exactly (xi = 2pi/3, q^2 = 1): is the particle-labelled kink R-matrix trivial? Compare nearby couplings."""
from ik_common import *
for d in (0.0, 1e-4, 1e-3, 1e-2):
    q = q_of_xi(2/3 + d); Rp = setup(q)
    off = max(np.abs(Rp(np.exp(l)) - np.diag(np.diag(Rp(np.exp(l))))).max() for l in (0.3, 1.5, 4.0))
    dev = max(np.abs(np.abs(np.linalg.eigvals(Rp(np.exp(l)))) - 1).max() for l in np.linspace(0.01, 12, 400))
    T3 = max(np.abs(np.abs(np.linalg.eigvals(transferN(Rp, 3, ls))) - 1).max() for ls in np.random.default_rng(0).uniform(-5, 5, (40, 2)))
    print(f"xi = (2/3 + {d:g}) pi: max offdiag of R {off:.1e}; 2-kink maxdev {dev:.1e}; 3-kink maxdev {T3:.1e}")
