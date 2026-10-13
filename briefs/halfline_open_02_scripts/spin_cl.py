"""d_5 uniform vacuum: linearized reflection in the spinor channel vs B -> 0 limit of the CAF+am3 spinor amplitude."""
import numpy as np
import mpmath as mp
from sym_vac2 import Prof
from dn_classical import masses, linear_K
from dn_particle import S_dn
n = 5
m, U = masses(n)
P = Prof(n)
for ths in (0.5, 0.35):
    Om = np.sqrt(2)*np.cosh(ths)
    K, op = linear_K(n, np.ones(n + 1), P, Om)
    print('open channels (mass):', [round(m[a], 4) for a in op])
    print('K =\n', np.round(K, 8))
    sub = [i for i, a in enumerate(op) if abs(m[a] - np.sqrt(2)) < 1e-6]
    Ks = K[np.ix_(sub, sub)]
    print('spinor block eig:', np.round(np.linalg.eigvals(Ks), 8))
    kcl = complex(mp.fprod([S_dn(5, j, n)(mp.mpc(ths), 1) for j in range(1, n + 1)]))
    kss = complex(S_dn(4, 5, n)(mp.mpc(ths), 1)); kss2 = complex(S_dn(4, 4, n)(mp.mpc(ths), 1))
    print('CAF+am3 K_s(B=0) = prod_j S_sj(th, 1) =', np.round(kcl, 8))
