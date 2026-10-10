"""The sign of q as operator identities (book: sec-crossing-sign, sec-sector-status).

1. Rcheck(x; -q) = (Pi x 1) Rcheck(x; q) (Pi x 1), Pi = (-1)^g, g = parity of the number of minus signs in the weight.
2. For N solitons: T_1(-q) = S (-1)^G Pi_1^N T_1(q) S^{-1}, with S the product of the parities of the even-numbered
   solitons and G the total grading. S commutes with every twist Omega_alpha, so -q maps the sector (Q, alpha) to
   (Q, alpha + pi N gamma), gamma = (1,...,1), with the spectrum multiplied by a phase that is constant on the sector.
"""
from spin_common import *
for n, alg in ((2, 'c'), (3, 'c'), (2, 'a'), (3, 'a')):
    D, Wn, setup = make(n, alg)
    g = np.array([bin(i).count('1') % 2 for i in range(D)]); Pi = (-1.0)**g
    P = flip(D, D)
    for om in (2.37, 1.61, 0.3013):
        q = np.exp(-1j*np.pi*om); Rp, Rm = setup(q), setup(-q)      # particle-labelled R = P Rcheck
        A = np.kron(np.diag(Pi), np.eye(D))
        r1 = max(np.abs(P @ Rm(x) - A @ (P @ Rp(x)) @ A).max() for x in (1.7 + 0.4j, 0.3 - 1.1j, 2.2))
        out = [f"Rcheck: {r1:.0e}"]
        for Nn in (2, 3, 4):
            if D**Nn > 600: continue
            labels = np.array(list(itertools.product(range(D), repeat=Nn)))
            G = g[labels].sum(1) % 2
            Sdiag = np.prod([Pi[labels[:, j]] for j in range(1, Nn, 2)], axis=0)    # solitons 2, 4, ...
            Pi1N = Pi[labels[:, 0]]**Nn
            ls = np.random.default_rng(3).uniform(-3, 3, Nn - 1)
            Tp, Tm = transferN(Rp, D, ls), transferN(Rm, D, ls)
            rhs = (Sdiag*(-1.0)**G*Pi1N)[:, None]*Tp*(1/Sdiag)[None, :]
            out.append(f"N={Nn}: {np.abs(Tm - rhs).max():.0e}")
        print(f"{alg}_{n} spinor, omega={om}: residuals  " + ";  ".join(out), flush=True)
