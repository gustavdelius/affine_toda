# Lemma: Rcheck(x;-q) = (Pi x 1) Rcheck(x;q) (Pi x 1), Pi = (-1)^{# minus signs}.
# Consequence (derived by pairing collisions): T_1(-q) = S [ Pi_tot Pi_1^N T_1(q) ] S^{-1}, S a product of one-soliton parities:
#   N odd: S = Pi_2 Pi_4 ... Pi_{N-1};  N even: S = Pi_2 Pi_4 ... Pi_N.
import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from spin_common import *
def kronlist(ms):
    out = np.ones((1, 1))
    for m in ms: out = np.kron(out, m)
    return out
for n, alg in ((2, 'c'), (3, 'c'), (2, 'a'), (3, 'a')):
    D, Wn, setup = make(n, alg); P = flip(D, D)
    Pi = np.diag([(-1.0)**sum(1 for w in Wn[i] if w < 0) for i in range(D)]); I = np.eye(D)
    for om in (0.3, 0.7, 1.61, 2.37):
        q = np.exp(-1j*np.pi*om); Rp, Rm = setup(q), setup(-q)
        lem = max(np.abs(P@Rm(x) - np.kron(Pi, I)@P@Rp(x)@np.kron(Pi, I)).max() for x in (0.37, 1.9, 7.3, 1.2+0.5j))
        line = f"{alg}_{n} omega={om}: lemma residual {lem:.1e}"
        for Nn in (2, 3, 4, 5):
            if D**Nn > 1100: continue
            ls = np.random.default_rng(Nn).uniform(-3, 3, Nn - 1)
            Tp, Tm = transferN(Rp, D, ls), transferN(Rm, D, ls)
            Pis = [kronlist([Pi if k == j else I for k in range(Nn)]) for j in range(Nn)]
            Ptot = kronlist([Pi]*Nn)
            S = np.eye(D**Nn)
            for j in range(1, Nn, 2): S = S @ Pis[j]   # Pi_2, Pi_4, ... (0-based 1,3,...)
            if Nn % 2 == 1:
                S = np.eye(D**Nn)
                for j in range(1, Nn-1, 2): S = S @ Pis[j]
            rhs = S @ Ptot @ np.linalg.matrix_power(Pis[0], Nn) @ Tp @ S
            line += f"; N={Nn}: {np.abs(Tm - rhs).max():.1e}"
        print(line, flush=True)
