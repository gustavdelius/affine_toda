# q -> 1/q (= conj q): spectra per charge sector compared with conj spectra at q (real rapidities)
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from spin_common import *
for n, alg, om in ((2, 'c', 0.3), (3, 'c', 1.61), (3, 'a', 2.37)):
    D, Wn, setup = make(n, alg); q = -np.exp(-1j*np.pi*om); R1, R2 = setup(q), setup(1/q)
    for Nn in (2, 3):
        tot, first = weights_multi(Wn, Nn); ls = np.random.default_rng(5).uniform(-2, 2, Nn-1)
        T1, T2 = transferN(R1, D, ls), transferN(R2, D, ls); w = 0; wm = 0
        for k in sorted({tuple(t) for t in tot}):
            idx = np.where(np.all(np.isclose(tot, k), axis=1))[0]
            e1 = np.linalg.eigvals(T1[np.ix_(idx, idx)]); e2 = np.linalg.eigvals(T2[np.ix_(idx, idx)])
            w = max(w, np.abs(np.sort_complex(np.round(np.conj(e1), 8)) - np.sort_complex(np.round(e2, 8))).max())
            wm = max(wm, np.abs(np.sort(np.abs(e1)) - np.sort(np.abs(e2))).max())
        print(f"{alg}_{n} N={Nn}: max|spec(1/q) - conj spec(q)| {w:.1e}; max moduli difference {wm:.1e}")
