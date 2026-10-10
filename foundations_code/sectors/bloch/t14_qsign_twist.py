"""The sign of q as a Bloch twist (book: sec-crossing-sign, prp-bloch-sectors).

For N spinor solitons, the Bethe-Yang transfer matrix of soliton 1 at -q, in every charge sector Q and at alpha = 0,
has the spectrum of the transfer matrix at q with the twist exp(i alpha.mu_1), alpha = pi N gamma, gamma = (1,...,1),
times a phase c_Q common to the sector. The twist is trivial for even N, in particular for Q = 0.
"""
from spin_common import *
for n, alg, Nn, om in ((2, 'c', 2, 0.3), (2, 'c', 3, 0.3), (2, 'c', 3, 2.37), (2, 'c', 4, 0.7),
                       (3, 'c', 2, 1.61), (3, 'c', 3, 0.3013), (2, 'a', 3, 1.61), (3, 'a', 2, 0.3013)):
    D, Wn, setup = make(n, alg)
    tot, first = weights_multi(Wn, Nn)
    q = np.exp(-1j*np.pi*om); Rp, Rm = setup(q), setup(-q)
    rng = np.random.default_rng(7); lss = [rng.uniform(-3, 3, Nn - 1) for _ in range(3)]
    Tp = [transferN(Rp, D, ls) for ls in lss]; Tm = [transferN(Rm, D, ls) for ls in lss]
    a = np.pi*Nn*np.ones(n); worst = 0; nsec = 0
    for kk in sorted({tuple(t) for t in tot}):
        idx = np.where(np.all(np.isclose(tot, kk), axis=1))[0]; nsec += 1
        best = 9e9
        for c in (1, -1, 1j, -1j):
            d = max(np.abs(np.sort_complex(np.round(c*np.linalg.eigvals(np.exp(1j*(first[idx] @ a))[:, None]*Tp[k][np.ix_(idx, idx)]), 9))
                           - np.sort_complex(np.round(np.linalg.eigvals(Tm[k][np.ix_(idx, idx)]), 9))).max() for k in range(len(lss)))
            best = min(best, d)
        worst = max(worst, best)
    print(f"{alg}_{n} spinor, N={Nn}, omega={om}: spec(-q, alpha=0) = c_Q spec(q, alpha=pi N gamma) in all {nsec} charge sectors, "
          f"max difference {worst:.1e}", flush=True)
