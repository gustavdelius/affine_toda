"""Spinor solitons, N=2: scan the Bloch torus exp(i alpha.m_1) for the neutral sector Q=0 (and Q=e_1): min over alpha of max ||s|-1|."""
from spin_common import *
import itertools, sys
n = int(sys.argv[1]); alg = sys.argv[2]; steps = int(sys.argv[3])
D, Wn, setup = make(n, alg)
tot, first = weights_multi(Wn, 2)
lx = np.concatenate([np.linspace(0.01, 5, 120), np.linspace(5, 25, 30)])
grid = [np.array(a)*2*np.pi/steps for a in itertools.product(range(steps), repeat=n)]
for om in (2.37, 1.61, 0.73, 0.3):
    q = np.exp(-1j*np.pi*om); Rp = setup(q); Rs = [Rp(np.exp(l)) for l in lx]
    for name, Qv in (("Q=0", np.zeros(n)), ("Q=e1", np.eye(n)[0])):
        idx = np.where(np.all(np.isclose(tot, Qv), axis=1))[0]
        blocks = [R[np.ix_(idx, idx)] for R in Rs]; best = (9e9, None); n_unbroken = 0
        for a in grid:
            Om = np.exp(1j*(first[idx] @ a))
            d = 0
            for B in blocks:
                d = max(d, np.abs(np.abs(np.linalg.eigvals(Om[:, None]*B)) - 1).max())
                if d > best[0] and d > 1e-6: break
            if d < 1e-6: n_unbroken += 1
            if d < best[0]: best = (d, a/np.pi)
        print(f"{alg}_{n} omega={om} {name} (dim {len(idx)}): alpha=0 dev {max(np.abs(np.abs(np.linalg.eigvals(B))-1).max() for B in blocks):.2e}; "
              f"best alpha/pi={np.round(best[1],3)} dev {best[0]:.2e}; unbroken grid points {n_unbroken}/{len(grid)}", flush=True)
