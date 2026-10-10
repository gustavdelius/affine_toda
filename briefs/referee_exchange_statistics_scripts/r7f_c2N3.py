# Three c_2^(1) spinor solitons at alpha=0: breaking by charge sector at the physical q and at -q; also N=2 sign c_Q
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from spin_common import *            # QSIGN unset: setup(q) uses q as given
D, Wn, setup = make(2, 'c')
rng = np.random.default_rng(3)
for om in (0.3, 0.7):
    qphys = -np.exp(-1j*np.pi*om)
    for lab, qq in (('q_phys', qphys), ('-q_phys', -qphys)):
        Rp = setup(qq); Nn = 3
        tot, first = weights_multi(Wn, Nn); keys = sorted({tuple(t) for t in tot})
        dev = {k: 0 for k in keys}
        lss = [rng.uniform(-4, 4, 2) for _ in range(300)] + [np.array([a, b]) for a in (-0.005, 0.003, 0.01) for b in (-0.007, 0.004, 2.0)]
        for ls in lss:
            T = transferN(Rp, D, ls)
            for k in keys:
                idx = np.where(np.all(np.isclose(tot, k), axis=1))[0]
                dev[k] = max(dev[k], np.abs(np.abs(np.linalg.eigvals(T[np.ix_(idx, idx)])) - 1).max())
        br = {k: v for k, v in dev.items() if v > 1e-6}
        print(f"omega={om} {lab}: broken sectors (max||s|-1|) " + (", ".join(f"Q={k}:{v:.1e}" for k, v in br.items()) or "none"), flush=True)
# N=2: eigenvalue sign c_Q between q and -q
for om in (0.3,):
    qphys = -np.exp(-1j*np.pi*om); Rp, Rm = setup(qphys), setup(-qphys)
    tot, first = weights_multi(Wn, 2)
    T1, T2 = transferN(Rp, D, [0.83]), transferN(Rm, D, [0.83])
    for k in sorted({tuple(t) for t in tot}):
        idx = np.where(np.all(np.isclose(tot, k), axis=1))[0]
        e1 = np.sort_complex(np.round(np.linalg.eigvals(T1[np.ix_(idx, idx)]), 10)); e2 = np.sort_complex(np.round(np.linalg.eigvals(T2[np.ix_(idx, idx)]), 10))
        e2m = np.sort_complex(np.round(-np.linalg.eigvals(T2[np.ix_(idx, idx)]), 10))
        print(f"  N=2 Q={k}: spec equal {np.allclose(e1, e2)}, spec equal up to sign -1 {np.allclose(e1, e2m)}")
