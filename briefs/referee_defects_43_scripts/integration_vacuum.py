"""Integration check, ch. 43 item 4: prp-dpt-vacuum, the sigma < 0 direction phi = psi = t v, and the total static energy
with linear layers of width e^{-A t/2} (A = max_j beta alpha_j.v).  Uses defect_common.defect_potential, an_data."""
import sys, numpy as np
sys.path.insert(0, 'code')
from defect_common import defect_potential, an_data, bulk_potential
rng = np.random.default_rng(5)
for r in (1, 2, 3):
    d = an_data(r); h = r + 1
    for typ in (1, -1):
        v = rng.normal(size=r); v /= np.linalg.norm(v)
        A = max(d['alpha'] @ v)
        sig = -0.7
        for t in (2.0, 6.0, 12.0):
            Bd = defect_potential(t*v, t*v, sig, r, 1.0, 1.0, 'real', typ)[0].real
            s = sig if typ == 1 else 1/sig
            closed = s*(np.exp(t*(d['alpha'] @ v)).sum() - h)
            # total static energy: linear layers on both sides, width l = e^{-A t/2}
            l = np.exp(-A*t/2); xs = np.linspace(0, l, 2001)
            Vs = np.array([bulk_potential(t*v*(1 - x/l), r)[0].real for x in xs])
            Ebulk = 2*(0.5*t**2/l + np.trapz(Vs, xs))
            print(f'r={r} typ={typ:+d} t={t:5.1f}: B_d={Bd: .3e} closed={closed: .3e} bulk={Ebulk: .3e} total={Bd+Ebulk: .3e}')
