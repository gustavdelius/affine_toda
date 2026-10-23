import numpy as np
from allS import *
import allS as AS
R0 = lambda theta: S.R(np.exp(mu*theta))                 # pure R-matrix factor, no scalar
AS.S0 = R0
def label(z):
    s = np.angle(z)/np.pi; best = None
    for L in range(-14, 15):
        for sig, nm in ((0, '+'), (1, '-'), (0.5, '+i'), (-0.5, '-i')):
            d_ = (s - (sig - L*omega)) % 2; d_ = min(d_, 2 - d_)
            if best is None or d_ < best[0]: best = (d_, L, nm)
    return f"{best[2]}q^{best[1]:+d}" + ("" if abs(abs(z) - 1) < 1e-6 and best[0] < 1e-6 else f"(?|z|={abs(z):.3f})")
res = {}
for a, b in [(3, 1), (3, 2), (1, 1), (1, 2), (2, 2)]:
    ph = 2*np.pi*(np.arange(40) + 0.31)/40
    xs = 1.27*np.exp(1j*ph); ths = np.log(xs)/mu
    ks = np.array([K_top(a, b, th)[0] for th in ths])
    for d in range(0, 13):
        V = np.vander(xs, d+1, increasing=True)
        Am = np.hstack([V, -(ks[:, None]*V[:, 1:])]); sol, *_ = np.linalg.lstsq(Am, ks, rcond=None)
        num, den = sol[:d+1], np.concatenate([[1], sol[d+1:]])
        err = np.abs(np.polyval(num[::-1], xs)/np.polyval(den[::-1], xs) - ks).max()/np.abs(ks).max()
        if err < 1e-9: break
    zs, ps = np.roots(num[::-1]), np.roots(den[::-1])
    # cancel common roots
    zs_l, ps_l = list(zs), list(ps)
    for z in list(zs_l):
        m = [p for p in ps_l if abs(p - z) < 1e-6]
        if m: zs_l.remove(z); ps_l.remove(m[0])
    lead = num[np.nonzero(np.abs(num) > 1e-12)[0][-1]]/den[np.nonzero(np.abs(den) > 1e-12)[0][-1]]
    print(f"k_{a}{b}(x): degree {d}, fit {err:.0e}")
    print(f"     zeros at x = {sorted(label(z) for z in zs_l)}")
    print(f"     poles at x = {sorted(label(p) for p in ps_l)}")
    res[(a, b)] = (zs_l, ps_l, num, den)
np.save('kfits.npy', np.array([res], dtype=object), allow_pickle=True)
