"""Breather amplitude of the d_n s=+-1 matrices vs candidate block forms."""
import mpmath as mp
from dn_scal import DN, blk
h = mp.mpf(1)/2
def cand(sig, mu, T, om):
    b = lambda x: blk(x, mu, T)
    if sig == 1:
        return (b(om/2 - T)*b(h - T)**2*b(-om/2 - h)*b(T/2 - om/2)*b(T/2 + om/2 + h)*b(T/2 - h)
                / (b(T/2)*b(h)*b(h - om/2)*b(T + om/2 - h)*b(T/2 - h + om/2)*b(T/2 + h - om/2)))
    return -b(om/2)*b(h - T)*b(-om/2 - h)*b(T/2 + om/2 + h)*b(T/2)/(b(T/2 + h)*b(T/2 + om/2))
for n, om, mu in ((4, '0.37', ('0.13', '0.21')), (4, '0.61', ('-0.17', '0.11')), (5, '0.53', ('0.21', '-0.07')),
                  (6, '0.29', ('0.08', '0.15'))):
    D = DN(n, mp.mpf(om)); mu = mp.mpc(*mu)
    for sig in (1, -1):
        kb = D.breather(mu, sig)
        print(f'n={n} om={om} sig={sig}: K_B = {mp.nstr(kb, 15)}, candidate ratio {mp.nstr(kb/cand(sig, mu, D.T, D.om), 15)}', flush=True)
