import sys
import mpmath as mp
from a5lib import *
mp.mp.dps = 40
n = 4; T = Theory(n); h = 5
chi1 = 9*mp.pi/5; chi2 = 18*mp.pi/5
def bg_of(d1, d2):
    c1 = d1*mp.exp(I*chi1); c2 = d2*mp.exp(I*chi2)
    return Background(T, [(1, c1), (4, mp.conj(c1)), (2, c2), (3, mp.conj(c2))])
for d1, d2 in ((0.2579, 0.5637),):
    bg = bg_of(mp.mpf(d1), mp.mpf(d2))
    print([mp.nstr(bg.tau[j](0), 6) for j in range(h)])
P = (1, 1, -1, 1, -1)
def d2_of(d1, guess=mp.mpf('0.5637')):
    return mp.findroot(lambda d2: bg_of(d1, d2).tau[2](0).real, guess)
def R(d1, eta=mp.mpf('1e-15')):
    d2 = d2_of(d1); bg = bg_of(d1, d2); x = -eta
    t = [(bg.tau[j](x), bg.tau[j](x, 1)) for j in range(h)]
    eh = [mp.sqrt(v) for v in bg.psi(x)]
    dphi = [t[(k-1) % h][1]/t[(k-1) % h][0] - t[k][1]/t[k][0] for k in range(h)]
    return [ (dphi[k] + sum(P[j]*eh[j]*((k == j) - (k == (j+1) % h)) for j in range(h))).real for k in range(h)], d2, [bg.tau[j](0).real for j in range(h)]
for d1 in ('0.2', '0.2579', '0.3'):
    r, d2, t0 = R(mp.mpf(d1)); print(d1, mp.nstr(d2, 8), [mp.nstr(v, 5) for v in r])
