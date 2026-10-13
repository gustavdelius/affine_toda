"""zeroth-order boundary-condition residual |phi'(0) + sum_j C_j alpha_j e^{psi_j(0)/2}| of the approach configurations."""
import mpmath as mp
from a5lib import *
mp.mp.dps = 50
def bcres(bg, C):
    h = bg.th.h; t = bg.taus(mp.mpf(0)); eh = [mp.sqrt(v) for v in bg.psi(0)]
    dphi = [t[(k-1) % h][1]/t[(k-1) % h][0] - t[k][1]/t[k][0] for k in range(h)]
    r = [dphi[k] + sum(C[j]*eh[j]*((k == j) - (k == (j+1) % h)) for j in range(h)) for k in range(h)]
    return max(abs(v) for v in r)
T = Theory(5); s = mp.sin(mp.pi/3); d0 = 1/(1+s); ce = 5*mp.pi/3
P = (1, 1, -1, -1, -1, -1)
for e in ('1e-2', '1e-4', '1e-6'):
    ee = mp.mpf(e)
    chi = Background(T, [(2, d0*mp.exp(I*(ce+ee))), (4, d0*mp.exp(-I*(ce+ee)))])
    dp = Background(T, [(2, d0*(1-ee)*mp.exp(I*ce)), (4, d0*(1-ee)*mp.exp(-I*ce))])
    print('a_5 ++---- e', e, 'chi-path', mp.nstr(bcres(chi, P), 3), ' d-path', mp.nstr(bcres(dp, P), 3))
T = Theory(2); s = mp.sin(mp.pi/3); d0 = 1/(1+s)
for e in ('1e-2', '1e-4', '1e-6'):
    ee = mp.mpf(e)
    chi = Background(T, [(1, d0*mp.exp(I*(mp.pi+ee))), (2, d0*mp.exp(-I*(mp.pi+ee)))])
    dp = Background(T, [(1, d0*(1-ee)*mp.exp(I*mp.pi)), (2, d0*(1-ee)*mp.exp(-I*mp.pi))])
    print('a_2 -++ e', e, 'chi-path', mp.nstr(bcres(chi, (-1, 1, 1)), 3), ' d-path', mp.nstr(bcres(dp, (-1, 1, 1)), 3))
exec(open('t4b.py').read().split('eta = ')[0].replace('mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60', 'mp.mp.dps = 50'))
d1s = mp.mpf('0.257789025012388'); d2s = mp.findroot(lambda d2: bg_of(d1s, d2).tau[2](0).real, mp.mpf('0.5637'))
for e in ('1e-2', '1e-4', '1e-6'):
    ee = mp.mpf(e)
    print('a_4 exc e', e, 'd2-path', mp.nstr(bcres(bg_of(d1s, d2s*(1-ee)), P), 3), ' d1-path', mp.nstr(bcres(bg_of(d1s*(1-ee), d2s), P), 3))
