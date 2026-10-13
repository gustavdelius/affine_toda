import mpmath as mp
from a5lib import *
for dps in (50, 90):
    mp.mp.dps = dps
    n = 5; T = Theory(n); h = 6; s = mp.sin(mp.pi/3); d0 = 1/(1+s); chi_end = 5*mp.pi/3
    Om = T.mass(1)*mp.cosh(mp.mpf('0.6'))
    for e in ('1e-5', '1e-7', '1e-9'):
        bg = Background(T, [(2, d0*(1-mp.mpf(e))*mp.exp(I*chi_end)), (4, d0*(1-mp.mpf(e))*mp.exp(-I*chi_end))])
        opn, K, Kf, res = kmatrix(bg, (1, 1, -1, -1, -1, -1), Om)
        print(dps, e, mp.nstr(K[(1,1)], 12), mp.nstr(abs(K[(5,1)]), 12), 'unit', mp.nstr(unitarity_defect(opn, Kf), 3), 'res', mp.nstr(res, 3),
              'pde', mp.nstr(pde_residual(bg, bg.wave(1, I*mp.sqrt(Om**2-1)), mp.mpf('-0.01'), Om), 3))
