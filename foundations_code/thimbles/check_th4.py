import numpy as np, time, sys
from lie import An
from model import Chain
from thimble_exact import thimble_integral
from exactZ import Z_grid
g = An(2); lam = g.fundamental_weights_rep(1)[(1,)]
ch = Chain(g, 1, 1.0, lam)
cps = list(np.load('cps_a2L1_k1.00.npy'))
Ua, Uc, Ud = cps[0], cps[2], cps[3]
mua, Va = ch.tangent(Ua); muc, Vc = ch.tangent(Uc); mud, Vd = ch.tangent(Ud)
nc, nd = 1, -1   # recomputed with these V in check_th2.py
for hbar in (0.15, 0.1):
    Zex = Z_grid(ch, hbar, h=0.02)
    Zc = ch.one_loop(Uc, hbar, Vc, muc)*(1+hbar*ch.two_loop_c1(Uc)); Zd = ch.one_loop(Ud, hbar, Vd, mud)*(1+hbar*ch.two_loop_c1(Ud))
    pred = (nc*Zc + nd*Zd)/Zex
    for rc, nr in ((640, 128), (2560, 256)):
        t = time.time()
        Za, dIm = thimble_integral(ch, Ua, hbar, Va, mua, nrho=nr, nphi=64, rho_cut=rc)
        res = (Zex + 2*Za)/Zex
        print("hbar=%.2f rho_cut=%d: residual (Z_exact - n_a Z_a - n_b Z_b)/Z = %.4e ; complex-pair 2-loop prediction = %.4e ; difference = %.2e  (%.0fs)" % (
            hbar, rc, res.real, pred.real, abs(res - pred), time.time()-t)); sys.stdout.flush()
