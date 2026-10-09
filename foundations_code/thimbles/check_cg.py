"""Independent check of Z_exact for L=1 via the Coulomb-gas (Fourier) series, in mpmath; plus a broader critical-point search."""
import numpy as np, mpmath as mp, itertools, time
from lie import An
from model import Chain
from exactZ import Z_grid
from newton_batch import find_cps
mp.mp.dps = 40
g = An(2); lam = g.fundamental_weights_rep(1)[(1,)]; kappa = 1.0
ch = Chain(g, 1, kappa, lam)
Gram = g.alpha @ g.alpha.T
al_lam = g.alpha @ lam  # integers
for hbar in (1.0, 0.5, 0.25):
    k = mp.mpf(kappa)/hbar; M = int(4*float(k) + 25)
    s = mp.mpf(0)
    fac = [mp.mpf(1)]
    for i in range(1, M+1): fac.append(fac[-1]*i)
    pw = [k**i for i in range(M+1)]
    for m in itertools.product(range(M+1), repeat=3):
        q2 = float(np.array(m) @ Gram @ np.array(m))
        ph = np.pi*float(np.dot(m, al_lam))
        s += pw[m[0]]*pw[m[1]]*pw[m[2]]/(fac[m[0]]*fac[m[1]]*fac[m[2]]) * mp.cos(ph) * mp.e**(-mp.mpf(hbar)*q2/4)
    lam2 = float(lam @ lam)
    Z = (2*mp.pi*hbar)**(-2) * mp.e**(-mp.pi**2*lam2/hbar) * mp.pi*hbar * mp.e**(-3*k) * s
    Zg = Z_grid(ch, hbar, h=0.02)
    print("hbar=%.2f  Coulomb-gas series Z=%s   grid Z=%.12e   rel diff=%.1e" % (hbar, mp.nstr(Z, 15), Zg.real, abs(float(Z)-Zg.real)/float(Z)))
t = time.time()
cps, cnt = find_cps(ch, nstart=60000, re_spread=7.0, im_spread=4.0, seed=7)
S = np.array([ch.action(u) for u in cps]); o = np.argsort(S.real)
print("broad search: %d critical points (%.0fs); those with |Im S|<=%.1f and Re S<30:" % (len(cps), time.time()-t, ch.Imax+0.5))
for i in o:
    if abs(S[i].imag) <= ch.Imax + 0.5 and S[i].real < 30:
        print("   S=%.6f%+.6fi  hits=%d  u=%s" % (S[i].real, S[i].imag, cnt[i], np.round(cps[i].ravel(), 4)))
