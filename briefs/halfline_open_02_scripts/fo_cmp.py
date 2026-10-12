"""Fateev-Onofri hep-th/0207152 (r1c) R_1 = (x-1)(h/2+x)(h-x)(h/2+1-x)(h/2)(h/2-1)/[(h-1)(h/2+w+)(h/2-w+)(h/2-w-)(h/2+w-)],
w+- = w_r +- w_0, blocks (z) = sinh(th/2 + i pi z/2h)/sinh(th/2 - i pi z/2h); B_r: w_0 = x/2.
Compare with the d_n eps = 0 amplitude K_1 at h = 2n - 2, x = B/2, w_r = x/2."""
from dn_particle import *
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from t6 import K1_eps0
def R1_FO(th, h, x, w0, wr):
    b = lambda z: blk(z, th, h)
    wp, wm = wr + w0, wr - w0
    return (b(x - 1)*b(h/2 + x)*b(h - x)*b(h/2 + 1 - x)*b(mp.mpf(h)/2)*b(mp.mpf(h)/2 - 1)
            / (b(h - 1)*b(h/2 + wp)*b(h/2 - wp)*b(h/2 - wm)*b(h/2 + wm)))
for n in (4, 5, 6, 7):
    h = 2*n - 2
    w = max(abs(R1_FO(th, h, B/2, B/4, B/4)/K1_eps0(n)(th, B) - 1) for th in THS for B in BS)
    print(n, 'FO (r1c) at h = 2n-2, x = B/2, w_0 = w_r = x/2  vs  d_n eps=0 K_1:', mp.nstr(w, 3))
