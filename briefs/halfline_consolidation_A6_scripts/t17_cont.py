"""Continuation omega -> -2/B of the d_n s = +-1 breather amplitudes; closed forms in book blocks; axioms with the book's
d_n S-matrix; classical limit; one loop against Kim-Cho hep-th/9505138 (K1), (K2) for d_4."""
import numpy as np
import mpmath as mp
import cmath
from dn_scal import blk

mp.mp.dps = 30
TC = '/home/gustav/Git/affine_toda/code/'
ns = {}
src = open(TC + 'smatrix_simply_laced.py').read().split('\nB = 0.37')[0]
src = src.replace("open('toda_classical.py')", "open('" + TC + "toda_classical.py')")
exec(src, ns)
setup, Sbook = ns['setup'], ns['S']
H = mp.mpf(1)/2


def KB_imag(sig, mu, T, om, cdd=()):
    b = lambda x: blk(x, mu, T)
    if sig == 1:
        k = (b(om/2 - T)*b(H - T)**2*b(-om/2 - H)*b(T/2 - om/2)*b(T/2 + om/2 + H)*b(T/2 - H)
             / (b(T/2)*b(H)*b(H - om/2)*b(T + om/2 - H)*b(T/2 - H + om/2)*b(T/2 + H - om/2)))
    else:
        k = -b(om/2)*b(H - T)*b(-om/2 - H)*b(T/2 + om/2 + H)*b(T/2)/(b(T/2 + H)*b(T/2 + om/2))
    c = T/2 - H
    for x in cdd:          # CDD factor P_x = (T/2 - x)_I (T/2 + x)_I of the scalar factor, at mu_a and mu_b
        for y in (T/2 - x, T/2 + x):
            k *= blk(y, mu - c, T)*blk(y, mu + c, T)
    return k


def K1(sig, th, B, n, cdd_k=()):
    """continued amplitude; cdd_k: list of k with x = k omega/2."""
    om = -2/B
    T = (n - 1)*om
    mu = -1j*T*th/mp.pi
    return KB_imag(sig, mu, T, om, [k*om/2 for k in cdd_k])


def bb(x, th, h):
    return mp.sinh(th/2 + 1j*mp.pi*x/(2*h))/mp.sinh(th/2 - 1j*mp.pi*x/(2*h))


def K1_closed(sig, th, B, n):
    h = 2*n - 2
    b = lambda x: bb(x, th, h)
    if sig == -1:
        return -b(1)*b(h/2)*b(h/2 + 1 - B/2)/(b(1 - B/2)*b(h + B/2)*b(h/2 - B/2)*b(h/2 + 1))
    return (b(h/2 - 1)*b(h/2 + 1 - B/2)*b(h/2 + B/2)*b(1 + B/2)*b(B/2)*b(h - 1 - B/2)
            / (b(h - 1)*b(h + B/2)**2*b(1 - B/2)*b(h/2)*b(h/2 + 1 + B/2)*b(h/2 - 1 - B/2)))


if __name__ == '__main__':
    # 1. closed forms, unitarity, crossing-unitarity, classical limit, n = 4..7
    dev, devax, devcl = 0, 0, 0
    for n in (4, 5, 6, 7):
        h, bl, mass, _ = setup('d', n)
        i1 = 0                              # node 0 of the chain = vector particle 1 (lightest)
        for sig in (1, -1):
            for Bv in (mp.mpf('0.37'), mp.mpf('1.3')):
                for th in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
                    k = K1(sig, th, Bv, n)
                    dev = max(dev, abs(k/K1_closed(sig, th, Bv, n) - 1))
                    devax = max(devax, abs(K1(sig, th, Bv, n)*K1(sig, -th, Bv, n) - 1),
                                abs(K1(sig, th, Bv, n)*K1(sig, th + 1j*mp.pi, Bv, n)
                                    / complex(Sbook(bl[(i1, i1)], h, float(Bv), complex(2*th))) - 1))
            th = mp.mpc('0.6', '0.2')
            devcl = max(devcl, abs(K1(sig, th, mp.mpf('1e-15'), n) - 1))
    print(f'continuation = closed forms [{mp.nstr(dev, 3)}]; unitarity and crossing-unitarity with book S_11 '
          f'[{mp.nstr(devax, 3)}]; B -> 0 limit is 1 [{mp.nstr(devcl, 3)}]  (n = 4..7, sig = +-1)')

    # 2. one loop, d_4: d/d(beta^2) at 0 with B = beta^2/(2 pi) + O(beta^4)
    h = 6
    kc1 = lambda t: 1j/24*(mp.sinh(t)/mp.cosh(t) - mp.sinh(t)/(mp.cosh(t) - 1) + 2*mp.sinh(t)/(2*mp.cosh(t) + mp.sqrt(3))
                          - 2*mp.sinh(t)/(2*mp.cosh(t) - 1))
    kc2 = lambda t: 1j/24*(mp.sinh(t)/mp.cosh(t) - mp.sinh(t)/(mp.cosh(t) - 1) - 4*mp.sinh(t)/(2*mp.cosh(t) - mp.sqrt(3))
                          + 2*mp.sinh(t)/(2*mp.cosh(t) + mp.sqrt(3)) + 2*mp.sinh(t)/(2*mp.cosh(t) + 1))
    h4, bl4, _, _ = setup('d', 4)
    def S11(t, Bv):
        r = mp.mpf(1)
        for x, kk in bl4[(0, 0)].items():
            r *= (bb(x - 1, 2*t, h4)*bb(x + 1, 2*t, h4)/(bb(x - 1 + Bv, 2*t, h4)*bb(x + 1 - Bv, 2*t, h4)))**kk
        return r
    e = mp.mpf('1e-12')
    def d1(f, t):
        return (f(e, t) - f(-e, t))/(2*e)/(2*mp.pi)
    ths = [mp.mpf('0.37'), mp.mpc('1.2', '0.3'), mp.mpc('-0.8', '0.1')]
    mp.mp.dps = 40
    for sig in (1, -1):
        for cdd in ((), (0, 1, 2)):
            f1 = lambda Bv, t: K1(sig, t, Bv, 4, cdd)
            f2 = lambda Bv, t: K1(sig, t + 1j*mp.pi/6, Bv, 4, cdd)*K1(sig, t - 1j*mp.pi/6, Bv, 4, cdd)*mp.mpc(S11(t, Bv))
            r1 = max(abs(d1(f1, t) - kc1(t)) for t in ths)
            r1s = max(abs(d1(f1, t) - kc1(t + 1j*mp.pi)) for t in ths)
            r2 = max(abs(d1(f2, t) - kc2(t)) for t in ths)
            r2s = max(abs(d1(f2, t) - kc2(t + 1j*mp.pi)) for t in ths)
            print(f'sig={sig:+d} CDD k={cdd}: |K_1 one loop - Kim-Cho (K1)| = {mp.nstr(r1, 3)}; vs (K1) at theta+i pi: '
                  f'{mp.nstr(r1s, 3)}; |K_2 = K_1 K_1 S_11 one loop - (K2)| = {mp.nstr(r2, 3)}; at theta+i pi {mp.nstr(r2s, 3)}')
