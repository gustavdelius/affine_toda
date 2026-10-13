"""Semiclassical (omega -> infinity) limit of the scalar factors: (beta^2/4h) d(-i log A)/d theta, beta^2 = 4 pi/(omega+1),
as in an_nondiagonal_scalar.py check 7, for a_n A^(+-) (book formulas) and d_n A^(sig) (this work); compared with the
Neumann time advance of a_n, (1/2) ln(1 - sin^2(pi/h)/cosh^2 theta), and with (1/2)(beta^2/4h) d/dtheta(-i log F(2 theta))."""
import mpmath as mp
from dn_scal import DN, logprod, g2, blk

mp.mp.dps = 30
TCns = {'mp': mp}
_src = open('/home/gustav/Git/affine_toda/textbook_code/an_nondiagonal_scalar.py').read()
exec(_src.split('# ---------------------------------------------------------------- 3. scalar factor and breather identity')[1].split('d3 = 0')[0], TCns)
exec(_src.split('# ---------------------------------------------------------------- 9. s = -1: scalar factor and breather identity')[1].split('mp.mp.dps = 30')[0], TCns)
scalar_data, scalar_minus = TCns['scalar_data'], TCns['scalar_minus']


def semi(Afun, T, om, h, th):
    f = lambda t: mp.log(Afun(-1j*T*t/mp.pi))
    d = (f(th + mp.mpf('1e-8')) - f(th - mp.mpf('1e-8')))/mp.mpf('2e-8')
    return mp.re(-1j*d)*4*mp.pi/(om + 1)/(4*h)


th = mp.mpf('0.6')
for omv in (40, 80, 160):
    om = mp.mpf(omv)
    out = []
    for n in (2, 3):
        h = n + 1
        T, Ap, F11 = scalar_data(n, om)
        Am = scalar_minus(n, om)
        neu = mp.log(1 - (mp.sin(mp.pi/h)/mp.cosh(th))**2)/2
        fF = semi(lambda m: F11(2*m), T, om, h, th)
        out.append(f'a_{n}: A+ {mp.nstr(semi(Ap, T, om, h, th), 8)}, A- {mp.nstr(semi(Am, T, om, h, th), 8)}, '
                   f'Neumann {mp.nstr(neu, 8)}, F11(2th)/2 {mp.nstr(fF/2, 8)}')
    for n in (4, 5):
        h = 2*n - 2
        D = DN(n, om)
        T = D.T
        r = [mp.nstr(semi(lambda m: D.A(m, s), T, om, h, th), 8) for s in (1, -1)]
        r0 = mp.nstr(semi(D.A0, T, om, h, th), 8)
        fF = semi(lambda m: D.F(2*m), T, om, h, th)
        out.append(f'd_{n}: A(sig=+1) {r[0]}, A(sig=-1) {r[1]}, A0 {r0}, F(2th)/2 {mp.nstr(fF/2, 8)}, '
                   f'(1/2)ln(1-sin^2(pi/h)/cosh^2) {mp.nstr(mp.log(1 - (mp.sin(mp.pi/h)/mp.cosh(th))**2)/2, 8)}')
    print(f'omega={omv}:')
    for o in out:
        print('   ' + o, flush=True)
