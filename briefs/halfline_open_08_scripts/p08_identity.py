"""Check: on beta = b_{-n,n}, alpha = b_{-n+1,n-1}, v = v_{n-1} = pi/2 - (2n-1) pi B/6, the pole of K_1^beta at
i u1, u1 = 2pi/3 - v, has residue (1/2) g^2 Gamma^2 K_2^alpha'(theta_R), theta_R = i(pi/3 - v), where
(i/2) g^2 = Res_{iv} K_1^alpha, i Gamma^2 = Res_{2 pi i/3} S_11, and K_2^alpha has a simple zero at theta_R."""
from fractions import Fraction as Fr
import mpmath as mp
from p08_common import *

def num_deriv(f, z, h=mp.mpf('1e-15')):
    return (f(z + h) - f(z - h))/(2*h)

S11 = S(1, 1, 2)
for n, B in ((2, Fr(9, 20)), (2, Fr(2, 5)), (2, Fr(3, 5)), (3, Fr(3, 10)), (4, Fr(1, 5)), (2, Fr(1, 4)), (3, Fr(3, 20))):
    Bn = fnum(B)
    pv = Fr(3, 2) - (2*n - 1)*B/2          # v = pi pv/3
    pu1 = Fr(1, 2) + (2*n - 1)*B/2         # u1 = pi pu1/3 = 2pi/3 - v
    pR = pu1 - 1                           # theta_R = i pi pR/3 = i(pi/3 - v)... check: pi/3 - v = pi(1 - pv)/3
    assert pR == 1 - pv
    alpha, beta = K1nm(-n + 1, n - 1), K1nm(-n, n)
    K2alpha = K1nm(n - 1, -n + 1)          # K_2 on b_{-n+1,n-1} = K_1 on b_{n-1,-n+1}
    o, Ra = alpha.laurent(pv, B)           # residue of K_1^alpha at iv
    assert o == 1
    g2 = Ra/(I/2)
    o, RS = S11.laurent(Fr(2), B)
    Gam2 = RS/I
    oR = K2alpha.order_at(pR, B) if pR != 0 else None
    thR = I*PI*fnum(pR)/3
    K2p = num_deriv(lambda z: K2alpha(z, Bn), thR)
    o1, Rex = beta.laurent(pu1, B)
    pred = g2*Gam2*K2p/2
    tag = 'window' if Fr(1, 2*n - 1) < B < Fr(2, 2*n - 1) else 'below window'
    print(f'n={n} B={B} ({tag}): order of K_2^alpha at theta_R: {oR};  g^2={mp.nstr(g2,6)}  Gamma^2={mp.nstr(Gam2,6)}')
    print(f'    exact Res = {mp.nstr(Rex, 10)},  (1/2) g^2 Gamma^2 K_2^alpha\'(theta_R) = {mp.nstr(pred, 10)},  ratio = {mp.nstr(Rex/pred, 12)}')
