"""Coordinator check: O(beta^2) of the d_4 continued amplitudes vs Kim-Cho (29), (40) read in the TeX."""
import mpmath as mp
mp.mp.dps = 40
h = 6
P = lambda x, th: mp.sinh(th/2 + 1j*mp.pi*x/(2*h))/mp.sinh(th/2 - 1j*mp.pi*x/(2*h))
def br(x, B, th): return P(x-1, th)*P(x+1, th)/(P(x-1+B, th)*P(x+1-B, th))
S11 = lambda B, th: br(1, B, th)*br(5, B, th)
K1 = lambda B, th: P(h-B/2, th)*P(h/2-1, th)*P(h/2+1-B/2, th)/(P(h-1, th)*P(1-B/2, th)*P(h/2, th)*P(h/2-B/2, th))
KN = lambda B, th: K1(2-B, th)
Kcl = lambda th: -P(h/2-1, th)*P(h/2+1, th)/(P(1, th)*P(h-1, th)*P(h/2, th)**2)
Km = lambda B, th: K1(B, th)/Kcl(th)
KN2 = lambda B, th: KN(B, th + 1j*mp.pi/6)*KN(B, th - 1j*mp.pi/6)*S11(B, 2*th)
s3 = mp.sqrt(3)
kc1 = lambda th: 1j/24*(mp.tanh(th) - mp.sinh(th)/(mp.cosh(th)-1) + 2*mp.sinh(th)/(2*mp.cosh(th)+s3) - 2*mp.sinh(th)/(2*mp.cosh(th)-1))
kc2 = lambda th: 1j/24*(mp.tanh(th) - mp.sinh(th)/(mp.cosh(th)-1) - 4*mp.sinh(th)/(2*mp.cosh(th)-s3)
                        + 2*mp.sinh(th)/(2*mp.cosh(th)+s3) + 2*mp.sinh(th)/(2*mp.cosh(th)+1))
# coefficient of beta^2: B = beta^2/2pi + O(beta^4), so d/d(beta^2) = (1/2pi) d/dB at B=0
for th in (mp.mpf('0.37'), mp.mpf('1.21') + 0.1j):
    for name, f, ref in (('K^N_1', KN, kc1), ('K^N_2', KN2, kc2), ('K^(-)_1(th+i pi)', lambda B, t: Km(B, t + 1j*mp.pi), kc1),
                         ('K^(-)_1(th)', Km, kc1)):
        f0 = f(mp.mpf(0), th); d = mp.diff(lambda B: f(B, th), 0)/(2*mp.pi)
        print(f'{name:18s} th={mp.nstr(th,4)}: K(B=0)={mp.nstr(f0,8)}  |coef - KimCho| = {mp.nstr(abs(d - ref(th)), 3)}')
# Kim-Cho's exact self-dual K_1 = [1/2][3/2][5/2]^2[7/2][9/2], [x] = (x-1/2)(x+1/2)/((x-1/2+B/2)(x+1/2-B/2))
sq = lambda x, B, th: P(x - 0.5, th)*P(x + 0.5, th)/(P(x - 0.5 + B/2, th)*P(x + 0.5 - B/2, th))
KC = lambda B, th: sq(0.5, B, th)*sq(1.5, B, th)*sq(2.5, B, th)**2*sq(3.5, B, th)*sq(4.5, B, th)
for th in (mp.mpf('0.37'), mp.mpf('1.1')):
    d1 = mp.diff(lambda B: KC(B, th), 0)/(2*mp.pi)
    print('Kim-Cho self-dual: one loop vs (29):', mp.nstr(abs(d1 - kc1(th)), 3), ' |K^N_1/KC - 1| at B=0.6:', mp.nstr(abs(KN(mp.mpf('0.6'), th)/KC(mp.mpf('0.6'), th) - 1), 3))
