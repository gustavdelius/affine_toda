"""Threshold mechanism: the pole of K_2 on alpha = b_{-n+1,n-1} at x_w = 1/2 - (2n-1)B/2 (theta = i pi x_w/3) creates
b_{-n+1,n}.  Its residue -> 0 as B -> 1/(2n-1)^-, and for B > 1/(2n-1) the strip contains instead a simple ZERO at
theta_R = -i pi x_w/3 = i((2n-1) pi B/6 - pi/6): the reflection point of particle 2 in the double-triangle diagram.
Also: energy of the would-be state b_{-n+1,n} and the threshold e(alpha) + m_1."""
from fractions import Fraction as Fr
import mpmath as mp
from p08_common import *
for n in (2, 3):
    print(f'n = {n}: threshold B = 1/{2*n-1}')
    K2a = K1nm(n - 1, -n + 1)          # K_2 on b_{-n+1,n-1}
    for B in [Fr(1, 2*n - 1) - Fr(k, 400) for k in (40, 10, 2)] + [Fr(1, 2*n - 1) + Fr(k, 400) for k in (2, 10, 40)]:
        xw = Fr(1, 2) - (2*n - 1)*B/2
        if xw > 0:
            o, R = K2a.laurent(xw, B)
            s = f'pole at x={xw}: order {o}, -i Res = {mp.nstr((R/I).real, 6)}'
        else:
            o = K2a.order_at(-xw, B)
            s = f'at x={-xw}: order {o} (negative = zero)'
        Bn = fnum(B)
        de = E(-n + 1, n, Bn) - E(-n + 1, n - 1, Bn)
        print(f'   B={str(B):7s} {s};  e(b_{-n+1},{n}) - e(b_{-n+1},{n-1}) = {mp.nstr(de, 6)} m_1')
