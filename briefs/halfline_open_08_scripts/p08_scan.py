"""Scan -i Res of the pole p1 = 1/2 + (2n-1)B/2 of K_1 on b_{-n,n} across B (a_2)."""
from fractions import Fraction as Fr
import mpmath as mp
from p08_common import *

for n in (1, 2, 3, 4):
    lo, hi = Fr(1, 2*n - 1), Fr(2, 2*n - 1)
    print(f'n = {n}: window 1/(2n-1) = {lo} < B < {hi};  p1 = 1/2 + (2n-1)B/2; p3 = 3/2 - (2n-1)B/2 (crossed)')
    for k in range(1, 40):
        B = Fr(k, 40)
        if B >= 1:
            break
        p1 = Fr(1, 2) + (2*n - 1)*B/2
        if not (0 < p1 < Fr(3, 2)):
            continue
        K = K1nm(-n, n)
        o = K.order_at(p1, B)
        if o <= 0:
            print(f'   B={B}: order {o}')
            continue
        order, c = K.laurent(p1, B)
        val = c/I if order == 1 else c
        tag = 'below window' if B < lo else ('WINDOW' if B < hi else '')
        print(f'   B={str(B):6s} p1={str(p1):6s} order {order}  -i Res = {mp.nstr(val.real, 6):>12s}  {tag}')
