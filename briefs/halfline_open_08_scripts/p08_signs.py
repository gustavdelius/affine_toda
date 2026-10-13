"""Residue signs -i Res of K_a on the a_2 boundary states (DG (6.7)), at B = 9/20 and other B."""
import sys
from fractions import Fraction as Fr
import mpmath as mp
from p08_common import *

def table(B, states):
    Bn = fnum(B)
    for st in states:
        for a in (1, 2):
            K = Ka2(a, *st)
            pl = K.strip_poles(B)
            out = []
            for p, o in pl.items():
                if o > 0:
                    order, c = K.laurent(p, B)
                    # check numerically
                    th0 = I*PI*fnum(p)/3
                    eps = mp.mpf('1e-12')
                    num = K(th0 + eps, Bn)*eps**order
                    assert abs(num/c - 1) < 1e-8, (st, a, p, num, c)
                    val = c/I if order == 1 else c
                    out.append(f'{p}(ord {order}): {mp.nstr(val, 6)}')
            print(f'  b_{st[0]},{st[1]}  K_{a}: ' + '; '.join(out))

B = Fr(9, 20)
states = [(0,0),(-1,1),(1,-1),(-2,2),(2,-2),(1,0),(0,1),(-3,3),(3,-3),(1,1)]
print('B = 9/20: poles p (theta = i pi p/3), and -i Res for simple poles (leading coeff for higher)')
table(B, states)
