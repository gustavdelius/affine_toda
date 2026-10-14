import numpy as np
from fold import folded, describe
def drev(r):
    p = [0]*(r+1); p[0], p[1], p[r-1], p[r] = r-1, r, 0, 1
    for j in range(2, r-1): p[j] = r - j
    return p
def drho(r):
    p = [0]*(r+1); p[0], p[r], p[1], p[r-1] = r, 1, r-1, 0
    for j in range(2, r-1): p[j] = r - j
    return p
def ds1(r):
    p = list(range(r+1)); p[0], p[1] = 1, 0; p[r-1], p[r] = r, r-1
    return p
import sys
which = sys.argv[1]
if which == 'a2n-1':
    for r in [4, 6, 8]:
        describe(folded('d', r, drev(r), f"a_{r-1}^(2) [d_{r} reflection]"))
if which == 'dn+1':
    for r in [5, 6, 7, 8]:
        describe(folded('d', r, ds1(r), f"d_{r-1}^(2) [d_{r} (01)(r-1 r)]"))
if which == 'a2n':
    for r in [4, 6, 8]:
        describe(folded('d', r, drho(r), f"a_{r-2}^(2) [d_{r} Z4]"))
if which == 'e':
    describe(folded('e', 6, [1, 3, 5, 2, 4, 6, 0][:0] or [None]*0 if False else None, "x") if False else None)
