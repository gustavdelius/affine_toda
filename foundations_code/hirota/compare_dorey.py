"""Compare Hirota transmission factors (data/X_<alg>.json from fluct_scan.py) with the block formula
X_ab(z) = prod_{x in A_ab} (z - cos(pi(x-1)/h)) / (z - cos(pi(x+1)/h)),  S_ab = prod_{x in A_ab} {x}
(A_ab from Dorey's root-system formula).  Also tabulates zeros/poles in the bound-state half-plane
(z = kappa/m_b in (0,1]), threshold zeros/poles (z = 0) and multiple zeros."""
import json, sys
import numpy as np
from lie import Algebra
from dorey import blocks, X_from_blocks


def X_dorey(A):
    bl, col = blocks(A, 1)
    out = {}
    for (a, b), d in bl.items():
        Xr = X_from_blocks(A, d)
        X = {}
        for q, m in Xr.items():
            assert m % 2 == 0
            X[q] = -m // 2
        out[(a, b)] = {q: m for q, m in X.items() if m}
        # block content of S_ab: {x}^(-m/2) with {2h-x} = {x}^{-1}
        S = {}
        for x, m in d.items():
            if x > A.h:
                x2, m2 = 2 * A.h - x, -m
            else:
                x2, m2 = x, m
            S[x2] = S.get(x2, 0) + m2
        out[(a, b, 'S')] = {x: -m // 2 for x, m in S.items() if m}
    return out


if __name__ == '__main__':
    for nm in sys.argv[1:]:
        A = Algebra(nm[0], int(nm[1:]))
        dat = json.load(open(f'data/X_{nm}.json'))
        Xd = X_dorey(A)
        bad = 0
        for key, v in dat['pairs'].items():
            a, b = map(int, key.split(','))
            Nh = {int(q): m for q, m in v['N'].items()}
            if Nh != Xd[(a, b)]:
                bad += 1
                print('MISMATCH', nm, a, b, Nh, Xd[(a, b)])
        print(nm, 'pairs', len(dat['pairs']), 'mismatches', bad, 'worst', dat['worst'])
