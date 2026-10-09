"""Dorey's root-system formula for the real-coupling S-matrix block content A_ab, and the
conjectured soliton transmission factor  X_ab(z) = prod_{x in A_ab} (z - cos(pi(x-1)/h))/(z - cos(pi(x+1)/h))."""
import numpy as np
from lie import Algebra


def colouring(A):
    r = A.r
    c = {1: 1}
    stack = [1]
    while stack:
        i = stack.pop()
        for j in range(1, r + 1):
            if j != i and abs(A.Cfin[i - 1, j - 1]) > 0.5 and j not in c:
                c[j] = -c[i]
                stack.append(j)
    return c


def blocks(A, variant=0):
    """returns dict (a,b) -> dict x -> multiplicity, x in 1..2h-1 (blocks {x})"""
    r, h = A.r, A.h
    al = A.alpha[1:]
    lam = A.Cinv @ al                    # rows: fundamental weights
    col = colouring(A)
    def refl(i):
        a = al[i - 1]
        return np.eye(r) - np.outer(a, a)
    wW = np.eye(r); wB = np.eye(r)
    for i in range(1, r + 1):
        if col[i] == 1:
            wW = refl(i) @ wW
        else:
            wB = refl(i) @ wB
    w = wB @ wW
    wi = np.linalg.inv(w)
    phi = {i: (np.eye(r) - wi) @ lam[i - 1] for i in range(1, r + 1)}
    res = {}
    for a in range(1, r + 1):
        for b in range(1, r + 1):
            d = {}
            for p in range(1, h + 1):
                if variant == 0:
                    N = lam[a - 1] @ np.linalg.matrix_power(wi, p) @ phi[b]
                    x = 2 * p - 1 + (col[a] - col[b]) // 2
                else:
                    N = lam[a - 1] @ np.linalg.matrix_power(w, p) @ phi[b]
                    x = 2 * p - 1 - (col[a] - col[b]) // 2
                N = int(round(N))
                if N:
                    xx = x % (2 * h)
                    d[xx] = d.get(xx, 0) + N
            res[(a, b)] = {x: m for x, m in d.items() if m}
    return res, col


def X_from_blocks(A, bl):
    """exponents N_q of (z - cos(pi q/h)), q in 0..h, from block content"""
    h = A.h
    N = {}
    def add(q, m):
        q = q % (2 * h)
        if q > h:
            q = 2 * h - q
        N[q] = N.get(q, 0) + m
    for x, m in bl.items():
        add(x - 1, m)
        add(x + 1, -m)
    return {q: m for q, m in N.items() if m}


if __name__ == '__main__':
    import sys
    for nm in sys.argv[1:]:
        A = Algebra(nm[0], int(nm[1:]))
        for v in (0, 1):
            bl, col = blocks(A, v)
            print(nm, 'variant', v, 'colours', col)
            for a in range(1, A.r + 1):
                for b in range(a, A.r + 1):
                    print(f"  S_{a}{b} blocks {dict(sorted(bl[(a,b)].items()))}  X: {dict(sorted(X_from_blocks(A, bl[(a,b)]).items()))}")
