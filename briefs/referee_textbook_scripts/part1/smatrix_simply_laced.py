"""Real-coupling S-matrices of the simply-laced theories (chapters 12-13).

S_ab = prod_x {x}^{N_x}, {x} = (x-1)(x+1)/((x-1+B)(x+1-B)),
(x) = sinh(theta/2 + i pi x/2h)/sinh(theta/2 - i pi x/2h), with Dorey's block content
S_ab = prod_{p=1}^{h} {2p-1+(c(b)-c(a))/2}^{-N_p/2}, N_p = lambda_a . w^p gamma_b, gamma_b=(1-w^{-1})lambda_b,
w = s_black s_white (eq-dorey-formula; the same convention as variant 1 of foundations_code/hirota/dorey.py).
Note: the variant with w^{-p} gives a different, also unitary, crossing-symmetric and bootstrap-consistent
product that is NOT the physical S-matrix (it fails the e8 and a_n checks below).
Checks: unitarity, crossing, B -> 2-B invariance, e8 S_11 = {1}{11}{19}{29}, the a_n closed
form (eq-an-smatrix), the bootstrap (eq-bootstrap) at every fusing of e7 and e8 with masses
from the Perron-Frobenius vector, and the sinh-Gordon tree-level amplitude (sec-tree-level).
"""
import itertools
import cmath
import numpy as np

exec(open('toda_classical.py').read().split('# ---------------------------------------------------------------- root data')[1]
     .split('ok_pf = True')[0])


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


SIGN = -1   # x = 2p-1+(c(b)-c(a))/2; the other sign does not even give even multiplicities


def setup(typ, n):
    C, al, roots, h, A, nj, M2 = data(typ, n)
    lamw = np.linalg.inv(C) @ al
    col = {0: 1}; stack = [0]
    while stack:
        i = stack.pop()
        for j in range(n):
            if j != i and C[i, j] < 0 and j not in col:
                col[j] = -col[i]; stack.append(j)
    refl = lambda v: np.eye(n) - np.outer(v, v)
    sW = np.eye(n); sB = np.eye(n)
    for i in range(n):
        if col[i] == 1: sW = refl(al[i]) @ sW
        else: sB = refl(al[i]) @ sB
    w = sB @ sW; wi = np.linalg.inv(w)
    gam = [(np.eye(n) - wi) @ lamw[a] for a in range(n)]
    blocks = {}
    for a in range(n):
        for b in range(n):
            d = {}
            for p in range(1, h+1):
                N = int(round(lamw[a] @ np.linalg.matrix_power(w, p) @ gam[b]))
                if N:
                    x = 2*p - 1 + SIGN*(col[a] - col[b])//2
                    x %= 2*h
                    if x > h:                      # {2h-x} = {x}^{-1}
                        x, N = 2*h - x, -N
                    d[x] = d.get(x, 0) + N
            assert all(k % 2 == 0 for k in d.values())
            blocks[(a, b)] = {x: -k//2 for x, k in d.items() if k}
    pf = np.abs(np.linalg.eigh(C)[1][:, 0])
    orbits = [[np.linalg.matrix_power(w, p) @ gam[a] for p in range(h)] for a in range(n)]
    return h, blocks, pf/pf.min(), orbits


def paren(x, h, th):
    return cmath.sinh(th/2 + 1j*cmath.pi*x/(2*h))/cmath.sinh(th/2 - 1j*cmath.pi*x/(2*h))


def brace(x, h, B, th):
    return paren(x-1, h, th)*paren(x+1, h, th)/(paren(x-1+B, h, th)*paren(x+1-B, h, th))


def S(bl, h, B, th):
    r = 1
    for x, k in bl.items():
        r *= brace(x, h, B, th)**k
    return r


B = 0.37; ths = [0.31 + 0.2j, -0.77 + 1.1j, 1.3 + 2.4j]
for typ, n in [('e', 8), ('e', 7), ('d', 4), ('d', 6), ('a', 4)]:
    h, bl, mass, orbits = setup(typ, n)
    uni = all(abs(S(bl[(a, b)], h, B, t)*S(bl[(a, b)], h, B, -t) - 1) < 1e-10 for a in range(n) for b in range(n) for t in ths)
    sym = all(bl[(a, b)] == bl[(b, a)] for a in range(n) for b in range(n))
    dual = all(abs(S(bl[(a, b)], h, B, t) - S(bl[(a, b)], h, 2-B, t)) < 1e-10 for a in range(n) for b in range(n) for t in ths)
    # crossing S_ab(i pi - th) = S_{a bbar}(th); bbar found as the node whose amplitude matches
    cross = True
    for a in range(n):
        for b in range(n):
            ok = any(all(abs(S(bl[(a, b)], h, B, 1j*np.pi - t) - S(bl[(a, c)], h, B, t)) < 1e-9 for t in ths)
                     and abs(mass[c] - mass[b]) < 1e-9 for c in range(n))
            cross &= ok
    report(f'{typ}{n}: unitarity, S_ab=S_ba, B->2-B invariance, crossing to a conjugate of equal mass',
           uni and sym and dual and cross)

h, bl, mass, orbits = setup('e', 8)
light = int(np.argmin(mass))
report('e8: S_11 = {1}{11}{19}{29}', bl[(light, light)] == {1: 1, 11: 1, 19: 1, 29: 1})

# a_n closed form S_ab = prod_{p=|a-b|+1, step 2}^{a+b-1} {p}; nodes are labelled so that mass_a = sin(pi a/h)
ok_an = True
for n in (2, 3, 4, 5):
    h = n+1
    def San(a, b, th):
        r = 1
        for p in range(abs(a-b)+1, a+b, 2):
            r *= brace(p, h, B, th)
        return r
    # bootstrap a + b -> a+b (a+b<=n): u_ab^{a+b} = pi (a+b)/h; ubar_{a,c}^b = pi b/h ... check with generic formula below
    for a in range(1, n+1):
        for b in range(1, n+1):
            c = a + b
            if c > n: continue
            ma, mb, mc = [np.sin(np.pi*k/h) for k in (a, b, c)]
            ub_ac = np.arccos((ma**2 + mc**2 - mb**2)/(2*ma*mc))   # interior angle between sides a and c
            ub_bc = np.arccos((mb**2 + mc**2 - ma**2)/(2*mb*mc))
            for d in range(1, n+1):
                for t in ths:
                    lhs = San(d, c, t); rhs = San(d, a, t - 1j*ub_ac)*San(d, b, t + 1j*ub_bc)
                    ok_an &= abs(lhs - rhs) < 1e-9
report('a_n (n=2..5): closed form satisfies the bootstrap for every fusing a+b -> a+b', ok_an)
ok_match = True
for n in (2, 3, 4, 5):
    h_, bl_, mass_, _ = setup('a', n)
    for a in range(n):
        for b in range(n):
            closed = {}
            for p in range(abs(a-b)+1, a+b+2, 2):
                closed[p] = closed.get(p, 0) + 1
            ok_match &= all(abs(S(bl_[(a, b)], h_, B, t) - S(closed, h_, B, t)) < 1e-9 for t in ths)
report('a_n (n=2..5): Dorey formula = closed form prod_{p=|a-b|+1, step 2}^{a+b-1} {p} for all a, b', ok_match)

# bootstrap for e7, e8 using Dorey's rule for the fusings
for typ, n in [('e', 7), ('e', 8)]:
    h, bl, mass, orbits = setup(typ, n)
    nf = 0; ok = True
    for a, b, c in itertools.product(range(n), repeat=3):
        if not any(np.allclose(x + y + z, 0) for x in orbits[a] for y in orbits[b] for z in orbits[c]):
            continue
        ma, mb, mc = mass[a], mass[b], mass[c]
        if not (ma + mb > mc + 1e-9 and mb + mc > ma + 1e-9 and ma + mc > mb + 1e-9):
            continue
        nf += 1
        ub_ac = np.arccos((ma**2 + mc**2 - mb**2)/(2*ma*mc))
        ub_bc = np.arccos((mb**2 + mc**2 - ma**2)/(2*mb*mc))
        for d in range(n):
            for t in ths[:2]:
                lhs = S(bl[(d, c)], h, B, t)
                rhs = S(bl[(d, a)], h, B, t - 1j*ub_ac)*S(bl[(d, b)], h, B, t + 1j*ub_bc)
                ok &= abs(lhs - rhs) < 1e-8*max(1, abs(lhs))
    report(f'{typ}{n}: bootstrap holds at all {nf} fusings a+b->c given by Dorey\'s rule', ok)

# sinh-Gordon tree level: S = 1 + i M/(4 m^2 sinh th), M = -lambda, lambda = m^2 b^2  (V = m^2/b^2 (cosh b phi -1))
b2 = 1e-4; Bsg = (b2/(4*np.pi))/(1 + b2/(8*np.pi)); th = 0.83
exact = (cmath.sinh(th) - 1j*np.sin(np.pi*Bsg/2))/(cmath.sinh(th) + 1j*np.sin(np.pi*Bsg/2))
tree = 1 - 1j*b2/(4*np.sinh(th))
report('sinh-Gordon: exact S agrees with tree level 1 - i b^2/(4 sinh theta) to O(b^4)', abs(exact - tree) < 10*b2**2)
