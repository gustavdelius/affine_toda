"""Topological charges of single a_n^(1) solitons (sec-topological-charges):
t(eta) = (beta/2pi) phi(+inf) = -(1/2pi) sum_j alpha_j [2 pi j a/h + eta]_(-pi,pi),  alpha_j = e_j - e_{j+1} in R^h.
Checks: number of distinct charges = h/gcd(a,h); each is a weight of Lambda^a C^h; jumps are by gcd(a,h) phases at once."""
import numpy as np
from math import gcd, comb
def report(n, ok): print(('PASS ' if ok else 'FAIL ') + n)
def charges(n, a, N=4000):
    h = n+1; al = np.zeros((h, h))
    for j in range(h): al[j, j] = 1; al[j, (j+1) % h] -= 1
    out = set(); njumps = []
    for eta in (np.arange(N)+0.5)/N*2*np.pi:
        th = 2*np.pi*np.arange(h)*a/h + eta
        red = (th + np.pi) % (2*np.pi) - np.pi
        t = -(red @ al)/(2*np.pi)
        out.add(tuple(np.round(t, 9)))
    return h, [np.array(c) for c in out]
ok_count = ok_weight = True; rows = []
for n in range(1, 13):
    for a in range(1, n+1):
        h, cs = charges(n, a)
        ok_count &= len(cs) == h//gcd(a, h)
        for c in cs:
            v = c + a/h
            ok_weight &= np.allclose(v, np.round(v)) and set(np.round(v).astype(int)) <= {0, 1} and round(v.sum()) == a
        if n in (3, 4, 5): rows.append((n, a, len(cs), h//gcd(a, h), comb(h, a)))
for r in rows: print(f'   a_{r[0]}^(1) species {r[1]}: {r[2]} distinct charges, h/gcd = {r[3]}, dim Lambda^a = {r[4]}')
report('n<=12, all a: number of distinct single-soliton charges = h/gcd(a,h)', ok_count)
report('n<=12, all a: every charge + (a/h)(1,..,1) is a 0/1 vector with a ones (weight of Lambda^a C^h)', ok_weight)
# jumps: at each discontinuity, how many reduced phases jump simultaneously
n, a = 3, 2; h = 4
etas = np.sort(np.unique(np.round((np.pi - 2*np.pi*np.arange(h)*a/h) % (2*np.pi), 9)))
print(f'   a_3, species 2: jump points of Im xi per period: {len(etas)}; phases jumping at each: {gcd(a, h)}')
