"""Exact check of Lemma 6.13 (nested scales with signed exponents and cutoff b).

F_j(t) = int_{s_i <= s_par(i), s_j <= t} prod_{i in subtree(j)} s_i * max(s_i,b)^{-e_i} ds
is computed exactly as a piecewise sum of powers.  Compared with the bound
t^{2h} hat t^{-C} prod_j 2/min(E_j, 2h_j).
"""
import numpy as np
from collections import defaultdict
rng = np.random.default_rng(1)

def random_tree(nint):
    # random rooted tree with nint internal nodes; children lists; leaves implicit
    parent = [-1]
    for i in range(1, nint):
        parent.append(rng.integers(0, i))
    children = defaultdict(list)
    for i, p in enumerate(parent):
        if p >= 0:
            children[p].append(i)
    return parent, children

def analyse(children, e, b, t):
    nint = len(e)
    h = {}; C = {}
    def rec(j):
        hj = 1; Cj = e[j]
        for c in children[j]:
            rec(c); hj += h[c]; Cj += C[c]
        h[j] = hj; C[j] = Cj
    rec(0)
    E = {j: 2*h[j]-C[j] for j in range(nint)}
    if min(E.values()) <= 0:
        return None
    # exact piecewise integration
    # For s >= b, F_j(s) = sum_k A_k s^{alpha_k} (dict); for s<=b F_j(s)= low_j s^{2h_j}
    low = {}; high = {}
    def integ(j):
        for c in children[j]:
            integ(c)
        # below b: integrand s * b^{-e_j} * prod low_c s^{2h_c}
        lc = b**(-e[j]) * np.prod([low[c] for c in children[j]]) if children[j] else b**(-e[j])
        low[j] = lc/(2*h[j])  # coefficient of s^{2h_j}
        Fb = low[j]*b**(2*h[j])
        # above b: integrand s^{1-e_j} prod_c F_c(s)
        poly = {1.0 - e[j]: 1.0}
        for c in children[j]:
            newp = defaultdict(float)
            for a1, c1 in poly.items():
                for a2, c2 in high[c].items():
                    newp[a1+a2] += c1*c2
            poly = dict(newp)
        res = defaultdict(float)
        res[0.0] += Fb
        for a, cc in poly.items():
            if abs(a+1) < 1e-9:
                raise ValueError('log term')
            res[a+1] += cc/(a+1)
            res[0.0] -= cc*b**(a+1)/(a+1)
        high[j] = dict(res)
    integ(0)
    if t <= b:
        val = low[0]*t**(2*h[0])
    else:
        val = sum(cc*t**a for a, cc in high[0].items())
    that = max(t, b)
    bound = t**(2*h[0])*that**(-C[0])*np.prod([2/min(E[j], 2*h[j]) for j in range(nint)])
    return val, bound, E

worst = 0; count = 0
for trial in range(20000):
    nint = rng.integers(1, 7)
    parent, children = random_tree(nint)
    e = list(rng.uniform(-6, 4, size=nint))
    b = 10**rng.uniform(-6, 0)
    t = 10**rng.uniform(-7, 0.5)
    try:
        r = analyse(children, e, b, t)
    except ValueError:
        continue
    if r is None:
        continue
    val, bound, E = r
    count += 1
    if val <= 0:
        print('nonpositive integral?', val)
    worst = max(worst, val/bound)
print('trees tested', count, 'max ratio integral/bound', worst)
