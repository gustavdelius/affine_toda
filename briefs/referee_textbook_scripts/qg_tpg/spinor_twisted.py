"""Spinor of U_q(d_{l+1}^(2)) in the conventions of ch. 7 (sec-tpg 'When it fails').

Finite subalgebra (remove node 0, both end nodes short): U_q(b_l).
1. builds e_i,f_i,k_i (i=0..l) on (C^2)^{(x)l}, checks all relations incl. q-Serre
   with q_i = q^{alpha_i^2/2} (long q, short q^{1/2});
2. dimension of the commutant of U_q(b_l) on S(x)S  (= l+1 iff multiplicity-free);
3. solves Jimbo's equations (homogeneous gradation, book coproduct): dimension of
   the solution space, eigenvalues of Rc(z) vs the twisted TPG of DGZ 1996 sec. 4.3
   in book units: rho_{Lambda^{l-j}} = prod_{i=1}^j <i/2>_{s_i},
   <a>_- = (1 - z q^{2a})/(z - q^{2a}),  <a>_+ = (1 + z q^{2a})/(z + q^{2a}),
   and vs the untwisted rule eq-tpg-rule (all signs -).
"""
import numpy as np
from functools import reduce
from qaff import check_relations, solve_jimbo, coprod


def commutant_sv(pairs, kd):
    import scipy.sparse as sps
    M = kd.shape[0]
    key = [tuple(np.round(r, 10)) for r in kd]
    allowed = [(a, b) for a in range(M) for b in range(M) if key[a] == key[b]]
    idx = {ab: n for n, ab in enumerate(allowed)}
    G = 0
    for a_, b_ in pairs:
        ri, ci, va = [], [], []
        for (u, c), n in idx.items():
            vs = np.nonzero(a_[c])[0]; ri.extend(u*M+vs); ci.extend([n]*len(vs)); va.extend(a_[c][vs])
        for (c, v), n in idx.items():
            us = np.nonzero(b_[:, c])[0]; ri.extend(us*M+v); ci.extend([n]*len(us)); va.extend(-b_[:, c][us])
        L = sps.csr_matrix((va, (ri, ci)), shape=(M*M, len(allowed)))
        G = G + (L.conj().T @ L)
    w = np.linalg.eigvalsh(G.toarray())
    return np.sqrt(np.abs(w))

rng = np.random.default_rng(7)
sp_ = np.array([[0, 1], [0, 0]], dtype=complex); sm_ = sp_.T.copy(); I2 = np.eye(2)


def site(op, i, l):
    return reduce(np.kron, [op if j == i else I2 for j in range(l)])


def spinor(l, t):
    q = t*t
    sz = [site(np.diag([0.5, -0.5]), i, l) for i in range(l)]   # s_i
    e = [site(sm_, 0, l)] + [site(sp_, i, l)@site(sm_, i+1, l) for i in range(l-1)] + [site(sp_, l-1, l)]
    f = [m.conj().T.copy() for m in e]
    # h_0 = -2 s_1 (alpha_0 = -eps_1 short), h_i = s_i - s_{i+1} (long), h_l = 2 s_l (short)
    h = [-2*sz[0]] + [sz[i]-sz[i+1] for i in range(l-1)] + [2*sz[l-1]]
    qi = [t] + [q]*(l-1) + [t]
    k = [np.diag(qi[i]**np.diag(h[i])) for i in range(l+1)]
    A = np.zeros((l+1, l+1), int)
    for i in range(l+1):
        A[i, i] = 2
    A[0, 1], A[1, 0] = -2, -1
    for i in range(1, l-1):
        A[i, i+1] = A[i+1, i] = -1
    A[l-1, l], A[l, l-1] = -1, -2
    if l == 1:
        raise ValueError
    return (e, f, k), A.tolist(), qi


def ang(a, z, q, s):
    return (1 + s*z*q**(2*a))/(z + s*q**(2*a))


def comb(n, k):
    from math import comb as c
    return c(n, k)


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def match(ev, want, tol=1e-6):
    ev = list(ev); ok = True
    for val, m in want:
        hits = [i for i, v in enumerate(ev) if abs(v - val) < tol*max(1, abs(val))]
        if len(hits) != m:
            ok = False
        for i in sorted(hits, reverse=True):
            ev.pop(i)
    return ok and not ev


t = complex(0.83, 0.41); q = t*t
for l in (2, 3, 4):
    rep, A, qi = spinor(l, t)
    err, bad = check_relations(*rep, A, qi, verbose=True)
    report(f'l={l}: spinor is a rep of U_q(d_{l+1}^(2)) (max residual {err:.1e})', not bad)
    # commutant of U_q(b_l) (nodes 1..l) on S(x)S
    e, f, k = rep
    sub = ([e[i] for i in range(1, l+1)], [f[i] for i in range(1, l+1)], [k[i] for i in range(1, l+1)])
    D = coprod(sub, sub)
    kd = np.array([np.diag(D[3*i+2]) for i in range(l)]).T
    pairs = [(D[3*m], D[3*m]) for m in range(l)] + [(D[3*m+1], D[3*m+1]) for m in range(l)]
    Gsv = commutant_sv(pairs, kd)
    dimcomm = int((Gsv < 1e-6).sum())
    report(f'l={l}: commutant of U_q(b_l) on S(x)S has dim {dimcomm} = l+1 (multiplicity-free); smallest sv {np.round(Gsv[:l+3],8)}', dimcomm == l+1)
    x = complex(rng.normal(), rng.normal()); y = complex(rng.normal(), rng.normal()); z = x/y
    R, sv2 = solve_jimbo(rep, x, y)
    report(f'l={l}: Jimbo solution space 1-dim (sing. values {abs(sv2[0]):.1e}, {abs(sv2[1]):.1e})', abs(sv2[0]) > 1e-4 and abs(sv2[1]) < 1e-6)
    R = R/R[0, 0]
    ev = np.linalg.eigvals(R)
    n = 2*l+1
    # twisted TPG (DGZ 1996 sec 4.3), translated to book units, s_i = (-1)^i
    for label, signs in (('twisted TPG, s_i=(-1)^i', [(-1)**i for i in range(1, l+1)]),
                         ('untwisted rule, all s_i=-1', [-1]*l),
                         ('twisted TPG with s_i=-(-1)^i', [-(-1)**i for i in range(1, l+1)])):
        want = []; rho = 1
        for j in range(0, l+1):
            if j > 0:
                rho = rho*ang(j/2, z, q, signs[j-1])
            want.append((rho, comb(n, l-j)))
        ok = match(ev, want)
        print(f'      l={l}: eigenvalues match {label}: {ok}')
