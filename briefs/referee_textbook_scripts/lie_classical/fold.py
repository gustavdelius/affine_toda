"""Ch. 6 tbl-foldings, ch. 19 table, appendix tbl-app-twisted: fold simply-laced affine diagrams.
Folded roots alpha_O = orbit average, labels |O| n_O (book sec-affine-folding)."""
import numpy as np, itertools, math
from functools import reduce
from lie import *

def affine(typ, n):
    R = RootSystem(typ, n)
    A = np.vstack([-R.theta, R.al]); nj = np.concatenate([[1], R.marks])
    return R, A, nj

def cartan(A):
    G = A @ A.T
    return np.array([[2*G[i, j]/G[i, i] for j in range(len(A))] for i in range(len(A))])

def automorphisms(C):
    """all permutations p with C[p][:,p] == C (backtracking)."""
    N = len(C); out = []
    def bt(p):
        k = len(p)
        if k == N: out.append(tuple(p)); return
        for c in range(N):
            if c in p: continue
            if all(abs(C[k, j] - C[c, p[j]]) < 1e-9 and abs(C[j, k] - C[p[j], c]) < 1e-9 for j in range(k)) and abs(C[k,k]-C[c,c])<1e-9:
                bt(p + [c])
    bt([]); return out

def orbits_of(p):
    seen = set(); orbs = []
    for i in range(len(p)):
        if i in seen: continue
        o = [i]; j = p[i]
        while j != i: o.append(j); j = p[j]
        seen |= set(o); orbs.append(sorted(o))
    return orbs

def isomorphic(C1, C2):
    if len(C1) != len(C2): return False
    for p in itertools.permutations(range(len(C1))):
        p = list(p)
        if np.allclose(C1[np.ix_(p, p)], C2): return True
    return False

def intvec(v):
    v = np.asarray(v, float); v = v/np.min(np.abs(v))
    for k in range(1, 13):
        if np.allclose(v*k, np.round(v*k), atol=1e-9):
            w = np.round(v*k).astype(int); g = reduce(math.gcd, w); return w//g
    raise ValueError(v)

def untwisted_cartan(typ, n):
    R, A, nj = affine(typ, n); return cartan(A), A, nj

def a2n2(n):   # a_{2n}^{(2)}: chain, lengths 1/2, 1,...,1, 2 (book normalisation, longest 2)
    if n == 1: return np.array([[2., -4.], [-1., 2.]])
    L = [0.5] + [1.0]*(n-1) + [2.0]
    G = np.diag(L)
    for i in range(n):
        G[i, i+1] = G[i+1, i] = -min(L[i], L[i+1]) if L[i] != L[i+1] else -L[i]/2
    return np.array([[2*G[i, j]/G[i, i] for j in range(n+1)] for i in range(n+1)])

def expected(name, n=None):
    if name == 'a2n-1(2)': return untwisted_cartan('B', n)[0].T
    if name == 'dn+1(2)':  return untwisted_cartan('C', n)[0].T if n >= 2 else None
    if name == 'e6(2)':    return untwisted_cartan('F', 4)[0].T
    if name == 'd4(3)':    return untwisted_cartan('G', 2)[0].T
    if name == 'a2n(2)':   return a2n2(n)
    if name == 'bn(1)':    return untwisted_cartan('B', n)[0]
    if name == 'cn(1)':    return untwisted_cartan('C', n)[0]
    if name == 'f4(1)':    return untwisted_cartan('F', 4)[0]
    if name == 'g2(1)':    return untwisted_cartan('G', 2)[0]

def fold(typ, N, want_orbits, fixes0):
    """pick an automorphism with the given sorted orbit-size pattern; fixes0: whether alpha_0 is fixed."""
    R, A, nj = affine(typ, N); C = cartan(A)
    for p in automorphisms(C):
        orbs = orbits_of(p)
        if sorted(len(o) for o in orbs) == sorted(want_orbits) and ((p[0] == 0) == fixes0):
            # linear map: T alpha_j = alpha_p(j); check consistency (orthogonal)
            T = np.linalg.solve(A[1:], A[list(p)][1:]).T
            assert np.allclose(T @ A.T, A[list(p)].T) and np.allclose(T @ T.T, np.eye(N))
            aO = np.array([A[o].mean(0) for o in orbs]); lab = np.array([len(o)*nj[o[0]] for o in orbs])
            assert np.allclose(lab @ aO, 0)
            return dict(perm=p, orbs=orbs, aO=aO, lab=lab, R=R, A=A, nj=nj)
    raise RuntimeError('no automorphism')

def masses(aO, lab):
    M2 = sum(l*np.outer(a, a) for l, a in zip(lab, aO))
    ev = np.linalg.eigvalsh(M2); ev = ev[ev > 1e-9]
    return np.sort(np.sqrt(ev))

cases = [  # (folded name, n, parent typ, parent N, orbit pattern, fixes alpha_0, book h, book hv, kind k)
 ('bn(1)', 3, 'D', 4, [1, 1, 1, 2], True, 6, 5, 1), ('bn(1)', 4, 'D', 5, [1, 1, 1, 1, 2], True, 8, 7, 1),
 ('cn(1)', 2, 'A', 3, [1, 1, 2], True, 4, 3, 1), ('cn(1)', 3, 'A', 5, [1, 1, 2, 2], True, 6, 4, 1),
 ('f4(1)', 4, 'E', 6, [1, 1, 1, 2, 2], True, 12, 9, 1), ('g2(1)', 2, 'D', 4, [1, 1, 3], True, 6, 4, 1),
 ('a2n-1(2)', 2, 'D', 4, [1, 2, 2], False, 3, 4, 2), ('a2n-1(2)', 3, 'D', 6, [1, 2, 2, 2], False, 5, 6, 2),
 ('a2n-1(2)', 4, 'D', 8, [1, 2, 2, 2, 2], False, 7, 8, 2),
 ('dn+1(2)', 2, 'D', 4, [1, 2, 2], False, 3, 4, 2), ('dn+1(2)', 3, 'D', 5, [1, 1, 2, 2], False, 4, 6, 2),
 ('dn+1(2)', 4, 'D', 6, [1, 1, 1, 2, 2], False, 5, 8, 2),
 ('a2n(2)', 1, 'D', 4, [1, 4], False, 3, 3, 2), ('a2n(2)', 2, 'D', 6, [1, 2, 4], False, 5, 5, 2),
 ('a2n(2)', 3, 'D', 8, [1, 2, 2, 4], False, 7, 7, 2),
 ('e6(2)', 4, 'E', 7, [1, 1, 2, 2, 2], False, 9, 12, 2), ('d4(3)', 2, 'E', 6, [1, 3, 3], False, 4, 6, 3),
]
for name, n, typ, N, pat, f0, hb, hvb, k in cases:
    F = fold(typ, N, pat, f0)
    Cf = np.array([[2*F['aO'][i] @ F['aO'][j]/(F['aO'][i] @ F['aO'][i]) for j in range(len(pat))] for i in range(len(pat))])
    E = expected(name, n)
    iso = isomorphic(Cf, E); isoT = isomorphic(Cf, E.T)
    L2 = np.array([a @ a for a in F['aO']])
    g = reduce(math.gcd, [int(x) for x in F['lab']])
    lab = F['lab']//g                                    # coprime labels
    dual_book = lab*L2/2                                 # book formula n_j alpha_j^2/2 (longest root 2)
    u, s, vt = np.linalg.svd(Cf.T); dual_kac = intvec(vt[-1])   # integer left null vector (Kac's a_j^vee)
    # sign/orientation: make positive
    dual_kac = np.abs(dual_kac)
    ms = masses(F['aO'], F['lab']); ms_std = masses(F['aO'], lab)
    print(f"{name:9s} n={n} <- {typ.lower()}_{N}^(1): orbits {F['orbs']}, |alpha_O|^2={np.round(L2,4).tolist()}, "
          f"Cartan==expected:{iso} (==transpose:{isoT}); labels |O|n_O={F['lab'].tolist()} gcd={g}; "
          f"h=sum coprime labels={lab.sum()} (book {hb}); hv(Kac,integer)={dual_kac.sum()} (book {hvb}); "
          f"hv(book formula n alpha^2/2)={dual_book.sum():g}; k*that={k*dual_book.sum():g}")
    print(f"           masses (labels |O|n_O): {np.round(ms,6).tolist()}  ratios {np.round(ms/ms[0],6).tolist()};"
          f" with coprime labels: {np.round(ms_std,6).tolist()}")
