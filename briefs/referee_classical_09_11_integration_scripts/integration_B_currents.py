"""a_2^(1) local conserved currents of spin s, built from the field equation alone.

psi = i beta phi has components c_0, c_1, c_2 (sum zero) in the basis of an_alpha;
d_+ d_- c_k = -(E_k - E_{k-1}),  E_j = exp(alpha_j . psi) = exp(c_j - c_{j+1})   (m = 1).
Jet variables q[k][c] = D^k c_c with D = d_- = d_t - d_x, c = 0, 1 (c_2 = -c_0 - c_1).
Solve d_+ T = D X with T a polynomial in the q of weight s+1 and X = sum_j E_j * (poly of weight s-1).
Then Q = int (T - X) dx is conserved, since d_+T - d_-X = d_t(T-X) + d_x(T+X).
"""
import itertools
import sympy as sp

K = 6
q = [[sp.Symbol(f'q{k}_{c}') for c in range(2)] for k in range(K + 2)]
E = sp.symbols('E0 E1 E2')


def comp(k, c):            # D^k c_c for c = 0, 1, 2
    return q[k][c] if c < 2 else -q[k][0] - q[k][1]


def adq(j):                # alpha_j . D psi = D(c_j - c_{j+1})
    return comp(1, j % 3) - comp(1, (j + 1) % 3)


def D(expr):               # total d_- derivative on the jet
    out = 0
    for k in range(1, K + 1):
        for c in range(2):
            out += sp.diff(expr, q[k][c]) * q[k + 1][c]
    for j in range(3):
        out += sp.diff(expr, E[j]) * E[j] * adq(j)
    return sp.expand(out)


F = [-(E[c] - E[(c - 1) % 3]) for c in range(2)]   # d_+ d_- c_c
DF = [[None] * (K + 1) for _ in range(2)]
for c in range(2):
    DF[c][0] = sp.expand(F[c])
    for k in range(1, K):
        DF[c][k] = D(DF[c][k - 1])


def Dplus_poly(T):         # d_+ of a polynomial in the q only
    out = 0
    for k in range(1, K + 1):
        for c in range(2):
            d = sp.diff(T, q[k][c])
            if d != 0:
                out += d * DF[c][k - 1]
    return sp.expand(out)


def monomials(w):
    """All monomials of weight w in q[k][c] (weight k)."""
    vars_w = [(q[k][c], k) for k in range(1, w + 1) for c in range(2)]
    res = set()

    def rec(start, left, cur):
        if left == 0:
            res.add(sp.Mul(*cur)); return
        for i in range(start, len(vars_w)):
            v, k = vars_w[i]
            if k <= left:
                rec(i, left - k, cur + [v])
    rec(0, w, [])
    return sorted(res, key=sp.default_sort_key)



from sympy.polys.matrices import DomainMatrix
GENS = None


def toQ(v):
    if isinstance(v, sp.Basic):
        return sp.QQ.from_sympy(sp.Rational(v))
    return sp.QQ(v) if not isinstance(v, type(sp.QQ(1))) else v


def pdict(expr):
    if expr == 0:
        return {}
    return dict(sp.Poly(expr, *GENS).terms())


def nullspace_QQ(cols):
    """cols: list of dicts monomial->coeff. Returns (basis of nullspace as lists, keys)."""
    keys = sorted({k for c in cols for k in c})
    idx = {k: i for i, k in enumerate(keys)}
    rows = [[sp.QQ(0)] * len(cols) for _ in keys]
    for j, c in enumerate(cols):
        for k, v in c.items():
            rows[idx[k]][j] = toQ(v)
    if not keys:
        return [[sp.QQ(1) if i == j else sp.QQ(0) for i in range(len(cols))] for j in range(len(cols))]
    M = DomainMatrix(rows, (len(keys), len(cols)), sp.QQ)
    N = M.nullspace()
    return [list(r) for r in N.to_Matrix().tolist()]


def rank_QQ(vecs):
    if not vecs:
        return 0
    M = DomainMatrix([[toQ(x) for x in v] for v in vecs], (len(vecs), len(vecs[0])), sp.QQ)
    return M.rank()


def currents(s):
    global GENS
    GENS = [x for row in q for x in row] + list(E)
    mT = monomials(s + 1)
    mX = [E[j] * m for j in range(3) for m in (monomials(s - 1) if s > 1 else [sp.Integer(1)])]
    cols = [pdict(Dplus_poly(m)) for m in mT] + [pdict(-D(m)) for m in mX]
    ns = nullspace_QQ(cols)
    triv = []
    for g in monomials(s):
        dT = pdict(D(g)); dX = pdict(Dplus_poly(g))
        vec = [dT.get(sp.Poly(m, *GENS).monoms()[0], 0) for m in mT] + \
              [dX.get(sp.Poly(m, *GENS).monoms()[0], 0) for m in mX]
        triv.append(vec)
    rt = rank_QQ(triv)
    rall = rank_QQ(triv + ns)
    pick = None
    for v in ns:
        if rank_QQ(triv + [v]) > rt:
            pick = v; break
    if pick is None:
        return len(ns), rt, rall, None, None
    T = sp.expand(sum(sp.QQ.to_sympy(toQ(c)) * m for c, m in zip(pick[:len(mT)], mT)))
    X = sp.expand(sum(sp.QQ.to_sympy(toQ(c)) * m for c, m in zip(pick[len(mT):], mX)))
    return len(ns), rt, rall, T, X


if __name__ == '__main__':
    import sys
    for s in map(int, sys.argv[1:] or ['1', '2', '3']):
        n, rt, rall, T, X = currents(s)
        print(f'spin {s}: nullspace {n}, trivial rank {rt}, total rank {rall}, nontrivial {rall - rt}')
        if T is not None:
            chk = sp.expand(Dplus_poly(T) - D(X))
            print('   residual', chk, '\n   T =', T, '\n   X =', X, flush=True)
