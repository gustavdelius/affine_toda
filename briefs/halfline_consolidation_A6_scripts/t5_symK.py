"""Exact d_n vector K-matrices (s = +-1, and eps = 0) from the book's coideal, symbolically in w = (-q)^(1/2), y.
q = -w^2, eps_* = w/(1 + w^2) (eps_*^2 = 1/((1-q)(1-q^-1))). Generators as boundary_common.dn_vector."""
import sys
import pickle
import sympy as sp
from sympy.polys.matrices import DomainMatrix

n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
w, y = sp.symbols('w y')
q = -w**2
N = 2*n
idx = lambda i: i - 1 if i > 0 else N + i


def unit(r, c):
    M = sp.zeros(N, N); M[idx(r), idx(c)] = 1; return M


Es = {}
for i in range(1, n):
    Es[i] = unit(i, i + 1) - unit(-(i + 1), -i)
Es[n] = unit(n - 1, -n) - unit(n, -(n - 1))
Es[0] = unit(-2, 1) - unit(-1, 2)


def gens(yy):
    g = {}
    for j, E in Es.items():
        F = E.T
        H = [(E*F - F*E)[k, k] for k in range(N)]
        k = sp.diag(*[q**int(hh) for hh in H])
        g[j] = (yy*E, F/yy, k)
    return g


def solve(eps):
    gin, gout = gens(y), gens(1/y)
    Ks = sp.symbols('K0:%d' % (N*N))
    K = sp.Matrix(N, N, Ks)
    eqs = []
    for j in range(n + 1):
        e, f, k = gin[j]
        bin_ = e + f*k/q + eps[j]*k
        e2, f2, k2 = gout[j]
        bout = e2 + f2*k2/q + eps[j]*k2
        Mq = K*bin_ - bout*K
        eqs += list(Mq)
    A = sp.Matrix([[sp.together(eq).as_numer_denom()[0].coeff(Ks[c]) if False else sp.diff(eq, Ks[c])
                    for c in range(N*N)] for eq in eqs])
    dom = sp.FractionField(sp.QQ, w, y) if False else None
    DM = DomainMatrix.from_Matrix(A).to_field()
    ns = DM.nullspace()
    return ns.to_Matrix()


es = w/(1 + w**2)
out = {}
for name, eps in (('s+', [es]*(n + 1)), ('s-', [-es] + [es]*n), ('0', [0]*(n + 1))):
    NS = solve(eps)
    print(name, 'nullity', NS.shape[0], flush=True)
    v = NS.row(0)
    Kmat = sp.Matrix(N, N, list(v))
    # clear denominators: make entries polynomials in w, y, y^-1
    dens = [sp.fraction(sp.together(x))[1] for x in Kmat if x != 0]
    L = sp.lcm(dens)
    Kmat = (Kmat*L).applyfunc(lambda x: sp.factor(sp.cancel(x)))
    g = sp.gcd([x for x in Kmat if x != 0])
    Kmat = Kmat.applyfunc(lambda x: sp.factor(sp.cancel(x/g)))
    out[name] = Kmat
    sp.pprint(Kmat, wrap_line=False, num_columns=250)
pickle.dump({k: sp.srepr(v) for k, v in out.items()}, open(f'symK_d{n}.pkl', 'wb'))
