"""A1: full Dirichlet, prp-no-dirichlet (ch. 30).  Independent exact check.

Compatibility (eq-boundary-compatibility) at x = 0 with phi(0) = phi_0 fixed and d_x phi free:
  linear sigma:  K rho(a_t(lam)) - rho'(a_t(1/lam)) K = 0,  rho' = rho (antipode) or rho^* (transposition),
  a_t(lam) = (beta/2) d_x phi . H + (m/2) sum_j c_j w_j (lam E_j - lam^-1 E_-j),  w_j = exp(beta alpha_j.phi_0/2).
Unknown: K (all N^2 entries).  Conditions: K H_a = rho'(H_a) K for every Cartan basis element (d_x phi free),
and the step part.  Generic lam = 3/7 and generic boundary values w_j with prod (w_j)^{n_j} = 1 (a real phi_0).
Expect: only K = 0 for a_2 (V, adjoint), a_3 (V), d_4 (V), b_2 (V), c_2 (V); a nonzero K for a_1 (V).
Antilinear sigma at imaginary coupling: K a_t(lam) + sigma(a_t(1/lam)) K = 0 with complex d_x phi free,
sigma = dagger or X -> -conj(X), a_t = (i beta/2) d_x phi.H + ...;  real-linear system in Re K, Im K.
"""
import sympy as sp
from itertools import product

lam = sp.Rational(3, 7)


def E(N, i, j):
    M = sp.zeros(N, N)
    M[i, j] = 1
    return M


def comm(A, B):
    return A*B - B*A


# ---------------------------------------------------------------- representations
def sl_vector(n):
    """sl(n+1) defining rep: E_j = e_{j,j+1}, j = 1..n (0-based j-1), E_0 = e_{N,1} = E_{-theta}."""
    N = n + 1
    Ep = [E(N, j, j+1) for j in range(n)]
    Em = [M.T for M in Ep]
    H = [comm(a, b) for a, b in zip(Ep, Em)]
    E0p = E(N, N-1, 0)          # E_{-theta}
    E0m = E0p.T
    return N, Ep, Em, H, E0p, E0m, [1]*(n+1)


def adjoint_of(rep):
    """Adjoint representation built from a faithful rep: basis = matrices spanning the algebra."""
    N, Ep, Em, H, E0p, E0m, marks = rep
    # generate the Lie algebra spanned by Ep, Em
    basis = []

    def add(M):
        vecs = [list(b) for b in basis] + [list(M)]
        if sp.Matrix(vecs).rank() > len(basis):
            basis.append(M)
            return True
        return False
    frontier = Ep + Em
    for M in frontier:
        add(M)
    changed = True
    while changed:
        changed = False
        for A in list(basis):
            for B in list(basis):
                if add(comm(A, B)):
                    changed = True
    d = len(basis)
    Bm = sp.Matrix([list(b) for b in basis]).T     # columns = basis vectors (flattened)

    def ad(X):
        cols = []
        for b in basis:
            c = Bm.solve_least_squares(sp.Matrix(list(comm(X, b))))
            cols.append(c)
        return sp.Matrix.hstack(*cols)
    return d, [ad(M) for M in Ep], [ad(M) for M in Em], [ad(M) for M in H], ad(E0p), ad(E0m), marks


def so_vector(n2):
    """so(2n) or so(2n+1) vector rep, form J = antidiag(1); index i <-> e_i, N-1-i <-> e_{-i}."""
    N = n2
    n = N // 2
    def Eij(i, j):          # E_{eps_i - eps_j}
        return E(N, i, j) - E(N, N-1-j, N-1-i)
    def Fij(i, j):          # E_{eps_i + eps_j}
        return E(N, i, N-1-j) - E(N, j, N-1-i)
    Ep = [Eij(i, i+1) for i in range(n-1)]
    if N % 2 == 0:
        Ep.append(Fij(n-2, n-1))                 # d_n: eps_{n-1} + eps_n
        marks = [1, 1] + [2]*(n-3) + [1, 1]       # d_n Kac labels (node 0, 1, ..., n)
    else:
        Ep.append(E(N, n-1, n) - E(N, n, N-n))    # b_n short root eps_n (middle index n)
        marks = [1, 1] + [2]*(n-1)                # b_n: n_0 = n_1 = 1, others 2
    Em = [M.T for M in Ep]
    H = [comm(a, b) for a, b in zip(Ep, Em)]
    E0p = Fij(0, 1).T                              # E_{-theta}, theta = eps_1 + eps_2
    return N, Ep, Em, H, E0p, E0p.T, marks


def sp4_vector():
    """c_2 = sp(4): basis e1, e2, e-2, e-1 with symplectic form; roots eps1-eps2 (short), 2 eps2 (long)."""
    N = 4
    Ep = [E(N, 0, 1) - E(N, 2, 3), E(N, 1, 2)]
    Em = [M.T for M in Ep]
    H = [comm(a, b) for a, b in zip(Ep, Em)]
    E0p = E(N, 3, 0)                               # E_{-theta}, theta = 2 eps1
    return N, Ep, Em, H, E0p, E0p.T, [1, 2, 1]


# ---------------------------------------------------------------- solving
def solve_linear(rep, prime, w):
    N, Ep, Em, H, E0p, E0m, marks = rep
    Pp = [E0p] + Ep            # node 0 first
    Pm = [E0m] + Em
    if prime == 'antipode':
        rp = lambda X: X
    else:
        rp = lambda X: -X.T
    K = sp.Matrix(N, N, sp.symbols(f'k0:{N*N}'))
    eqs = []
    for Ha in H:
        eqs += list(K*Ha - rp(Ha)*K)
    step = sp.zeros(N, N)
    for j in range(len(Pp)):
        step += w[j]*(lam*K*Pp[j] - K*Pm[j]/lam - rp(Pp[j])*K/lam + lam*rp(Pm[j])*K)
    eqs += list(step)
    A, _ = sp.linear_eq_to_matrix(eqs, list(K))
    return N*N - A.rank()


def solve_antilinear(rep, kind, w):
    """Imaginary coupling, sigma antilinear; w_j complex (boundary values), lam real."""
    N, Ep, Em, H, E0p, E0m, marks = rep
    Pp = [E0p] + Ep
    Pm = [E0m] + Em
    sig = (lambda X: X.H) if kind == 'dagger' else (lambda X: -X.conjugate())
    a = sp.symbols(f'a0:{N*N}', real=True)
    b = sp.symbols(f'b0:{N*N}', real=True)
    K = sp.Matrix(N, N, [a[i] + sp.I*b[i] for i in range(N*N)])
    eqs = []
    # Cartan part: (i beta/2)[K (z.H) + sigma(z.H) ...];  z = u + i v, u, v real and free.
    # sigma((i/2) z.H) = sigma applied to the matrix (i/2)(u+iv).H.
    for Ha in H:
        for zc in (1, sp.I):
            X = sp.I*zc*Ha/2
            eqs += list(K*X + sig(X)*K)
    step = sp.zeros(N, N)
    for j in range(len(Pp)):
        X = w[j]*(lam*Pp[j] - Pm[j]/lam)          # a_t(lam) step term of node j
        Y = w[j]*(Pp[j]/lam - lam*Pm[j])          # a_t(1/lam) step term of node j
        step += K*X + sig(Y)*K
    eqs += list(step)
    real_eqs = []
    for e in eqs:
        e = sp.expand(e)
        real_eqs += [sp.re(e), sp.im(e)]
    A, _ = sp.linear_eq_to_matrix(real_eqs, list(a) + list(b))
    return 2*N*N - A.rank()


def generic_w(marks):
    """Generic boundary values w_j > 0 with prod w_j^{n_j} = 1 (a real phi_0); c_j absorbed."""
    vals = [sp.Rational(2, 1), sp.Rational(3, 5), sp.Rational(7, 4), sp.Rational(5, 9), sp.Rational(11, 6),
            sp.Rational(4, 13)]
    w = vals[1:len(marks)]
    prod = sp.Integer(1)
    for wj, nj in zip(w, marks[1:]):
        prod *= wj**nj
    w0 = 1/prod                     # n_0 = 1
    return [w0] + w


def generic_wc(marks):
    vals = [sp.Rational(2, 3) + sp.I/2, sp.Rational(3, 5) - sp.I, sp.Rational(-7, 4) + sp.I/3, sp.Rational(5, 9) + 2*sp.I]
    w = vals[:len(marks)-1]
    prod = sp.Integer(1)
    for wj, nj in zip(w, marks[1:]):
        prod *= wj**nj
    return [sp.simplify(1/prod)] + w


if __name__ == '__main__':
    reps = {
        'a1 V': sl_vector(1),
        'a2 V': sl_vector(2),
        'a3 V': sl_vector(3),
        'd4 V': so_vector(8),
        'b2 V': so_vector(5),
        'c2 V': sp4_vector(),
    }
    reps['a2 adj'] = adjoint_of(reps['a2 V'])
    for name, rep in reps.items():
        w = generic_w(rep[6])
        dims = [solve_linear(rep, p, w) for p in ('antipode', 'transpose')]
        print(f'{name:8s} dim of solution space for K: antipode {dims[0]}, transposition {dims[1]}')
    for name in ('a1 V', 'a2 V', 'a3 V'):
        rep = reps[name]
        wc = generic_wc(rep[6])
        dims = [solve_antilinear(rep, k, wc) for k in ('dagger', 'minusconj')]
        print(f'{name:8s} antilinear sigma, imaginary coupling: dagger {dims[0]}, -conj {dims[1]} (real dimension)')
