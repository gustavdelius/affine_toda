"""A2-A3: half-coweight boundaries (prp-soliton-preserving-half, ch. 30) and missing single-soliton classes.

Written from the book's definitions (ch. 10 sec-topological-charges, eq-an-soliton, eq-n-soliton,
eq-interaction-coefficient; ch. 30 eq-boundary-lax-at, display before prp-soliton-preserving).

1. Exact rational enumeration of the charges of single solitons of a_n (theta_j = 2 pi j a/h + Im xi; the
   charge -(1/2pi) sum_j alpha_j [theta_j] changes only where some theta_j crosses pi), reduced mod 2 Lambda_W,
   union over species; count the nontrivial classes without a single-soliton charge (expect 0, 2, 5, 21).
2. Compatibility at imaginary coupling with sigma = dagger in the defining rep of a_2, a_3: the general K
   (all N^2 entries, real-linear system) for Re phi(0) = (pi/beta) nu (nu a weight) and for generic Re phi(0).
3. Classical reflection equation for constant diagonal K, sigma = dagger, real spectral parameters, a_2 and a_3.
4. Half-soliton: static soliton with Re xi = 0 obeys phi(x) + conj phi(-x) = 2 Re phi(0) and
   (beta/pi) Re phi(0) = its charge.  Three-soliton image with Re xi_3 = -ln A_ab(theta),
   xi_2 = -conj(xi_1) - ln A_aa(2 theta) - ln A_ab(theta) satisfies phi(x,t) + conj phi(-x,t) = const.
Units beta = m = 1 where irrelevant.
"""
from fractions import Fraction as Fr
from itertools import product, combinations
from math import gcd
import numpy as np
import sympy as sp


# ---------------------------------------------------------------- 1. charges
def frac_red(x):
    """Reduce a rational x to (-1/2, 1/2]."""
    y = x - (x.numerator // x.denominator)          # in [0,1)
    return y - 1 if y > Fr(1, 2) else y


def species_charges(n, a):
    h = n + 1
    # breakpoints r where (j a/h + r) = 1/2 mod 1
    cuts = sorted({(Fr(1, 2) - Fr(j*a, h)) % 1 for j in range(h)})
    pts = [(cuts[i] + (cuts[i+1] if i+1 < len(cuts) else cuts[0] + 1))/2 for i in range(len(cuts))]
    charges = set()
    for r in pts:
        t = [frac_red(Fr(j*a, h) + r) for j in range(h)]     # [theta_j]/2pi, j = 0..n
        # charge = -sum_j alpha_j t_j, alpha_j = e_j - e_{j+1} (j=1..n), alpha_0 = e_h - e_1  (e_0 == e_h)
        c = [Fr(0)]*h                                         # components on e_1..e_h (index 0..h-1)
        for j in range(h):
            # alpha_j = e_j - e_{j+1} with e_0 = e_h: positions (j-1 mod h) and (j mod h) in 0-based e_1..e_h
            c[(j - 1) % h] -= t[j]
            c[j % h] += t[j]
        shift = c[0] - (c[0].numerator // c[0].denominator)   # common fractional part (= -a/h mod 1)
        v = tuple(int(x - shift) for x in c)
        assert all((x - shift).denominator == 1 for x in c)
        # normalise modulo the all-ones vector to a 0/1 vector
        mn = min(v)
        v = tuple(x - mn for x in v)
        assert set(v) <= {0, 1} and sum(v) in (a, h - a) or True
        charges.add(v)
    return charges


def missing_classes(n):
    h = n + 1
    reached = set()
    nch = {}
    for a in range(1, h):
        ch = species_charges(n, a)
        nch[a] = len(ch)
        for v in ch:
            assert sorted(set(v)) == [0, 1] and sum(v) == a, (n, a, v)
            par = tuple(x % 2 for x in v)
            comp = tuple(1 - x for x in par)
            reached.add(min(par, comp))
    classes = set()
    for S in product((0, 1), repeat=h):
        if 0 < sum(S) < h:
            comp = tuple(1 - x for x in S)
            classes.add(min(S, comp))
    missing = classes - reached
    return len(classes), sorted(missing), nch


for n in (1, 2, 3, 4, 5, 6):
    tot, miss, nch = missing_classes(n)
    print(f'a_{n}: nontrivial classes {tot}, missing {len(miss)}; charges per species {nch}')
    if n == 3:
        print('   a_3 missing classes (0/1 vectors on e_1..e_4, up to complement):', miss)

# ---------------------------------------------------------------- 2. compatibility, sigma = dagger
def a_n_rep(n):
    N = n + 1
    def E(i, j):
        M = sp.zeros(N, N); M[i, j] = 1; return M
    Ep = [E(N-1, 0)] + [E(j, j+1) for j in range(n)]       # node 0 first: E_{alpha_0} = e_{N,1}
    Em = [M.T for M in Ep]
    # simple roots in R^N (e_1..e_N): alpha_0 = e_N - e_1, alpha_j = e_j - e_{j+1}
    roots = []
    for j in range(N):
        v = [0]*N
        if j == 0:
            v[N-1] += 1; v[0] -= 1
        else:
            v[j-1] += 1; v[j] -= 1
        roots.append(v)
    H = [E(i, i) for i in range(N)]                        # weights e_i: H = diag basis (gl, fine)
    return N, Ep, Em, H, roots


def compat_dim(n, rephi, lam=sp.Rational(5, 3), imphi=None, check_K=None):
    """Real dimension of {K : K a_t(lam) + a_t(1/lam)^dagger K = 0} with d_x Re phi free, d_x Im phi = 0.
    a_t(lam) = (i/2) d_x phi . H + (1/2) sum_j w_j (lam E_j - E_-j/lam), w_j = exp(i alpha_j . phi/2)."""
    N, Ep, Em, H, roots = a_n_rep(n)
    if imphi is None:
        imphi = [sp.Rational(k+1, 7) for k in range(N)]
    phi = [rephi[i] + sp.I*imphi[i] for i in range(N)]
    w = [sp.exp(sp.I*sum(r[i]*phi[i] for i in range(N))/2) for r in roots]
    a = sp.symbols(f'a0:{N*N}', real=True)
    b = sp.symbols(f'b0:{N*N}', real=True)
    K = sp.Matrix(N, N, [a[i] + sp.I*b[i] for i in range(N*N)]) if check_K is None else check_K
    eqs = []
    for Ha in H:                       # d_x Re phi free (traceless combos suffice; gl basis is harmless here)
        X = sp.I*Ha/2
        eqs += list(K*X + X.H*K)
    step = sp.zeros(N, N)
    for j in range(N):
        X = w[j]*(lam*Ep[j] - Em[j]/lam)/2
        Y = w[j]*(Ep[j]/lam - lam*Em[j])/2
        step += K*X + Y.H*K
    eqs += list(step)
    if check_K is not None:
        return max(abs(complex(sp.N(e, 30))) for e in eqs)
    real_eqs = []
    for e in eqs:
        e = sp.expand(e)
        real_eqs += [sp.re(e), sp.im(e)]
    A, _ = sp.linear_eq_to_matrix(real_eqs, list(a) + list(b))
    return 2*N*N - A.rank()


pi = sp.pi
for n in (2, 3):
    N = n + 1
    # weights nu = e_1, e_1 + e_2 (beta = 1): Re phi(0) = pi nu (+ 2 pi Lambda)
    for nu in ([1] + [0]*n, [1, 1] + [0]*(n-1), [2] + [0]*n):
        rephi = [pi*x for x in nu]
        d = compat_dim(n, rephi)
        Kexp = sp.diag(*[sp.exp(sp.I*pi*x) for x in nu])
        res = compat_dim(n, rephi, check_K=Kexp)
        print(f'a_{n}: Re phi(0) = pi*{nu}: solution space real dim {d}; residual of K = exp(i pi nu.H): {res:.1e}')
    rephi = [sp.Rational(3, 10)*pi, sp.Rational(-1, 7)*pi] + [sp.Rational(2, 9)*pi]*(N-2)
    print(f'a_{n}: generic Re phi(0): solution space real dim {compat_dim(n, rephi)}')


# ---------------------------------------------------------------- 3. reflection equation, diagonal K
def r_matrix(n, z):
    N = n + 1
    h = N
    # orthonormal Cartan basis of sl(N) in the defining rep
    Hs = []
    for k in range(1, N):
        d = np.zeros(N); d[:k] = 1; d[k] = -k
        Hs.append(np.diag(d/np.sqrt(k*(k+1))))
    I = np.eye(N)
    r = (z**h + 1)*sum(np.kron(Ha, Ha) for Ha in Hs)
    for i in range(N):
        for k in range(i+1, N):
            Ea = np.outer(I[i], I[k]); Ema = Ea.T
            ht = k - i
            r = r + 2*(z**ht*np.kron(Ea, Ema) + z**(h-ht)*np.kron(Ema, Ea))
    return r/(4*(z**h - 1))


def pt(M, N, which):
    T = M.reshape(N, N, N, N)
    if which == 1:
        T = T.transpose(2, 1, 0, 3)
    elif which == 2:
        T = T.transpose(0, 3, 2, 1)
    elif which == 12:
        T = T.transpose(2, 3, 0, 1)
    return T.reshape(N*N, N*N)


rng = np.random.default_rng(1)
for n in (2, 3):
    N = n + 1
    worst = 0
    for trial in range(3):
        l1, l2 = rng.uniform(0.3, 2.5, 2)
        K = np.diag(rng.normal(size=N) + 1j*rng.normal(size=N))
        K1 = np.kron(K, np.eye(N)); K2 = np.kron(np.eye(N), K)
        rd = lambda z: pt(r_matrix(n, z), N, 12).conj()          # r^{sigma1 sigma2}, sigma = dagger
        r1 = lambda z: pt(r_matrix(n, z), N, 1).conj()
        r2 = lambda z: pt(r_matrix(n, z), N, 2).conj()
        lhs = (K1@K2@r_matrix(n, l1/l2) + rd(l2/l1)@K1@K2 + K1@r2(l1*l2)@K2 + K2@r1(1/(l1*l2))@K1)
        worst = max(worst, np.abs(lhs).max())
        # identity check r^{t1t2}(z) = -r(1/z)
        worst = max(worst, np.abs(pt(r_matrix(n, l1), N, 12) + r_matrix(n, 1/l1)).max())
    print(f'a_{n}: reflection equation with constant diagonal K, sigma = dagger: residual {worst:.1e}')


# ---------------------------------------------------------------- 4. half-solitons and three-soliton images
def tau_multi(n, sols, x, t):
    """eq-n-soliton for a_n: sols = list of (a, theta, xi); returns tau_0..tau_n."""
    h = n + 1
    om = np.exp(2j*np.pi/h)
    taus = []
    for j in range(h):
        tot = 0
        for S in product((0, 1), repeat=len(sols)):
            term = 1
            idx = [k for k in range(len(sols)) if S[k]]
            for k in idx:
                a, th, xi = sols[k]
                ma = 2*np.sin(np.pi*a/h)
                term *= om**(j*a)*np.exp(ma*(x*np.cosh(th) - t*np.sinh(th)) + xi)
            for k, l in combinations(idx, 2):
                a, b = sols[k][0], sols[l][0]
                d = sols[k][1] - sols[l][1]
                term *= (np.cosh(d) - np.cos(np.pi*(a - b)/h))/(np.cosh(d) - np.cos(np.pi*(a + b)/h))
            tot += term
        taus.append(tot)
    return taus


def field(n, sols, x, t, ref=None):
    """phi = i sum_j alpha_j ln tau_j (beta = 1, alpha^2 = 2), logs continued along x from far left."""
    h = n + 1
    xs = np.linspace(-25, x, 4001)
    taus = np.array([tau_multi(n, sols, xx, t) for xx in xs])     # (len, h)
    logs = np.log(taus[0]) + np.concatenate([[0*taus[0]], np.cumsum(np.log(taus[1:]/taus[:-1]), axis=0)])[-1]
    # roots in R^h
    phi = np.zeros(h, complex)
    for j in range(h):
        v = np.zeros(h)
        if j == 0:
            v[h-1] += 1; v[0] -= 1
        else:
            v[j-1] += 1; v[j] -= 1
        phi += 1j*v*logs[j]
    return phi


worst_half = 0
for n, a, eta in ((2, 1, 0.3), (3, 2, 1.1), (3, 1, -0.7), (4, 2, 2.0)):
    sols = [(a, 0.0, 1j*eta)]
    p0 = field(n, sols, 0.0, 0.0)
    for x in (0.4, 1.3):
        s = field(n, sols, x, 0.0) + field(n, sols, -x, 0.0).conj()
        worst_half = max(worst_half, np.abs(s - 2*p0.real).max())
    pinf = field(n, sols, 25.0, 0.0)
    worst_half = max(worst_half, np.abs(pinf/(2*np.pi) - p0.real/np.pi).max())
print(f'half-soliton symmetry and Re phi(0) = (pi/beta) charge: worst {worst_half:.1e}')

worst3 = 0
for n, a, b, th, xi1, eta3 in ((2, 1, 1, 0.8, 0.4 + 0.3j, 0.5), (3, 1, 2, 0.6, -0.2 + 1.0j, 1.3), (3, 2, 1, 1.1, 0.3 - 0.4j, -0.9)):
    h = n + 1
    A = lambda p, q, d: (np.cosh(d) - np.cos(np.pi*(p - q)/h))/(np.cosh(d) - np.cos(np.pi*(p + q)/h))
    xi3 = -np.log(A(a, b, th)) + 1j*eta3
    xi2 = -np.conj(xi1) - np.log(A(a, a, 2*th)) - np.log(A(a, b, th))
    sols = [(a, th, xi1), (a, -th, xi2), (b, 0.0, xi3)]
    for t in (-1.5, 0.7):
        c0 = field(n, sols, 0.0, t) + field(n, sols, 0.0, t).conj()
        for x in (0.5, 2.0):
            s = field(n, sols, x, t) + field(n, sols, -x, t).conj()
            worst3 = max(worst3, np.abs(s - c0).max())
print(f'three-soliton image phi(x)+conj phi(-x) = const: worst {worst3:.1e}')
