"""Soliton reflection matrices (Part VII, chapter 33; sec-soliton-reflection).

Part VI conventions throughout (boundary_common.py): Delta(e)=e(x)1+k(x)e, Delta(f)=f(x)k^-1+1(x)f, coideal
generators b_j = e_j + q^-1 f_j k_j + eps_j k_j, a K-matrix is an intertwiner K b_j = b_j K.

Checks:
 1. XXZ reflection matrices of any spin (sec-xxz-spin-j).  The Delius-Nepomechie closed form hep-th/0204076
    (3.8)-(3.13), with the square root over the first product only, solves their linear equations (3.6)-(3.7)
    for j <= 2 at real and complex eta, but not at the physical imaginary eta for j >= 1 (branches of the square
    roots); their linear equations are the book's intertwining equations: the solutions agree up to a diagonal
    gauge under y = e^u, q = e^eta, eps_1 = e^xi/(2 kappa sinh eta), eps_0 = -e^-xi/(2 kappa sinh eta).
 2. Book spin-j K (physical q, j <= 2) solves the reflection equation with R^(j,j) and with R^(1/2,j)
    and R^(1,2); the same at -q, where K equals K(q) up to gauge at y (half-integer j) or iy (integer j).
 3. Boundary fusion (sec-boundary-fusion): the spin-1 K from the coideal equals the fused
    (1 (x) K(y1)) R(y1 y2) (1 (x) K(y2)) restricted to V_1 inside V(y q^(-1/2)) (x) V(y q^(1/2)).
 4. d_n^(1) vector solitons, n = 4, 5 (sec-dn-vector-k): R-check from the intertwiner = tensor-product-graph
    eigenvalues (eq-tpg-rule) in the homogeneous gradation; the Delius-George existence conditions (math/0208043
    (22)) in their conventions, plus the eps = 0 solution they omit; in book conventions K: V(y) -> V(1/y)
    exists iff eps = 0 (antidiagonal K) or all eps_j^2 = eps_*^2 = 1/((1-q)(1-q^-1)), and solves the reflection
    equation; conjugation by the Cartan torus gives the continuous families.
 5. Delius-Gandenberger scalar factor hep-th/9904002 (3.9) (sec-an-diagonal-k): unitarity and the crossing
    relation (3.8), with the a_n soliton scalar factor F_11 (h-generalization of hep-th/9806003 (2.6)) itself
    checked against bulk unitarity and crossed-transmission unitarity; the CDD factors (3.24) solve (3.22)-(3.23).
 6. Semiclassical limit (sec-boundary-semiclassics): (beta^2/4h) d(-i log A_1)/dtheta -> ln(1 + sin(pi/h)/cosh theta),
    the (+...+) time delay of eq-boundary-time-delay.
Runtime: about a minute.
"""
import itertools
import numpy as np
import mpmath as mp
from boundary_common import (report, q_phys, qnum, sl2_spin, dn_vector, rcheck, k_intertwiner, intertwiners,
                             cop, coideal_gen, dual, kron, prop_residual, nodes_of, cartan_d, check_relations, _mk)

rng = np.random.default_rng(11)
OMEGA = 0.37
q = q_phys(OMEGA)
EPS_STAR = 1/np.sqrt((1 - q)*(1 - 1/q))


def gauge_equiv(A, B):
    """A = D1 B D2 with D1, D2 diagonal  <=>  the elementwise ratio A/B has rank 1 (all entries nonzero)."""
    s = np.linalg.svd(A/B, compute_uv=False)
    return s[1]/s[0]


def rcheck_w(A, B):
    """R-check A (x) B -> B (x) A (Jimbo's equations, Part VI coproduct), normalized to 1 on hw (x) hw, solved only
    for the entries that preserve the k-eigenvalues (weight filter); fast for the d_n vector representation."""
    import scipy.sparse as sps
    dA, dB = A['dim'], B['dim']
    js = nodes_of(A)
    win = list(zip(*[np.round(np.diag(cop(A, B, 'k', j)), 9) for j in js]))
    wout = list(zip(*[np.round(np.diag(cop(B, A, 'k', j)), 9) for j in js]))
    cols = [(r, c) for r in range(dA*dB) for c in range(dA*dB) if wout[r] == win[c]]
    blocks = []
    for j in js:
        for kind in 'ef':
            D, Dp = cop(A, B, kind, j), cop(B, A, kind, j)
            rows, cc, vals = [], [], []
            for i, (r, c) in enumerate(cols):
                for t in np.nonzero(D[c, :])[0]:        # E_rc D
                    rows.append(r*dA*dB + t); cc.append(i); vals.append(D[c, t])
                for t in np.nonzero(Dp[:, r])[0]:       # - D' E_rc
                    rows.append(t*dA*dB + c); cc.append(i); vals.append(-Dp[t, r])
            blocks.append(sps.csr_matrix((vals, (rows, cc)), shape=((dA*dB)**2, len(cols))))
    M = sps.vstack(blocks)
    ev, U = np.linalg.eigh((M.conj().T @ M).toarray())
    assert ev[1] > 1e-8*ev[-1], 'intertwiner not unique'
    R = np.zeros((dA*dB, dA*dB), complex)
    for i, (r, c) in enumerate(cols):
        R[r, c] = U[i, 0]
    i_in = A['hw']*dB + B['hw']; i_out = B['hw']*dA + A['hw']
    return R/R[i_out, i_in]


# ---------------------------------------------------------------- 1. Delius-Nepomechie spin-j K-matrices
def dn_K(j2, u, eta, xi, ka):
    """Delius-Nepomechie (3.8)-(3.13), K^(j)_{mn}, m,n = 1..2j+1, symmetric; square root over the first product."""
    j = j2/2; d = j2 + 1; sh = np.sinh

    def F(l, n):
        return sh(2*u-(n+l)*eta)*sh((2*j-n-l)*eta)/(sh(xi+u+(j+.5-n-l)*eta)*sh(xi+u+(j-.5-n-l)*eta))
    cache = {}

    def J(m, n, k):
        if k == 0:
            return 1.0 if m >= 1 else 0.0
        if k < 0 or m < 1:
            return 0.0
        if (m, n, k) in cache:
            return cache[(m, n, k)]
        if m == 1:      # (3.10)
            val = sum(np.prod([F(l, n) for l in ls]) for ls in itertools.combinations(range(int(round(2*j-1-n)) + 1), k)
                      if all(ls[i+1] - ls[i] >= 2 for i in range(k - 1)))
        else:           # (3.12)-(3.13)
            a = sh(xi+u+(j-n+1.5)*eta)*sh(2*u+eta)/(sh((n-m+1)*eta)*sh(xi-u+(j-m+1.5)*eta))
            b = -sh(xi+u+(j-m+2.5)*eta)*sh(2*u-(n-m)*eta)/(sh((n-m+1)*eta)*sh(xi-u+(j-m+1.5)*eta))
            c = (-sh((m-2)*eta)*sh((2*j-m+3)*eta)/(sh((n-m+2)*eta)*sh((n-m+1)*eta)**2)
                 * sh(2*u+(n-m+2)*eta)*sh(2*u-(n-m+1)*eta)*sh(2*u-(n-m)*eta)
                 / (sh(xi-u+(j-m+2.5)*eta)*sh(xi-u+(j-m+1.5)*eta)))
            val = a*J(m-1, n-1, k) + b*J(m-1, n, k) + c*J(m-2, n, k-1)
        cache[(m, n, k)] = val
        return val
    K = np.zeros((d, d), complex)
    for m in range(1, d + 1):
        for n in range(m, d + 1):
            r = np.prod([sh((2*j-m+1-l)*eta)/sh((n-1-l)*eta) for l in range(n - m)])
            pre = (ka**(n-m)*np.sqrt(r + 0j)*np.prod([sh((n-1-l)*eta)/sh((m-1-l)*eta) for l in range(m - 1)])
                   * np.prod([sh(2*u-l*eta) for l in range(n - m)])
                   * np.prod([sh(xi+u+(l-j+.5)*eta) for l in range(int(round(2*j-n)) + 1)])
                   * np.prod([sh(xi-u-(l-j+.5)*eta) for l in range(m - 1)]))
            K[m-1, n-1] = K[n-1, m-1] = pre*sum(ka**(2*k)*J(m, n, k) for k in range((j2 + m - n)//2 + 1))
    return K


def dn_equations(K, j2, u, eta, xi, ka):
    """Residual vector of Delius-Nepomechie (3.6) and (3.7) (linear in K)."""
    j = j2/2; d = j2 + 1
    om = lambda n: np.sqrt(np.sinh(n*eta)*np.sinh((2*j+1-n)*eta)/np.sinh(eta)**2 + 0j) if 1 <= n <= 2*j else 0
    Kf = lambda m, n: K[m-1, n-1] if (1 <= m <= d and 1 <= n <= d) else 0
    e1 = np.exp(-xi)/(2*ka*np.sinh(eta)); e2 = np.exp(xi)/(2*ka*np.sinh(eta)); E = np.exp
    v = []
    for m in range(1, d + 1):
        for n in range(1, d + 1):
            v.append(E(-u-eta*(j+1.5-m))*om(m-1)*Kf(m-1, n) + E(u-eta*(j+.5-m))*om(m)*Kf(m+1, n)
                     - e1*E(-2*eta*(j+1-m))*Kf(m, n)
                     - E(u-eta*(j+.5-n))*om(n)*Kf(m, n+1) - E(-u-eta*(j+1.5-n))*om(n-1)*Kf(m, n-1)
                     + e1*E(-2*eta*(j+1-n))*Kf(m, n))
            v.append(E(-u+eta*(j+.5-m))*om(m)*Kf(m+1, n) + E(u+eta*(j+1.5-m))*om(m-1)*Kf(m-1, n)
                     + e2*E(2*eta*(j+1-m))*Kf(m, n)
                     - E(u+eta*(j+1.5-n))*om(n-1)*Kf(m, n-1) - E(-u+eta*(j+.5-n))*om(n)*Kf(m, n+1)
                     - e2*E(2*eta*(j+1-n))*Kf(m, n))
    return np.array(v)


u0, xi0, ka0 = 0.37 + 0.1j, 0.52 - 0.2j, 0.8 + 0.3j
def dn_res(K, j2, eta):
    return np.abs(dn_equations(K, j2, u0, eta, xi0, ka0)).max()/np.abs(K).max()


def dn_null(j2, eta):
    """Unique solution of DN (3.6)-(3.7), from the null space of the linear system."""
    d = j2 + 1
    M = np.array([dn_equations(np.eye(d*d)[i].reshape(d, d), j2, u0, eta, xi0, ka0) for i in range(d*d)]).T
    _, sv, Vh = np.linalg.svd(M)
    return Vh[-1].conj().reshape(d, d), sv[-2]/sv[0], sv[-1]/sv[0]


ETA_PHYS = 1j*np.pi*(1 - OMEGA)      # e^eta = -e^{-i pi omega}, the physical q
for label, eta in (('real eta', 0.41), ('complex eta', 0.41 + 0.2j)):
    res = max(dn_res(dn_K(j2, u0, eta, xi0, ka0), j2, eta) for j2 in (1, 2, 3, 4))
    report(f'Delius-Nepomechie (3.8)-(3.13) solve their (3.6)-(3.7), j = 1/2..2, {label} [{res:.0e}]', res < 1e-10)
res_half = dn_res(dn_K(1, u0, ETA_PHYS, xi0, ka0), 1, ETA_PHYS)
res_int = max(dn_res(dn_K(j2, u0, ETA_PHYS, xi0, ka0), j2, ETA_PHYS) for j2 in (2, 3, 4))
report(f'... at the physical (imaginary) eta, principal-branch square roots: j = 1/2 solves [{res_half:.0e}], '
       f'j >= 1 does not [{res_int:.1f}] (the phase caveat of DN footnote 1)', res_half < 1e-10 and res_int > 1e-3)
ok = True
for eta in (0.41, 0.41 + 0.2j, ETA_PHYS):
    qq = np.exp(eta)
    e1b = np.exp(xi0)/(2*ka0*np.sinh(eta)); e0b = -np.exp(-xi0)/(2*ka0*np.sinh(eta))
    for j2 in (1, 2, 3, 4):
        Kl, gap, nul = dn_null(j2, eta)
        null, _ = k_intertwiner(sl2_spin(j2, qq, np.exp(u0)), sl2_spin(j2, qq, np.exp(-u0)), [1/qq, 1/qq], [e0b, e1b])
        ok &= len(null) == 1 and nul < 1e-12 and gap > 1e-6 and gauge_equiv(Kl, null[0]) < 1e-10
        if eta != ETA_PHYS:
            ok &= gauge_equiv(dn_K(j2, u0, eta, xi0, ka0), null[0]) < 1e-10
report('DN (3.6)-(3.7) have a unique solution = book intertwiner up to diagonal gauge (real, complex and physical eta): '
       'y = e^u, q = e^eta, eps_1 = e^xi/(2 kappa sinh eta), eps_0 = -e^-xi/(2 kappa sinh eta), j <= 2', ok)
# the reading of (3.8): with the second product also under the square root the entries m >= 2 are wrong
eta = 0.41
Kwrong = dn_K(2, u0, eta, xi0, ka0)
Kwrong[1, 2] /= np.sqrt(np.sinh(2*eta)/np.sinh(eta)); Kwrong[2, 1] = Kwrong[1, 2]
report('... the square root in DN (3.8) must not extend over the second product (control fails)',
       np.abs(dn_equations(Kwrong, 2, u0, eta, xi0, ka0)).max()/np.abs(Kwrong).max() > 1e-3)
# K(0) is proportional to 1 and kappa -> 0 gives a diagonal K (DN remarks after (3.14))
K0 = dn_K(4, 1e-10, 0.41, xi0, ka0); Kd = dn_K(4, u0, 0.41, xi0, 1e-12)
report('DN: K^(j)(0) is proportional to 1; kappa = 0 gives a diagonal K',
       np.abs(K0 - K0[0, 0]*np.eye(5)).max()/abs(K0[0, 0]) < 1e-8 and np.abs(Kd - np.diag(np.diag(Kd))).max()/np.abs(Kd).max() < 1e-10)


# ---------------------------------------------------------------- 2. spin-j reflection equations, book conventions
def K_spin(j2, y, eps, qq=q):
    null, _ = k_intertwiner(sl2_spin(j2, qq, y), sl2_spin(j2, qq, 1/y), [1/qq, 1/qq], eps)
    assert len(null) == 1
    return null[0]


def re_mixed(Rf, Kf, mu, nu, y1, y2):
    """R_{nu mu}(y1/y2) (1 (x) K_mu(y1)) R_{mu nu}(y1 y2) (1 (x) K_nu(y2))
       = (1 (x) K_nu(y2)) R_{nu mu}(y1 y2) (1 (x) K_mu(y1)) R_{mu nu}(y1/y2)   on V_mu(y1) (x) V_nu(y2)
    (eq-reflection-twisted with mubar = mu, nubar = nu: for a_1^(1) every multiplet is self-conjugate)."""
    Km, Kn = Kf(mu, y1), Kf(nu, y2)
    dm, dn = Km.shape[0], Kn.shape[0]
    L = Rf(nu, mu, y1/y2) @ kron(np.eye(dn), Km) @ Rf(mu, nu, y1*y2) @ kron(np.eye(dm), Kn)
    R = kron(np.eye(dm), Kn) @ Rf(nu, mu, y1*y2) @ kron(np.eye(dn), Km) @ Rf(mu, nu, y1/y2)
    return prop_residual(L, R)


eps_gen = [0.31 - 0.52j, -0.7 + 0.2j]
for qq, eps, lab in ((q, eps_gen, 'physical q'), (-q, [-e for e in eps_gen], '-q with eps -> -eps')):
    Rf = lambda a, b, z, qq=qq: rcheck(sl2_spin(a, qq, z), sl2_spin(b, qq, 1.0))[0]
    Kf = lambda a, y, qq=qq, eps=eps: K_spin(a, y, eps, qq)
    worst = 0
    for mu, nu in ((1, 1), (2, 2), (3, 3), (4, 4), (1, 2), (1, 3), (1, 4), (2, 4)):
        r, _ = re_mixed(Rf, Kf, mu, nu, 1.3 + 0.2j, 0.7 - 0.4j)
        worst = max(worst, r)
    report(f'spin-j K (coideal) solves the reflection equation with R^(j,j) and R^(1/2,j), R^(1,2), {lab} [{worst:.0e}]',
           worst < 1e-9)
# K at -q is K at q up to diagonal gauge, at y (half-integer spin) or i y (integer spin)
ok = True
for j2 in (1, 2, 3, 4):
    yv = 0.9 + 0.3j
    Km_ = K_spin(j2, yv*(1 if j2 % 2 else 1j), eps_gen, -q)
    ok &= gauge_equiv(K_spin(j2, yv, eps_gen, q) + 1e-300, Km_ + 1e-300) < 1e-10
report('spin j: K(theta) at -q equals K at q up to diagonal gauge, with y -> y (j half-integer) or y -> iy (j integer)', ok)

# ---------------------------------------------------------------- 3. boundary fusion, spin 1
y = 0.83 + 0.27j
ok = False
for a in (0.5, -0.5):
    V1 = sl2_spin(2, q, y)
    A_, B_ = sl2_spin(1, q, y*q**a), sl2_spin(1, q, y*q**(-a))
    pairs = [(V1[(k, j)], cop(A_, B_, k, j)) for j in (0, 1) for k in 'efk']
    emb_in, _ = intertwiners(pairs)
    if len(emb_in) != 1:
        continue
    yo = 1/y
    V1o = sl2_spin(2, q, yo)
    Ao, Bo = sl2_spin(1, q, yo*q**a), sl2_spin(1, q, yo*q**(-a))      # = V(1/y2) (x) V(1/y1)
    emb_out, _ = intertwiners([(V1o[(k, j)], cop(Ao, Bo, k, j)) for j in (0, 1) for k in 'efk'])
    y1, y2 = y*q**a, y*q**(-a)
    K1 = lambda yy: K_spin(1, yy, eps_gen)
    M = kron(np.eye(2), K1(y1)) @ rcheck(sl2_spin(1, q, y1*y2), sl2_spin(1, q, 1.0))[0] @ kron(np.eye(2), K1(y2))
    lhs = M @ emb_in[0]
    rhs = emb_out[0] @ K_spin(2, y, eps_gen)
    ok = len(emb_out) == 1 and prop_residual(lhs, rhs)[0] < 1e-10
    shift = a
report(f'boundary fusion: (1(x)K(y1)) R(y1 y2) (1(x)K(y2)) restricted to V_1(y) in V(y q^{shift:+.1f}) (x) V(y q^{-shift:+.1f}) '
       '= spin-1 intertwiner K_1(y)', ok)


# ---------------------------------------------------------------- 4. d_n^(1) vector solitons
def dn_hom(n, qq, z):
    """Vector representation in the homogeneous gradation: only e_0, f_0 carry the spectral parameter."""
    V = dn_vector(n, qq, 1.0)
    g = {k: v for k, v in V.items() if isinstance(k, tuple)}
    g[('e', 0)] = z*g[('e', 0)]; g[('f', 0)] = g[('f', 0)]/z
    return _mk(g, 0)


def chevalley(V, y):
    """Vbar(y): V composed with e <-> f, k -> k^-1 (as an_vector(conj=True)); V given at y = 1."""
    g = {}
    for j in nodes_of(V):
        g[('e', j)] = y*V[('f', j)]; g[('f', j)] = V[('e', j)]/y; g[('k', j)] = np.linalg.inv(V[('k', j)])
    return _mk(g, None)


br = lambda a, z, qq: (1 - z*qq**(2*a))/(z - qq**(2*a))
for n in (4, 5):
    z = 0.71 + 0.33j
    R = rcheck_w(dn_hom(n, q, z), dn_hom(n, q, 1.0))
    ev = np.linalg.eigvals(R)
    dims = {'2w1': n*(2*n + 1) - 1, 'w2': n*(2*n - 1), '0': 1}
    # Part VI coproduct: q -> q^-1 relative to ch. 8 (eq-tpg-rule); rho_{2w1} = 1
    rho = {'2w1': 1, 'w2': br(1, z, 1/q), '0': br(1, z, 1/q)*br(n - 1, z, 1/q)}
    ok = all(np.sum(np.abs(ev - rho[k]) < 1e-8) == dims[k] for k in dims)
    report(f'd_{n}: R-check(z) (homogeneous) has eigenvalues 1, <1>, <1><{n-1}> on 2w1, w2, 0 (eq-tpg-rule, q -> 1/q, '
           f'C = 4n, 4n-4, 0)', ok)

for n in (4, 5):
    N = 2*n; h = 2*n - 2
    nodes = range(n + 1)
    # Delius-George (8), Qhat_j = q^(h_j/2)(x_j^+ + x_j^-) + eps_j (q^(h_j) - 1), in their representation (12)-(14)
    # (x^+ and x^- with entries +-1, i.e. x^+ = E_j, x^- = F_j of dn_vector): with k = q^h this is
    # q^(1/2) E_j + q^(-3/2) F_j k_j + eps_j (k_j - 1); homogeneous gradation; K: V(eta x) -> V(eta/x)
    def dg_exists(eta, eps, x=0.77 + 0.2j):
        Vin, Vout = dn_hom(n, q, eta*x), dn_hom(n, q, eta/x)
        gen = lambda V, j: (np.sqrt(q)*V[('e', j)] + q**-1.5*V[('f', j)] @ V[('k', j)]
                            + eps[j]*(V[('k', j)] - np.eye(N)))
        null, _ = intertwiners([(gen(Vin, j), gen(Vout, j)) for j in nodes])
        return len(null)
    edg = 1j/(np.sqrt(q) - 1/np.sqrt(q))
    signs = rng.choice([-1, 1], size=n + 1)
    ok = (dg_exists((-1)**n*q**(1 - n), edg*np.ones(n + 1)) == 1 and dg_exists(-(-1)**n*q**(1 - n), edg*signs) == 1
          and dg_exists(1.0, edg*np.ones(n + 1)) == 0 and dg_exists(q**(2 - 2*n), edg*np.ones(n + 1)) == 0
          and dg_exists((-1)**n*q**(1 - n), 1.3*edg*np.ones(n + 1)) == 0)
    report(f'd_{n}: Delius-George (22): K: V(eta x) -> V(eta/x) exists for eta = +-(-1)^n q^(1-n), '
           'eps_j = +-i/(q^1/2 - q^-1/2) (any signs), not for other eta or |eps_j|', ok)
    report(f'd_{n}: ... and also for eps = 0 at the same eta, a solution that (22) omits',
           dg_exists((-1)**n*q**(1 - n), np.zeros(n + 1)) == 1)
    report(f'd_{n}: ... eps_DG^2 = eps_*^2 = 1/((1-q)(1-q^-1)), the a_n value',
           abs(edg**2 - EPS_STAR**2) < 1e-12)

    # book conventions, principal gradation: the antisoliton multiplet, the intertwiner, the reflection equation
    V1 = dn_vector(n, q, 1.0)
    yy = 0.8 + 0.1j
    iso = intertwiners([(dual(dn_vector(n, q, yy))[(k, j)], chevalley(V1, -yy/q)[(k, j)]) for j in nodes for k in 'efk'])[0]
    iso2 = intertwiners([(chevalley(V1, yy)[(k, j)], dn_vector(n, q, yy)[(k, j)]) for j in nodes for k in 'efk'])[0]
    report(f'd_{n}: V(z)* = Vbar(-z/q) (antisoliton at the same rapidity) and Vbar(y) = V(y): the vector multiplet '
           'is self-conjugate at equal y', len(iso) == 1 and len(iso2) == 1)

    def kb(yv, eps):
        return k_intertwiner(dn_vector(n, q, yv), dn_vector(n, q, 1/yv), [1/q]*(n + 1), eps)[0]
    cases = {'eps = 0': (np.zeros(n + 1), 1), 'all eps_j = eps_*': (EPS_STAR*np.ones(n + 1), 1),
             'random signs, eps_j = +-eps_*': (EPS_STAR*rng.choice([-1, 1], size=n + 1), 1),
             'generic eps': (rng.normal(size=n + 1) + 1j*rng.normal(size=n + 1), 0),
             '|eps_j| = 1': (np.ones(n + 1), 0), 'i eps_*': (1j*EPS_STAR*np.ones(n + 1), 0),
             'one eps_j = 0': (np.r_[0, EPS_STAR*np.ones(n)], 0)}
    ok = all(len(kb(0.9 + 0.4j, e)) == want for e, want in cases.values())
    report(f'd_{n}: K: V(y) -> V(1/y) (book coideal, c_j = q^-1) exists (uniquely) iff eps = 0 or all eps_j^2 = eps_*^2', ok)
    K0 = kb(0.9 + 0.4j, np.zeros(n + 1))[0]
    anti = K0[np.arange(N), N - 1 - np.arange(N)]
    off = K0.copy(); off[np.arange(N), N - 1 - np.arange(N)] = 0
    report(f'd_{n}: eps = 0 gives an antidiagonal K: v_i -> +-v_(-i) with equal moduli (the analogue of K = 1 for a_n)',
           np.abs(off).max() < 1e-12*np.abs(anti).max() and np.ptp(np.abs(anti)) < 1e-12*np.abs(anti).max())
    cache = {}

    def Rv(z):
        key = complex(z)
        if key not in cache:
            cache[key] = rcheck_w(dn_vector(n, q, z), dn_vector(n, q, 1.0))
        return cache[key]
    worst = 0
    for e in (np.zeros(n + 1), EPS_STAR*np.ones(n + 1), EPS_STAR*np.r_[-1, np.ones(n)]):
        Kf = lambda yv, e=e: kb(yv, e)[0]
        y1, y2 = 1.3 + 0.2j, 0.7 - 0.4j
        L = Rv(y1/y2) @ kron(np.eye(N), Kf(y1)) @ Rv(y1*y2) @ kron(np.eye(N), Kf(y2))
        Rr = kron(np.eye(N), Kf(y2)) @ Rv(y1*y2) @ kron(np.eye(N), Kf(y1)) @ Rv(y1/y2)
        worst = max(worst, prop_residual(L, Rr)[0])
    report(f'd_{n}: these K solve the reflection equation (braid form, V self-conjugate) [{worst:.0e}]', worst < 1e-9)
    # Cartan-torus conjugation: Sigma (diagonal, commuting with R (x) ) gives a family of solutions
    # a torus element exp(i alpha.H): phases by the weights eps_1..eps_n, -eps_n..-eps_1 of the basis vectors
    wts = np.array(list(range(1, n + 1)) + [-k for k in range(n, 0, -1)], float)
    Sig = np.diag(np.exp(0.37j*wts))
    comm = np.abs(kron(Sig, Sig) @ Rv(0.6 + 0.2j) - Rv(0.6 + 0.2j) @ kron(Sig, Sig)).max()
    Ks = lambda yv: Sig @ kb(yv, EPS_STAR*np.ones(n + 1))[0] @ np.linalg.inv(Sig)
    L = Rv(y1/y2) @ kron(np.eye(N), Ks(y1)) @ Rv(y1*y2) @ kron(np.eye(N), Ks(y2))
    Rr = kron(np.eye(N), Ks(y2)) @ Rv(y1*y2) @ kron(np.eye(N), Ks(y1)) @ Rv(y1/y2)
    report(f'd_{n}: Sigma K Sigma^-1 with Sigma in the Cartan torus also solves the RE (the continuous families) '
           f'[{comm:.0e}, {prop_residual(L, Rr)[0]:.0e}]', comm < 1e-10 and prop_residual(L, Rr)[0] < 1e-9)

# ---------------------------------------------------------------- 5. Delius-Gandenberger scalar factor, CDD factors
mp.mp.dps = 30


def blk(a, mu, h, lam):
    """DG (3.4): (a)_I = [a]/[-a], [a] = sin(pi (mu + a)/(h lam))."""
    return mp.sin(mp.pi*(mu + a)/(h*lam))/mp.sin(mp.pi*(mu - a)/(h*lam))


def _rich(f, N):
    return 2*f(2*N) - f(N)          # Richardson step for products whose log-terms are O(1/j^2)


def F11(mu, h, lam, N=400):
    """a_n soliton scalar factor (amplitude of identical solitons): G98 (2.6) with 3 -> h."""
    def lg(M):
        s = 0
        for j in range(1, M + 1):
            a = h*j*lam
            s += (mp.loggamma(mu + a - (h-1)*lam + 1) + mp.loggamma(mu + a - lam) + mp.loggamma(-mu + a - h*lam + 1)
                  + mp.loggamma(-mu + a) - mp.loggamma(-mu + a - (h-1)*lam + 1) - mp.loggamma(-mu + a - lam)
                  - mp.loggamma(mu + a - h*lam + 1) - mp.loggamma(mu + a))
        return s
    return -mp.exp(_rich(lg, N))


def A1(mu, n, lam, N=300):
    """DG (3.9): A_1 = -((n+1) lam/4)_I prod_k G_2k(2mu, -3/2 h lam + 1) G_2k(2mu, -1/2 h lam)
                        / [G_2k(2mu, -3/2 n lam - lam/2 + 1) G_2k(2mu, -1/2 n lam - 3/2 lam)]."""
    h = n + 1
    G = lambda j, m, a: mp.loggamma(m + j*h*lam + a) - mp.loggamma(-m + j*h*lam + a)

    def lg(M):
        return sum(G(2*k, 2*mu, -1.5*h*lam + 1) + G(2*k, 2*mu, -0.5*h*lam) - G(2*k, 2*mu, -1.5*n*lam - 0.5*lam + 1)
                   - G(2*k, 2*mu, -0.5*n*lam - 1.5*lam) for k in range(1, M + 1))
    return -blk(h*lam/4, mu, h, lam)*mp.exp(_rich(lg, N))


lam = mp.mpf(OMEGA)*3
mu = mp.mpc('0.13', '0.21')
for n in (2, 3, 4):
    h = n + 1
    ST = lambda m: mp.sin(mp.pi*m)/mp.sin(mp.pi*(lam - m))*F11(m, h, lam)
    r1 = abs(F11(mu, h, lam)*F11(-mu, h, lam) - 1); r2 = abs(ST(h*lam/2 - mu)*ST(h*lam/2 + mu) - 1)
    report(f'a_{n}: F_11 is unitary and the crossed transmission amplitude is unitary [{mp.nstr(max(r1, r2), 2)}]',
           max(r1, r2) < 1e-8)
    u1 = abs(A1(mu, n, lam)*A1(-mu, n, lam) - 1)
    c1 = abs(A1(-mu + h*lam/4, n, lam)/A1(mu + h*lam/4, n, lam)/F11(2*mu, h, lam) - 1)
    report(f'a_{n}: DG (3.9) satisfies A_1(mu)A_1(-mu) = 1 and the crossing relation (3.8) '
           f'[{mp.nstr(u1, 2)}, {mp.nstr(c1, 2)}]', u1 < 1e-12 and c1 < 1e-8)
ok = True
for n in (2, 3, 4, 5):
    h = n + 1
    for c in range(1, n + 1):
        sig = lambda m: -blk(c*lam/2, m, h, lam)*blk(h*lam/2 - c*lam/2, m, h, lam)
        ok &= abs(sig(mu)*sig(-mu) - 1) < 1e-20 and abs(sig(-mu + h*lam/4) - sig(mu + h*lam/4)) < 1e-20
        ok &= abs(sig(mu) - mp.fprod([sig(mu + h*lam/2 - k*lam) for k in range(1, n + 1)])) < 1e-20
report('DG (3.24): sigma_c = -(c lam/2)_I (h lam/2 - c lam/2)_I solves (3.22) and (3.23), n = 2..5, all c', ok)
sig_bad = lambda m: blk(lam/2, m, 5, lam)*blk(5*lam/2 - lam/2, m, 5, lam)
report('... without the minus sign (3.23) fails for even n (control)',
       abs(sig_bad(mu) - mp.fprod([sig_bad(mu + 5*lam/2 - k*lam) for k in range(1, 5)])) > 1e-3)

# ---------------------------------------------------------------- 6. semiclassical limit of A_1
# beta^2 = 4 pi/(lam + 1); M_1 = 2 h m_1/beta^2; time delay for C = +1 (eq-boundary-time-delay):
# Delta t = (2/(m_1 sinh theta)) ln(1 + s), s = sin(pi/h)/cosh theta;  d delta/d theta = M_1 sinh theta Delta t
# => (beta^2/(4h)) d(-i log A_1)/d theta -> ln(1 + sin(pi/h)/cosh theta) as lam -> infinity.
mp.mp.dps = 25
for n in (2, 3):
    h = n + 1
    th = mp.mpf('0.6')
    vals = []
    for L_ in (40, 80, 160):
        lamL = mp.mpf(L_)
        f = lambda t: mp.log(A1(-1j*h*lamL*t/(2*mp.pi), n, lamL, N=60))
        d = (f(th + mp.mpf('1e-6')) - f(th - mp.mpf('1e-6')))/mp.mpf('2e-6')
        beta2 = 4*mp.pi/(lamL + 1)
        vals.append(mp.re(-1j*d)*beta2/(4*h))
    target = mp.log(1 + mp.sin(mp.pi/h)/mp.cosh(th))
    errs = [abs(v - target) for v in vals]
    report(f'a_{n}: semiclassical limit of A_1 = (+...+) time delay: errors {[mp.nstr(e, 2) for e in errs]} '
           f'at lam = 40, 80, 160 (decreasing like 1/lam)', errs[2] < errs[1] < errs[0] and errs[2] < 0.05*target)
