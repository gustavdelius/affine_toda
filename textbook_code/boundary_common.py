"""Shared helpers for the half-line scripts (Part VII); prints nothing.

Conventions are those of Part VI (sec-qg-coproduct, sec-amplitude-and-gradation):

  coproduct   Delta(e_j) = e_j (x) 1 + k_j (x) e_j,  Delta(f_j) = f_j (x) k_j^-1 + 1 (x) f_j,  Delta(k_j) = k_j (x) k_j
  antipode    S(e_j) = -k_j^-1 e_j,  S(f_j) = -f_j k_j,  S(k_j) = k_j^-1
  physical q  q = -exp(-i pi omega) = exp(-4 pi^2 i/beta^2)          (eq-q-physical)
  rapidity    x(theta) = exp(2 T theta), principal gradation.

All representations are written in the principal gradation with the variable
y = x^(1/h_s), h_s = sum_j n_j s_j / s, so that every e_j carries one power of y and
every f_j one power of 1/y:  y(theta) = exp(lambda theta) for sine-Gordon (h_s = 2),
y(theta) = exp(omega theta) for the a_n vector solitons (h_s = h = n+1).  Using y
avoids the branch of x^(1/h).

A representation is a dict {('e', j): M, ('f', j): M, ('k', j): M, 'hw': index of the
highest-weight basis vector, 'dim': d}.

Coideal generators (sec-coideal, to be written):  b_j = e_j + c_j f_j k_j + eps_j k_j,
with Delta(b_j) = (b_j - eps_j k_j) (x) 1 + k_j (x) b_j  (a left coideal).
The parity-symmetric normalization is c_j = q_j^-1 (then q_j^-1 f_j k_j = k_j^(1/2) f_j k_j^(1/2)).

Intertwiners are found by a dense SVD null-space computation, with no filtering by weight.
"""
import numpy as np

kron = np.kron


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def qnum(n, q):
    return (q**n - q**(-n))/(q - 1/q)


def q_phys(omega):
    """eq-q-physical: q = -exp(-i pi omega)."""
    return -np.exp(-1j*np.pi*omega)


# ---------------------------------------------------------------- representations
def _mk(gens, hw):
    d = gens[('k', 0)].shape[0]
    rep = dict(gens)
    rep['hw'] = hw
    rep['dim'] = d
    return rep


def sl2_spin(j2, q, y):
    """Spin j = j2/2 evaluation representation of U_q(a_1^(1)), eq-spin-j and eq-evaluation-rep,
    in the principal gradation: e_1 = yE, f_1 = F/y, k_1 = K, e_0 = yF, f_0 = E/y, k_0 = K^-1.
    Basis v_0 (highest weight) ... v_{2j}.  For spin 1/2, v_0 is the soliton (T=+1)."""
    d = j2 + 1
    E = np.zeros((d, d), complex); F = np.zeros((d, d), complex)
    for k in range(d - 1):
        F[k+1, k] = 1
    for k in range(1, d):
        E[k-1, k] = qnum(k, q)*qnum(d - k, q)
    K = np.diag([q**(j2 - 2*k) for k in range(d)]).astype(complex)
    g = {('e', 1): y*E, ('f', 1): F/y, ('k', 1): K,
         ('e', 0): y*F, ('f', 0): E/y, ('k', 0): np.linalg.inv(K)}
    return _mk(g, 0)


def an_vector(n, q, y, conj=False):
    """Vector representation of U_q(a_n^(1)) in the principal gradation.
    Basis v_1..v_{n+1} (indices 0..n), weights eps_i; alpha_j = eps_j - eps_{j+1} (j=1..n),
    alpha_0 = eps_{n+1} - eps_1.  E_j = e_{j,j+1}, F_j = E_j^T, H_j = e_jj - e_{j+1,j+1}.
    V(y):    e_j = y E_j, f_j = F_j/y, k_j = q^{H_j}.
    Vbar(y): e_j = y F_j, f_j = E_j/y, k_j = q^{-H_j}  (V composed with the Chevalley
             involution e<->f, k->k^-1), weights -eps_i.  At the physical q this is the
             antisoliton at the same rapidity: V(z)* = Vbar(-z/q) and -1/q = exp(i pi omega)."""
    N = n + 1
    g = {}
    for j in range(N):
        a, b = (j - 1) % N, j % N          # node j>=1: rows a=j-1 -> b=j  (eps_j - eps_{j+1} in 1-based)
        if j == 0:
            a, b = n, 0                    # alpha_0 = eps_{n+1} - eps_1
        E = np.zeros((N, N), complex); E[a, b] = 1
        F = E.T.copy()
        H = np.zeros(N); H[a] += 1; H[b] -= 1
        if not conj:
            g[('e', j)] = y*E; g[('f', j)] = F/y; g[('k', j)] = np.diag(q**H).astype(complex)
        else:
            g[('e', j)] = y*F; g[('f', j)] = E/y; g[('k', j)] = np.diag(q**(-H)).astype(complex)
    return _mk(g, n if conj else 0)


def dn_vector(n, q, y):
    """Vector representation (2n-dim) of U_q(d_n^(1)), n>=4, principal gradation (h = 2n-2).
    Basis v_1..v_n, v_{-n}..v_{-1} (indices 0..n-1, n..2n-1; v_{-i} at index 2n-i).
    alpha_i = eps_i - eps_{i+1} (i=1..n-1), alpha_n = eps_{n-1} + eps_n, alpha_0 = -eps_1 - eps_2.
    For later use; its relations are checked by check_relations."""
    N = 2*n
    idx = lambda i: i - 1 if i > 0 else N + i      # i in +-1..+-n
    def unit(r, c):
        M = np.zeros((N, N), complex); M[idx(r), idx(c)] = 1; return M
    Es = {}
    for i in range(1, n):
        Es[i] = unit(i, i+1) - unit(-(i+1), -i)
    Es[n] = unit(n-1, -n) - unit(n, -(n-1))
    Es[0] = unit(-2, 1) - unit(-1, 2)
    g = {}
    for j, E in Es.items():
        F = E.conj().T
        H = np.real(np.diag(E @ F - F @ E))
        g[('e', j)] = y*E; g[('f', j)] = F/y; g[('k', j)] = np.diag(q**H).astype(complex)
    return _mk(g, 0)


def check_relations(rep, cartan, q):
    """Max residual of the Drinfeld-Jimbo relations (all q_i = q: simply laced) including q-Serre."""
    nodes = range(len(cartan))
    e = lambda j: rep[('e', j)]; f = lambda j: rep[('f', j)]; k = lambda j: rep[('k', j)]
    res = 0.0
    for i in nodes:
        ki = np.linalg.inv(k(i))
        for j in nodes:
            a = cartan[i][j]
            res = max(res, np.abs(k(i) @ e(j) @ ki - q**a*e(j)).max(),
                      np.abs(k(i) @ f(j) @ ki - q**(-a)*f(j)).max(),
                      np.abs(k(i) @ k(j) - k(j) @ k(i)).max())
            c = e(i) @ f(j) - f(j) @ e(i)
            want = (k(i) - ki)/(q - 1/q) if i == j else 0*c
            res = max(res, np.abs(c - want).max())
            if i != j:
                m = 1 - a
                for X in (e, f):
                    s = sum((-1)**t*_qbinom(m, t, q)*np.linalg.matrix_power(X(i), m - t) @ X(j) @
                            np.linalg.matrix_power(X(i), t) for t in range(m + 1))
                    res = max(res, np.abs(s).max())
    return res


def _qbinom(m, t, q):
    num = np.prod([qnum(m - s, q) for s in range(t)]) if t else 1
    den = np.prod([qnum(s + 1, q) for s in range(t)]) if t else 1
    return num/den


def cartan_a(n):
    N = n + 1
    if n == 1:
        return [[2, -2], [-2, 2]]
    C = [[0]*N for _ in range(N)]
    for i in range(N):
        C[i][i] = 2; C[i][(i+1) % N] = -1; C[i][(i-1) % N] = -1
    return C


def cartan_d(n):
    # nodes 0..n; 0 attached to 2, n attached to n-2 (with alpha_n = eps_{n-1}+eps_n)
    C = [[0]*(n+1) for _ in range(n+1)]
    for i in range(n+1):
        C[i][i] = 2
    edges = [(0, 2)] + [(i, i+1) for i in range(1, n-1)] + [(n-2, n)]
    for a, b in edges:
        C[a][b] = C[b][a] = -1
    return C


# ---------------------------------------------------------------- Hopf structure
def nodes_of(rep):
    return sorted(key[1] for key in rep if isinstance(key, tuple) and key[0] == 'e')


def cop(A, B, kind, j):
    """Part VI coproduct on A (x) B."""
    IA, IB = np.eye(A['dim']), np.eye(B['dim'])
    if kind == 'e':
        return kron(A[('e', j)], IB) + kron(A[('k', j)], B[('e', j)])
    if kind == 'f':
        return kron(A[('f', j)], np.linalg.inv(B[('k', j)])) + kron(IA, B[('f', j)])
    return kron(A[('k', j)], B[('k', j)])


def dual(rep):
    """pi*(a) = pi(S(a))^T with the Part VI antipode."""
    g = {}
    for j in nodes_of(rep):
        k = rep[('k', j)]; ki = np.linalg.inv(k)
        g[('e', j)] = (-ki @ rep[('e', j)]).T
        g[('f', j)] = (-rep[('f', j)] @ k).T
        g[('k', j)] = ki.T
    return _mk(g, None)


def coideal_gen(rep, j, c, eps):
    """b_j = e_j + c f_j k_j + eps k_j on a single representation."""
    return rep[('e', j)] + c*rep[('f', j)] @ rep[('k', j)] + eps*rep[('k', j)]


def coideal_gen_tensor(A, B, j, c, eps):
    """Delta(b_j) = (e_j + c f_j k_j) (x) 1 + k_j (x) b_j on A (x) B."""
    return kron(coideal_gen(A, j, c, 0), np.eye(B['dim'])) + kron(A[('k', j)], coideal_gen(B, j, c, eps))


# ---------------------------------------------------------------- dense intertwiner solver
def intertwiners(pairs, tol=1e-9):
    """All X (dim_W x dim_U) with X A = B X for every (A, B) in pairs (A on U, B on W).
    Null space of the stacked linear system  (A^T (x) 1 - 1 (x) B) vec(X) = 0  (column-major vec),
    with no weight filtering: dense SVD for up to 400 unknowns, otherwise the eigenvectors of the
    sparse Gram matrix M^H M (eigenvalues = squared singular values).
    Returns (list of X, relative singular values sorted ascending)."""
    import scipy.sparse as sps
    A0, B0 = pairs[0]
    dU, dW = A0.shape[0], B0.shape[0]
    if dU*dW <= 400:
        IU, IW = np.eye(dU), np.eye(dW)
        M = np.vstack([kron(A.T, IW) - kron(IU, B) for A, B in pairs])
        if M.shape[0] < M.shape[1]:
            M = np.vstack([M, np.zeros((M.shape[1] - M.shape[0], M.shape[1]))])
        _, s, Vh = np.linalg.svd(M, full_matrices=False)
        rel = s/s[0]
        vecs = [Vh[i].conj() for i in range(len(s)) if rel[i] < tol]
    else:
        IU, IW = sps.identity(dU, format='csr'), sps.identity(dW, format='csr')
        G = None
        for A, B in pairs:
            Mg = sps.kron(sps.csr_matrix(A.T), IW) - sps.kron(IU, sps.csr_matrix(B))
            G = Mg.conj().T @ Mg if G is None else G + Mg.conj().T @ Mg
        ev, U = np.linalg.eigh(G.toarray())
        ev = np.clip(ev, 0, None)
        rel = np.sqrt(ev/ev[-1])
        vecs = [U[:, i] for i in range(len(ev)) if rel[i] < max(tol, 1e-6)]   # Gram: null sigma ~ sqrt(eps)
    null = [v.reshape(dU, dW).T for v in vecs]
    return null, np.sort(rel)


def rcheck(A, B, tol=1e-9):
    """R-check: A (x) B -> B (x) A intertwining the Part VI coproduct (Jimbo's equations),
    normalized to 1 on hw_A (x) hw_B -> hw_B (x) hw_A.  Returns (R, nullity)."""
    pairs = []
    for j in nodes_of(A):
        for kind in 'efk':
            pairs.append((cop(A, B, kind, j), cop(B, A, kind, j)))
    null, rel = intertwiners(pairs, tol)
    assert len(null) >= 1, 'no intertwiner'
    R = null[0]
    i_in = A['hw']*B['dim'] + B['hw']
    i_out = B['hw']*A['dim'] + A['hw']
    return R/R[i_out, i_in], len(null)


def k_intertwiner(Vin, Vout, c, eps, tol=1e-9):
    """K: Vin -> Vout with K b_j = b_j K for the coideal generators b_j (all nodes)."""
    pairs = [(coideal_gen(Vin, j, c[j], eps[j]), coideal_gen(Vout, j, c[j], eps[j]))
             for j in nodes_of(Vin)]
    return intertwiners(pairs, tol)


# ---------------------------------------------------------------- reflection equations
def one_k(K, dim_left):
    return kron(np.eye(dim_left), K)


def prop_residual(L, R):
    """min_c |L - c R| / |L| and the optimal c."""
    c = np.vdot(R.ravel(), L.ravel())/np.vdot(R.ravel(), R.ravel())
    return np.linalg.norm(L - c*R)/np.linalg.norm(L), c


def re_preserving(Rf, Kf, y1, y2):
    """Soliton-preserving reflection equation in braid form, K(y): V(y) -> V(1/y):
    R(y1/y2) (1 (x) K(y1)) R(y1 y2) (1 (x) K(y2)) = (1 (x) K(y2)) R(y1 y2) (1 (x) K(y1)) R(y1/y2).
    Rf(z) is R-check on V (x) V at ratio z of the y-variables."""
    K1, K2 = Kf(y1), Kf(y2); d = K1.shape[0]
    L = Rf(y1/y2) @ one_k(K1, d) @ Rf(y1*y2) @ one_k(K2, d)
    R = one_k(K2, d) @ Rf(y1*y2) @ one_k(K1, d) @ Rf(y1/y2)
    return prop_residual(L, R)


def re_twisted(Rf, Kf, mu, nu, y1, y2):
    """Soliton-conjugating reflection equation for multiplets mu (rapidity theta_1, y1) and
    nu (theta_2, y2), from factorization (DM (2.8)); K_mu(y): mu(y) -> mubar(1/y):
    R_{nubar mubar}(y1/y2) K_mu(y1)_2 R_{mu nubar}(y1 y2) K_nu(y2)_2
        = K_nu(y2)_2 R_{nu mubar}(y1 y2) K_mu(y1)_2 R_{mu nu}(y1/y2).
    Rf(a, b, z): R-check on a (x) b at ratio z; Kf(a, y): K_a(y); bar(a) = the conjugate label."""
    bar = {'V': 'W', 'W': 'V'}
    Kmu, Knu = Kf(mu, y1), Kf(nu, y2); d = Kmu.shape[0]
    L = Rf(bar[nu], bar[mu], y1/y2) @ one_k(Kmu, d) @ Rf(mu, bar[nu], y1*y2) @ one_k(Knu, d)
    R = one_k(Knu, d) @ Rf(nu, bar[mu], y1*y2) @ one_k(Kmu, d) @ Rf(mu, nu, y1/y2)
    return prop_residual(L, R)
