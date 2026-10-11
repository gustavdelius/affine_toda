"""Checks for sec-qsp (chapter 34, quantum symmetric pairs and universal K-matrices).

Conventions are those of Part VI and of boundary_common.py: Delta(e) = e (x) 1 + k (x) e,
Delta(f) = f (x) k^-1 + 1 (x) f, physical q = -exp(-i pi omega), book coideal generators
b_j = e_j + c_j f_j k_j + eps_j k_j with c_j = q_j^-1 (eq-coideal-generators), a left coideal.
Kolb (arXiv:1207.6036) uses the same coproduct (his (3.3)-(3.4)) and right coideal generators
B_i = F_i - c_i E_i K_i^-1 + s_i K_i^-1 for the Satake diagram (X, tau) = (emptyset, id).

 1. q-Onsager (q-Dolan-Grady) relations of the sine-Gordon coideal (sec-q-onsager,
    eq-q-dolan-grady): F_ij(b_i, b_j) = q c_i [2]^2 (b_i b_j - b_j b_i) for (i,j) = (0,1), (1,0),
    on spin 1/2 (x) 1/2, spin 1 and 1/2 (x) 1 (x) boundary singlet, for generic eps, generic c,
    physical and generic q; independent of eps.  Onsager limit q = 1:
    [b0,[b0,[b0,b1]]] = 4 [b0,b1] (Dolan-Grady with A = 2 b_0).  Existence and uniqueness of K.
 2. Kolb's generators and the map omega (sec-letzter-kolb, eq-omega-map, eq-kolb-dictionary):
    Kolb's B_i is a right coideal; omega (e<->f, k->k^-1) is an algebra automorphism with
    Delta omega = (omega (x) omega) Delta^op; omega(B_i) = e_i - c_i f_i k_i + s_i k_i, so
    c_book = -c_Kolb, eps = s.  Kolb's q-Onsager relations (7.19) in his conventions; his
    (7.10) for a_ij = -1 at s = 0, and the obstruction of his Remark 5.12 for s != 0:
    F_ij(B_i,B_j) + q c_i B_j = q s_j (c_i - (q^-1 - 1)^2 s_i^2) K_i^-2 K_j^-1.
 3. Adjacent nodes, book conventions (prp-serre-obstruction): for a_ij = -1
    F_ij(b_i,b_j) = q c_i b_j + eps_j((1-q)(1-q^-1) eps_i^2 - q c_i) k_i^2 k_j,
    on a_2, a_3 (V (x) Vbar (x) V, V (x) Vbar) and d_4 (V (x) V) with random eps and c; b_i, b_j
    commute when a_ij = 0.  Admissibility table: for a_2, a_3, d_4 and sign/zero patterns of eps,
    "all obstructions vanish" <=> "K exists" (prp-an-k-existence, prp-dn-k).
 4. One-dimensional B-modules (prp-one-dim-modules): scalar relations force beta in {0} u {+-eps_*}^N;
    in the a_n vector representation V(y) (n = 2, 3) there is no one-dimensional B-submodule at
    eps = 0, and for s = +-1 exactly one, C_eps itself, at y^(n+1) = s(-q)^((n+1)/2) (theta = -i pi/2
    for s = +1).
 5. Table rows of tbl-physics-qsp that are marked "checked": sine-Gordon (q-Onsager, s in S),
    a_n eps = 0 (standard, K = 1), a_n s = +-1 (generalized q-Onsager, s not in S, K exists),
    d_n eps = 0 (standard, K anti-diagonal), d_n +-eps_* (generalized q-Onsager); classical limit
    of B_0: e_j + f_j fixed by theta(X)(y) = -J X(1/y)^T J^-1 (x loop rotation for n even).
"""
import numpy as np
from boundary_common import (report, q_phys, qnum, sl2_spin, an_vector, dn_vector, cartan_a, cartan_d,
                             cop, coideal_gen, k_intertwiner, nodes_of, _mk)

rng = np.random.default_rng(34)
kron = np.kron


def tensor(A, B):
    """Tensor product representation A (x) B of U_q with the Part VI coproduct."""
    g = {}
    for j in nodes_of(A):
        for kind in 'efk':
            g[(kind, j)] = cop(A, B, kind, j)
    return _mk(g, A['hw']*B['dim'] + B['hw'])


def omega_rep(A):
    """A composed with omega: e -> f, f -> e, k -> k^-1."""
    g = {}
    for j in nodes_of(A):
        g[('e', j)] = A[('f', j)]; g[('f', j)] = A[('e', j)]; g[('k', j)] = np.linalg.inv(A[('k', j)])
    return _mk(g, A['hw'])


def serre(x, y, a, q):
    """Kolb's F_ij(x, y) = sum_n (-1)^n [1-a choose n]_q x^(1-a-n) y x^n (palindromic q-Serre)."""
    m = 1 - a
    def qb(m, t):
        num = np.prod([qnum(m - s, q) for s in range(t)]) if t else 1
        den = np.prod([qnum(s + 1, q) for s in range(t)]) if t else 1
        return num/den
    return sum((-1)**t*qb(m, t)*np.linalg.matrix_power(x, m - t) @ y @ np.linalg.matrix_power(x, t)
               for t in range(m + 1))


def rel(L, R):
    return np.linalg.norm(L - R)/max(np.linalg.norm(L), np.linalg.norm(R), 1e-300)


def crand(n=1, s=1.0):
    return s*(rng.normal(size=n) + 1j*rng.normal(size=n))


def eps_star(q):
    return 1/np.sqrt((1 - q)*(1 - 1/q))


OMEGA = 0.37
QP = q_phys(OMEGA)                       # physical q (sine-Gordon: lambda = omega)
QG = 0.83*np.exp(0.41j)                  # a generic complex q

# ---------------------------------------------------------------- 1. q-Onsager relations, sine-Gordon
worst = 0.0; worst_bad = np.inf
for q in (QP, QG):
    for trial in range(3):
        c = [1/q, 1/q] if trial == 0 else list(crand(2))
        eps = list(crand(2, 1.5))
        y1, y2, y3 = np.exp(crand(3, 0.4))
        reps = [tensor(sl2_spin(1, q, y1), sl2_spin(1, q, y2)),                      # 1/2 (x) 1/2
                sl2_spin(2, q, y3),                                                    # spin 1
                tensor(sl2_spin(1, q, y1), sl2_spin(2, q, y2))]                        # 1/2 (x) 1
        for W in reps:
            # boundary singlet C_eps at the right end: b_j acts as e_j + c f_j k_j + eps_j k_j on W
            b = [coideal_gen(W, j, c[j], eps[j]) for j in (0, 1)]
            for i, j in ((0, 1), (1, 0)):
                rho = q*c[i]*(q + 1/q)**2
                L = serre(b[i], b[j], -2, q)
                worst = max(worst, rel(L, rho*(b[i] @ b[j] - b[j] @ b[i])))
                worst_bad = min(worst_bad, rel(L, 1.1*rho*(b[i] @ b[j] - b[j] @ b[i])))
report('q-Dolan-Grady F_ij(b_i,b_j) = q c_i [2]^2 [b_i,b_j], sine-Gordon, generic eps and c, '
       'physical and generic q, 1/2x1/2, 1, 1/2x1 [%.1e]' % worst, worst < 1e-10)
report('  (control: rho -> 1.1 rho fails [%.2f])' % worst_bad, worst_bad > 1e-3)

# with c = q^-1 the constant is [2]^2 = 4 cos^2(pi lambda) at the physical q
report('rho = (q+q^-1)^2 = 4 cos^2(pi lambda) at q = -exp(-i pi lambda)',
       abs((QP + 1/QP)**2 - 4*np.cos(np.pi*OMEGA)**2) < 1e-12)


# Onsager limit q = 1: loop algebra of sl_2, b_j = e_j + f_j + eps_j
def sl2_spin_q1(j2, y):
    d = j2 + 1
    E = np.zeros((d, d), complex); F = np.zeros((d, d), complex)
    for k in range(d - 1):
        F[k+1, k] = 1
    for k in range(1, d):
        E[k-1, k] = k*(d - k)
    I = np.eye(d)
    return _mk({('e', 1): y*E, ('f', 1): F/y, ('k', 1): I, ('e', 0): y*F, ('f', 0): E/y, ('k', 0): I.copy()}, 0)


worst = 0.0
for trial in range(3):
    eps = crand(2, 1.5); y1, y2 = np.exp(crand(2, 0.4))
    W = tensor(sl2_spin_q1(1, y1), sl2_spin_q1(2, y2))
    b = [coideal_gen(W, j, 1.0, eps[j]) for j in (0, 1)]
    for i, j in ((0, 1), (1, 0)):
        comm = lambda X, Y: X @ Y - Y @ X
        A, As = 2*b[i], 2*b[j]
        worst = max(worst, rel(comm(A, comm(A, comm(A, As))), 16*comm(A, As)))
report('Onsager limit q=1: Dolan-Grady [A,[A,[A,A*]]] = 16[A,A*] for A = 2b_0, A* = 2b_1 [%.1e]' % worst,
       worst < 1e-12)

# sine-Gordon K exists and is unique for generic eps (eq-sg-k-coideal), physical q
nulls = []
for trial in range(3):
    eps = list(crand(2, 1.5)); y = np.exp(crand(1, 0.5)[0])
    null, _ = k_intertwiner(sl2_spin(1, QP, y), sl2_spin(1, QP, 1/y), [1/QP, 1/QP], eps)
    nulls.append(len(null))
report('sine-Gordon: K: V(y) -> V(1/y) exists and is unique for generic (eps_0, eps_1) %s' % nulls,
       all(n == 1 for n in nulls))

# ---------------------------------------------------------------- 2. Kolb's generators and omega
def kolb_B(A, j, c, s):
    """Kolb's B_j = F_j - c E_j K_j^-1 + s K_j^-1 on the representation A."""
    Ki = np.linalg.inv(A[('k', j)])
    return A[('f', j)] - c*A[('e', j)] @ Ki + s*Ki


worst_rc = 0.0; worst_om = 0.0; worst_img = 0.0
q = QG
for n in (1, 2):
    for trial in range(2):
        ys = np.exp(crand(2, 0.4))
        A = sl2_spin(1, q, ys[0]) if n == 1 else an_vector(n, q, ys[0])
        B = sl2_spin(2, q, ys[1]) if n == 1 else an_vector(n, q, ys[1], conj=True)
        W = tensor(A, B)
        dB, dA = B['dim'], A['dim']
        P = np.zeros((dA*dB, dA*dB))
        for a in range(dA):
            for bb in range(dB):
                P[bb*dA + a, a*dB + bb] = 1           # flip A (x) B -> B (x) A
        for j in nodes_of(A):
            c, s = crand(2)
            # right coideal: Delta(B_j) = B_j (x) K_j^-1 + 1 (x) (F_j - c E_j K_j^-1)
            lhs = kolb_B(W, j, c, s)
            rhs = kron(kolb_B(A, j, c, s), np.linalg.inv(B[('k', j)])) + kron(np.eye(dA), kolb_B(B, j, c, 0))
            worst_rc = max(worst_rc, rel(lhs, rhs))
            # Delta(omega x) = (omega (x) omega) Delta^op(x): on A (x) B,
            # (A (x) B)(Delta(omega x)) = P^T (wB (x) wA)(Delta x) P with wX = X o omega
            wA, wB = omega_rep(A), omega_rep(B)
            for kind, wkind in (('e', 'f'), ('f', 'e')):
                worst_om = max(worst_om, rel(cop(A, B, wkind, j), P.T @ cop(wB, wA, kind, j) @ P))
            worst_om = max(worst_om, rel(np.linalg.inv(cop(A, B, 'k', j)), P.T @ cop(wB, wA, 'k', j) @ P))
            # omega(B_j) = e_j - c f_j k_j + s k_j = book generator with c_book = -c, eps = s
            worst_img = max(worst_img, rel(kolb_B(omega_rep(W), j, c, s), coideal_gen(W, j, -c, s)))
report("Kolb's B_i = F_i - c_i E_i K_i^-1 + s_i K_i^-1 is a right coideal (a_1, a_2 tensor products) [%.1e]"
       % worst_rc, worst_rc < 1e-12)
report('omega (e<->f, k->k^-1) satisfies Delta o omega = (omega x omega) o Delta^op [%.1e]' % worst_om,
       worst_om < 1e-12)
report('omega(B_i) = e_i - c_i f_i k_i + s_i k_i: c_book = -c_Kolb, eps = s [%.1e]' % worst_img,
       worst_img < 1e-12)

# Kolb (7.19): q-Onsager relations in his conventions, any s
worst = 0.0
for q in (QP, QG):
    c = crand(2); s = crand(2, 1.5); y1, y2 = np.exp(crand(2, 0.4))
    W = tensor(sl2_spin(1, q, y1), sl2_spin(2, q, y2))
    Bm = [kolb_B(W, j, c[j], s[j]) for j in (0, 1)]
    for i, j in ((0, 1), (1, 0)):
        worst = max(worst, rel(serre(Bm[i], Bm[j], -2, q), q*(q + 1/q)**2*c[i]*(Bm[j] @ Bm[i] - Bm[i] @ Bm[j])))
report("Kolb (7.19), q-Onsager algebra in Kolb's conventions, generic c, s [%.1e]" % worst, worst < 1e-10)

# Kolb (7.10) for a_ij = -1 with tau = id, Z_i = -1: F_ij(B_i,B_j) = -q c_i B_j at s = 0;
# for s != 0 the obstruction of his Remark 5.12
worst0 = 0.0; worst_obs = 0.0; obs_size = np.inf
for n in (2, 3):
    q = QG
    C = cartan_a(n)
    y = np.exp(crand(3, 0.4))
    W = tensor(tensor(an_vector(n, q, y[0]), an_vector(n, q, y[1], conj=True)), an_vector(n, q, y[2]))
    c = crand(n + 1)
    for s in (np.zeros(n + 1), crand(n + 1, 2.0)):
        Bm = [kolb_B(W, j, c[j], s[j]) for j in range(n + 1)]
        Kinv = [np.linalg.inv(W[('k', j)]) for j in range(n + 1)]
        for i in range(n + 1):
            for j in range(n + 1):
                if C[i][j] != -1:
                    continue
                L = serre(Bm[i], Bm[j], -1, q) + q*c[i]*Bm[j]
                obs = q*s[j]*(c[i] - (1/q - 1)**2*s[i]**2)*Kinv[i] @ Kinv[i] @ Kinv[j]
                if not s.any():
                    worst0 = max(worst0, np.linalg.norm(L)/np.linalg.norm(Bm[j]))
                else:
                    worst_obs = max(worst_obs, rel(L, obs))
                    obs_size = min(obs_size, np.linalg.norm(obs)/np.linalg.norm(Bm[j]))
report("Kolb (7.10), a_ij=-1, s=0: F_ij(B_i,B_j) = -q c_i B_j (a_2, a_3) [%.1e]" % worst0, worst0 < 1e-10)
report("Kolb Remark 5.12 made explicit: F_ij(B_i,B_j) + q c_i B_j = q s_j (c_i - (q^-1-1)^2 s_i^2) "
       "K_i^-2 K_j^-1 (a_2, a_3, generic s) [%.1e; obstruction/|B| >= %.2f]" % (worst_obs, obs_size),
       worst_obs < 1e-10 and obs_size > 1e-3)

# ---------------------------------------------------------------- 3. adjacent nodes, book conventions
def book_obstruction(W, i, j, c, eps, q):
    k = W[('k', i)]; kj = W[('k', j)]
    return eps[j]*((1 - q)*(1 - 1/q)*eps[i]**2 - q*c[i])*k @ k @ kj


worst = 0.0; worst_comm = 0.0
for q in (QP, QG):
    cases = []
    for n in (2, 3):
        y = np.exp(crand(3, 0.4))
        cases.append((cartan_a(n), tensor(tensor(an_vector(n, q, y[0]), an_vector(n, q, y[1], conj=True)),
                                          an_vector(n, q, y[2])) if n == 2 else
                      tensor(an_vector(n, q, y[0]), an_vector(n, q, y[1], conj=True))))
    y = np.exp(crand(2, 0.4))
    cases.append((cartan_d(4), tensor(dn_vector(4, q, y[0]), dn_vector(4, q, y[1]))))
    for C, W in cases:
        nn = len(C)
        c = crand(nn); eps = crand(nn, 1.5)
        b = [coideal_gen(W, j, c[j], eps[j]) for j in range(nn)]
        for i in range(nn):
            for j in range(nn):
                if i == j:
                    continue
                if C[i][j] == -1:
                    L = serre(b[i], b[j], -1, q)
                    R = q*c[i]*b[j] + book_obstruction(W, i, j, c, eps, q)
                    worst = max(worst, rel(L, R))
                elif C[i][j] == 0:
                    worst_comm = max(worst_comm, rel(b[i] @ b[j], b[j] @ b[i]))
report('book: F_ij(b_i,b_j) = q c_i b_j + eps_j((1-q)(1-q^-1)eps_i^2 - q c_i) k_i^2 k_j for a_ij=-1 '
       '(a_2, a_3, d_4; random eps, c; physical and generic q) [%.1e]' % worst, worst < 1e-10)
report('book: b_i b_j = b_j b_i for a_ij = 0 (a_3, d_4) [%.1e]' % worst_comm, worst_comm < 1e-10)


def obstructions_vanish(C, eps, q, tol=1e-9):
    es2 = eps_star(q)**2
    for i in range(len(C)):
        for j in range(len(C)):
            if C[i][j] == -1 and abs(eps[j]*(eps[i]**2/es2 - 1)) > tol:
                return False
    return True


def k_exists(kind, n, eps, q):
    y = np.exp(0.31 + 0.17j)
    if kind == 'a':
        null, _ = k_intertwiner(an_vector(n, q, y), an_vector(n, q, 1/y, conj=True), [1/q]*(n + 1), eps)
    else:
        null, _ = k_intertwiner(dn_vector(n, q, y), dn_vector(n, q, 1/y), [1/q]*(n + 1), eps)
    return len(null)


q = QP
es = eps_star(q)
agree = True; summary = []
for kind, n in (('a', 2), ('a', 3), ('d', 4)):
    C = cartan_a(n) if kind == 'a' else cartan_d(n)
    N = len(C)
    pats = {'zero': [0]*N, 'all +': [es]*N, 'one -': [-es] + [es]*(N - 1), 'two -': [-es, -es] + [es]*(N - 2),
            'mixed 0': [0] + [es]*(N - 1), 'one off': [1.3*es] + [es]*(N - 1), 'one nonzero': [es] + [0]*(N - 1),
            'random': list(crand(N, 1.5))}
    for name, eps in pats.items():
        ov = obstructions_vanish(C, eps, q)
        nk = k_exists(kind, n, eps, q)
        agree &= (ov == (nk >= 1))
        summary.append('%s%d %s: %s/%d' % (kind, n, name, 'ok' if ov else 'obs', nk))
report('admissibility: "all Serre obstructions vanish" <=> "K exists" (a_2, a_3, d_4; 8 eps patterns each)',
       agree)
print('   ', '; '.join(summary))

# ---------------------------------------------------------------- 4. one-dimensional B-submodules
# Gauge y onto node 0 (D = diag(y^i)): the finite-node generators b_1..b_n are y-independent and node 0
# carries z = y^(n+1).  A one-dimensional B-submodule of V(y) is a common eigenvector of b_1..b_n (which lies
# in the joint kernel of the commutators [b_i, b_j], i, j >= 1) on which b_0(z) also acts by a scalar:
# b_0(z) v = z E_0 v + z^-1 q^-1 F_0 k_0 v + eps_0 k_0 v = beta_0 v is linear in (z, 1/z, beta_0).
q = QP; es = eps_star(q); mq = np.exp(-1j*np.pi*OMEGA)          # -q = exp(-i pi omega)
info = []; ok = True
for n in (2, 3):
    N = n + 1
    V = an_vector(n, q, 1.0)
    for name, eps, s in (('eps=0', [0]*N, 0), ('s=+1', [es]*N, 1), ('s=-1', [-es] + [es]*n, -1)):
        b = [coideal_gen(V, j, 1/q, eps[j]) for j in range(N)]
        M = np.vstack([b[i] @ b[j] - b[j] @ b[i] for i in range(1, N) for j in range(i + 1, N)])
        _, sv, Vh = np.linalg.svd(M)
        ker = [Vh[k].conj() for k in range(len(sv)) if sv[k]/sv[0] < 1e-10]
        found = []
        for v in ker:
            if max(np.linalg.norm(b[j] @ v - np.vdot(v, b[j] @ v)*v) for j in range(1, N)) > 1e-8:
                continue
            beta = [None] + [np.vdot(v, b[j] @ v) for j in range(1, N)]
            A = np.column_stack([V[('e', 0)] @ v, (1/q)*V[('f', 0)] @ V[('k', 0)] @ v, -v])
            sol = np.linalg.lstsq(A, -eps[0]*V[('k', 0)] @ v, rcond=None)[0]
            if np.linalg.norm(A @ sol + eps[0]*V[('k', 0)] @ v) < 1e-10 and abs(sol[0]*sol[1] - 1) < 1e-10:
                beta[0] = sol[2]; found.append((sol[0], beta))
        if s == 0:
            ok &= (len(found) == 0)
            info.append('%s %s: none' % ('a%d' % n, name))
        else:
            ok &= len(found) == 1 and np.allclose(found[0][1], eps, atol=1e-10) and \
                abs(found[0][0] - s*mq**(N/2)) < 1e-10
            info.append('a%d %s: beta = eps at y^N = s(-q)^(N/2)' % (n, name))
report('one-dimensional B-submodules of V(y): none for eps = 0; for s = +-1 exactly C_eps (beta = eps) '
       'at y^(n+1) = s(-q)^((n+1)/2) (n = 2, 3)', ok)
print('   ', '; '.join(info))

# scalar version of prp-one-dim-modules: beta_j((1-q)(1-q^-1)beta_i^2 - 1) = 0 on adjacent nodes
ok = True
for beta, allowed in (([0, 0, 0], True), ([es, es, es], True), ([-es, es, es], True), ([es, es, 0], False),
                      ([es, 0, 0], False), ([1.2*es, es, es], False)):
    lhs = [beta[j]*((1 - q)*(1 - 1/q)*beta[i]**2 - 1) for i in range(3) for j in range(3) if i != j]
    ok &= (max(abs(np.array(lhs))) < 1e-12) == allowed
report('one-dimensional modules of the a_2 generalized q-Onsager algebra: beta in {0} u {+-eps_*}^3 only', ok)

# ---------------------------------------------------------------- 5. table rows
q = QP; es = eps_star(q)
# a_n eps = 0: K = 1 (standard QSP, s = 0); a_n s = +-1: K exists (generalized q-Onsager, s not in S)
rows = []
for n in (2, 3):
    y = np.exp(0.4 + 0.2j)
    for eps, name in (([0]*(n + 1), 'eps=0'), ([es]*(n + 1), 's=+1'), ([-es] + [es]*n, 's=-1')):
        null, _ = k_intertwiner(an_vector(n, q, y), an_vector(n, q, 1/y, conj=True), [1/q]*(n + 1), eps)
        Kk = null[0]/null[0][0, 0] if null else None
        isone = Kk is not None and np.allclose(Kk, np.eye(n + 1), atol=1e-9)
        rows.append((n, name, len(null), isone))
report('table: a_n eps=0 -> K = 1; s = +-1 -> unique non-diagonal K (n = 2, 3) %s'
       % [(r[0], r[1], r[2]) for r in rows],
       all(r[2] == 1 for r in rows) and all(r[3] == (r[1] == 'eps=0') for r in rows))
# d_4 eps = 0: anti-diagonal; +-eps_*: exists
y = np.exp(0.4 + 0.2j)
null0, _ = k_intertwiner(dn_vector(4, q, y), dn_vector(4, q, 1/y), [1/q]*5, [0]*5)
K0 = null0[0]
anti = np.linalg.norm(K0*(1 - np.fliplr(np.eye(8)))) < 1e-9*np.linalg.norm(K0)
nullp, _ = k_intertwiner(dn_vector(4, q, y), dn_vector(4, q, 1/y), [1/q]*5, [es, -es, es, es, -es])
report('table: d_4 eps=0 -> unique anti-diagonal K; eps_j = +-eps_* -> unique K (%d, %d)'
       % (len(null0), len(nullp)), len(null0) == 1 and anti and len(nullp) == 1)

# classical limit of B_0: e_j + f_j is fixed by theta(X)(y) = -J X(1/y)^T J^-1, J = diag((-1)^i),
# for the a_n vector at q = 1 (for n even combined with y -> -y on node 0)
ok = True
for n in (2, 3, 4):
    N = n + 1
    J = np.diag([(-1)**i for i in range(N)]).astype(complex)
    y = np.exp(0.3 + 0.2j)
    V1 = an_vector(n, 1.0, y); Vinv = an_vector(n, 1.0, 1/y); Vm = an_vector(n, 1.0, -1/y)
    for j in range(N):
        X = V1[('e', j)] + V1[('f', j)]
        tgt = Vinv if (j > 0 or n % 2 == 1) else Vm
        th = -J @ (tgt[('e', j)] + tgt[('f', j)]).T @ np.linalg.inv(J)
        ok &= np.allclose(th, X)
report('classical limit: theta fixes e_j + f_j (a_2, a_3, a_4; loop rotation y -> -y needed on node 0 for n even)',
       ok)
