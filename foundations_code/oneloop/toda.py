"""Simply-laced affine Toda data, Hirota single solitons, and Hirota fluctuation solutions.

Static Hirota equations (m = 1, imaginary coupling, phi = (i/beta) sum_j alpha_j log tau_j):
    tau_j tau_j'' - tau_j'^2 = n_j ( tau_j^2 - prod_{i~j} tau_i^{|C_ij|} )
Linearized (fluctuation e^{kappa x - i w t} P_j(E), w^2 = m_b^2 - kappa^2):
    sum_{k+l=K} c_{j,l} p_{j,k} [ (kappa + k mu - l mu)^2 + m_b^2 - kappa^2 ]
        = n_j [ 2 sum c_{j,l} p_{j,k} - sum_{i~j} (p_i prod_{l~j, l!=i} tau_l)_K ]
"""
import numpy as np
import itertools

def affine(alg, r):
    """affine Cartan matrix C (nodes 0..r), Kac labels n, list of diagram automorphisms (permutations)."""
    edges = []
    if alg == 'a':
        N = r + 1
        edges = [(j, (j+1) % N) for j in range(N)]
        n = [1]*N
        rot = [(j+1) % N for j in range(N)]
        autos = [rot]
    elif alg == 'd':
        # nodes 0,1 -> 2 ; chain 2..r-2 ; r-1, r -> r-2
        edges = [(0, 2), (1, 2)] + [(j, j+1) for j in range(2, r-2)] + [(r-2, r-1), (r-2, r)]
        if r == 4:
            edges = [(0, 2), (1, 2), (3, 2), (4, 2)]
        n = [1, 1] + [2]*(r-3) + [1, 1]
        s1 = list(range(r+1)); s1[0], s1[1] = 1, 0; s1[r-1], s1[r] = r, r-1          # (01)(r-1 r)
        rev = [0]*(r+1)
        rev[0], rev[1], rev[r-1], rev[r] = r-1, r, 0, 1
        for j in range(2, r-1): rev[j] = r - j
        if r % 2 == 0:
            autos = [s1, rev]                       # special automorphisms Z2 x Z2
        else:
            rho = [0]*(r+1)                         # Z4: 0->r->1->r-1->0, j->r-j
            rho[0], rho[r], rho[1], rho[r-1] = r, 1, r-1, 0
            for j in range(2, r-1): rho[j] = r - j
            autos = [rho]
    elif alg == 'e':
        if r == 6:
            # Bourbaki: 1-3-4-5-6, 2-4, 0-2
            edges = [(1, 3), (3, 4), (4, 5), (5, 6), (2, 4), (0, 2)]
            n = [1, 1, 2, 2, 3, 2, 1]
            rot = [6, 0, 5, 1, 4, 3, 2]   # 0->6? define 1->6->0->1, 3->5->2->3
            rot = [None]*7
            rot[1], rot[6], rot[0] = 6, 0, 1
            rot[3], rot[5], rot[2] = 5, 2, 3
            rot[4] = 4
            autos = [rot]
        elif r == 7:
            edges = [(0, 1), (1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (2, 4)]
            n = [1, 2, 2, 3, 4, 3, 2, 1]
            ref = [7, 6, 2, 5, 4, 3, 1, 0]
            autos = [ref]
        elif r == 8:
            edges = [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 0), (2, 4)]
            n = [1, 2, 3, 4, 6, 5, 4, 3, 2]
            autos = []
    N = r + 1
    C = 2*np.eye(N)
    for i, j in edges:
        C[i, j] -= 1; C[j, i] -= 1
    if alg == 'a' and r == 1:
        C = np.array([[2., -2.], [-2., 2.]])
    n = np.array(n, float)
    assert np.allclose(C @ n, 0), "Kac labels wrong"
    for p in autos:
        P = np.zeros((N, N)); P[p, range(N)] = 1
        assert np.allclose(P @ C @ P.T, C), "automorphism wrong"
    return C, n, autos

def perm_matrix(p):
    N = len(p); P = np.zeros((N, N)); P[p, range(N)] = 1   # (P v)_{p(j)} = v_j
    return P

def species(alg, r, tol=1e-8):
    """eigenvectors of N C (m=1), split degenerate eigenspaces by the diagram automorphisms.
    Returns list of (mass, vector) for nonzero eigenvalues, and C, n."""
    C, n, autos = affine(alg, r)
    NC = np.diag(n) @ C
    w, V = np.linalg.eig(NC)
    order = np.argsort(w.real)
    w, V = w[order].real, V[:, order]
    out = []
    i = 0
    while i < len(w):
        j = i
        while j+1 < len(w) and abs(w[j+1] - w[i]) < 1e-7: j += 1
        if w[i] > 1e-9:
            S = V[:, i:j+1]
            Q, _ = np.linalg.qr(S)
            if Q.shape[1] > 1:
                # diagonalize random combination of automorphisms restricted to the eigenspace
                rng = np.random.default_rng(1)
                Mx = sum(rng.normal()*perm_matrix(p) for p in autos) if autos else None
                Pm = [perm_matrix(p) for p in autos]
                # build commuting combination: use products of the automorphisms (group generated)
                G = [np.eye(len(n))]
                for _ in range(6):
                    G = G + [g @ P for g in G for P in Pm]
                    # dedupe
                    U = []
                    for g in G:
                        if not any(np.allclose(g, u) for u in U): U.append(g)
                    G = U
                Qp = np.linalg.pinv(Q)
                # find a combination that is diagonalizable with distinct eigenvalues on the space
                coeffs = rng.normal(size=len(G)) + 1j*rng.normal(size=len(G))
                Mr = sum(c*(Qp @ g @ Q) for c, g in zip(coeffs, G))
                ww, U = np.linalg.eig(Mr)
                for k in range(U.shape[1]):
                    out.append((np.sqrt(w[i]), Q @ U[:, k]))
            else:
                out.append((np.sqrt(w[i]), Q[:, 0]))
        i = j + 1
    return out, C, n, autos

def nbrs(C):
    N = len(C)
    return [[(i, int(round(-C[i, j]))) for i in range(N) if i != j and C[i, j] < -0.5] for j in range(N)]

def polymul(a, b, K):
    out = np.zeros(K+1, complex)
    for i in range(min(len(a), K+1)):
        for j in range(min(len(b), K+1-i)):
            out[i+j] += a[i]*b[j]
    return out

def soliton(C, n, mu, delta, Kmax=None):
    """coefficients c[j,k] of tau_j = sum_k c_{jk} E^k, E = e^{mu x}."""
    if Kmax is None: Kmax = int(2*max(n)) + 2
    N = len(n); NC = np.diag(n) @ C
    nb = nbrs(C)
    c = np.zeros((N, Kmax+1), complex)
    c[:, 0] = 1; c[:, 1] = delta
    for K in range(2, Kmax+1):
        rhs = np.zeros(N, complex)
        for j in range(N):
            lhs_nl = 0.5*mu*mu*sum(c[j, p]*c[j, K-p]*(2*p-K)**2 for p in range(1, K))
            sq = sum(c[j, p]*c[j, K-p] for p in range(1, K))
            prod = np.zeros(K+1, complex); prod[0] = 1
            for i, mult in nb[j]:
                for _ in range(mult):
                    t = c[i, :K+1].copy(); t[K] = 0      # exclude linear term at order K
                    prod = polymul(prod, c[i, :K+1], K)
            # prod includes linear terms c_{i,K}; remove them: coefficient of E^K linear in c_{.,K} is sum_i mult c_{iK}
            prodK = prod[K] - sum(mult*c[i, K] for i, mult in nb[j])   # c[:,K] currently 0 anyway
            rhs[j] = n[j]*(sq - prodK) - lhs_nl
        Mk = mu*mu*K*K*np.eye(N) - NC
        # resonance K mu = m_b makes Mk singular; take the minimal-norm solution (no extra soliton added)
        c[:, K] = np.linalg.lstsq(Mk, rhs, rcond=1e-10)[0]
    return c

def degree(c, tol=1e-9):
    s = np.max(np.abs(c[:, :c.shape[1]//2 + 1]), axis=1)
    return [max([k for k in range(c.shape[1]) if abs(c[j, k]) > tol*max(1, s[j])] or [-1]) for j in range(c.shape[0])]

def fluct(C, n, c, mu, mb, p0, kappa, Kmax=None):
    """coefficients p[j,K] of the Hirota fluctuation sigma_j = e^{kappa x} sum_K p_{jK} E^K."""
    if Kmax is None: Kmax = c.shape[1] - 1
    N = len(n); NC = np.diag(n) @ C
    nb = nbrs(C)
    p = np.zeros((N, Kmax+1), complex)
    p[:, 0] = p0
    lam = mb*mb - kappa*kappa
    def br(k, l): return (kappa + k*mu - l*mu)**2 + lam
    for K in range(1, Kmax+1):
        rhs = np.zeros(N, complex)
        for j in range(N):
            # LHS terms other than (k=K, l=0)
            lhs = sum(c[j, l]*p[j, K-l]*br(K-l, l) for l in range(1, K+1))
            # RHS: n_j [2 sum_{l>=1} c_l p_{K-l} - sum_i (p_i prod others)_K, excluding p_{i,K} * 1 terms]
            two = 2*sum(c[j, l]*p[j, K-l] for l in range(1, K+1))
            tot = 0
            plist = []
            for i, mult in nb[j]:
                plist += [i]*mult
            for idx, i in enumerate(plist):
                prod = np.zeros(K+1, complex); prod[0] = 1
                for idx2, i2 in enumerate(plist):
                    if idx2 != idx: prod = polymul(prod, c[i2, :K+1], K)
                # (p_i * prod)_K excluding p_{i,K}*prod_0
                tot += sum(p[i, k]*prod[K-k] for k in range(0, K))
            rhs[j] = n[j]*(two - tot) - lhs
        Mk = (br(K, 0))*np.eye(N) - NC
        p[:, K] = np.linalg.solve(Mk, rhs)
    return p

def Xfactor(C, n, c, mu, mb, p0, z, Kmax=None):
    kappa = z*mb
    p = fluct(C, n, c, mu, mb, p0, kappa, Kmax)
    D = degree(c)
    Dp = degree(p, tol=1e-8)
    ratio = np.array([p[j, D[j]]/c[j, D[j]] for j in range(len(n))])
    # least squares X: ratio = X p0
    X = np.vdot(p0, ratio)/np.vdot(p0, p0)
    resid = np.linalg.norm(ratio - X*p0)/max(1e-300, np.linalg.norm(ratio))
    return X, resid, D, Dp, p
