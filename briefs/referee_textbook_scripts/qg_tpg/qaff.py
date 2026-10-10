"""Generic numerical tools: U_q(g^) relations and Jimbo's equations.

Conventions of ch. 7 (eq-coproduct, sec-uqghat):
  k_i e_j k_i^-1 = q_i^{a_ij} e_j, [e_i,f_j] = delta_ij (k_i-k_i^-1)/(q_i-q_i^-1),
  q_i = q^{alpha_i^2/2}, long roots alpha^2 = 2,
  Delta(k)=k(x)k, Delta(e)=e(x)k + 1(x)e, Delta(f)=f(x)1 + k^-1(x)f,
  homogeneous gradation: e_0 -> x e_0, f_0 -> f_0/x.
"""
import numpy as np
from math import comb


def qnum(n, q):
    return (q**n - q**(-n))/(q - 1/q)


def qfact(n, q):
    r = 1
    for j in range(1, n+1):
        r *= qnum(j, q)
    return r


def qbinom(n, s, q):
    return qfact(n, q)/(qfact(s, q)*qfact(n-s, q))


def check_relations(e, f, k, A, qi, tol=1e-9, verbose=False):
    """Return max residual over all U_q(g^) relations (incl. q-Serre with q_i)."""
    n = len(e); N = e[0].shape[0]; I = np.eye(N); res = {}
    for i in range(n):
        ki_inv = np.linalg.inv(k[i])
        for j in range(n):
            res[f'kk{i}{j}'] = np.abs(k[i]@k[j]-k[j]@k[i]).max()
            res[f'ke{i}{j}'] = np.abs(k[i]@e[j]@ki_inv - qi[i]**A[i][j]*e[j]).max()
            res[f'kf{i}{j}'] = np.abs(k[i]@f[j]@ki_inv - qi[i]**(-A[i][j])*f[j]).max()
            rhs = (k[i]-ki_inv)/(qi[i]-1/qi[i]) if i == j else 0*I
            res[f'ef{i}{j}'] = np.abs(e[i]@f[j]-f[j]@e[i]-rhs).max()
            if i != j:
                m = 1 - A[i][j]
                for X, nm in ((e, 'e'), (f, 'f')):
                    S = 0*I
                    for s in range(m+1):
                        S = S + (-1)**s*qbinom(m, s, qi[i])*np.linalg.matrix_power(X[i], m-s)@X[j]@np.linalg.matrix_power(X[i], s)
                    res[f'serre{nm}{i}{j}'] = np.abs(S).max()
    bad = {kk: v for kk, v in res.items() if v > tol}
    if verbose and bad:
        print('violated:', bad)
    return max(res.values()), bad


def grade(e, f, k, x, s):
    """Apply gradation s: e_i -> x^{s_i} e_i, f_i -> x^{-s_i} f_i."""
    return [x**s[i]*e[i] for i in range(len(e))], [x**(-s[i])*f[i] for i in range(len(e))], k


def coprod(gx, gy):
    (ex, fx, kx), (ey, fy, ky) = gx, gy
    N1 = ex[0].shape[0]; N2 = ey[0].shape[0]
    I1, I2 = np.eye(N1), np.eye(N2); out = []
    for i in range(len(ex)):
        out.append(np.kron(ex[i], kx[i]) + np.kron(I1, ey[i]))
        out.append(np.kron(fx[i], I2) + np.kron(np.linalg.inv(kx[i]), fy[i]))
        out.append(np.kron(kx[i], ky[i]))
    return out


def solve_jimbo(rep, x, y, s=None, rep2=None, return_dim=False):
    """Solve Rc Delta_{x,y}(a) = Delta_{y,x}(a) Rc on V(x)(x)V(y) (V=W)."""
    e, f, k = rep
    if s is None:
        s = [1] + [0]*(len(e)-1)
    gx = grade(e, f, k, x, s); gy = grade(e, f, k, y, s)
    A = coprod(gx, gy); B = coprod(gy, gx)
    M = A[0].shape[0]
    pairs = [(A[3*m], B[3*m]) for m in range(len(e))] + [(A[3*m+1], B[3*m+1]) for m in range(len(e))]
    kd = np.array([np.diag(A[3*i+2]) for i in range(len(e))]).T   # M x (r+1)
    null, sv, idx = intertwiner_space(pairs, kd, kd)
    null = null[:, 0]
    R = np.zeros((M, M), dtype=complex)
    for (u, v), n in idx.items():
        R[u, v] = null[n]
    if return_dim:
        return R, sv[-3:]
    return R, sv[-2:]


def P_flip(N):
    P = np.zeros((N*N, N*N))
    for i in range(N):
        for j in range(N):
            P[j*N+i, i*N+j] = 1
    return P


def intertwiner_space(pairs, kdL, kdR, nvec=1):
    """Solve R a = b R for all (a,b) in pairs, R restricted to entries (u,v) with
    equal k-weights kdL[u] == kdR[v] (k's diagonal). Sparse normal equations.
    Returns (lowest nvec null vectors as columns, singular values descending (last 3), index)."""
    import scipy.sparse as sps
    M = kdL.shape[0]
    keyL = [tuple(np.round(r, 10)) for r in kdL]; keyR = [tuple(np.round(r, 10)) for r in kdR]
    allowed = [(a, b) for a in range(M) for b in range(M) if keyL[a] == keyR[b]]
    idx = {ab: n for n, ab in enumerate(allowed)}
    G = None
    for a_, b_ in pairs:
        arow = [(np.nonzero(a_[c])[0], a_[c][np.nonzero(a_[c])[0]]) for c in range(M)]
        bcol = [(np.nonzero(b_[:, c])[0], b_[:, c][np.nonzero(b_[:, c])[0]]) for c in range(M)]
        ri, ci, va = [], [], []
        for (u, c), n in idx.items():
            vs, vals = arow[c]
            ri.extend(u*M+vs); ci.extend([n]*len(vs)); va.extend(vals)
        for (c, v), n in idx.items():
            us, vals = bcol[c]
            ri.extend(us*M+v); ci.extend([n]*len(us)); va.extend(-vals)
        L = sps.csr_matrix((va, (ri, ci)), shape=(M*M, len(allowed)))
        g = (L.conj().T @ L)
        G = g if G is None else G + g
    w, V = np.linalg.eigh(G.toarray())
    sv = np.sqrt(np.abs(w))
    return V[:, :nvec], sv[:3][::-1], idx
