import numpy as np, itertools
# U_q(d_{n+1}^{(2)}) "vector" representation, dimension N = 2n+2.
# basis: w_1..w_n (0..n-1), w_0 (n), w_{-n}..w_{-1} (n+1..2n), s (2n+1)
def idx(n):
    I = {}
    for i in range(1, n+1): I[i] = i-1; I[-i] = 2*n+1-i
    I[0] = n; I['s'] = 2*n+1
    return I
def weights(n):
    I = idx(n); N = 2*n+2; W = np.zeros((N, n))
    for i in range(1, n+1): W[I[i], i-1] = 1; W[I[-i], i-1] = -1
    return W
def E(N, a, b):
    M = np.zeros((N, N), complex); M[a, b] = 1; return M
def rep(n, q, x):
    I = idx(n); N = 2*n+2; W = weights(n); qn = lambda k: (q**k - q**-k)/(q - 1/q)
    e, f, k, qi = {}, {}, {}, {}
    for i in range(1, n):                       # long roots, q_i = q^2
        e[i] = E(N, I[i], I[i+1]) + E(N, I[-(i+1)], I[-i]); f[i] = e[i].T.copy()
        h = W[:, i-1] - W[:, i]; k[i] = np.diag(q**(2*h)).astype(complex); qi[i] = q**2
    c = np.sqrt(qn(2))                           # short roots, q_i = q
    e[n] = c*(E(N, I[n], I[0]) + E(N, I[0], I[-n])); f[n] = e[n].T.copy()
    k[n] = np.diag(q**(2*W[:, n-1])).astype(complex); qi[n] = q
    e[0] = x*c*(E(N, I['s'], I[1]) + E(N, I[-1], I['s'])); f[0] = (c/x)*(E(N, I[1], I['s']) + E(N, I['s'], I[-1]))
    k[0] = np.diag(q**(-2*W[:, 0])).astype(complex); qi[0] = q
    return e, f, k, qi
def check_relations(n, q, x=1.7):
    e, f, k, qi = rep(n, q, x); err = 0
    nodes = range(n+1)
    for i in nodes:
        for j in nodes:
            comm = e[i]@f[j] - f[j]@e[i]
            tgt = (k[i] - np.linalg.inv(k[i]))/(qi[i] - 1/qi[i]) if i == j else 0*comm
            err = max(err, np.abs(comm - tgt).max())
            # k_i e_j k_i^-1 = q_i^{a_ij} e_j : extract a_ij
    # Cartan matrix from k-conjugation
    A = np.zeros((n+1, n+1))
    for i in nodes:
        for j in nodes:
            M = k[i]@e[j]@np.linalg.inv(k[i]); nz = np.abs(e[j]) > 1e-12
            r = (M[nz]/e[j][nz]); A[i, j] = np.round(np.log(r[0].real)/np.log(qi[i]), 6)
            err = max(err, np.abs(r - r[0]).max())
    # quantum Serre relations
    def qbin(m, r, qq):
        qn_ = lambda k: (qq**k - qq**-k)/(qq - 1/qq); fact = lambda k: np.prod([qn_(t) for t in range(1, k+1)]) if k else 1.0
        return fact(m)/(fact(r)*fact(m-r))
    for i in nodes:
        for j in nodes:
            if i == j: continue
            m = int(1 - A[i, j]); S = 0
            for r in range(m+1):
                S = S + (-1)**r*qbin(m, r, qi[i])*np.linalg.matrix_power(e[i], m-r)@e[j]@np.linalg.matrix_power(e[i], r)
            err = max(err, np.abs(S).max())
    return A, err
if __name__ == '__main__':
  for n in [2, 3, 4]:
    A, err = check_relations(n, 0.7)
    print(f"n={n}: Cartan matrix of the algebra acting on C^{2*n+2}:\n{A.astype(int)}\n   max violation of [e,f], k e k^-1 and quantum Serre relations: {err:.1e}")
