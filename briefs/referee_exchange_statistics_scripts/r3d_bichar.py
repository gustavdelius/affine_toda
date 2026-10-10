# Coproduct change by a bicharacter F = exp(i pi B(mu_a, mu_b)), B real bilinear: Rc -> F Rc F^{-1}.
# Claim (derived): T_1' = e^{i pi Phi} e^{i pi A(Q, mu_1^out)} T_1 e^{-i pi Phi}, Phi = sum_{i<k} B(w_i, w_k), A = B - B^T.
import numpy as np, itertools
rng = np.random.default_rng(2)
W = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [-1, 1]]); d = len(W)
def rand_cons():
    R = np.zeros((d*d, d*d), complex)
    for a, b, c, e in itertools.product(range(d), repeat=4):
        if np.array_equal(W[a]+W[b], W[c]+W[e]): R[c*d+e, a*d+b] = rng.normal()+1j*rng.normal()
    return R
P = np.zeros((d*d, d*d))
for a in range(d):
    for b in range(d): P[b*d+a, a*d+b] = 1
B = rng.normal(size=(2, 2)); A = B - B.T
F = np.diag([np.exp(1j*np.pi*W[a]@B@W[b]) for a in range(d) for b in range(d)])
def embed(M, j, N):   # M on factors (0, j) of N, ordinary
    out = np.zeros((d**N, d**N), complex); basis = list(itertools.product(range(d), repeat=N)); ix = {s: i for i, s in enumerate(basis)}
    for s in basis:
        for c, e in itertools.product(range(d), repeat=2):
            amp = M[c*d+e, s[0]*d+s[j]]
            if amp != 0:
                t = list(s); t[0], t[j] = c, e; out[ix[tuple(t)], ix[s]] += amp
    return out, basis
for N in (2, 3, 4):
    Rcs = [rand_cons() for _ in range(N)]
    T = np.eye(d**N, dtype=complex); Tp = np.eye(d**N, dtype=complex)
    for j in range(1, N):
        Mj, basis = embed(P@Rcs[j], j, N); Mjp, _ = embed(P@F@Rcs[j]@np.linalg.inv(F), j, N)
        T = Mj@T; Tp = Mjp@Tp
    Phi = np.diag([np.exp(1j*np.pi*sum(W[s[i]]@B@W[s[k]] for i in range(N) for k in range(i+1, N))) for s in basis])
    Q = [sum(W[x] for x in s) for s in basis]
    Tw = np.diag([np.exp(1j*np.pi*(Q[i]@A@W[s[0]])) for i, s in enumerate(basis)])
    print(f"N={N}: |T' - Phi Tw T Phi^-1| = {np.abs(Tp - Phi@Tw@T@np.linalg.inv(Phi)).max():.1e}")
