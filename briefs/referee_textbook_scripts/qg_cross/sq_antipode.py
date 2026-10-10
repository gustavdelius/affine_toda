"""Item B: the S^2 shift V(x)** = V(x t) and the crossing point V(x)* = V(x x_c) in the homogeneous
gradation, for U_q(sl2^), U_q(sl3^) (vector), U_q(c_2^(1)) (vector, untwisted non-simply-laced) and
U_q(d_3^(2)) (spinor, twisted).  Generic q.  Two coproducts:
  ch.7 : Delta(e)=e(x)k+1(x)e, Delta(f)=f(x)1+k^-1(x)f ; S(e)=-e k^-1, S(f)=-k f, S(k)=k^-1
  VI   : Delta(e)=e(x)1+k(x)e, Delta(f)=f(x)k^-1+1(x)f ; S(e)=-k^-1 e, S(f)=-f k
Normalisation: q_i = Q^{(alpha_i,alpha_i)/2}; we report t as a power of Q for the stated normalisation."""
import numpy as np, itertools
Qs = 0.83*np.exp(0.61j)          # generic complex Q (not a root of unity)
def Emat(N, a, b):
    M = np.zeros((N, N), complex); M[a, b] = 1; return M
def check(rep, cartan, qi, x=1.3+0.2j):
    e, f, k = rep(x); err = 0; nodes = list(e)
    for i in nodes:
        for j in nodes:
            c = e[i]@f[j] - f[j]@e[i]
            if i == j: c = c - (k[i] - np.linalg.inv(k[i]))/(qi[i] - 1/qi[i])
            err = max(err, abs(c).max(), abs(k[i]@e[j]@np.linalg.inv(k[i]) - qi[i]**cartan[i][j]*e[j]).max())
            if i != j:
                m = 1 - cartan[i][j]
                qn = lambda a, qq: (qq**a - qq**-a)/(qq - 1/qq)
                fac = lambda a, qq: np.prod([qn(t, qq) for t in range(1, a+1)]) if a else 1.0
                for X in (e, f):
                    S = sum((-1)**s*fac(m, qi[i])/(fac(s, qi[i])*fac(m-s, qi[i]))*np.linalg.matrix_power(X[i], m-s)@X[j]@np.linalg.matrix_power(X[i], s) for s in range(m+1))
                    err = max(err, abs(S).max())
    return err
def gens_dual(rep, x, conv, power):
    """power=1: dual rep a -> pi(S(a))^T ; power=2: double dual a -> pi(S^2(a)).
    ch.7: S(e)=-e k^-1, S(f)=-k f  => S^2(e)=k e k^-1, S^2(f)=k f k^-1
    VI  : S(e)=-k^-1 e, S(f)=-f k  => S^2(e)=k^-1 e k, S^2(f)=k^-1 f k"""
    e, f, k = rep(x); out = []
    for i in e:
        ki = np.linalg.inv(k[i])
        if power == 1:
            Se, Sf = (-e[i]@ki, -k[i]@f[i]) if conv == '7' else (-ki@e[i], -f[i]@k[i])
            out += [Se.T, Sf.T, ki.T]
        else:
            A, Ai = (k[i], ki) if conv == '7' else (ki, k[i])
            out += [A@e[i]@Ai, A@f[i]@Ai, k[i]]
    return out
def plain(rep, x):
    e, f, k = rep(x); out = []
    for i in e: out += [e[i], f[i], k[i]]
    return out
def equivalent(A, B):
    D = A[0].shape[0]
    M = np.vstack([np.kron(np.eye(D), a.T) - np.kron(b, np.eye(D)) for a, b in zip(A, B)])   # C A = B C  (row-major vec)
    s = np.linalg.svd(M, compute_uv=False); return s[-1] < 1e-9*s[0]
def scan(rep, conv, power, roots=(1,), lmax=16, Q=Qs, x0=0.77+0.31j):
    hits = []
    for r in roots:
        for l in range(-lmax, lmax+1):
            t = r*Q**l
            if equivalent(gens_dual(rep, x0, conv, power), plain(rep, x0*t)): hits.append((r, l))
    return hits
def lab(h):
    out = []
    for r, l in h:
        ph = np.angle(r)/(2*np.pi)
        out.append(('' if abs(ph) < 1e-9 else ('-' if abs(abs(ph)-0.5) < 1e-9 else f'e^(2pi i {ph:.3f})')) + f'Q^{l}')
    return out
# ---------------- U_q(sl2^), spin 1/2 (eq-evaluation-rep), q = Q
def sl2(x, Q=Qs):
    E = Emat(2, 0, 1); F = Emat(2, 1, 0); K = np.diag([Q, 1/Q])
    return {1: E, 0: x*F}, {1: F, 0: E/x}, {1: K, 0: np.linalg.inv(K)}
cart_sl2 = {1: {1: 2, 0: -2}, 0: {1: -2, 0: 2}}
# ---------------- U_q(sl3^) vector, q = Q, h^vee = 3
def sl3(x, Q=Qs):
    e = {1: Emat(3, 0, 1), 2: Emat(3, 1, 2), 0: x*Emat(3, 2, 0)}
    f = {0: Emat(3, 0, 2)/x, 1: Emat(3, 1, 0), 2: Emat(3, 2, 1)}
    k = {1: np.diag([Q, 1/Q, 1]), 2: np.diag([1, Q, 1/Q]), 0: np.diag([1/Q, 1, Q])}
    return e, f, k
cart_sl3 = {0: {0: 2, 1: -1, 2: -1}, 1: {0: -1, 1: 2, 2: -1}, 2: {0: -1, 1: -1, 2: 2}}
# ---------------- U_q(c_2^(1)) vector; book/Kac normalisation long^2=2: q_long = Q, q_short = Q^(1/2); use Q = P^2
def c2(x, P=np.sqrt(Qs)):
    Q = P**2
    e = {1: Emat(4, 0, 1) + Emat(4, 2, 3), 2: Emat(4, 1, 2), 0: x*Emat(4, 3, 0)}
    f = {1: Emat(4, 1, 0) + Emat(4, 3, 2), 2: Emat(4, 2, 1), 0: Emat(4, 0, 3)/x}
    k = {1: np.diag([P, 1/P, P, 1/P]), 2: np.diag([1, Q, 1/Q, 1]), 0: np.diag([1/Q, 1, 1, Q])}
    return e, f, k
cart_c2 = {0: {0: 2, 1: -1, 2: 0}, 1: {0: -2, 1: 2, 2: -2}, 2: {0: 0, 1: -1, 2: 2}}
# ---------------- U_q(d_3^(2)) spinor (Part VI's c_2 soliton rep): Kac normalisation, short^2=2: q_short=Q, q_long=Q^2
sp_, sm_ = Emat(2, 0, 1), Emat(2, 1, 0)
def d32(x, Q=Qs):
    n = 2; W = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n)])
    site = lambda op, i: np.kron(op, np.eye(2)) if i == 0 else np.kron(np.eye(2), op)
    e = {1: site(sp_, 0)@site(sm_, 1), 2: site(sp_, 1), 0: x*site(sm_, 0)}
    f = {1: e[1].T.copy(), 2: e[2].T.copy(), 0: site(sp_, 0)/x}
    k = {1: np.diag(Q**(2*(W[:, 0]-W[:, 1]))), 2: np.diag(Q**(2*W[:, 1])), 0: np.diag(Q**(-2*W[:, 0]))}
    return e, f, k
cart_d32 = {0: {0: 2, 1: -2, 2: 0}, 1: {0: -1, 1: 2, 2: -1}, 2: {0: 0, 1: -2, 2: 2}}
cases = [('U_q(sl2^) spin 1/2       [h^vee=2, q_i=Q]', sl2, cart_sl2, {0: Qs, 1: Qs}, 2),
         ('U_q(sl3^) vector         [h^vee=3, q_i=Q]', sl3, cart_sl3, {0: Qs, 1: Qs, 2: Qs}, 3),
         ('U_q(c_2^(1)) vector      [h^vee=3, q_long=Q, q_short=Q^1/2]', c2, cart_c2, {0: Qs, 1: np.sqrt(Qs), 2: Qs}, 3),
         ('U_q(d_3^(2)) spinor      [h^vee=4 (Kac), q_short=Q, q_long=Q^2]', d32, cart_d32, {0: Qs, 1: Qs**2, 2: Qs}, 4)]
for name, rep, cart, qi, hv in cases:
    print(name, f': relations residual {check(rep, cart, qi):.1e}')
    for conv in ('7', 'VI'):
        dd = scan(rep, conv, 2, roots=(1, -1, 1j, -1j))
        d1 = scan(rep, conv, 1, roots=(1, -1, 1j, -1j)) if name.split()[0] != 'U_q(sl3^)' else []
        print(f'   coproduct {conv:2s}:  V(x)** = V(x t), t = {lab(dd)}   (prediction Q^{"+" if conv=="7" else "-"}{2*hv});'
              f'  V(x)* = V(x t): t = {lab(d1) if d1 else "none (V* not a shift of V)"}')
