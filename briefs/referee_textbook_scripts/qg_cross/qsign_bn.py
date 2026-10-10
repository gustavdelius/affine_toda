"""Same crossing-sign test for the a_{2n-1}^(2) solitons: U_q(b_n^(1)) spinor (Part VI, sec-representation-r-matrix):
node 0 long, e_0 = x s1^- s2^-, f_0 = s1^+ s2^+/x, k_0 = q^{-2(w1+w2)}, q_0 = q^2;
T = (2n-1) w + n - 1, zeros a_j = (2j-1) w + (j-1), q = s e^{-i pi w}."""
import numpy as np, itertools, sys
from rsolve import intertwiner
omega = 2.37
sp_, sm_, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def make(n, q):
    def site(op, i):
        mats = [I2]*n; mats[i] = op; out = mats[0]
        for m in mats[1:]: out = np.kron(out, m)
        return out
    W = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n)])
    def rep(x):
        e, f, k, qi = {}, {}, {}, {}
        for i in range(1, n):
            e[i] = site(sp_, i-1) @ site(sm_, i); f[i] = e[i].T.copy(); k[i] = np.diag(q**(2*(W[:, i-1] - W[:, i]))).astype(complex); qi[i] = q**2
        e[n] = site(sp_, n-1); f[n] = e[n].T.copy(); k[n] = np.diag(q**(2*W[:, n-1])).astype(complex); qi[n] = q
        e[0] = x*site(sm_, 0)@site(sm_, 1); f[0] = site(sp_, 0)@site(sp_, 1)/x; k[0] = np.diag(q**(-2*(W[:, 0] + W[:, 1]))).astype(complex); qi[0] = q**2
        return e, f, k, qi
    return rep, W
def analyse(n, s):
    q = s*np.exp(-1j*np.pi*omega); D = 2**n; rep, W = make(n, q)
    P = np.zeros((D*D, D*D))
    for i in range(D):
        for j in range(D): P[j*D+i, i*D+j] = 1
    def Rh(x):
        Rc = intertwiner(rep(x), rep(1.0), range(n+1), W, W)[0]; return Rc/Rc[0, 0]
    def dual(x):
        e, f, k, _ = rep(x); out = []
        for i in range(n+1):
            ki = np.linalg.inv(k[i]); out += [(-ki@e[i]).T, (-f[i]@k[i]).T, ki.T]
        return out
    def plain(x):
        e, f, k, _ = rep(x); out = []
        for i in range(n+1): out += [e[i], f[i], k[i]]
        return out
    x0 = 0.83+0.29j; found = []
    for ph in (1, -1):
        for l in range(-4*n, 4*n+1):
            t = ph*q**l
            M = np.vstack([np.kron(np.eye(D), A.T) - np.kron(B, np.eye(D)) for A, B in zip(dual(x0), plain(x0*t))])
            sv = np.linalg.svd(M, compute_uv=False)
            if sv[-1] < 1e-9*sv[0]: found.append((ph, l, t, np.linalg.svd(M)[2][-1].conj().reshape(D, D)))
    ph, l, tc, C = found[0]; Ci = np.linalg.inv(C)
    def pt1(M): return M.reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(D*D, D*D)
    def cu(x):
        R = P @ Rh(x); lhs = pt1(np.linalg.inv(R)); rhs = np.kron(Ci, np.eye(D)) @ (P @ Rh(x*tc)) @ np.kron(C, np.eye(D))
        c = np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel()); return c, np.abs(lhs - c*rhs).max()/np.abs(lhs).max()
    T = (2*n-1)*omega + n - 1; A = [(2*j-1)*omega + (j-1) for j in range(1, n+1)]
    c_of = lambda t: np.prod([np.sin(np.pi*(t - a)) for a in A])
    ks = []
    for t in (0.37+0.21j, 1.3-0.4j):
        x = np.exp(2j*np.pi*t); c, res = cu(1/x); ks.append(((c_of(T - t)/c_of(t))/c, res))
    return [(p_, l_) for p_, l_, _, _ in found], abs(tc - np.exp(2j*np.pi*T)) < 1e-9, ks
for n in [int(a) for a in sys.argv[1:]] or [3]:
    for s in (1, -1):
        found, xc_ok, ks = analyse(n, s)
        print(f"b_{n}^(1) spinor, q={'+' if s>0 else '-'}e^(-i pi w): x_c = {['%sq^%d' % ('+' if p>0 else '-', l) for p, l in found]}; "
              f"x_c = e^(2 pi i T): {xc_ok}; kappa = {[complex(np.round(k, 9)) for k, _ in ks]} (residual {max(r for _, r in ks):.1e})")
