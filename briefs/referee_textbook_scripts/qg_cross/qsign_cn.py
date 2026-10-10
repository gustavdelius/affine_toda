"""Crossing sign of the c_n^(1) spinor S-matrix (U_q(d_{n+1}^(2)) spinor) as a function
of the SIGN of q.  Part VI conventions: Delta(e)=e(x)1+k(x)e, Delta(f)=f(x)k^-1+1(x)f,
homogeneous gradation e_0 = x sigma_1^-, q = s*exp(-i pi omega), T = n omega + (n-1)/2,
zeros a_j = j omega + (j-1)/2.  kappa = [c(T-t)/c(t)] / c_u(1/x) is the number the
book calls the crossing-test sign ((-1)^n for s=+1)."""
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
        e[0] = x*site(sm_, 0); f[0] = site(sp_, 0)/x; k[0] = np.diag(q**(-2*W[:, 0])).astype(complex); qi[0] = q
        return e, f, k, qi
    return rep, W
def analyse(n, s):
    q = s*np.exp(-1j*np.pi*omega); D = 2**n
    rep, W = make(n, q)
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
        for l in range(-2*n-2, 2*n+3):
            t = ph*q**l
            M = np.vstack([np.kron(np.eye(D), A.T) - np.kron(B, np.eye(D)) for A, B in zip(dual(x0), plain(x0*t))])
            sv = np.linalg.svd(M, compute_uv=False)
            if sv[-1] < 1e-9*sv[0]:
                found.append((ph, l, t, np.linalg.svd(M)[2][-1].conj().reshape(D, D)))
    ph, l, tc, C = found[0]
    Ci = np.linalg.inv(C)
    def pt1(M): return M.reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(D*D, D*D)
    def cu(x):
        R = P @ Rh(x); lhs = pt1(np.linalg.inv(R)); rhs = np.kron(Ci, np.eye(D)) @ (P @ Rh(x*tc)) @ np.kron(C, np.eye(D))
        c = np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel())
        return c, np.abs(lhs - c*rhs).max()/np.abs(lhs).max()
    T = n*omega + (n-1)/2; A = [j*omega + (j-1)/2 for j in range(1, n+1)]
    c_of = lambda t: np.prod([np.sin(np.pi*(t - a)) for a in A])
    ks = []
    for t in (0.37+0.21j, 1.3-0.4j):
        x = np.exp(2j*np.pi*t); c, res = cu(1/x); ks.append(((c_of(T - t)/c_of(t))/c, res))
    # crossing point must equal x(i pi) = e^{2 pi i T}
    xc_ok = abs(tc - np.exp(2j*np.pi*T)) < 1e-9
    return q, [(p_, l_) for p_, l_, _, _ in found], xc_ok, ks, C, Rh
if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [1, 2, 3]
    for n in ns:
        for s in (+1, -1):
            q, found, xc_ok, ks, C, Rh = analyse(n, s)
            lab = ['%sq^%d' % ('+' if p_ > 0 else '-', l_) for p_, l_ in found]
            print(f"n={n} q={'+' if s>0 else '-'}e^(-i pi w): V*(x)=V(x t) for t in {lab}; x_c = e^(2 pi i T): {xc_ok}; "
                  f"kappa = {[complex(np.round(k_, 9)) for k_, _ in ks]} (universal-crossing residual {max(r for _, r in ks):.1e})")
