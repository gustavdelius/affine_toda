import numpy as np, itertools
from rsolve import intertwiner
# --- n = 2 spinor representation of U_q(d_3^(2)) at unimodular q, with the same conventions
n2 = 2; omega = 2.37; q = np.exp(-1j*np.pi*omega)
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def site(op, i):
    mats = [I2]*n2; mats[i] = op; out = mats[0]
    for m in mats[1:]: out = np.kron(out, m)
    return out
W2 = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n2)])
def rep(x):
    e, f, k, qi = {}, {}, {}, {}
    e[1] = site(sp, 0) @ site(sm, 1); f[1] = e[1].T.copy(); k[1] = np.diag(q**(2*(W2[:, 0] - W2[:, 1]))).astype(complex); qi[1] = q**2
    e[2] = site(sp, 1); f[2] = e[2].T.copy(); k[2] = np.diag(q**(2*W2[:, 1])).astype(complex); qi[2] = q
    e[0] = x*site(sm, 0); f[0] = site(sp, 0)/x; k[0] = np.diag(q**(-2*W2[:, 0])).astype(complex); qi[0] = q
    return e, f, k, qi
D = 4; P = np.zeros((16, 16))
for i in range(D):
    for j in range(D): P[j*D+i, i*D+j] = 1
def Rh(x):
    Rc = intertwiner(rep(x), rep(1.0), range(n2+1), W2, W2)[0]; return Rc/Rc[0, 0]
def dual(x):
    e, f, k, _ = rep(x); out = []
    for i in range(n2+1):
        ki = np.linalg.inv(k[i]); out += [(-ki@e[i]).T, (-f[i]@k[i]).T, ki.T]
    return out
def plain(x):
    e, f, k, _ = rep(x); out = []
    for i in range(n2+1): out += [e[i], f[i], k[i]]
    return out
x0 = 0.83+0.29j
for ph in (1, -1):
    for l in range(-6, 7):
        t = ph*q**l
        M = np.vstack([np.kron(np.eye(D), A.T) - np.kron(B, np.eye(D)) for A, B in zip(dual(x0), plain(x0*t))])
        s = np.linalg.svd(M, compute_uv=False)
        if s[-1] < 1e-9*s[0]:
            C = np.linalg.svd(M)[2][-1].conj().reshape(D, D); tc = t; print(f"n=2: crossing point x(i pi) = {'+' if ph>0 else '-'}q^{l}")
Ci = np.linalg.inv(C)
def pt1(M): return M.reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(D*D, D*D)
def cu(x):
    R = P @ Rh(x); lhs = pt1(np.linalg.inv(R)); rhs = np.kron(Ci, np.eye(D)) @ (P @ Rh(x*tc)) @ np.kron(C, np.eye(D))
    return np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel())
# n=2 calibration: mu = 4 omega + 1, T = 2 omega + 1/2 ; c(t) = sin pi(t - w) sin pi(t - 2w - 1/2)
T = 2*omega + 0.5; A = [omega, 2*omega + 0.5]
c_of = lambda t: np.prod([np.sin(np.pi*(t - a)) for a in A])
for t in (0.37+0.21j, 1.3-0.4j):
    x = np.exp(2j*np.pi*t)
    ratio = (c_of(T - t)/c_of(t))/cu(1/x)
    print(f"n=2: [c(T-t)/c(t)] / c_u(1/x) = {ratio.real:+.6f}{ratio.imag:+.6f}i   at t = {t}")
