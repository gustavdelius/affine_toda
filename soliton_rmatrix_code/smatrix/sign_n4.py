import numpy as np, itertools, time
from rsolve import intertwiner
n4 = 4; omega = 2.37; q = np.exp(-1j*np.pi*omega); D = 16
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def site(op, i):
    mats = [I2]*n4; mats[i] = op; out = mats[0]
    for m in mats[1:]: out = np.kron(out, m)
    return out
W4 = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n4)])
def rep(x):
    e, f, k, qi = {}, {}, {}, {}
    for i in range(1, n4):
        e[i] = site(sp, i-1) @ site(sm, i); f[i] = e[i].T.copy(); k[i] = np.diag(q**(2*(W4[:, i-1] - W4[:, i]))).astype(complex); qi[i] = q**2
    e[n4] = site(sp, n4-1); f[n4] = e[n4].T.copy(); k[n4] = np.diag(q**(2*W4[:, n4-1])).astype(complex); qi[n4] = q
    e[0] = x*site(sm, 0); f[0] = site(sp, 0)/x; k[0] = np.diag(q**(-2*W4[:, 0])).astype(complex); qi[0] = q
    return e, f, k, qi
br = lambda x, l, s: (x - s*q**l)/(1 - s*x*q**l)
def rho(x): 
    out = [1.0]
    for j in range(1, n4+1): out.append(out[-1]*br(x, 2*j, (-1)**(j+1)))
    return out
t0 = time.time(); x0 = 1.37+0.41j
Rc = intertwiner(rep(x0), rep(1.0), range(n4+1), W4, W4)[0]; Rc = Rc/Rc[0, 0]; lam = rho(x0)
Pk = []
for k in range(n4+1):
    M = np.eye(D*D, dtype=complex)
    for j in range(n4+1):
        if j != k: M = M @ (Rc - lam[j]*np.eye(D*D))/(lam[k] - lam[j])
    Pk.append(M)
Rx = lambda x: sum(r*Pm for r, Pm in zip(rho(x), Pk))
print(f"n=4 closed form check {np.abs(Rx(x0) - Rc).max():.1e}  [{time.time()-t0:.0f}s]")
P = np.zeros((D*D, D*D))
for i in range(D):
    for j in range(D): P[j*D+i, i*D+j] = 1
def dual(x):
    e, f, k, _ = rep(x); out = []
    for i in range(n4+1):
        ki = np.linalg.inv(k[i]); out += [(-ki@e[i]).T, (-f[i]@k[i]).T, ki.T]
    return out
def plain(x):
    e, f, k, _ = rep(x); out = []
    for i in range(n4+1): out += [e[i], f[i], k[i]]
    return out
xa = 0.83+0.29j
for ph in (1, -1):
    for l in (-8, 8):
        tt = ph*q**l
        M = np.vstack([np.kron(np.eye(D), A.T) - np.kron(B, np.eye(D)) for A, B in zip(dual(xa), plain(xa*tt))])
        s = np.linalg.svd(M, compute_uv=False)
        if s[-1] < 1e-9*s[0]:
            C = np.linalg.svd(M)[2][-1].conj().reshape(D, D); tc = tt; print(f"n=4 crossing point x(i pi) = {'+' if ph>0 else '-'}q^{l}")
Ci = np.linalg.inv(C)
def pt1(M): return M.reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(D*D, D*D)
def cu(x):
    R = P @ Rx(x); lhs = pt1(np.linalg.inv(R)); rhs = np.kron(Ci, np.eye(D)) @ (P @ Rx(x*tc)) @ np.kron(C, np.eye(D))
    c = np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel()); return c, np.abs(lhs - c*rhs).max()/np.abs(lhs).max()
T = n4*omega + (n4 - 1)/2; A = [j*omega + ((j+1) % 2)/2 for j in range(1, n4+1)]
c_of = lambda t: np.prod([np.sin(np.pi*(t - a)) for a in A])
for t in (0.37+0.21j, 1.3-0.4j):
    x = np.exp(2j*np.pi*t); c, res = cu(1/x)
    r = (c_of(T - t)/c_of(t))/c
    print(f"n=4: universal-crossing residual {res:.1e};  [c(T-t)/c(t)] / c_u(1/x) = {r.real:+.6f}{r.imag:+.6f}i")
