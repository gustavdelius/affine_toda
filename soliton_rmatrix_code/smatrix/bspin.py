"""spinor representation of U_q(b_n^(1)) (symmetry of imaginary a_{2n-1}^(2) Toda solitons)"""
import numpy as np, itertools, os, math
from rsolve import intertwiner
N = int(os.environ.get('NSPIN', '3')); omega = float(os.environ.get('OMEGA', '2.37'))
q = float(__import__("os").environ.get("QSIGN", "1"))*np.exp(-1j*np.pi*omega); D = 2**N
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def site(op, i):
    mats = [I2]*N; mats[i] = op; out = mats[0]
    for m in mats[1:]: out = np.kron(out, m)
    return out
Wn = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=N)])
def rep(x):
    e, f, k, qi = {}, {}, {}, {}
    for i in range(1, N):
        e[i] = site(sp, i-1) @ site(sm, i); f[i] = e[i].T.copy(); k[i] = np.diag(q**(2*(Wn[:, i-1] - Wn[:, i]))).astype(complex); qi[i] = q**2
    e[N] = site(sp, N-1); f[N] = e[N].T.copy(); k[N] = np.diag(q**(2*Wn[:, N-1])).astype(complex); qi[N] = q
    e[0] = x*site(sm, 0) @ site(sm, 1); f[0] = site(sp, 0) @ site(sp, 1)/x
    k[0] = np.diag(q**(-2*(Wn[:, 0] + Wn[:, 1]))).astype(complex); qi[0] = q**2
    return e, f, k, qi
# Cartan matrix of b_n^(1): node 0 attached to node 2 (n >= 3) or node 2 = short (n = 2)
roots = {0: -np.eye(N)[0] - np.eye(N)[1]}
for i in range(1, N): roots[i] = np.eye(N)[i-1] - np.eye(N)[i]
roots[N] = np.eye(N)[N-1]
def a_ij(i, j): return 2*np.dot(roots[i], roots[j])/np.dot(roots[i], roots[i])
def check(x=0.7+0.3j):
    e, f, k, qi = rep(x); err = 0
    for i in range(N+1):
        for j in range(N+1):
            comm = e[i] @ f[j] - f[j] @ e[i]
            tgt = (k[i] - np.linalg.inv(k[i]))/(qi[i] - 1/qi[i]) if i == j else 0*comm
            err = max(err, np.abs(comm - tgt).max())
            if i != j:
                a = int(round(a_ij(i, j))); m = 1 - a; Q = qi[i]
                def qb(n_): return (Q**n_ - Q**(-n_))/(Q - 1/Q)
                def qf(n_): return np.prod([qb(r) for r in range(1, n_+1)]) if n_ > 0 else 1.0
                ser = sum((-1)**r*qf(m)/(qf(r)*qf(m - r))*np.linalg.matrix_power(e[i], m - r) @ e[j] @ np.linalg.matrix_power(e[i], r) for r in range(m+1))
                err = max(err, np.abs(ser).max())
            kk = k[i] @ e[j] @ np.linalg.inv(k[i]); err = max(err, np.abs(kk - qi[i]**a_ij(i, j)*e[j]).max())
    return err
if __name__ == "__main__":
    print(f"U_q(b_{N}^(1)) spinor ({D}-dim): Cartan row 0 = {[int(round(a_ij(0, j))) for j in range(N+1)]};  relations + q-Serre max error {check():.1e}")
    x0 = 1.37+0.41j
    Rc, null, _, _ = intertwiner(rep(x0), rep(1.0), range(N+1), Ws := Wn, Wn); Rc = Rc/Rc[0, 0]
    ev = np.linalg.eigvals(Rc); groups = []
    for v in ev:
        for g in groups:
            if abs(g[0] - v) < 1e-7*max(1, abs(v)): g[1] += 1; break
        else: groups.append([v, 1])
    print("eigenvalue multiplicities:", sorted([m for _, m in groups], reverse=True), " (so(2n+1) Lambda^k dims:", [int(round(math.comb(2*N+1, k))) for k in range(N, -1, -1)], ")")
