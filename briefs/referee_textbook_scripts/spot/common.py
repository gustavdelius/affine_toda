import numpy as np, itertools
from rsolve import intertwiner
n = 3
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def site(op, i):
    mats = [I2]*n; mats[i] = op; out = mats[0]
    for m in mats[1:]: out = np.kron(out, m)
    return out
Ws = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n)])
def spin_rep(x, q):
    e, f, k, qi = {}, {}, {}, {}
    for i in range(1, n):
        e[i] = site(sp, i-1) @ site(sm, i); f[i] = e[i].T.copy()
        k[i] = np.diag(q**(2*(Ws[:, i-1] - Ws[:, i]))).astype(complex); qi[i] = q**2
    e[n] = site(sp, n-1); f[n] = e[n].T.copy(); k[n] = np.diag(q**(2*Ws[:, n-1])).astype(complex); qi[n] = q
    e[0] = x*site(sm, 0); f[0] = site(sp, 0)/x; k[0] = np.diag(q**(-2*Ws[:, 0])).astype(complex); qi[0] = q
    return e, f, k, qi
P = np.zeros((64, 64))
for i in range(8):
    for j in range(8): P[j*8+i, i*8+j] = 1
def bracket(x, l, s, q): return (x - s*q**l)/(1 - s*x*q**l)
class Spinor33:
    """exact closed form R_33(x) = sum_k rho_k P_k ; projectors extracted once per q"""
    def __init__(self, q):
        self.q = q
        x0 = 1.37 + 0.41j
        Rc, null, _, _ = intertwiner(spin_rep(x0, q), spin_rep(1.0, q), range(n+1), Ws, Ws)
        Rc = Rc/Rc[0, 0]
        lam = self.rho(x0)                      # eigenvalues for 35, 21, 7, 1
        self.Pk = []
        for k in range(4):
            M = np.eye(64, dtype=complex)
            for j in range(4):
                if j != k: M = M @ (Rc - lam[j]*np.eye(64))/(lam[k] - lam[j])
            self.Pk.append(M)
        self.check = np.abs(sum(lam[k]*self.Pk[k] for k in range(4)) - Rc).max()
        self.ranks = [int(round(np.trace(Pm).real)) for Pm in self.Pk]
    def rho(self, x):
        q = self.q; b2, b4, b6 = bracket(x, 2, 1, q), bracket(x, 4, -1, q), bracket(x, 6, 1, q)
        return [1.0, b2, b2*b4, b2*b4*b6]
    def R(self, x):
        r = self.rho(x); return sum(r[k]*self.Pk[k] for k in range(4))
