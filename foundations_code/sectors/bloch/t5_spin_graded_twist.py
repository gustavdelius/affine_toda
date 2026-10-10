"""Spinor: (i) graded transfer matrix = Bloch twist exp(i pi (G-1) g_1) times bosonic; (ii) alpha=0 breaking by total-weight sector, N=2."""
from spin_common import *
for n, alg in ((2, 'c'), (3, 'c'), (3, 'a')):
    D, Wn, setup = make(n, alg)
    g = np.array([bin(i).count('1') % 2 for i in range(D)])
    Pg = np.zeros((D*D, D*D))
    for i in range(D):
        for j in range(D): Pg[j*D+i, i*D+j] = (-1)**(g[i]*g[j])
    def transfer_graded(Rg, ls):
        Nn = len(ls) + 1; T = np.eye(D**Nn, dtype=complex)
        for j, l in enumerate(ls, start=1):
            v = T
            for k in range(j, 1, -1): v = apply_pair(Pg, v, k-1, k, D, Nn)
            v = apply_pair(Rg(np.exp(l)), v, 0, 1, D, Nn)
            for k in range(2, j+1): v = apply_pair(Pg, v, k-1, k, D, Nn)
            T = v
        return T
    om = 1.61; q = np.exp(-1j*np.pi*om); Rp = setup(q)
    Qd = np.diag([(-1.0)**(g[i]*g[j]) for i in range(D) for j in range(D)]); Rg = lambda x: Qd @ Rp(x)
    for Nn in (2, 3):
        if D**Nn > 600: continue
        gs = np.array([[g[i] for i in c] for c in itertools.product(range(D), repeat=Nn)])
        Om = np.diag((-1.0)**(gs[:, 0]*(gs.sum(1) - 1)))
        ls = np.random.default_rng(0).uniform(-3, 3, Nn - 1)
        err = np.abs(transfer_graded(Rg, ls) - Om @ transferN(Rp, D, ls)).max()
        print(f"{alg}_{n} spinor, N={Nn}: |T_graded - (-1)^((G-1) g_1) T_bos| = {err:.1e}")
    # alpha = 0, N = 2: breaking by total weight sector
    tot, first = weights_multi(Wn, 2)
    keys = sorted({tuple(t) for t in tot}, key=lambda t: (sum(abs(a) for a in t), t))
    lx = np.concatenate([np.linspace(0.005, 6, 300), np.linspace(6, 30, 60)])
    for om in (2.37, 1.61, 0.3013):   # 0.3 is a root of unity (q^20 = 1): rank-3 projectors degenerate
        q = np.exp(-1j*np.pi*om); Rp = setup(q); Rs = [Rp(np.exp(l)) for l in lx]; res = []
        for kk in keys:
            idx = np.where(np.all(np.isclose(tot, kk), axis=1))[0]
            d = max(np.abs(np.abs(np.linalg.eigvals(R[np.ix_(idx, idx)])) - 1).max() for R in Rs)
            res.append(f"Q={tuple(float(a) for a in kk)}(dim {len(idx)}): {d:.1e}")
        print(f"   omega={om}, alpha=0, N=2: " + "; ".join(res))
