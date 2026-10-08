import numpy as np, itertools
from dtw import rep, weights
from rsolve import intertwiner
q, n = 0.7, 3
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)   # basis |1>=s=+1/2, |0>=s=-1/2? use index0=+1/2
# convention: index 0 = spin +1/2, index 1 = spin -1/2 ; sp raises (1->0)
def site(op, i):
    mats = [I2]*n; mats[i] = op; out = mats[0]
    for m in mats[1:]: out = np.kron(out, m)
    return out
def spin_weights():
    return np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n)])
def spin_rep(x):
    Wt = spin_weights(); e, f, k, qi = {}, {}, {}, {}
    for i in range(1, n):
        e[i] = site(sp, i-1) @ site(sm, i); f[i] = e[i].conj().T
        k[i] = np.diag(q**(2*(Wt[:, i-1] - Wt[:, i]))).astype(complex); qi[i] = q**2
    e[n] = site(sp, n-1); f[n] = e[n].conj().T; k[n] = np.diag(q**(2*Wt[:, n-1])).astype(complex); qi[n] = q
    e[0] = x*site(sm, 0); f[0] = site(sp, 0)/x; k[0] = np.diag(q**(-2*Wt[:, 0])).astype(complex); qi[0] = q
    return e, f, k, qi
def Rgen(Arep, Brep, wA, wB):
    Rc, null, _, _ = intertwiner(Arep, Brep, range(n+1), wA, wB); return Rc/np.abs(Rc).max(), null
Ws, Wv = spin_weights(), weights(n)
if __name__ == '__main__':
    # verify algebra relations (same checks as for the vector)
    import dtw
    dtw_rep_backup = dtw.rep
    dtw.rep = lambda nn, qq, x: spin_rep(x)
    A, err = dtw.check_relations(n, q); dtw.rep = dtw_rep_backup
    print(f"spinor rep: Cartan matrix rows {A.astype(int).tolist()}, relation violation {err:.1e}")
    Ws, Wv = spin_weights(), weights(n)
    def special(fun, D):
        found = []
        for l in np.arange(-7, 7.5, 0.5):
            for sgn in (+1, -1):
                dims = []
                for d in (1e-5, 1e-7):
                    R, _ = fun(sgn*q**l*(1+d)); R = R/np.linalg.norm(R); s = np.linalg.svd(R, compute_uv=False)
                    dims.append(int(np.sum(s > 1e-3*s[0])))
                if dims[0] == dims[1] < D: found.append((sgn, l, dims[1]))
        return found
    for name, fun, D in [("spinor x spinor (R_33)", lambda x: Rgen(spin_rep(x), spin_rep(1.0), Ws, Ws), 64),
                         ("vector x spinor (R_13)", lambda x: Rgen(rep(n, q, x), spin_rep(1.0), Wv, Ws), 64)]:
        _, null = fun(1.3+0.4j)
        print(f"{name}: intertwiner unique: {null == 1}")
        for sgn, l, d in special(fun, D):
            if sgn > 0: print(f"     x* = +-q^{l:+.1f}: image dimension {d}")
