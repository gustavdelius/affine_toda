import numpy as np, itertools, importlib
import dtw, rsolve
from rsolve import intertwiner
q = 0.7
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def make(n):
    def site(op, i):
        mats = [I2]*n; mats[i] = op; out = mats[0]
        for m in mats[1:]: out = np.kron(out, m)
        return out
    Wt = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n)])
    def spin_rep(x):
        e, f, k, qi = {}, {}, {}, {}
        for i in range(1, n):
            e[i] = site(sp, i-1) @ site(sm, i); f[i] = e[i].conj().T
            k[i] = np.diag(q**(2*(Wt[:, i-1] - Wt[:, i]))).astype(complex); qi[i] = q**2
        e[n] = site(sp, n-1); f[n] = e[n].conj().T; k[n] = np.diag(q**(2*Wt[:, n-1])).astype(complex); qi[n] = q
        e[0] = x*site(sm, 0); f[0] = site(sp, 0)/x; k[0] = np.diag(q**(-2*Wt[:, 0])).astype(complex); qi[0] = q
        return e, f, k, qi
    return spin_rep, Wt
br = lambda x, l, s: (x - s*q**l)/(1 - s*x*q**l)        # <l>_s
def test(n, xs):
    spin_rep, Wt = make(n); worst = 0
    from math import comb
    dims = [comb(2*n+1, j) for j in range(n+1)]          # Lambda^j of the (2n+1)-dim vector, j = 0..n
    for x in xs:
        R, null, _, _ = intertwiner(spin_rep(x), spin_rep(1.0), range(n+1), Wt, Wt)
        ev = np.linalg.eigvals(R)
        # predicted eigenvalue ratios relative to Lambda^n:  rho_{Lambda^{n-k}} / rho_{Lambda^n} = prod_{j=1..k} <2j>_{(-1)^{j+1}}
        pred = {}
        for kk in range(n+1):
            pred[dims[n-kk]] = np.prod([br(x, 2*j, (-1)**(j+1)) for j in range(1, kk+1)]) if kk else 1.0
        # match: find rho_top as eigenvalue with multiplicity dims[n]
        vals, cnt = [], []
        for v in ev:
            for t, u in enumerate(vals):
                if abs(u - v) < 1e-7*max(1, abs(v)): cnt[t] += 1; break
            else: vals.append(v); cnt.append(1)
        byd = {c: v for v, c in zip(vals, cnt)}
        top = byd[dims[n]]
        for dd, p in pred.items(): worst = max(worst, abs(byd[dd]/top - p)/max(1, abs(p)))
    return worst, null, dims
xs = [1.3+0.4j, -0.6+0.9j, 2.2-0.7j]
for n in [2, 3, 4]:
    worst, null, dims = test(n, xs)
    print(f"n={n}: spinor (x) spinor = " + " + ".join(str(d) for d in dims[::-1]) +
          f"; R unique: {null == 1}; max deviation from  R = sum_k [prod_(j<=k) <2j>_((-1)^(j+1))] P_(Lambda^(n-k)) : {worst:.1e}")
