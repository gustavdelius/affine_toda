# Lemma 6.11 step 4: the two terms of Delta_kl on Z^2 (no walls), -Laplacian normalisation.
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
def ball(c, rho):
    R = int(np.ceil(rho))+1
    return [(c[0]+i, c[1]+j) for i in range(-R, R+1) for j in range(-R, R+1) if i*i+j*j <= rho*rho]
def green_row(A, x):
    idx = {p: i for i, p in enumerate(A)}; n = len(A)
    rows, cols, vals = [], [], []
    for p, i in idx.items():
        rows.append(i); cols.append(i); vals.append(4.0)
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            q = (p[0]+d[0], p[1]+d[1])
            if q in idx: rows.append(i); cols.append(idx[q]); vals.append(-1.0)
    L = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    e = np.zeros(n); e[idx[x]] = 1.0
    g = spl.spsolve(L, e)          # G_A(x, .), symmetric
    return {p: g[i] for p, i in idx.items()}
rho = 40.0
xl = (0, 0); Bl = ball(xl, rho); Gl = green_row(Bl, xl); Blset = set(Bl)
print(" r   term1   term2   G_tildeB(xl,xk)  Delta   (1/2pi)ln(rho/r)")
for r in [1, 2, 4, 8, 16, 30, 40, 45, 60]:
    xk = (r, 0); Bk = ball(xk, rho); Gk = green_row(Bk, xk); Bkset = set(Bk)
    mu = {}
    for z, g in Gk.items():
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            y = (z[0]+d[0], z[1]+d[1])
            if y not in Bkset: mu[y] = mu.get(y, 0.0)+g
    t1 = Gk[xl] if xl in Bkset else 0.0
    t2 = sum(mu[y]*Gl[y] for y in mu if y in Blset)
    Bt = ball(xl, r+rho+2); Gt = green_row(Bt, xl)
    print(f"{r:3d} {t1:7.3f} {t2:7.3f} {Gt[xk]:10.3f} {t1+t2:8.3f} {np.log(rho/r)/(2*np.pi):8.3f}")
