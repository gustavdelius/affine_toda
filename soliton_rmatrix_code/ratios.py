import numpy as np
from spinor import spin_rep, Rgen, Ws, q, n
def eig_ratios(x):
    R, _ = Rgen(spin_rep(x), spin_rep(1.0), Ws, Ws); ev = np.linalg.eigvals(R)
    groups = []
    for v in ev:
        for g in groups:
            if abs(g[0] - v) < 1e-7*max(1, abs(v)): g[1] += 1; break
        else: groups.append([v, 1])
    d = {m: v for v, m in groups}
    if sorted(d) != [1, 7, 21, 35]: return None
    return d[21]/d[35], d[7]/d[35], d[1]/d[35]
xs = 1.37*np.exp(2j*np.pi*(np.arange(40)+0.31)/40)
data = [eig_ratios(x) for x in xs]
ok = [i for i, d in enumerate(data) if d is not None]; xs = xs[ok]; data = np.array([data[i] for i in ok])
print(f"{len(ok)} sample points with clean (35,21,7,1) eigenvalue grouping")
def fit_rational(x, r, dmax=6):
    for d in range(1, dmax+1):
        # r * (1 + b1 x + ... + bd x^d) = a0 + ... + ad x^d
        A = np.hstack([np.vander(x, d+1, increasing=True), -(r[:, None]*np.vander(x, d+1, increasing=True)[:, 1:])])
        sol, *_ = np.linalg.lstsq(A, r, rcond=None)
        a, b = sol[:d+1], np.concatenate([[1], sol[d+1:]])
        res = np.abs(np.polyval(a[::-1], x)/np.polyval(b[::-1], x) - r).max()/np.abs(r).max()
        if res < 1e-9: return d, a, b, res
    return None
def qpow(z):
    return f"{np.angle(z)/np.pi:+.3f}pi phase, |z| = q^{np.log(abs(z))/np.log(q):+.3f}"
for lab, col in [("rho_21/rho_35", 0), ("rho_7/rho_35", 1), ("rho_1/rho_35", 2)]:
    d, a, b, res = fit_rational(xs, data[:, col])
    zeros, poles = np.roots(a[::-1]), np.roots(b[::-1])
    # cancel common roots
    print(f"{lab}: degree {d}, fit residual {res:.1e}")
    for z in zeros: print(f"      zero at x = {qpow(z)}")
    for p in poles: print(f"      pole at x = {qpow(p)}")
