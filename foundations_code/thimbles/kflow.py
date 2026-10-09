"""Generic-n intersection numbers n_sigma = <R, K_sigma> by downward flows from sigma.
K_sigma chart: (t, c) -> sigma + delta_t, delta_0 = i V e^{-D tau} c, c in S^{n-1} (uniform in the asymptotic
linear chart), flowed downward in deviation coordinates (expm1, full relative precision near sigma).
Orientation (d_t, d_c) ~ (radial, angular). Intersection points Im q = 0 solved by Newton in (y, t)."""
import numpy as np
from flows import St, tangent_basis
from thimble_exact import rk4 as rk4dev, field as fielddev

def start_points(ch, V, mu, C, amp=0.02):
    tau = np.log(1.0/amp)/mu.min()
    return 1j*((C*np.exp(-mu*tau)[None, :]) @ V.T)

def scan(ch, sig, V, mu, dirs, floor, dt=0.01, tmax=50.0, gcap=5.0, ymax=12.0):
    s0 = sig.reshape(-1)
    D = start_points(ch, V, mu, dirs)
    B = D.shape[0]
    best = np.full(B, np.inf); tbest = np.zeros(B)
    alive = np.ones(B, bool); t = 0.0
    with np.errstate(all='ignore'):
        for k in range(int(tmax/dt)):
            idx = np.where(alive)[0]
            if len(idx) == 0: break
            Da = rk4dev(ch, sig, D[idx], dt, -1, gcap); D[idx] = Da; t += dt
            Qa = s0[None, :] + Da
            ReS = St(ch, Qa).real
            dead = (~np.all(np.isfinite(Qa), axis=1)) | ~np.isfinite(ReS) | (ReS < floor) | (np.max(np.abs(Qa.imag), axis=1) > ymax)
            m = np.linalg.norm(Qa.imag, axis=1)
            upd = (m < best[idx]) & ~dead
            best[idx[upd]] = m[upd]; tbest[idx[upd]] = t
            alive[idx[dead]] = False
    return best, tbest

def newton_hit(ch, sig, V, mu, c0, t0, floor, dt=0.004, gcap=5.0, maxit=30, tol=1e-10):
    n = ch.n; s0 = sig.reshape(-1)
    E = tangent_basis(c0)
    y = np.zeros(n-1); t = t0; hy = 1e-6
    for it in range(maxit):
        if t <= 0 or t > 100 or np.linalg.norm(y) > 1.5: return None
        Ys = np.vstack([y] + [y + hy*np.eye(n-1)[i] for i in range(n-1)] + [y - hy*np.eye(n-1)[i] for i in range(n-1)])
        Cs = c0[None, :] + Ys @ E.T
        Cs /= np.linalg.norm(Cs, axis=1, keepdims=True)
        D = start_points(ch, V, mu, Cs)
        nst = max(1, int(np.ceil(t/dt))); h = t/nst
        with np.errstate(all='ignore'):
            for k in range(nst):
                D = rk4dev(ch, sig, D, h, -1, gcap)
        if not np.all(np.isfinite(D)): return None
        q = s0 + D[0]
        Ft = fielddev(ch, sig, D[:1], -1, gcap)[0]
        dY = [(D[1+i] - D[n+i])/(2*hy) for i in range(n-1)]
        T = np.column_stack([Ft] + dY)
        r = q.imag
        if np.linalg.norm(r) < tol:
            if St(ch, q[None, :])[0].real < floor: return None
            return dict(p=q.real.copy(), q=q, t=t, sign=int(np.sign(np.linalg.det(T.imag))), ReS=St(ch, q[None, :])[0].real,
                        cond=np.linalg.cond(T.imag))
        try:
            d = np.linalg.solve(T.imag, -r)
        except np.linalg.LinAlgError:
            return None
        nd = np.linalg.norm(d)
        if nd > 0.3: d *= 0.3/nd
        t += d[0]; y += d[1:]
    return None

def intersection_number(ch, sig, floor, ndir=3000, ncand=80, seed=0, thr=0.6):
    mu, V = ch.tangent(sig)
    n = ch.n
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(ndir, n)); dirs = X/np.linalg.norm(X, axis=1, keepdims=True)
    best, tbest = scan(ch, sig, V, mu, dirs, floor)
    order = np.argsort(best)
    hits = []; tried = 0
    for i in order:
        if best[i] > thr or tried >= ncand: break
        tried += 1
        r = newton_hit(ch, sig, V, mu, dirs[i], tbest[i], floor)
        if r is None: continue
        if not any(np.max(np.abs(r['p'] - h['p'])) < 1e-6 for h in hits):
            hits.append(r)
    return sum(h['sign'] for h in hits), hits, mu, V, tried, int(np.sum(best < thr))
