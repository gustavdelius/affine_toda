"""Robust n=2 intersection numbers: winding number of (Im q1, Im q2) over cells of the (t, phi) chart
of K_sigma, with adaptive refinement in phi. Real points have Re St >= floor = cos(eps) min_R Re S - |sin eps| Imax."""
import numpy as np
from scipy.optimize import minimize
from flows import rk4_step, St

def realmin_ReS(ch, nstart=200, seed=0):
    rng = np.random.default_rng(seed)
    base = None
    from model import lin_interp
    b = lin_interp(ch).real.reshape(-1)
    def f(x):
        U = x.reshape(1, ch.L, ch.r).astype(complex)
        return ch.action(U)[0].real
    best = np.inf
    for k in range(nstart):
        x0 = b + rng.normal(0, 1.0, b.shape)
        r = minimize(f, x0, method='BFGS')
        best = min(best, r.fun)
    return best

def floor_of(ch, rmin):
    return np.cos(ch.eps)*rmin - abs(np.sin(ch.eps))*ch.Imax - 1e-9

def trajectories(ch, sig, V, phis, floor, delta=0.02, dt=0.01, tmax=40.0, gcap=5.0, ymax=12.0):
    s0 = sig.reshape(-1)
    C = np.column_stack([np.cos(phis), np.sin(phis)])
    Q = s0[None, :] + 1j*delta*(C @ V.T)
    nst = int(tmax/dt)
    traj = np.full((nst+1, len(phis), 2), np.nan+0j)
    traj[0] = Q
    alive = np.ones(len(phis), bool)
    with np.errstate(all='ignore'):
        for k in range(1, nst+1):
            idx = np.where(alive)[0]
            if len(idx) == 0:
                traj = traj[:k]; break
            Qa = rk4_step(ch, Q[idx], dt, -1, gcap)
            Q[idx] = Qa
            ReS = St(ch, Qa).real
            dead = (~np.all(np.isfinite(Qa), axis=1)) | ~np.isfinite(ReS) | (ReS < floor) | (np.max(np.abs(Qa.imag), axis=1) > ymax)
            good = idx[~dead]
            traj[k, good] = Qa[~dead]
            alive[idx[dead]] = False
    return traj

def winding(traj):
    Y = traj.imag
    ang = np.arctan2(Y[..., 1], Y[..., 0])
    def d(a, b):
        return (b - a + np.pi) % (2*np.pi) - np.pi
    a00 = ang[:-1]; a10 = ang[1:]; a01 = np.roll(ang[:-1], -1, axis=1); a11 = np.roll(ang[1:], -1, axis=1)
    ok = np.isfinite(a00)&np.isfinite(a10)&np.isfinite(a01)&np.isfinite(a11)
    tot = d(a00, a10) + d(a10, a11) + d(a11, a01) + d(a01, a00)
    w = np.where(ok, np.round(tot/(2*np.pi)), 0).astype(int)
    # partial cells at the death boundary: count cells where some corner is nan but finite ones straddle a zero
    return w, ok

def intersection_number_n2(ch, sig, floor, nphi=512, delta=0.02, dt=0.01, tmax=40.0, jump_tol=0.05, maxlevels=8, maxphi=60000):
    mu, V = ch.tangent(sig)
    phis = np.linspace(0, 2*np.pi, nphi, endpoint=False)
    traj = trajectories(ch, sig, V, phis, floor, delta, dt, tmax)
    for lev in range(maxlevels):
        J = np.nanmax(np.abs(traj - np.roll(traj, -1, axis=1)).max(axis=2), axis=0)  # per phi-gap
        J = np.nan_to_num(J, nan=0.0)
        bad = np.where(J > jump_tol)[0]
        if len(bad) == 0 or len(phis) > maxphi:
            break
        nxt = (phis[(bad+1) % len(phis)] + 2*np.pi*((bad+1)//len(phis)))
        mids = 0.5*(phis[bad] + np.where(bad+1 < len(phis), phis[(bad+1) % len(phis)], phis[0]+2*np.pi))
        tn = trajectories(ch, sig, V, mids, floor, delta, dt, tmax)
        nt = max(traj.shape[0], tn.shape[0])
        def padt(T):
            if T.shape[0] < nt:
                T = np.concatenate([T, np.full((nt-T.shape[0],)+T.shape[1:], np.nan+0j)])
            return T
        traj = np.concatenate([padt(traj), padt(tn)], axis=1)
        phis = np.concatenate([phis, mids % (2*np.pi)])
        o = np.argsort(phis); phis = phis[o]; traj = traj[:, o]
    J = np.nan_to_num(np.nanmax(np.abs(traj - np.roll(traj, -1, axis=1)).max(axis=2), axis=0), nan=0.0)
    w, ok = winding(traj)
    hits = []
    for k, i in zip(*np.where(w != 0)):
        q = traj[k, i]
        hits.append(dict(sign=int(w[k, i]), k=int(k), phi=phis[i], q=q, ReS=St(ch, q[None, :])[0].real))
    # incomplete cells near death boundary with sign changes of both Im components: flag
    return sum(h['sign'] for h in hits), hits, mu, V, J.max(), len(phis)

def newton_hit_n2(ch, sig, V, phi0, t0, delta=0.02, dt=0.002, gcap=5.0, maxit=30, tol=1e-11):
    """Solve Im q(phi, t) = 0 exactly (RK4, small dt). Returns (phi, t, q, sign) or None."""
    from flows import flow_field
    s0 = sig.reshape(-1)
    hph = 1e-6
    phi, t = phi0, t0
    for it in range(maxit):
        if t <= 0 or t > 200:
            return None
        ph3 = np.array([phi, phi + hph, phi - hph])
        C = np.column_stack([np.cos(ph3), np.sin(ph3)])
        Q = s0[None, :] + 1j*delta*(C @ V.T)
        nst = max(1, int(np.ceil(t/dt))); h = t/nst
        with np.errstate(all='ignore'):
            for k in range(nst):
                Q = rk4_step(ch, Q, h, -1, gcap)
        if not np.all(np.isfinite(Q)):
            return None
        q = Q[0]
        Ft = flow_field(ch, q[None, :], -1, gcap)[0]
        dphi = (Q[1] - Q[2])/(2*hph)
        Jm = np.column_stack([Ft.imag, dphi.imag])  # d Im q / d(t, phi)
        r = q.imag
        if np.linalg.norm(r) < tol:
            sgn = int(np.sign(np.linalg.det(Jm)))
            return phi % (2*np.pi), t, q, sgn
        try:
            d = np.linalg.solve(Jm, -r)
        except np.linalg.LinAlgError:
            return None
        lim = 0.2
        if abs(d[0]) > lim: d = d*lim/abs(d[0])
        t += d[0]; phi += d[1]
    return None

def intersection_number_n2_refined(ch, sig, floor, nphi=512, delta=0.02, **kw):
    ns0, hits0, mu, V, mj, nph = intersection_number_n2(ch, sig, floor, nphi=nphi, delta=delta, **kw)
    sols = []
    dt = 0.01
    for hh in hits0:
        r = newton_hit_n2(ch, sig, V, hh['phi'], hh['k']*dt, delta=delta)
        if r is None:
            continue
        phi, t, q, sgn = r
        if St(ch, q[None, :])[0].real < floor:
            continue
        if not any(np.max(np.abs(q - s[2])) < 1e-6 for s in sols):
            sols.append((phi, t, q, sgn))
    n = sum(s[3] for s in sols)
    return n, sols, mu, V, hits0
