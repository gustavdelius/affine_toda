"""Gradient flows, intersection numbers n_sigma = <R, K_sigma>, and exact thimble integrals (n=2).

Conventions (weight exp(-St/hbar), St = e^{-i eps} S):
  thimble J_sigma: upward flow dq/dt = +conj(grad St) from sigma (Re St increases), oriented by V;
  dual K_sigma:   downward flow dq/dt = -conj(grad St) from sigma, oriented by iV.
  <J_sigma, K_sigma> = +1, R = sum_sigma n_sigma J_sigma with n_sigma = <R, K_sigma>,
  where R = R^n with its standard orientation; sign at p in R cap K = sign det Im(T_K).
"""
import numpy as np
from scipy.optimize import least_squares


def flow_field(ch, Q, sgn, gcap=5.0):
    """Q (B,n) complex. sgn=-1 downward, +1 upward, speed capped at gcap."""
    U = Q.reshape(Q.shape[0], ch.L, ch.r)
    G = (ch.phase * ch.grad(U)).reshape(Q.shape[0], -1)
    F = sgn * np.conj(G)
    nrm = np.linalg.norm(F, axis=1)
    fac = np.where(nrm > gcap, gcap / np.maximum(nrm, 1e-300), 1.0)
    return F * fac[:, None]


def rk4_step(ch, Q, dt, sgn, gcap):
    k1 = flow_field(ch, Q, sgn, gcap)
    k2 = flow_field(ch, Q + 0.5 * dt * k1, sgn, gcap)
    k3 = flow_field(ch, Q + 0.5 * dt * k2, sgn, gcap)
    k4 = flow_field(ch, Q + dt * k3, sgn, gcap)
    return Q + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def St(ch, Q):
    return ch.phase * ch.action(Q.reshape(Q.shape[0], ch.L, ch.r))


def sphere_dirs(n, N, rng):
    X = rng.normal(size=(N, n))
    return X / np.linalg.norm(X, axis=1, keepdims=True)


def scan_K(ch, sig, V, dirs, delta=0.02, dt=0.01, tmax=60.0, gcap=5.0, ymax=12.0):
    """Integrate downward flows from sigma + i delta V c; track min |Im q| along each."""
    n = ch.n
    s0 = sig.reshape(-1)
    Q = s0[None, :] + 1j * delta * (dirs @ V.T)
    B = Q.shape[0]
    best = np.full(B, np.inf)
    tbest = np.zeros(B)
    alive = np.ones(B, bool)
    floor = np.cos(ch.eps) * ch.Kmin - abs(np.sin(ch.eps)) * ch.Imax - 1e-9
    t = 0.0
    nsteps = int(tmax / dt)
    for k in range(nsteps):
        idx = np.where(alive)[0]
        if len(idx) == 0:
            break
        Qa = rk4_step(ch, Q[idx], dt, -1, gcap)
        Q[idx] = Qa
        t += dt
        m = np.linalg.norm(Qa.imag, axis=1)
        upd = m < best[idx]
        best[idx[upd]] = m[upd]
        tbest[idx[upd]] = t
        ReS = St(ch, Qa).real
        dead = (~np.isfinite(m)) | (ReS < floor) | (np.max(np.abs(Qa.imag), axis=1) > ymax)
        alive[idx[dead]] = False
    return best, tbest


def integrate_to(ch, Q0, t, dt=0.01, gcap=5.0, sgn=-1):
    """Integrate batch Q0 for time t (same t for all)."""
    nst = max(1, int(np.ceil(t / dt)))
    h = t / nst
    Q = Q0.copy()
    for k in range(nst):
        Q = rk4_step(ch, Q, h, sgn, gcap)
    return Q


def tangent_basis(c):
    n = len(c)
    M = np.eye(n) - np.outer(c, c)
    U, s, Vt = np.linalg.svd(M)
    E = U[:, : n - 1]
    if np.linalg.det(np.column_stack([c, E])) < 0:
        E[:, 0] *= -1
    return E


def refine_K_hit(ch, sig, V, c0, t0, delta=0.02, dt=0.005, gcap=5.0):
    """Solve Im q(c(y), t) = 0 for (y, t). Returns dict or None."""
    n = ch.n
    s0 = sig.reshape(-1)
    E = tangent_basis(c0)

    def qof(params):
        P = np.atleast_2d(params)
        cs = c0[None, :] + P[:, : n - 1] @ E.T
        cs = cs / np.linalg.norm(cs, axis=1, keepdims=True)
        ts = P[:, n - 1]
        out = []
        for c, t in zip(cs, ts):
            if t <= 0:
                out.append(np.full(n, np.nan + 0j))
                continue
            Q0 = (s0 + 1j * delta * (V @ c))[None, :]
            out.append(integrate_to(ch, Q0, t, dt, gcap)[0])
        return np.array(out)

    def res(p):
        q = qof(p)[0]
        r = q.imag
        if not np.all(np.isfinite(r)):
            return np.full(n, 1e3)
        return r

    p0 = np.concatenate([np.zeros(n - 1), [t0]])
    try:
        sol = least_squares(res, p0, method='lm', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=400, diff_step=1e-7)
    except Exception:
        return None
    if not np.all(np.isfinite(sol.fun)) or np.linalg.norm(sol.fun) > 1e-9:
        return None
    p = sol.x
    q = qof(p)[0]
    # tangent frame for orientation
    h = 1e-6
    cols = [flow_field(ch, q[None, :], -1, gcap)[0]]
    for i in range(n - 1):
        dp = np.zeros(n); dp[i] = h
        qp = qof(p + dp)[0]; qm = qof(p - dp)[0]
        cols.append((qp - qm) / (2 * h))
    T = np.column_stack(cols)
    c = c0 + E @ p[: n - 1]
    c = c / np.linalg.norm(c)
    Ec = tangent_basis(c)
    # orientation of (c, d c/dy) relative to (c, E): use projected derivative directions
    D = np.column_stack([(lambda cc: cc / np.linalg.norm(cc))(c0 + E @ (p[: n - 1] + h * np.eye(n - 1)[i])) - c for i in range(n - 1)]) / h
    orient = np.sign(np.linalg.det(np.column_stack([c, D])))
    sgn = np.sign(np.linalg.det(T.imag)) * orient
    return dict(p=q.real.copy(), t=p[n - 1], c=c, sign=int(sgn), detImT=np.linalg.det(T.imag),
                ReS=St(ch, q[None, :])[0].real, ImS=St(ch, q[None, :])[0].imag)


def intersection_number(ch, sig, ndir=400, seed=0, delta=0.02, ncand=40, dt_scan=0.01, tmax=60.0, verbose=False,
                        cand_thresh=0.3):
    mu, V = ch.tangent(sig)
    rng = np.random.default_rng(seed)
    n = ch.n
    if n == 1:
        dirs = np.array([[1.0], [-1.0]])
    elif n == 2:
        ph = np.linspace(0, 2 * np.pi, ndir, endpoint=False)
        dirs = np.column_stack([np.cos(ph), np.sin(ph)])
    else:
        dirs = sphere_dirs(n, ndir, rng)
    best, tbest = scan_K(ch, sig, V, dirs.copy(), delta=delta, dt=dt_scan, tmax=tmax)
    order = np.argsort(best)
    hits = []
    tried = 0
    for i in order:
        if best[i] > cand_thresh or tried >= ncand:
            break
        tried += 1
        if n == 1:
            # 1D: K is a curve; hit iff Im q passes through zero; refine t only
            pass
        r = refine_K_hit(ch, sig, V, dirs[i], tbest[i], delta=delta)
        if r is None:
            continue
        dup = any(np.max(np.abs(r['p'] - h['p'])) < 1e-5 for h in hits)
        if not dup:
            hits.append(r)
    nsig = sum(h['sign'] for h in hits)
    return nsig, hits, mu, V


def thimble_integral_n2(ch, sig, hbar, V=None, mu=None, nphi=512, dt=0.004, delta=1e-3, gcap=5.0, cut=50.0,
                        tmax=80.0):
    """Exact integral over J_sigma for n=2 (with normalization N), orientation V.
    Chart (t, phi) -> Phi^up_t(sigma + delta V (cos phi, sin phi)); each ray is frozen once Re(St-St_sigma)/hbar > cut."""
    assert ch.n == 2
    if V is None:
        mu, V = ch.tangent(sig)
    s0 = sig.reshape(-1)
    Ss = ch.phase * ch.action(sig)
    N = (2 * np.pi * hbar) ** (-ch.r * (ch.L + 1) / 2)
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    hph = 1e-6

    def start(phv):
        C = np.column_stack([np.cos(phv), np.sin(phv)])
        return s0[None, :] + delta * (C @ V.T)

    Q = np.concatenate([start(ph), start(ph + hph), start(ph - hph)])
    nst = int(tmax / dt)
    active = np.ones(nphi, bool)
    total = 0.0 + 0j
    fprev = None
    with np.errstate(all='ignore'):
        for k in range(nst + 1):
            idx = np.where(active)[0]
            if len(idx) == 0:
                break
            if k > 0:
                sel = np.concatenate([idx, idx + nphi, idx + 2 * nphi])
                Q[sel] = rk4_step(ch, Q[sel], dt, +1, gcap)
            q0, qp, qm = Q[:nphi], Q[nphi:2 * nphi], Q[2 * nphi:]
            Fq = flow_field(ch, q0, +1, gcap)
            dq = (qp - qm) / (2 * hph)
            Jc = Fq[:, 0] * dq[:, 1] - Fq[:, 1] * dq[:, 0]
            d = St(ch, q0) - Ss
            f = np.where(active, np.exp(-d / hbar) * Jc, 0)
            f[~np.isfinite(f)] = 0
            w = 0.5 if (k == 0) else 1.0
            total += w * np.sum(f) * dt
            active &= ~((d.real / hbar > cut) | ~np.isfinite(d.real))
    I_outer = total * (2 * np.pi / nphi)
    a = mu[0] * np.cos(ph) ** 2 + mu[1] * np.sin(ph) ** 2
    disk = np.sum(hbar / a * (1 - np.exp(-a * delta ** 2 / (2 * hbar)))) * (2 * np.pi / nphi)
    I_inner = np.linalg.det(V) * disk
    return N * np.exp(-Ss / hbar) * (I_outer + I_inner)


def thimble_integral_n2_chart(ch, sig, hbar, V=None, mu=None, nrho=48, nphi=128, dt=0.002, gcap=50.0, start_eps=1e-7,
                              rho_cut=40.0):
    """Exact integral over J_sigma (n=2) in asymptotic linear coordinates c in R^2:
    q(c) = Phi^up_T(sigma + V e^{-D T} c), smooth chart of J_sigma, oriented by V.
    Polar quadrature: Gauss-Legendre in rho, periodic trapezoid in phi."""
    assert ch.n == 2
    if V is None:
        mu, V = ch.tangent(sig)
    s0 = sig.reshape(-1)
    Ss = ch.phase * ch.action(sig)
    N = (2 * np.pi * hbar) ** (-ch.r * (ch.L + 1) / 2)
    rho_max = np.sqrt(2 * rho_cut * hbar / mu.min())
    T = np.log(rho_max / start_eps) / mu.min()
    x, w = np.polynomial.legendre.leggauss(nrho)
    rho = 0.5 * rho_max * (x + 1); wr = 0.5 * rho_max * w
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    R, P = np.meshgrid(rho, ph, indexing='ij')
    C = np.stack([R * np.cos(P), R * np.sin(P)], -1).reshape(-1, 2)
    hc = 1e-6 * max(1.0, rho_max)
    Cs = np.concatenate([C, C + [hc, 0], C - [hc, 0], C + [0, hc], C - [0, hc]])
    scale = np.exp(-mu * T)
    Q = s0[None, :] + (Cs * scale[None, :]) @ V.T
    nst = int(np.ceil(T / dt)); h = T / nst
    with np.errstate(all='ignore'):
        for k in range(nst):
            Q = rk4_step(ch, Q, h, +1, gcap)
    M = C.shape[0]
    q0 = Q[:M]
    d1 = (Q[M:2 * M] - Q[2 * M:3 * M]) / (2 * hc)
    d2 = (Q[3 * M:4 * M] - Q[4 * M:]) / (2 * hc)
    Jc = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    d = St(ch, q0) - Ss
    f = np.exp(-d / hbar) * Jc
    f[~np.isfinite(f)] = 0
    f = f.reshape(nrho, nphi)
    I = np.sum(f * (wr * rho)[:, None]) * (2 * np.pi / nphi)
    return N * np.exp(-Ss / hbar) * I, np.max(np.abs(d.imag))
