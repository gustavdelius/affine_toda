"""Claim 4, dynamical side, with the book's exact soliton solutions (boundary_images.Sol) and Lax pair (boundary_classical.a_t).

(i)  static self-mirror soliton (b, theta=0, xi = i y): phi(x) + conj phi(-x) = 2v, v = Re phi(0) in (pi/beta) x coweight class of q_b
(ii) three-soliton image solution: (a, th, xi1), (a, -th, xi2), (b, 0, xi3) with
        Re xi3 = -ln A_ab(th),  xi2 = -conj(xi1) - ln A_aa(2th) - ln A_ab(th)
(iii) Sklyanin's tau(lam) = tr(K_+ T(lam) T(1/lam)^dagger) on [-L, 0] (K_- = 1, vacuum at -L) is conserved along (ii)
     with K_+ = e^{i pi nu.H}, nu = beta v/pi, and not with K_+ = 1; control: the book's Theta P pair (trivial class) with K_+ = 1.
"""
import numpy as np
import scipy.linalg as sla
from booklib import load_classical, load_images

C = load_classical(); I = load_images()
Sol, interaction, roots_an, vac_round = I['Sol'], I['interaction'], I['roots_an'], I['vac_round']
beta, m = I['beta'], I['m']
C['beta'], C['m'] = beta, m                     # use the soliton conventions in the book's a_t
algebra, a_t = C['algebra'], C['a_t']


def taus_t(sol, x, t):
    """d tau_j / dt, same expansion as Sol.taus."""
    x = np.atleast_1d(np.asarray(x, float))
    k = sol.ma[:, None] * np.cosh(sol.th)[:, None]
    w = sol.ma[:, None] * np.sinh(sol.th)[:, None]
    expo = k * x[None, :] - w * t + sol.xi[:, None]
    dT = np.zeros((sol.h, x.size), complex)
    for S, c in zip(sol.subsets, sol.coef):
        if not S:
            continue
        e = c * np.exp(expo[list(S)].sum(0)); dw = -w[list(S)].sum(0)
        for j in range(sol.h):
            dT[j] += sol.om ** (j * sum(sol.sp[i] for i in S)) * dw * e
    return dT


def fields(sol, x, t):
    """phi, d_x phi, d_t phi on the increasing grid x (book conventions, continuous branches from the left)."""
    ph, dph, T = sol.phi(x, t)
    dT = taus_t(sol, x, t)
    al = roots_an(sol.n)
    dtph = (1j / beta) * sum(al[j][:, None] * (dT[j] / T[j])[None, :] for j in range(sol.h))
    return ph, dph, dtph, T


def a_x(lam, phi, dtphi, roots, Es, Hs, nj):
    """x component of the book's Lax pair (eq-boundary-lax-at), imaginary coupling beta -> i beta."""
    b = 1j * beta
    out = (b / 2) * sum(dtphi[k] * Hs[k] for k in range(len(Hs)))
    for a, E, n in zip(roots, Es, nj):
        c = np.sqrt(n * (a @ a) / 2)
        out = out + (m / 2) * c * np.exp(b * (a @ phi) / 2) * (lam * E + E.T / lam)
    return out


def transport(sol, t, lam, L, N):
    """T(lam) from x=-L to x=0, dT/dx = -a_x T, 4th-order Magnus with two Gauss points per step."""
    n = sol.n; roots, Es, Hs, nj, sub = algebra(f'a{n}(1)')
    hstep = L / N; xs = -L + hstep * np.arange(N)
    g = hstep / (2 * np.sqrt(3))
    pts = np.empty(2 * N); pts[0::2] = xs + hstep / 2 - g; pts[1::2] = xs + hstep / 2 + g
    ph, dph, dtph, _ = fields(sol, pts, t)
    T = np.eye(sol.h, dtype=complex)
    for k in range(N):
        A1 = -a_x(lam, ph[:, 2 * k], dtph[:, 2 * k], roots, Es, Hs, nj)
        A2 = -a_x(lam, ph[:, 2 * k + 1], dtph[:, 2 * k + 1], roots, Es, Hs, nj)
        Om = (hstep / 2) * (A1 + A2) + (np.sqrt(3) * hstep**2 / 12) * (A2 @ A1 - A1 @ A2)
        T = sla.expm(Om) @ T
    return T


def tau(sol, t, lam, Kp, L=16.0, N=1600):
    T1 = transport(sol, t, lam, L, N); T2 = transport(sol, t, 1 / lam, L, N)
    return np.trace(Kp @ T1 @ T2.conj().T)


def compat(sol, t, lam, Kp):
    """|K a_t(lam) + a_t(1/lam)^dagger K| at x = 0 (book a_t, imaginary coupling)."""
    roots, Es, Hs, nj, sub = algebra(f'a{sol.n}(1)')
    xg = np.linspace(-25, 0, 5001)
    ph, dph, _, _ = fields(sol, xg, t)
    A = a_t(lam, ph[:, -1], dph[:, -1], roots, Es, Hs, nj, coupling='imag')
    B = a_t(1 / lam, ph[:, -1], dph[:, -1], roots, Es, Hs, nj, coupling='imag').conj().T
    return np.linalg.norm(Kp @ A + B @ Kp) / (np.linalg.norm(A) + np.linalg.norm(B))


def class_of(v, n):
    """v = (pi/beta) nu with nu a coweight: return nu, its distance from the coweight lattice, and whether nu in 2 Lambda."""
    al = roots_an(n)
    c = np.array([beta * (al[k] @ v) / np.pi for k in range(n + 1)])      # alpha_k . nu
    off = np.max(np.abs(c - np.round(c)))
    triv = np.all(np.round(c).astype(int) % 2 == 0)
    return beta * v / np.pi, off, triv


# zero-curvature sanity check of a_x against the book's a_t on a two-soliton solution
sol = Sol(2, [1, 2], [0.4, -0.7], [0.3 + 0.2j, -0.5 + 1.1j])
roots, Es, Hs, nj, sub = algebra('a2(1)')
xg = np.linspace(-20, 3.0, 23001); i0 = 20000; dx = xg[1] - xg[0]; dt = 1e-4; t0 = 0.37; lam = 0.8
def at_ax(t):
    ph, dph, dtph, _ = fields(sol, xg, t)
    return ph, dph, dtph
P0 = at_ax(t0); Pp = at_ax(t0 + dt); Pm = at_ax(t0 - dt)
At = lambda P, i: a_t(lam, P[0][:, i], P[1][:, i], roots, Es, Hs, nj, coupling='imag')
Ax = lambda P, i: a_x(lam, P[0][:, i], P[2][:, i], roots, Es, Hs, nj)
dxat = (At(P0, i0 + 1) - At(P0, i0 - 1)) / (2 * dx)
dtax = (Ax(Pp, i0) - Ax(Pm, i0)) / (2 * dt)
zc = dxat - dtax + Ax(P0, i0) @ At(P0, i0) - At(P0, i0) @ Ax(P0, i0)
print(f'zero curvature d_x a_t - d_t a_x + [a_x, a_t] on a 2-soliton solution: {np.linalg.norm(zc):.1e} (|d_x a_t| = {np.linalg.norm(dxat):.1e})')

# (i) static self-mirror solitons
xs_full = np.linspace(-25, 25, 10001); i0 = 5000
for n, b, y in [(1, 1, 0.7), (2, 1, 0.9), (2, 2, -1.3), (3, 1, 0.4), (3, 2, 0.8), (3, 2, 2.0)]:
    s = Sol(n, [b], [0.0], [1j * y])
    ph, dph, T = s.phi(xs_full, 0.0)
    c = ph + np.conj(ph[:, ::-1])
    v = ph[:, i0].real
    nu, off, triv = class_of(v, n)
    D = np.diag(np.exp(1j * np.pi * nu))
    print(f'a_{n} static b={b} y={y}: |phi(x)+conj phi(-x) - 2v| = {np.max(np.abs(c - c[:, [i0]])):.1e}, '
          f'|2v - 2Re phi(0)| = {np.max(np.abs(c[:, i0] - 2 * v)):.1e}, |d_x Im phi(0)| = {np.max(np.abs(dph[:, i0].imag)):.1e}, '
          f'nu = beta v/pi = {np.round(nu, 4)}, dist to coweights {off:.1e}, in 2Lambda: {triv}; '
          f'compat K=D: {compat(s, 0.0, 0.7, D):.1e}, K=1: {compat(s, 0.0, 0.7, np.eye(n + 1)):.1e}; min|tau| {np.min(np.abs(T)):.2f}')

# (ii) three-soliton image solutions and (iii) tau conservation
def three(n, a, b, th, xi1, y3):
    h = n + 1
    Aab = interaction(a, b, th, 0.0, h).real; Aaa = interaction(a, a, th, -th, h).real
    xi3 = -np.log(Aab) + 1j * y3
    xi2 = -np.conj(xi1) - np.log(Aaa) - np.log(Aab)
    return Sol(n, [a, a, b], [th, -th, 0.0], [xi1, xi2, xi3])

for (n, a, b, th, xi1, y3) in [(2, 1, 1, 0.5, 1.2 + 0.3j, 0.9), (2, 2, 1, 0.6, 0.8 - 0.4j, 0.5), (3, 1, 2, 0.5, 1.0 + 0.2j, 0.8),
                               (3, 2, 1, 0.45, 0.9 + 0.6j, 0.3)]:
    s = three(n, a, b, th, xi1, y3)
    worst_img = 0; vs = []; dIm = 0; mint = np.inf
    for t in np.linspace(-8, 8, 9):
        ph, dph, T = s.phi(xs_full, t)
        c = ph + np.conj(ph[:, ::-1])
        worst_img = max(worst_img, np.max(np.abs(c - c[:, [i0]])))
        vs.append(ph[:, i0].real); dIm = max(dIm, np.max(np.abs(dph[:, i0].imag))); mint = min(mint, np.min(np.abs(T)))
    vs = np.array(vs); v = vs[0]
    nu, off, triv = class_of(v, n)
    D = np.diag(np.exp(1j * np.pi * nu))
    cD = max(compat(s, t, lam, D) for t in (-3.0, 0.0, 2.5) for lam in (0.6, 1.5))
    c1 = max(compat(s, t, lam, np.eye(n + 1)) for t in (-3.0, 0.0, 2.5) for lam in (0.6, 1.5))
    print(f'a_{n} three-soliton (a={a}, b={b}, th={th}): image residual {worst_img:.1e}, Re phi(0,t) spread {np.ptp(vs, axis=0).max():.1e}, '
          f'|d_x Im phi(0)| {dIm:.1e}, nu = {np.round(nu, 4)} (dist {off:.1e}, in 2Lambda: {triv}), min|tau_j| {mint:.2f}; compat K=D {cD:.1e}, K=1 {c1:.1e}')
    for lam in (0.6, 1.5):
        tD = [tau(s, t, lam, D) for t in (-6.0, -2.0, 0.0, 2.0, 6.0)]
        t1 = [tau(s, t, lam, np.eye(n + 1)) for t in (-6.0, -2.0, 0.0, 2.0, 6.0)]
        rel = lambda z: np.ptp(np.abs(np.array(z) - z[0])) / abs(z[0])
        print(f'    lam={lam}: tau with K_+=D relative variation {rel(tD):.1e}; with K_+=1 {rel(t1):.1e}   (|tau| ~ {abs(tD[0]):.2e})')

# control: the book's Theta P pair (trivial class), K_+ = 1
for (n, a, xi1) in [(2, 1, 0.4 + 0.25j), (3, 2, 0.1 + 0.6j)]:
    th = 0.6; h = n + 1
    Aaa = interaction(a, a, th, -th, h).real
    s = Sol(n, [a, a], [th, -th], [xi1, -np.conj(xi1) - np.log(Aaa)])
    for lam in (0.6, 1.5):
        t1 = [tau(s, t, lam, np.eye(h)) for t in (-6.0, -2.0, 0.0, 2.0, 6.0)]
        print(f'control a_{n} Theta P pair, K_+=1, lam={lam}: relative variation {np.ptp(np.abs(np.array(t1) - t1[0])) / abs(t1[0]):.1e}')
