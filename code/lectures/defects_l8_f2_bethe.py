"""Lecture notes on defects, Lecture 8, figure 2 (fig-l8-bethe): the zeros of the defect operator are Bethe roots.

Six solitons of rapidities theta_j on a circle, coupling lambda = sqrt(5) + 0.3, q = -e^{-i pi lambda}. In each sector
of charge m = 0 and m = 2 the defect operator D(y_d) = <0| L_6^-1 ... L_1^-1 |0> commutes with the transfer matrix
t_omega, omega = q^{-m/2}, so it is diagonal on the eigenvectors of t_omega. Its eigenvalues, multiplied by
prod_j (1 + Y_j^2), are polynomials Q(y_d) in 1/y_d^2 of degree (6 - m)/2.
(a) Their zeros in the plane of theta_d, reduced modulo i pi/lambda (Q depends on y_d^2 = e^{2 lambda theta_d}), in the
strip pi/2 - pi/(2 lambda) <= Im theta_d < pi/2 + pi/(2 lambda): circles for m = 0 (20 states, 3 zeros each), crosses
for m = 2 (15 states, 2 zeros each). Most lie on the line Im theta_d = pi/2 or on the edges of the strip, which are one
line; the soliton rapidities are marked by ticks.
(b) |Q| along that line for three states of the neutral sector, divided by sum_j |c_j| |w|^j, where
Q = sum_j c_j w^j and w = 1/y_d^2; the dips are the zeros on the line.

Checks: the eigenvalues of D are diagonal in the eigenbasis of t_omega; Q has the stated degree; the plotted zeros
satisfy the Bethe equations prod_j (1 + z_j^2)/(1 + q^-2 z_j^2) Q(q y_k) = -Q(y_k/q), z_j = y_j/y_k; the eigenvalues
of t_omega at a test rapidity follow from Baxter's formula; the zeros plotted in (b) are zeros of the plotted curves.
Writes lecture-notes/defects/figures/l8_bethe.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"
LAM = math.sqrt(5) + 0.3
q = -cmath.exp(-1j * math.pi * LAM)
THS = np.array([-1.0, -0.55, -0.15, 0.2, 0.6, 1.05])
N = len(THS)
YS = np.exp(LAM * THS)


def site_Linv(Y, K):
    """L^-1(Y) on V (x) W at s = 0, as A[b, a] (d x d), n = -K..K; 0 = v+, 1 = v-."""
    ns = np.arange(-K, K + 1)
    d = 2 * K + 1
    A = np.zeros((2, 2, d, d), dtype=complex)
    den = 1 + Y ** 2
    A[0, 0] = np.diag(q ** ns.astype(complex)) / den
    A[1, 1] = np.diag(q ** (-ns).astype(complex)) / den
    for i, n in enumerate(ns):
        if i + 1 < d:
            A[1, 0][i + 1, i] = 1j * Y * q ** (-n - 1) / den
        if i >= 1:
            A[0, 1][i - 1, i] = 1j * Y * q ** (n - 1) / den
    return A


def mpo(sites, d):
    X = np.eye(d, dtype=complex)[None, None]
    for A in sites:
        nb = X.shape[0]
        Xn = np.zeros((2 * nb, 2 * nb, d, d), dtype=complex)
        for b2 in (0, 1):
            for a2 in (0, 1):
                if A[b2, a2].any():
                    Xn[b2::2, a2::2] = np.einsum('mn,bank->bamk', A[b2, a2], X)
        X = Xn
    return X


def Q_op(yd):
    K = N + 1
    X = mpo([site_Linv(yj / yd, K) for yj in YS], 2 * K + 1)
    return np.prod(1 + (YS / yd) ** 2) * X[:, :, K, K]


def R_site(y, yj):
    z = y / yj
    al = q * (z ** 2 - 1) / (q ** 2 * z ** 2 - 1)
    be = z * (q ** 2 - 1) / (q ** 2 * z ** 2 - 1)
    B = np.zeros((2, 2, 2, 2), dtype=complex)
    B[0, 0][0, 0] = B[1, 1][1, 1] = 1
    B[1, 1][0, 0] = B[0, 0][1, 1] = al
    B[0, 1][1, 0] = B[1, 0][0, 1] = be
    return B


def t_op(y, om):
    X = mpo([R_site(y, yj) for yj in YS], 2)
    return om * X[:, :, 0, 0] + X[:, :, 1, 1] / om


chg = np.array([sum(1 if ((v >> (N - 1 - k)) & 1) == 0 else -1 for k in range(N)) for v in range(2 ** N)])
data = {}
ok_diag, ok_deg, bethe, baxter = True, True, 0.0, 0.0
for m in (0, 2):
    idx = np.where(chg == m)[0]
    om = q ** (-m / 2)
    ev, V = np.linalg.eig(t_op(cmath.exp(LAM * (0.13 + 0.07j)), om)[np.ix_(idx, idx)])
    Vi = np.linalg.inv(V)
    M = (N - m) // 2
    npts = M + 4
    w2s = 0.6 * np.exp(1j * (np.arange(npts) * 2 * np.pi / npts + 0.3))
    vals = []
    for w2 in w2s:
        dd = Vi @ Q_op(1 / np.sqrt(w2))[np.ix_(idx, idx)] @ V
        ok_diag = ok_diag and np.abs(dd - np.diag(np.diag(dd))).max() < 1e-9 * np.abs(dd).max()
        vals.append(np.diag(dd))
    vals = np.array(vals)
    polys, roots = [], []
    for st in range(len(idx)):
        cf = np.linalg.solve(np.vander(w2s, npts, increasing=True), vals[:, st])
        ok_deg = ok_deg and np.abs(cf[M + 1:]).max() < 1e-9 * np.abs(cf).max() and abs(cf[M]) > 1e-8 * np.abs(cf).max()
        cM = cf[:M + 1]
        polys.append(cM)
        Qf = lambda y: np.polyval(cM[::-1], 1 / y ** 2)  # noqa: E731
        for w2k in np.roots(cM[::-1]):
            yk = 1 / np.sqrt(w2k)
            z2 = (YS / yk) ** 2
            lhs = np.prod((1 + z2) / (1 + z2 / q ** 2)) * Qf(q * yk)
            bethe = max(bethe, abs(lhs + Qf(yk / q)) / abs(lhs))
            th = cmath.log(yk) / LAM
            # reduce Im theta into [pi/2 - pi/(2 lam), pi/2 + pi/(2 lam))
            per = math.pi / LAM
            im = (th.imag - (math.pi / 2 - per / 2)) % per + (math.pi / 2 - per / 2)
            roots.append((th.real, im, st))
        yv = cmath.exp(LAM * (0.31 - 0.2j))
        ydb = -1j * yv
        tval = (Vi @ t_op(yv, om)[np.ix_(idx, idx)] @ V)[st, st]
        z2 = (YS / ydb) ** 2
        bax = q ** (-N / 2) * (np.prod((1 + z2) / (1 + z2 / q ** 2)) * Qf(q * ydb) + Qf(ydb / q)) / Qf(ydb)
        baxter = max(baxter, abs(tval - bax) / abs(tval))
    data[m] = (polys, roots)

report("D is diagonal on the eigenvectors of t_omega in the sectors m = 0, 2", ok_diag)
report("the eigenvalues of Q are polynomials in 1/y_d^2 of degree (N - m)/2", ok_deg)
report(f"the plotted zeros satisfy the Bethe equations (relative {bethe:.0e})", bethe < 1e-10)
report(f"Baxter's formula gives the eigenvalues of t_omega (relative {baxter:.0e})", baxter < 1e-10)

# (b): three neutral states with zeros on the line Im theta_d = pi/2
polys0, roots0 = data[0]
on_line = {}
for (re, im, st) in roots0:
    if abs(im - math.pi / 2) < 1e-6:
        on_line.setdefault(st, []).append(re)
chosen = sorted(on_line, key=lambda st: -len(on_line[st]))[:3]
xs = np.linspace(-2.5, 2.5, 4001)
xs = np.sort(np.concatenate([xs] + [np.array(on_line[st]) + dx for st in chosen for dx in (-2e-5, 2e-5)]))
curves = []
ok_line = True
for st in chosen:
    cM = polys0[st]
    wl = np.exp(-2 * LAM * (xs + 0.5j * math.pi))
    Qline = np.polyval(cM[::-1], wl)
    norm = sum(abs(c) * np.abs(wl) ** j for j, c in enumerate(cM))  # sum_j |c_j| |w|^j >= |Q|
    curves.append(np.abs(Qline) / norm)
    for x0 in on_line[st]:
        w0 = cmath.exp(-2 * LAM * (x0 + 0.5j * math.pi))
        val = abs(np.polyval(cM[::-1], w0)) / sum(abs(c) * abs(w0) ** j for j, c in enumerate(cM))
        ok_line = ok_line and val < 1e-9
report("the zeros marked in (b) are zeros of the plotted |Q| on Im theta_d = pi/2", ok_line)

fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw=dict(width_ratios=[1.25, 1.0]))
per = math.pi / LAM
lo, hi = math.pi / 2 - per / 2, math.pi / 2 + per / 2
axa.axhspan(lo, hi, color=style.LIGHTGREY, lw=0, zorder=0)
for yv in (lo, math.pi / 2, hi):
    axa.axhline(yv, color=style.GREY, lw=0.6, ls=':' if yv != math.pi / 2 else '--', zorder=1)
r0 = np.array([(r, i) for r, i, _ in data[0][1]])
r2 = np.array([(r, i) for r, i, _ in data[2][1]])
axa.scatter(r0[:, 0], r0[:, 1], s=16, facecolors='none', edgecolors=style.BLUE, lw=0.9, label='$m=0$', zorder=3)
axa.scatter(r2[:, 0], r2[:, 1], s=16, marker='x', color=style.RED, lw=0.9, label='$m=2$', zorder=3)
for th in THS:
    axa.plot([th, th], [math.pi / 2 - 0.05, math.pi / 2 + 0.05], color=style.GREEN, lw=1.6, zorder=2)
axa.set_xlim(-2.6, 2.6)
axa.set_ylim(lo - 0.08, hi + 0.08)
axa.set_yticks([lo, math.pi / 2, hi])
axa.set_yticklabels([r'$\frac{\pi}{2}-\frac{\pi}{2\lambda}$', r'$\frac{\pi}{2}$', r'$\frac{\pi}{2}+\frac{\pi}{2\lambda}$'])
axa.set_xlabel(r'$\mathrm{Re}\,\theta_d$')
axa.set_ylabel(r'$\mathrm{Im}\,\theta_d$')
axa.legend(loc='center right', fontsize=7.5, handletextpad=0.3)
axa.set_title('(a) zeros of the defect operator', fontsize=8.5, color=style.GREY)
cols = (style.BLUE, style.ORANGE, style.PURPLE)
for c, cur, st in zip(cols, curves, chosen):
    axb.semilogy(xs, cur, color=c, lw=1.2)
    for x0 in on_line[st]:
        axb.plot([x0], [1e-3], marker='^', color=c, ms=4, clip_on=False)
for th in THS:
    axb.axvline(th, color=style.GREEN, lw=0.6, ls=':')
axb.set_ylim(1e-3, 2)
axb.set_xlim(-2.5, 2.5)
axb.set_xlabel(r'$x$')
axb.set_ylabel(r'$|Q|\,/\,\Sigma_j|c_jw^j|$')
axb.set_title(r'(b) along $\mathrm{Im}\,\theta_d=\frac{\pi}{2}$', fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=2.0)
style.save(fig, "l8_bethe", OUT)
raise SystemExit(0 if style.all_ok() else 1)
