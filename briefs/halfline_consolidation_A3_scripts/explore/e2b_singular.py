"""Diagnose the three-soliton configurations whose boundary value Re phi(0,t) jumps: locate zeros of tau_j on x <= 0."""
import numpy as np
from booklib import load_images

I = load_images()
Sol, interaction, roots_an = I['Sol'], I['interaction'], I['roots_an']
beta = I['beta']


def three(n, a, b, th, xi1, y3):
    h = n + 1
    Aab = interaction(a, b, th, 0.0, h).real; Aaa = interaction(a, a, th, -th, h).real
    xi3 = -np.log(Aab) + 1j * y3
    xi2 = -np.conj(xi1) - np.log(Aaa) - np.log(Aab)
    return Sol(n, [a, a, b], [th, -th, 0.0], [xi1, xi2, xi3])


def term_sizes(sol, x, t):
    x = np.atleast_1d(np.asarray(x, float))
    k = sol.ma[:, None] * np.cosh(sol.th)[:, None]; w = sol.ma[:, None] * np.sinh(sol.th)[:, None]
    expo = k * x[None, :] - w * t + sol.xi[:, None]
    tot = np.zeros(x.size)
    for S, c in zip(sol.subsets, sol.coef):
        tot += 1.0 if not S else np.abs(c * np.exp(expo[list(S)].sum(0)))
    return tot


def regularity(sol, xs, ts):
    worst = np.inf; where = None
    for t in ts:
        T, _ = sol.taus(xs, t)
        r = np.min(np.abs(T), axis=0) / term_sizes(sol, xs, t)
        i = np.argmin(r)
        if r[i] < worst:
            worst, where = r[i], (xs[i], t)
    return worst, where


xs = np.linspace(-20, 20, 4001); ts = np.linspace(-12, 12, 2401)
for (n, a, b, th, xi1, y3) in [(2, 1, 1, 0.5, 1.2 + 0.3j, 0.9), (2, 2, 1, 0.6, 0.8 - 0.4j, 0.5), (3, 1, 2, 0.5, 1.0 + 0.2j, 0.8)]:
    s = three(n, a, b, th, xi1, y3)
    r, wh = regularity(s, xs, ts)
    print(f'a_{n} (a={a}, b={b}, th={th}, xi1={xi1}, y3={y3}): min_j |tau_j|/sum|terms| = {r:.2e} at (x,t) = {wh}')
    # boundary value over time
    xg = np.linspace(-25, 0, 5001)
    vals = [(t, np.round(s.phi(xg, t)[0][:, -1].real * beta / np.pi, 3)) for t in np.linspace(-6, 6, 13)]
    print('   beta Re phi(0,t)/pi:', '; '.join(f't={t:+.0f}: {v}' for t, v in vals[::2]))

# scan phases for regular configurations
print('--- phase scans (min |tau|/sum|terms| over x in [-20,20], t in [-12,12])')
xs = np.linspace(-20, 20, 801); ts = np.linspace(-12, 12, 481)
for (n, a, b, th) in [(2, 2, 1, 0.6), (3, 1, 2, 0.5), (3, 2, 1, 0.45), (2, 1, 1, 0.5), (3, 1, 1, 0.5), (3, 2, 2, 0.5)]:
    best = (-1, None)
    for p1 in np.linspace(0, 2 * np.pi, 13)[:-1]:
        for y3 in np.linspace(0, 2 * np.pi, 13)[:-1]:
            s = three(n, a, b, th, 1.0 + 1j * p1, y3)
            r, _ = regularity(s, xs, ts)
            if r > best[0]:
                best = (r, (p1, y3))
    print(f'a_{n} (a={a}, b={b}, th={th}): best regularity {best[0]:.3f} at (Im xi1, y3) = {np.round(best[1], 3)}')
