"""Figure F34 (chapter 48, prp-krein-fundamental, sec-sector-status): eigenvalues of the Bethe-Yang transfer matrix
T_1 of three fundamental solitons of a_2^(1) (Jimbo R-matrix, q = e^{i mu}, real rapidities), and the breaking
region of two kinks of a_2^(2).

Left: eigenvalues s of T_1 in the complex plane, for the sector with all three colours (6 states) and for a sector with two
colours (3 states).  The two-colour eigenvalues are phases; the three-colour ones are pure phases or pairs (s, 1/conj s).
Right: the (psi, alpha) plane of the Theta-broken region of two a_2^(2) kinks in the sector Q = 0,
|alpha| < max(0, min(2|psi| - pi, pi - |psi|)) (table of sec-sector-status).

Checks: T_1^dagger = C T_1^{-1} C (part 2); spectrum invariant under s -> 1/conj s; two-colour sector eigenvalues
have modulus 1; three-colour sector has max||s|-1| > 0.05; R~ Hermiticity and colour-reversal identities (part 1);
the broken region is empty for |psi| <= pi/2, its largest alpha is pi/3 at |psi| = 2 pi/3, and it vanishes at |psi| = pi.
Writes part9-foundations/figures/ch48_f34_krein_pairing.svg."""
import itertools

import numpy as np
import style
from style import plt, report

style.setup()
n, N, mu = 3, 3, 1.0
q = np.exp(1j * mu)


def Rt(x):
    """normalised Jimbo R-matrix R(x)/(xq - 1/(xq)) on C^n (x) C^n"""
    M = np.zeros((n * n, n * n), complex)
    for i in range(n):
        for j in range(n):
            if i == j:
                M[i * n + i, i * n + i] = x * q - 1 / (x * q)
            else:
                M[i * n + j, i * n + j] = x - 1 / x
                M[i * n + j, j * n + i] = (q - 1 / q) * (x if i < j else 1 / x)
    return M / (x * q - 1 / (x * q))


def embed(M, i, j):
    D = n ** N
    out = np.zeros((D, D), complex)
    for idx in itertools.product(range(n), repeat=N):
        col = np.ravel_multi_index(idx, [n] * N)
        a, b = idx[i], idx[j]
        for a2 in range(n):
            for b2 in range(n):
                v = M[a2 * n + b2, a * n + b]
                if v != 0:
                    l = list(idx); l[i] = a2; l[j] = b2
                    out[np.ravel_multi_index(l, [n] * N), col] += v
    return out


theta = np.array([-0.85, -0.68, -0.32])
xs = np.exp(4 * theta)
T = np.eye(n ** N, dtype=complex)
for j in range(1, N):
    T = embed(Rt(xs[0] / xs[j]), 0, j) @ T

C = np.eye(n)[::-1]
CN = np.kron(np.kron(C, C), C)
P = np.zeros((n * n, n * n))
for i in range(n):
    for j in range(n):
        P[i * n + j, j * n + i] = 1
CC = np.kron(C, C)
y = 0.7 + 0.4j
report("part 1: Rt(y)^dagger = Rt(1/conj y) for complex y",
       np.abs(Rt(y).conj().T - Rt(1 / np.conj(y))).max() < 1e-12)
report("part 1: (C x C) Rt_12 (C x C) = Rt_21", np.abs(CC @ Rt(1.7) @ CC - P @ Rt(1.7) @ P).max() < 1e-12)
report("part 2: T_1^dagger = C T_1^{-1} C", np.abs(T.conj().T - CN @ np.linalg.inv(T) @ CN).max() < 1e-11)


def sector(content):
    idx = [k for k, c in enumerate(itertools.product(range(n), repeat=N)) if tuple(sorted(c)) == content]
    return np.linalg.eigvals(T[np.ix_(idx, idx)])


ev3 = sector((0, 1, 2))
ev2 = sector((0, 0, 1))
ev2b = sector((1, 2, 2))
pair = lambda ev: max(min(abs(e - 1 / np.conj(f)) for f in ev) for e in ev)
report("three-colour sector (6 states): spectrum invariant under s -> 1/conj s", pair(ev3) < 1e-9)
report("three-colour sector: max ||s|-1| > 0.05 (breaking)", np.abs(np.abs(ev3) - 1).max() > 0.05)
report("two-colour sectors (3 states): all eigenvalues have modulus 1",
       np.abs(np.abs(ev2) - 1).max() < 1e-9 and np.abs(np.abs(ev2b) - 1).max() < 1e-9)

# right panel
psi = np.linspace(-np.pi, np.pi, 2001)
amax = np.maximum(0, np.minimum(2 * np.abs(psi) - np.pi, np.pi - np.abs(psi)))
report("a_2^(2): broken region empty for |psi| <= pi/2", np.all(amax[np.abs(psi) <= np.pi / 2] == 0))
report("a_2^(2): largest alpha is pi/3, attained at |psi| = 2 pi/3",
       abs(amax.max() - np.pi / 3) < 1e-2 and abs(abs(psi[np.argmax(amax)]) - 2 * np.pi / 3) < 1e-2)
report("a_2^(2): region closes at |psi| = pi", abs(amax[-1]) < 1e-12)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [1, 1.05]})
t = np.linspace(0, 2 * np.pi, 400)
ax1.plot(np.cos(t), np.sin(t), color=style.GREY, lw=0.8)
ax1.axhline(0, color=style.LIGHTGREY, lw=0.6); ax1.axvline(0, color=style.LIGHTGREY, lw=0.6)
ax1.scatter(ev2.real, ev2.imag, s=26, color=style.BLUE, zorder=3, label="two colours")
ax1.scatter(ev3.real, ev3.imag, s=26, color=style.RED, zorder=3, label="three colours")
used = set()
for i, e in enumerate(ev3):
    if abs(abs(e) - 1) > 0.05:
        j = int(np.argmin([abs(e - 1 / np.conj(f)) for f in ev3]))
        if (j, i) not in used:
            used.add((i, j))
            ax1.plot([e.real, ev3[j].real], [e.imag, ev3[j].imag], color=style.RED, lw=0.8, ls="--")
ax1.set_aspect("equal")
lim = 1.15 * max(1, np.abs(ev3).max(), 1 / np.abs(ev3).min())
ax1.set_xlim(-lim, lim); ax1.set_ylim(-lim, lim)
ax1.set_xlabel(r"$\mathrm{Re}\,s$"); ax1.set_ylabel(r"$\mathrm{Im}\,s$")
ax1.set_title(r"$a_2^{(1)}$: eigenvalues of $T_1$, three solitons", fontsize=8.5)
ax1.legend(loc="upper right", fontsize=7, handletextpad=0.2, borderaxespad=0.1)

ax2.fill_between(psi / np.pi, -amax / np.pi, amax / np.pi, color=style.RED, alpha=0.55, lw=0)
edge = np.where(amax > 0, amax / np.pi, np.nan)
ax2.plot(psi / np.pi, edge, color=style.RED, lw=1); ax2.plot(psi / np.pi, -edge, color=style.RED, lw=1)
ax2.text(0.0, 0.0, "unbroken", ha="center", va="center", fontsize=8, color=style.BLUE)
ax2.text(0.0, 0.5, r"$Q\neq0$: unbroken for every $\alpha$", ha="center", va="center", fontsize=7, color=style.GREY)
ax2.text(0.76, 0.0, "broken", ha="center", va="center", fontsize=8, color=style.RED)
ax2.set_xlim(-1, 1); ax2.set_ylim(-0.6, 0.6)
ax2.set_xlabel(r"$\psi/\pi$"); ax2.set_ylabel(r"Bloch angle $\alpha/\pi$")
ax2.set_title(r"$a_2^{(2)}$: two kinks, sector $Q=0$", fontsize=8.5)
fig.tight_layout()
style.save(fig, "ch48_f34_krein_pairing", style.OUT_PART9)
raise SystemExit(0 if style.all_ok() else 1)
