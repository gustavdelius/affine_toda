"""Figure F30 (chapter 35, sec-boundary-spectrum-a2, eq-a2-energies, tbl-a2-boundary-states): the boundary states
b_{n,m} of a_2^(1) with a solitonic boundary.

Each bootstrap step in which a particle binds at theta = i pi p/3 raises the energy by m_1 cos(pi p/3):
  b_{n,m} -> b_{n+1,m}   (particle 1, p = 1/2 - (2n+1)B/2),
  b_{n,m} -> b_{n,m+1}   (particle 2, p = 1/2 - (2m+1)B/2),
  b_{-n,n} -> b_{-n-1,n+1} (particle 1, p = 3/2 - (2n+1)B/2).
The closed form is e(b_{n,m}) - e(b_{0,0}) = m_1 [F(n) + F(m)], F(k) = sin(k d) cos(pi/6 - k d)/sin d, d = pi B/6.
Existing states: n + m > 0, |n|, |m| < 1/(2B) + 1/2, and the diagonal tower b_{n,-n}, b_{-n,n}, 0 <= n < 3/(2B) + 1/2.

Checks: path sums of the steps reproduce the closed form; the table of tbl-a2-boundary-states at B = 9/20 (ten states);
e(b_{n,-n}) = m_1 sin^2(n d)/sin d, tends to the rotor spectrum m_1 (pi B/6) n^2 as B -> 0.

Writes part7-half-line/figures/ch36_f30_a2_boundary_states.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()


def F(k, B):
    d = np.pi * B / 6
    return np.sin(k * d) * np.cos(np.pi / 6 - k * d) / np.sin(d)


def E(n, m, B):
    return F(n, B) + F(m, B)


def exists(n, m, B):
    if n + m == 0:
        return abs(n) < 3 / (2 * B) + 0.5
    return n + m > 0 and abs(n) < 1 / (2 * B) + 0.5 and abs(m) < 1 / (2 * B) + 0.5


# ---- checks ----------------------------------------------------------------
B = 0.45
step1 = lambda n: np.cos(np.pi * (0.5 - (2 * n + 1) * B / 2) / 3)
step2 = step1
stepd = lambda n: np.cos(np.pi * (1.5 - (2 * n + 1) * B / 2) / 3)
ok = True
for n in range(0, 4):
    for m in range(0, 4):
        path = sum(step1(k) for k in range(n)) + sum(step2(k) for k in range(m))
        ok &= abs(path - E(n, m, B)) < 1e-12
    ok &= abs(sum(stepd(k) for k in range(n)) - E(-n, n, B)) < 1e-12
report("sums of the bootstrap steps (particle 1, particle 2, diagonal step) equal F(n) + F(m) (n, m = 0..3, B = 9/20)", ok)
table = {(0, 0): 0.0, (-1, 1): 0.2334, (1, -1): 0.2334, (-2, 2): 0.8829, (2, -2): 0.8829, (1, 0): 0.9588, (0, 1): 0.9588,
         (-3, 3): 1.807, (3, -3): 1.807, (1, 1): 1.918}
states = [(n, m) for n in range(-6, 7) for m in range(-6, 7) if exists(n, m, 9 / 20)]
report("at B = 9/20 the existence rule gives exactly the ten states of tbl-a2-boundary-states", set(states) == set(table))
report("their energies agree with the table to its precision", all(abs(E(n, m, 9 / 20) - e) < 6e-4 for (n, m), e in table.items()))
report("e(b_{n,-n}) = m_1 sin^2(n d)/sin d for n = 1..5 and several B",
       all(abs(E(n, -n, b) - np.sin(n * np.pi * b / 6) ** 2 / np.sin(np.pi * b / 6)) < 1e-12 for n in range(1, 6) for b in (0.1, 0.45, 0.9)))
report("rotor limit: e(b_{n,-n}) / [m_1 (pi B/6) n^2] -> 1 as B -> 0",
       all(abs(E(n, -n, 1e-4) / (np.pi * 1e-4 / 6 * n**2) - 1) < 1e-3 for n in (1, 2, 3)))
report("e(b_{-1,0}) = -m_1 cos(pi/6 + d) < 0 (below the ground state: excluded by n + m >= 0)",
       all(abs(E(-1, 0, b) + np.cos(np.pi / 6 + np.pi * b / 6)) < 1e-12 and E(-1, 0, b) < 0 for b in (0.1, 0.45, 0.9)))

# ---- figure ----------------------------------------------------------------
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [1, 1.15]})
cm = plt.get_cmap("viridis")
Bf = 9 / 20
OFF = {(1, 0): (16, -2), (0, 1): (0, 9), (1, 1): (0, 9), (0, 0): (14, -8)}
for (n, m) in states:
    diag = n + m == 0
    e = E(n, m, Bf)
    ax.scatter([n], [m], s=70, color=style.RED if diag else style.BLUE, zorder=3, alpha=.9)
    ax.annotate(f"{e:.2f}", (n, m), textcoords="offset points", xytext=OFF.get((n, m), (9, 5)), ha="center", fontsize=6.8)
for (n, m), (n2, m2) in (((0, 0), (1, 0)), ((0, 0), (0, 1)), ((1, 0), (1, 1)), ((0, 1), (1, 1)), ((0, 0), (-1, 1)), ((0, 0), (1, -1)),
                          ((-1, 1), (-2, 2)), ((1, -1), (2, -2)), ((-2, 2), (-3, 3)), ((2, -2), (3, -3))):
    ax.annotate("", xy=(n2, m2), xytext=(n, m), arrowprops=dict(arrowstyle="->", color=style.GREY, lw=.9, shrinkA=5, shrinkB=5))
ax.set_xlim(-3.7, 3.7)
ax.set_ylim(-3.7, 3.7)
ax.set_aspect("equal")
ax.set_xlabel(r"$n$ (particles of type 1)")
ax.set_ylabel(r"$m$ (type 2)")
ax.set_title(r"states $b_{n,m}$ and $(e-e_{0,0})/m_1$ at $B=9/20$", fontsize=8.2)
ax.axhline(0, color=style.GREY, lw=.5)
ax.axvline(0, color=style.GREY, lw=.5)

Bs = np.linspace(0.2, 1.0, 1200)
for n in range(-7, 8):
    for m in range(-7, 8):
        if n + m < 0:
            continue
        ys = np.array([E(n, m, b) if exists(n, m, b) else np.nan for b in Bs])
        if np.all(np.isnan(ys)):
            continue
        ax2.plot(Bs, ys, color=style.RED if n + m == 0 else style.BLUE, lw=1.1 if n + m == 0 else .8, alpha=.95 if n + m == 0 else .7)
ax2.axvline(Bf, color=style.GREY, lw=.8, ls="--")
ax2.set_xlim(0.2, 1.0)
ax2.set_ylim(0, 2.6)
ax2.set_xlabel("$B$")
ax2.set_ylabel(r"$(e-e_{0,0})/m_1$")
ax2.set_title(r"spectrum against $B$", fontsize=8.2)
ax2.plot([], [], color=style.RED, label=r"diagonal tower $b_{n,-n}$ (rotor)")
ax2.plot([], [], color=style.BLUE, label=r"$b_{n,m}$, $n+m>0$")
ax2.legend(fontsize=7, loc="lower right", bbox_to_anchor=(1.0, 0.1))
fig.tight_layout()
style.save(fig, "ch36_f30_a2_boundary_states", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
