"""Figure F13 (chapter 9, sec-dorey-rule): Dorey's rule for d4^(1).

The 24 roots of D4 projected onto the Coxeter plane of the bicoloured Coxeter element w = s_black s_white (centre
node white, the three outer nodes black; h = 6).  The four orbits Omega_a = {w^p gamma_a} are regular hexagons;
their radii are proportional to the masses m_a (outer nodes 1, centre sqrt 3).  One closed triangle
gamma_1 + w^p gamma_1 + w^q gamma_c = 0 is drawn: its sides are m_1, m_1, m_c, its angles multiples of pi/h.

Checks: gamma_a is alpha_a for white a and s_white alpha_a for black a; the four orbits have h = 6 elements, are
disjoint and cover the 24 roots; w has order 6 and rotates the projection by 2 pi/6; the radii are proportional to
the Perron-Frobenius vector of the Cartan matrix; the triple sums to zero in R^4; the projected sides equal the
masses and the interior angles are multiples of pi/6.
Writes part3-classical/figures/ch9_f13_dorey_d4.svg.
"""
import itertools

import numpy as np

import style
from style import plt, report

style.setup()

# nodes 0,1,2 = outer (black), 3 = centre (white)
C = 2 * np.eye(4)
for i in range(3):
    C[i, 3] = C[3, i] = -1
S = np.linalg.cholesky(C)                                   # rows: simple roots in R^4


def refl(a):
    return np.eye(len(a)) - 2 * np.outer(a, a) / (a @ a)


def roots_from_cartan(A):
    r = len(A)
    simple = [tuple(int(i == j) for j in range(r)) for i in range(r)]
    roots, todo = set(simple), list(simple)
    while todo:
        b = todo.pop()
        for i in range(r):
            nb = list(b)
            nb[i] -= int(round(sum(b[j] * A[i][j] for j in range(r))))
            nb = tuple(nb)
            if nb not in roots:
                roots.add(nb)
                todo.append(nb)
    return sorted(roots)


R = np.array(roots_from_cartan(C.tolist())) @ S
report("D4: 24 roots of length^2 2", len(R) == 24 and np.allclose((R**2).sum(1), 2))
s_white = refl(S[3])
W = np.eye(4)
for i in range(3):
    W = refl(S[i]) @ W                                      # black reflections commute
W = W @ s_white                                             # w = s_black s_white
report("D4: Coxeter element has order h = 6",
       np.allclose(np.linalg.matrix_power(W, 6), np.eye(4)) and
       not any(np.allclose(np.linalg.matrix_power(W, k), np.eye(4)) for k in range(1, 6)))
gam = [s_white @ S[i] for i in range(3)] + [S[3]]           # black: s_white alpha_a ; white: alpha_a
report("gamma_a = alpha_c + alpha_a for the black nodes", all(np.allclose(gam[i], S[3] + S[i]) for i in range(3)))
key = lambda x: tuple(np.round(x, 6) + 0.0)
orbits = [[np.linalg.matrix_power(W, p) @ g for p in range(6)] for g in gam]
keys = [{key(x) for x in o} for o in orbits]
report("four orbits of 6 distinct roots, disjoint, covering all 24 roots",
       all(len(k) == 6 for k in keys) and len(set().union(*keys)) == 24 and
       {key(x) for x in R} == set().union(*keys))
ev, evec = np.linalg.eig(W)
v = evec[:, np.argmin(np.abs(ev - np.exp(2j * np.pi / 6)))]
v /= np.linalg.norm(v)
proj = lambda x: x @ v.conj()
z = [np.array([proj(x) for x in o]) for o in orbits]
report("w rotates the projection by 2 pi/6",
       all(np.allclose(zz[1] / zz[0], np.exp(2j * np.pi / 6)) for zz in z))
rad = np.array([abs(zz[0]) for zz in z])
pf = np.abs(np.linalg.eigh(C)[1][:, 0])
report("orbit radii proportional to the Perron-Frobenius vector (1, 1, 1, sqrt 3)",
       np.allclose(rad / pf, (rad / pf)[0]) and np.isclose(rad[3] / rad[0], np.sqrt(3)))
report("sum rule sum_a m_a^2 = 2 h m^2 = 12 with m_outer^2 = 2, m_centre^2 = 6",
       np.isclose(((rad / rad[0]) ** 2 * 2).sum(), 12))

# a closed triangle: outer orbit 0 twice, and the centre orbit (1 + 1 -> sqrt 3)
trip = None
for p, q in itertools.product(range(6), repeat=2):
    a, b, c = orbits[0][0], orbits[0][p], orbits[3][q]
    if np.allclose(a + b + c, 0):
        trip = (a, b, c)
        break
report("triple gamma_1 + w^p gamma_1 + w^q gamma_c = 0 exists", trip is not None)
za, zb, zc = (proj(x) for x in trip)
sides = np.array([abs(za), abs(zb), abs(zc)])
report("projected sides are m_1 : m_1 : m_c = 1 : 1 : sqrt 3",
       np.allclose(sides / sides[0], [1, 1, np.sqrt(3)]))
ang = lambda x, y: np.arccos(np.clip((-x * np.conj(y)).real / abs(x) / abs(y), -1, 1))
angles = [ang(za, zb), ang(zb, zc), ang(zc, za)]
report("interior angles are multiples of pi/h = pi/6 and sum to pi",
       np.allclose(np.array(angles) / (np.pi / 6), np.round(np.array(angles) / (np.pi / 6))) and
       np.isclose(sum(angles), np.pi))

# ---- figure ----
fig, ax = plt.subplots(figsize=(4.6, 4.6))
cols = [style.BLUE, style.BLUE, style.ORANGE, style.RED]
names = [None, None, r"$\Omega_1,\Omega_3,\Omega_4$ (coincide, $m=1$)", r"$\Omega_c$ ($m=\sqrt{3}$)"]
for k in range(4):
    pts = np.concatenate([z[k], z[k][:1]])
    ax.plot(pts.real, pts.imag, color=cols[k], lw=0.9, alpha=0.6) if k >= 2 else None
    (ax.scatter(z[k].real, z[k].imag, s=26, color=cols[k], zorder=3, label=names[k]) if k >= 2 else None)
tri = np.array([0, za, za + zb, 0])
ax.fill(tri.real, tri.imag, color=style.GREY, alpha=0.25, lw=0)
for s, t, col in ((0, za, style.TEAL), (za, za + zb, style.GREEN), (za + zb, 0, cols[3])):
    ax.annotate("", xy=(t.real, t.imag), xytext=(s.real, s.imag),
                arrowprops=dict(arrowstyle="->", color=col, lw=1.6, shrinkA=0, shrinkB=0))
ax.set_aspect("equal")
ax.axis("off")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.06), fontsize=8, ncol=2, handletextpad=0.2, columnspacing=0.8)
ax.set_title(r"$d_4^{(1)}$: orbits $\Omega_a$ in the Coxeter plane; $\gamma_1+\gamma_1^\prime+\gamma_c=0$", fontsize=8.5)
fig.tight_layout()
style.save(fig, "ch9_f13_dorey_d4", style.OUT_PART3)
raise SystemExit(0 if style.all_ok() else 1)
