"""Figure F10 (chapter 5, sec-roots, sec-weyl-coxeter): root systems and the Coxeter plane.

(a) The G2 roots with the simple roots, the highest root and the fundamental Weyl chamber; the 12 chambers
    are cut out by the lines orthogonal to the six positive roots.
(b) The roots of E8 projected onto the Coxeter plane of the bicoloured Coxeter element: eight regular
    30-gons, one per orbit of w.  The radii are proportional to the Perron-Frobenius eigenvector of the
    Cartan matrix, i.e. to the classical masses of the e8^(1) Toda theory (sec-perron-frobenius).

Checks: G2 has 12 roots (6 long of length^2 2, 6 short of length^2 2/3), highest root 3a1+2a2, h=6, h^vee=4,
Coxeter element = rotation by 2 pi/6 with two orbits of 6; E8 has 240 roots, w has order 30, the projections
lie on 8 circles with 30 roots each, w rotates them by 2 pi/30, and the radii of the simple roots are
proportional to the Perron-Frobenius eigenvector.
Writes part2-algebra/figures/ch5_f10_root_systems.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()


def roots_from_cartan(A):
    """All roots as integer vectors in the simple-root basis, by closing under the simple reflections.
    A[i][j] = alpha_i^vee . alpha_j, and s_i(beta) = beta - (alpha_i^vee . beta) alpha_i."""
    r = len(A)
    simple = [tuple(int(i == j) for j in range(r)) for i in range(r)]
    roots, todo = set(simple), list(simple)
    while todo:
        b = todo.pop()
        for i in range(r):
            nb = list(b)
            nb[i] -= sum(b[j] * A[i][j] for j in range(r))
            nb = tuple(nb)
            if nb not in roots:
                roots.add(nb)
                todo.append(nb)
    return sorted(roots)


def refl(a):
    a = np.asarray(a, float)
    return np.eye(len(a)) - 2 * np.outer(a, a) / (a @ a)


# ---- G2 -------------------------------------------------------------------------------------------------
a1 = np.sqrt(2 / 3) * np.array([1.0, 0.0])                                   # short
a2 = np.sqrt(2) * np.array([np.cos(5 * np.pi / 6), np.sin(5 * np.pi / 6)])  # long
simple = np.array([a1, a2])
cartan = np.array([[2 * (x @ y) / (x @ x) for y in simple] for x in simple])  # a_ij = alpha_i^vee . alpha_j
report("G2: Cartan matrix [[2,-3],[-1,2]]", np.allclose(cartan, [[2, -3], [-1, 2]]))
cf = np.array(roots_from_cartan(np.rint(cartan).astype(int).tolist()))
R = cf @ simple
n2 = (R**2).sum(1)
report("G2: 12 roots, 6 of length^2 2 and 6 of length^2 2/3",
       len(R) == 12 and np.isclose(n2, 2).sum() == 6 and np.isclose(n2, 2 / 3).sum() == 6)
top = cf[cf.sum(1).argmax()]
theta = top @ simple
report("G2: highest root is 3 a1 + 2 a2, of length^2 2", list(top) == [3, 2] and np.isclose(theta @ theta, 2))
marks_dual = top * np.array([a @ a / 2 for a in simple])
report("G2: h = 1 + 3 + 2 = 6 and h^vee = 1 + 1 + 2 = 4", 1 + top.sum() == 6 and np.isclose(1 + marks_dual.sum(), 4))
w = refl(a1) @ refl(a2)
report("G2: s1 s2 is the rotation by 2 pi/6", np.isclose(np.linalg.det(w), 1) and
       np.isclose(abs(np.arctan2(w[1, 0], w[0, 0])), 2 * np.pi / 6))
key = lambda x: tuple(np.round(x, 6) + 0.0)
orbits, seen = [], set()
for r in R:
    if key(r) in seen:
        continue
    o, x = set(), r.copy()
    while key(x) not in o:
        o.add(key(x))
        x = w @ x
    seen |= o
    orbits.append(len(o))
report("G2: w has two orbits on the roots, each of h = 6 elements", orbits == [6, 6])

# ---- E8 -------------------------------------------------------------------------------------------------
# Bourbaki labelling: chain 1-3-4-5-6-7-8 with node 2 attached to node 4.
edges = [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (2, 4)]
C = 2 * np.eye(8, dtype=int)
for i, j in edges:
    C[i - 1, j - 1] = C[j - 1, i - 1] = -1
S = np.linalg.cholesky(C.astype(float))                    # rows: simple roots in R^8 with Gram matrix C
E8 = np.array(roots_from_cartan(C.tolist())) @ S
report("E8: 240 roots of length^2 2", len(E8) == 240 and np.allclose((E8**2).sum(1), 2))
colour, stack = {1: 0}, [1]
while stack:
    u = stack.pop()
    for i, j in edges:
        for a, b in ((i, j), (j, i)):
            if a == u and b not in colour:
                colour[b] = 1 - colour[u]
                stack.append(b)
Wc = np.eye(8)
for c in (1, 0):                                           # w = s_black s_white
    for i in range(1, 9):
        if colour[i] == c:
            Wc = refl(S[i - 1]) @ Wc
order_ok = np.allclose(np.linalg.matrix_power(Wc, 30), np.eye(8)) and \
    not any(np.allclose(np.linalg.matrix_power(Wc, k), np.eye(8)) for k in range(1, 30))
report("E8: bicoloured Coxeter element has order h = 30", order_ok)
ev, evec = np.linalg.eig(Wc)
v = evec[:, np.argmin(np.abs(ev - np.exp(2j * np.pi / 30)))]
v = v / np.linalg.norm(v)
z = E8 @ v.conj()                                          # complex coordinate in the Coxeter plane
z_rot = (E8 @ Wc.T) @ v.conj()                             # projection of w(alpha)
nz = np.abs(z) > 1e-9
rot = z_rot[nz] / z[nz]
report("E8: w rotates the projection by 2 pi/30",
       np.allclose(rot, rot[0]) and np.isclose(abs(np.angle(rot[0])), 2 * np.pi / 30))
rad = np.abs(z)
radii = np.unique(np.round(rad, 6))
counts = [int((np.abs(rad - r) < 1e-5).sum()) for r in radii]
report("E8: the projections lie on 8 circles", len(radii) == 8)
report("E8: each circle carries 30 roots (a regular 30-gon)", counts == [30] * 8)
lam, vec = np.linalg.eigh(C.astype(float))
pf = np.abs(vec[:, 0])
rs = np.array([rad[np.argmin(np.abs(E8 - S[i]).sum(1))] for i in range(8)])
report("E8: radii of the simple roots are proportional to the Perron-Frobenius eigenvector",
       np.allclose(rs / pf, (rs / pf)[0]))
report("E8: the 8 circles are those of the 8 simple roots", np.allclose(np.sort(rs), radii, atol=1e-5))

# ---- figure ---------------------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.5))

lim = 1.75


def perp(a):
    return np.array([-a[1], a[0]]) / np.linalg.norm(a)


for r in R[(cf >= 0).all(1)]:                              # mirrors: lines orthogonal to the positive roots
    p = perp(r)
    ax1.plot([-lim * p[0], lim * p[0]], [-lim * p[1], lim * p[1]], color=style.GREY, lw=0.7, ls=(0, (3, 3)))
d1 = perp(a2)
d1 = d1 if d1 @ a1 > 0 else -d1
d2 = perp(a1)
d2 = d2 if d2 @ a2 > 0 else -d2
t1, t2 = np.arctan2(d1[1], d1[0]), np.arctan2(d2[1], d2[0])
t2 = t2 + 2 * np.pi * np.round((t1 - t2) / (2 * np.pi))
arc = np.linspace(t1, t2, 30)
wedge = np.vstack([[0, 0], lim * np.c_[np.cos(arc), np.sin(arc)]])
ax1.fill(wedge[:, 0], wedge[:, 1], color=style.BLUE, alpha=0.18, lw=0)
for r in R:
    thick = any(np.allclose(r, s) for s in (a1, a2, theta))
    col = style.RED if np.isclose(r @ r, 2) else style.GREEN
    ax1.annotate("", xy=r, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=2.0 if thick else 0.9,
                                                          shrinkA=0, shrinkB=0, mutation_scale=7))
ax1.text(a1[0], a1[1] - 0.12, r"$\vec\alpha_1$", color=style.GREEN, ha="center", va="top")
ax1.text(a2[0] - 0.08, a2[1] - 0.02, r"$\vec\alpha_2$", color=style.RED, ha="right", va="center")
ax1.text(theta[0] + 0.06, theta[1], r"$\vec\theta$", color=style.RED, ha="left", va="center")
mid = (t1 + t2) / 2
ax1.annotate("fundamental\nchamber", xy=(1.1 * np.cos(mid), 1.1 * np.sin(mid)), xytext=(1.2, 1.5), color=style.BLUE,
             ha="left", va="center", fontsize=7.5,
             arrowprops=dict(arrowstyle="-", color=style.BLUE, lw=0.6))
ax1.set_xlim(-1.9, 2.0)
ax1.set_ylim(-1.9, 1.9)
ax1.set_aspect("equal")
ax1.axis("off")
ax1.set_title(r"(a) $G_2$: long and short roots", fontsize=9, color=style.GREY)

cmap = plt.get_cmap("viridis")
for idx, r in enumerate(radii):
    m = np.abs(rad - r) < 1e-5
    t = np.linspace(0, 2 * np.pi, 200)
    ax2.plot(r * np.cos(t), r * np.sin(t), color=style.GREY, lw=0.4, alpha=0.5, zorder=1)
    ax2.scatter(z[m].real, z[m].imag, s=6, color=cmap(0.1 + 0.8 * idx / 7), lw=0, zorder=3)
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_title(r"(b) $E_8$: eight 30-gons, radii $\propto m_a$", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch5_f10_root_systems", style.OUT_PART2)
raise SystemExit(0 if style.all_ok() else 1)
