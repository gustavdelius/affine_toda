"""Figure F22 (chapter 17, sec-picard-lefschetz): thimbles of the zero-dimensional i phi^3 integral.

S(z) = z^2/2 + i z^3/3. Critical points z = 0 (S = 0) and z = i (S = -1/6). Both have Im S = 0, so the
curves Im S = 0 (the imaginary axis and the curves x^2 = 3y(y-1)) contain all the flow lines.
Drawn: the real line, the thimble J_0 of the saddle at 0 (lower branch y = (1 - sqrt(1 + 4x^2/3))/2),
the dual thimble K_i of the saddle at i (upper branch, y = (1 + sqrt(1 + 4x^2/3))/2), and the segment of
the imaginary axis between the two saddles, a flow line that joins them (both saddles have the same Im S).
The real line meets K_0 once and K_i never, so nu_0 = 1 and nu_i = 0.

Checks: critical points and S(i) = -1/6; the closed forms of J_0 and K_i solve the flow equation
(Im S = 0, Re S increasing along J_0 and decreasing along K_i); a numerically integrated downward
flow from i stays in Im z >= 1; the integral over J_0 equals the real-line integral at lambda = 5, 10;
I = sqrt(2 pi/lambda) (1 - 5/(6 lambda) + O(lambda^-2)).
Writes part5-semiclassics/figures/ch18_f22_thimbles.svg.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp

import style
from style import plt, report

style.setup()
S = lambda z: z**2 / 2 + 1j * z**3 / 3
dS = lambda z: z + 1j * z**2
yJ = lambda x: (1 - np.sqrt(1 + 4 * x**2 / 3)) / 2     # J_0
yK = lambda x: (1 + np.sqrt(1 + 4 * x**2 / 3)) / 2     # K_i

# ---- checks ---------------------------------------------------------------
report("critical points 0 and i, S(i) = -1/6", abs(dS(0)) < 1e-14 and abs(dS(1j)) < 1e-14 and abs(S(1j) + 1 / 6) < 1e-14)
x = np.linspace(-6, 6, 4001)
zJ, zK = x + 1j * yJ(x), x + 1j * yK(x)
report("J_0 and K_i lie on Im S = 0", np.max(np.abs(S(zJ).imag)) < 1e-9 and np.max(np.abs(S(zK).imag)) < 1e-9)
report("Re S increases away from 0 along J_0", np.all(np.diff(S(zJ[x >= 0]).real) > 0))
report("Re S decreases away from i along K_i", np.all(np.diff(S(zK[x >= 0]).real) < 0))
# numerical downward flow dz/dt = -conj(S'(z)) from i + eps: tangent to K_i
sol = solve_ivp(lambda t, u: [-(np.conj(dS(u[0] + 1j * u[1]))).real, -(np.conj(dS(u[0] + 1j * u[1]))).imag],
                [0, 6], [0.01, 1.0], rtol=1e-10, atol=1e-12)
zf = sol.y[0] + 1j * sol.y[1]
report("downward flow from i never comes near the real line (Im z >= 1)", np.min(zf.imag) > 0.999)
for lam in (5.0, 10.0):
    real = quad(lambda t: np.exp(-lam * S(t)).real, -np.inf, np.inf, limit=400)[0]
    dyJ = lambda t: -(2 * t / 3) / np.sqrt(1 + 4 * t**2 / 3)          # y_J'(x)
    J = quad(lambda t: (np.exp(-lam * S(t + 1j * yJ(t))) * (1 + 1j * dyJ(t))).real, -np.inf, np.inf, limit=400)[0]
    report(f"integral over J_0 equals the real-line integral, lambda = {lam:g}", abs(J - real) < 1e-8)
lam = 60.0
real = quad(lambda t: np.exp(-lam * S(t)).real, -np.inf, np.inf, limit=400)[0]
report("I = sqrt(2 pi/lambda)(1 - 5/(6 lambda) + O(lambda^-2)), lambda = 60",
       abs(real / np.sqrt(2 * np.pi / lam) - (1 - 5 / (6 * lam))) < 6 / lam**2)   # remainder -> 385/(72 lam^2)

# ---- figure ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.2, 3.6))
xs = np.linspace(-2.6, 2.6, 400)
ax.axhline(0, color=style.GREY, lw=1.0)
ax.text(2.55, 0.08, r"real line $\mathbb{R}$", ha="right", va="bottom", fontsize=7.5, color=style.GREY)
ax.plot(xs, yJ(xs), color=style.BLUE, lw=2.2)
ax.annotate(r"thimble $\mathcal{J}_0$", (1.6, yJ(1.6)), (1.0, -1.45), fontsize=8, color=style.BLUE,
            arrowprops=dict(arrowstyle="-", color=style.BLUE, lw=0.7))
ax.plot(xs, yK(xs), color=style.RED, lw=1.6, ls="--")
ax.annotate(r"dual thimble $\mathcal{K}_i$", (1.7, yK(1.7)), (0.9, 3.15), fontsize=8, color=style.RED,
            arrowprops=dict(arrowstyle="-", color=style.RED, lw=0.7))
ax.plot([0, 0], [0, 1], color=style.GREY, lw=1.0, ls=":")
ax.text(0.07, 0.5, "flow line joining\nthe two saddles", fontsize=7, color=style.GREY, va="center")
ax.plot([0], [0], "o", color=style.BLUE, ms=6)
ax.plot([0], [1], "o", color=style.RED, ms=6)
ax.text(0, -0.42, r"$z=0$, $S=0$", ha="center", va="top", fontsize=8, color=style.BLUE)
ax.text(0, 1.38, r"$z=i$, $S=-\frac{1}{6}$", ha="center", va="bottom", fontsize=8, color=style.RED)
ax.set_xlim(-2.6, 2.6)
ax.set_ylim(-1.9, 3.4)
ax.set_aspect("equal")
ax.set_xlabel(r"$\mathrm{Re}\,z$")
ax.set_ylabel(r"$\mathrm{Im}\,z$")
fig.tight_layout()
style.save(fig, "ch18_f22_thimbles", style.OUT_PART5)
raise SystemExit(0 if style.all_ok() else 1)
