"""Figure F15 (chapter 10, sec-multi-solitons): the interaction coefficient and the shift in a two-soliton collision.

(a) World lines of two a_5^(1) solitons (species 1 and 2, rapidities 0.6 and -0.6, m = 1) before and after the
    collision, from the asymptotics of the exact two-soliton tau-function: soliton 1 is displaced by
    -ln A/(m_1 cosh th_1) > 0 in x and soliton 2 by +ln A/(m_2 cosh th_2) < 0, so both are shifted forward, along
    their own direction of motion (a time advance for each).
(b) A_ab(theta) for h = 6: 0 < A_ab < 1, A -> 1 as theta -> infinity.

Checks: the two-soliton tau_j = 1 + w^{j a} E_1 + w^{j b} E_2 + w^{j(a+b)} A E_1 E_2 solves the bilinear equation
(1/2)(D_x^2 - D_t^2) tau_j.tau_j = tau_j^2 - tau_{j-1} tau_{j+1} (m = 1) at random complex points; near soliton 1 the
tau-function tends to 1 + w^{ja} E_1 for t -> -infinity and to w^{jb} E_2 (1 + w^{ja} A E_1) for t -> +infinity,
and symmetrically for soliton 2; 0 < A_ab < 1 on the real line for every a, b of a_5^(1); A_ab has a pole at
theta = i pi (a+b)/h when a + b < h.
Writes part3-classical/figures/ch10_f15_two_soliton_shift.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
rng = np.random.default_rng(3)
h = 6
om = np.exp(2j * np.pi / h)
mass = lambda a: 2 * np.sin(np.pi * a / h)
A = lambda a, b, th: (np.cosh(th) - np.cos(np.pi * (a - b) / h)) / (np.cosh(th) - np.cos(np.pi * (a + b) / h))


def tau(j, a, b, th1, th2, xi1, xi2, x, t):
    """tau_j and its x, t derivatives up to second order (sum of exponentials)."""
    k = [(mass(a) * np.cosh(th1), -mass(a) * np.sinh(th1)), (mass(b) * np.cosh(th2), -mass(b) * np.sinh(th2))]
    terms = [(1.0, (0, 0), 0), (om ** (j * a), k[0], xi1), (om ** (j * b), k[1], xi2),
             (om ** (j * (a + b)) * A(a, b, th1 - th2), (k[0][0] + k[1][0], k[0][1] + k[1][1]), xi1 + xi2)]
    val = dx = dt = dxx = dtt = 0
    for c, (kx, kt), xi in terms:
        e = c * np.exp(kx * x + kt * t + xi)
        val += e; dx += kx * e; dt += kt * e; dxx += kx**2 * e; dtt += kt**2 * e
    return val, dx, dt, dxx, dtt


# ---- checks ----
worst = 0.0
for _ in range(40):
    a, b = rng.integers(1, h, 2)
    th1, th2 = rng.normal(size=2)
    xi1, xi2 = rng.normal(size=2) + 1j * rng.normal(size=2)
    x, t = rng.normal(size=2)
    T = [tau(j, a, b, th1, th2, xi1, xi2, x, t) for j in range(h)]
    for j in range(h):
        v, dx, dt, dxx, dtt = T[j]
        lhs = v * (dxx - dtt) - (dx**2 - dt**2)
        rhs = v**2 - T[(j - 1) % h][0] * T[(j + 1) % h][0]
        worst = max(worst, abs(lhs - rhs) / (1 + abs(rhs)))
report("two-soliton tau solves (1/2)(D_x^2 - D_t^2) tau.tau = tau_j^2 - tau_{j-1} tau_{j+1} (residual %.1e)" % worst,
       worst < 1e-9)

a, b, th1, th2, xi1, xi2 = 1, 2, 0.6, -0.6, 0.3 + 0.4j, -0.2 + 0.1j
m1, m2, c1, c2 = mass(a), mass(b), np.cosh(th1), np.cosh(th2)
v1, v2 = np.tanh(th1), np.tanh(th2)
ok1 = ok2 = True
for sgn in (-1, 1):
    T_ = 60.0 * sgn
    s = 0.7
    # soliton 1: x = v1 t + s ; soliton 2: x = v2 t + s
    for j in range(h):
        E1 = lambda x, t: np.exp(m1 * (x * c1 - t * np.sinh(th1)) + xi1)
        E2 = lambda x, t: np.exp(m2 * (x * c2 - t * np.sinh(th2)) + xi2)
        x, t = v1 * T_ + s, T_
        exact = tau(j, a, b, th1, th2, xi1, xi2, x, t)[0]
        asym = (1 + om ** (j * a) * E1(x, t)) if sgn < 0 else om ** (j * b) * E2(x, t) * (1 + om ** (j * a) * A(a, b, th1 - th2) * E1(x, t))
        ok1 &= abs(exact / asym - 1) < 1e-8
        x = v2 * T_ + s
        exact = tau(j, a, b, th1, th2, xi1, xi2, x, t)[0]
        asym = (1 + om ** (j * b) * E2(x, t)) if sgn > 0 else om ** (j * a) * E1(x, t) * (1 + om ** (j * b) * A(a, b, th1 - th2) * E2(x, t))
        ok2 &= abs(exact / asym - 1) < 1e-8
report("soliton 1: 1 + w^{ja}E1 before, w^{jb}E2 (1 + w^{ja} A E1) after the collision", ok1)
report("soliton 2: w^{ja}E1 (1 + w^{jb} A E2) before, 1 + w^{jb}E2 after the collision", ok2)
th = np.linspace(-6, 6, 2001)
pos = all(0 < A(a_, b_, th[th != 0]).min() and A(a_, b_, th).max() < 1 + 1e-12
          for a_ in range(1, h) for b_ in range(1, h))
report("0 < A_ab(theta) <= 1 for real theta, all species of a_5^(1)", pos)
report("A_ab(theta) -> 1 as theta -> infinity", all(abs(A(a_, b_, 15.0) - 1) < 1e-5 for a_ in range(1, h) for b_ in range(1, h)))
pole = all(abs(np.cosh(1j * np.pi * (a_ + b_) / h) - np.cos(np.pi * (a_ + b_) / h)) < 1e-12
           for a_ in range(1, h) for b_ in range(1, h) if a_ + b_ < h)
report("denominator of A_ab vanishes at theta = i pi (a+b)/h", pole)
shift1 = -np.log(A(a, b, th1 - th2)) / (m1 * c1)
shift2 = np.log(A(a, b, th1 - th2)) / (m2 * c2)
report("shifts: soliton 1 moves right and is displaced to the right, soliton 2 moves left and is displaced to the left", shift1 > 0 > shift2 and v1 > 0 > v2)

# ---- figure ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw=dict(width_ratios=[1, 1.15]))
tt = np.array([-3.0, 3.0])
# symmetric offsets: in-lines at -D/2 (soliton 1) and +D/2 (soliton 2) at t = 0, out-lines at +D/2 and -D/2 (shifts as above)
d1, d2 = shift1, -shift2
for v, off_in, off_out, col, lab in ((v1, -d1 / 2, d1 / 2, style.RED, r"soliton 1 (species 1, $\theta_1=0.6$)"),
                                     (v2, d2 / 2, -d2 / 2, style.BLUE, r"soliton 2 (species 2, $\theta_2=-0.6$)")):
    tin, tout = np.array([-3.0, 0.0]), np.array([0.0, 3.0])
    ax1.plot(off_in + v * tin, tin, color=col, lw=1.8, label=lab)
    ax1.plot(off_out + v * tout, tout, color=col, lw=1.8)
    ax1.plot(off_in + v * np.array([0, 3.0]), [0, 3.0], color=col, lw=0.9, ls=":")      # in-line continued
    ax1.plot(off_out + v * np.array([-3.0, 0]), [-3.0, 0], color=col, lw=0.9, ls=":")  # out-line continued
ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$t$")
ax1.set_xlim(-3.2, 3.2)
ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=7.5)
ax1.set_title(r"world lines: each soliton is shifted forward", fontsize=8.5)
for (a_, b_), col in (((1, 1), style.BLUE), ((1, 2), style.RED), ((1, 3), style.GREEN), ((2, 2), style.ORANGE),
                      ((2, 3), style.PURPLE), ((3, 3), style.TEAL)):
    ax2.plot(th[th >= 0], A(a_, b_, th[th >= 0]), color=col, label=r"$A_{%d%d}$" % (a_, b_))
ax2.axhline(1, color=style.GREY, lw=0.6, ls="--")
ax2.set_xlabel(r"$\theta$")
ax2.set_ylim(0, 1.05)
ax2.set_xlim(0, 4)
ax2.legend(ncol=3, fontsize=7.5, loc="lower right")
ax2.set_title(r"$A_{ab}(\theta)$, $a_5^{(1)}$ ($h=6$)", fontsize=8.5)
fig.tight_layout()
style.save(fig, "ch10_f15_two_soliton_shift", style.OUT_PART3)
raise SystemExit(0 if style.all_ok() else 1)
