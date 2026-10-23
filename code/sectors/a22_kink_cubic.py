"""The two-kink eigenvalue cubic of a_2^(2) (Takacs-Watts 1999, eq. 5.5), for the paragraph "Fundamental kinks of
a_2^(2)" after prp-scalar-breathers (chapter 48).

Conventions of Takacs-Watts: q = i e^{i pi^2/(3 xi)}, so q^2 = -e^{i psi} with psi = 2 pi^2/(3 xi) reduced to
(-pi, pi], and x = e^{2 pi theta/xi}.  The cubic is the one of a22cubic.py (checked against the eigenvalues of the
particle-labelled R-matrix in a22tw.py).  One PASS/FAIL line per check; exit status non-zero if one fails.
 1. at xi = 0.9 pi and theta = 0.55 two roots have moduli 2.118 and 0.472 = 1/2.118;
 2. for psi/pi on a grid in (-1, 1), avoiding |cos psi| < 0.05, and x > 0 on a log grid up to e^40, the roots are
    off the unit circle (max ||l|-1| > 1e-3) exactly when cos psi < 0, and on it (< 1e-6) otherwise;
 3. at large x two roots approach q^4 = e^{2 i psi} as lambda = q^4 (1 + delta), with
    x delta^2 -> -16 sin^2 psi cos psi (to 1% at x = 1e8, for six values of psi).
Runtime: a few seconds."""
import sys

import numpy as np


def cubic(x, q):
    c3 = (x - q**4) * (x + q**6)
    c2 = q**6 * (2 + q**2) + x * (q**2 - 1) * (1 - 3 * q**4 + q**8) - q**2 * (1 + 2 * q**2) * x**2
    c1 = q**6 * (2 + q**2) * x**2 + x * (q**2 - 1) * (1 - 3 * q**4 + q**8) - q**2 * (1 + 2 * q**2)
    c0 = (1 - q**4 * x) * (1 + q**6 * x)
    return np.roots([c3, c2, c1, c0])


ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


# 1. the example
xi = 0.9 * np.pi
q = 1j * np.exp(1j * np.pi**2 / (3 * xi))
mods = np.sort(np.abs(cubic(np.exp(2 * np.pi * 0.55 / xi), q)))
check(f"xi = 0.9 pi, theta = 0.55: root moduli {np.round(mods, 4)}; 2.118 and 0.472 = 1/2.118",
      abs(mods[2] - 2.118) < 5e-4 and abs(mods[0] - 0.472) < 5e-4 and abs(mods[0] * mods[2] - 1) < 1e-10)

# 2. off the unit circle exactly when cos psi < 0
lxs = np.concatenate([np.linspace(1e-3, 6, 600), np.linspace(6, 40, 200)])
good, n = True, 0
for psi_pi in np.arange(-0.99, 1.0, 0.02):
    psi = np.pi * psi_pi
    if abs(np.cos(psi)) < 0.05:
        continue
    qq = 1j * np.exp(1j * psi / 2)                     # q^2 = -e^{i psi}
    dev = max(np.abs(np.abs(cubic(np.exp(l), qq)) - 1).max() for l in lxs)
    good &= (dev > 1e-3) if np.cos(psi) < 0 else (dev < 1e-6)
    n += 1
check(f"{n} values of psi: roots off the unit circle for some x > 0 exactly when cos psi < 0", good)

# 3. large-x asymptotics
good, worst = True, 0.0
for psi_pi in (-0.8, -0.6, -0.3, 0.2, 0.55, 0.9):
    psi = np.pi * psi_pi
    qq = 1j * np.exp(1j * psi / 2)
    x = 1e8
    r = cubic(x, qq)
    near = r[np.argsort(np.abs(r - qq**4))[:2]]
    d2 = x * (near / qq**4 - 1) ** 2
    pred = -16 * np.sin(psi) ** 2 * np.cos(psi)
    err = np.abs(d2 - pred).max() / abs(pred)
    worst = max(worst, err)
    good &= err < 1e-2
check(f"x delta^2 -> -16 sin^2 psi cos psi at x = 1e8 (worst relative error {worst:.1e})", good)
sys.exit(0 if ok else 1)
