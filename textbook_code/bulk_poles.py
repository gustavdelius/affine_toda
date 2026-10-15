"""Bulk vacuum-energy density of simply-laced ATFT in terms of the bare coupling, continued to imaginary coupling
(sec-collapse-thresholds, chapter 46).

Fateev's mass-mu relation and the bulk energy f = mbar^2 sin(pi/h)/(8 sin(pi B/h) sin(pi(1-B)/h))
(as quoted in Ahn-Fateev-Kim-Rim-Yang, hep-th/0002213, eqs. (mmu),(fg), with H=h), with b^2 -> -kappa,
kappa = beta^2/4pi in the book's normalization.  At fixed bare coupling z,
   f  propto  (pi z gamma(1-kappa))^{1/(1-kappa)} * G(kappa),
   G = Gamma(1-1/(h(1-k))) Gamma(1+k/(h(1-k))) / ( B Gamma(1/(h(1-k))) Gamma(1-k/(h(1-k))) ),  B=-k/(1-k).

Checks, for h = 2, 3, 4, 5, 6, 12 (one PASS/FAIL line each):
  - on kappa in (0.02, kappa_max), with kappa_max between the third and fourth predicted pole, the poles of G
    (zeros of 1/G, located by bisecting every sign change of 1/G on a grid of 2e5 points) are exactly the neutral
    thresholds kappa = 1 - 1/(k h), k = 1, 2, 3;
  - each of these is a simple pole: (kappa - kappa_k) G(kappa) tends to the same non-zero limit for
    kappa - kappa_k = +-1e-7, +-1e-9;
  - G is finite at kappa = 1/2 (beta^2 = 2 pi) for h >= 3;
  - G < 0 (finite and negative) on (0, 1 - 1/h), below the first threshold; since gamma(1-kappa) > 0 there, the
    bulk energy at fixed z has the sign of G (used in the remark after lem-stability, chapter 47);
  - for h = 2 the poles sit at odd xi = kappa/(1-kappa) (the sine-Gordon poles).
Exit status non-zero if a check fails.  Runtime: a few seconds."""
import sys

import numpy as np
from scipy.special import gamma as Gam

np.seterr(all="ignore")


def Gfun(k, h):
    a = 1 / (h * (1 - k))
    c = k / (h * (1 - k))
    B = -k / (1 - k)
    return Gam(1 - a) * Gam(1 + c) / (B * Gam(a) * Gam(1 - c))


ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


for h in (2, 3, 4, 5, 6, 12):
    poles = np.array([1 - 1 / (kk * h) for kk in (1, 2, 3, 4)])
    kmax = 0.5 * (poles[2] + poles[3])
    kap = np.linspace(0.02, kmax, 200001) + np.sqrt(2) * 1e-7   # offset: no grid point on a pole
    inv = 1 / Gfun(kap, h)
    found = []
    for i in np.flatnonzero(np.sign(inv[:-1]) != np.sign(inv[1:])):
        lo, hi = kap[i], kap[i + 1]                 # bisect the sign change of 1/G
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if np.sign(1 / Gfun(mid, h)) == np.sign(1 / Gfun(lo, h)):
                lo = mid
            else:
                hi = mid
        if min(abs(Gfun(lo, h)), abs(Gfun(hi, h))) > 1e6:   # |G| large on both sides: a pole (else a zero of G)
            found.append(0.5 * (lo + hi))
    found = np.array(found)
    match = len(found) == 3 and np.allclose(found, poles[:3], atol=1e-10)
    check(f"h={h}: G has poles at kappa = {np.round(found, 6)} only; predicted 1-1/(kh) = {np.round(poles[:3], 6)}",
          match)
    res = []
    for p in poles[:3]:
        r = [eps * Gfun(p + eps, h) for eps in (1e-7, -1e-7, 1e-9, -1e-9)]
        res.append(abs(r[0] - r[2]) < 1e-3 * abs(r[2]) and abs(r[1] - r[3]) < 1e-3 * abs(r[3])
                   and abs(r[0] - r[1]) < 1e-3 * abs(r[0]) and abs(r[0]) > 1e-12)
    check(f"h={h}: the three poles are simple (residue stable under 1/eps scaling)", all(res))
    kb = np.linspace(1e-4, poles[0] - 1e-4, 20001)
    Gb = Gfun(kb, h)
    check(f"h={h}: G finite and negative below the first threshold (max G = {Gb.max():.3g})",
          np.all(np.isfinite(Gb)) and np.all(Gb < 0))
    if h >= 3:
        check(f"h={h}: G finite at beta^2 = 2 pi (kappa = 1/2), |G| = {abs(Gfun(0.5, h)):.4g}",
              np.isfinite(Gfun(0.5, h)) and abs(Gfun(0.5, h)) < 1e3)
    else:
        xi = poles[:3] / (1 - poles[:3])
        check(f"h=2: poles at xi = {np.round(xi, 10)}, the odd integers 1, 3, 5", np.allclose(xi, [1, 3, 5]))
sys.exit(0 if ok else 1)
