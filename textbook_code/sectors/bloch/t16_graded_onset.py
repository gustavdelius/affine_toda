"""Onset of the breaking for N kinks of a_2^(2) with graded exchange (the row "graded exchange" of the table in
sec-sector-status, chapter 48): even total charge Q at Bloch angle alpha = pi, odd Q at alpha = 0.

Conventions of ik_common.py (Takacs-Watts q, psi/pi in (-1, 1], couplings offset by 0.0013 pi in xi).  The
deviation of a sector is max ||s|-1| over the eigenvalues of the twisted transfer matrix Omega_alpha T_1 on it,
maximized over the rapidities sampled: the clustered line log x_{1j} = t (j = 2..N), t in [-2.5, 2.5] with fine
steps near 0, with a small transverse offset, plus random sets (seed 9).  "Unbroken" means < 1e-8, "broken" > 1e-4.
One PASS/FAIL line per check; exit status non-zero if one fails.
 1. N = 4: unbroken at |psi| = 0.497 pi and broken at 0.502 pi (the threshold 2 pi/N = 0.5 pi);
 2. N = 5: unbroken at |psi| = 0.4097 pi and broken at 0.4122 pi, so the onset lies near 0.41 pi, above 2 pi/5.
Output: out/out_t16_graded_onset.txt.  Runtime: about 10 min."""
import sys

import numpy as np

from ik_common import *

ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name, flush=True)


def max_dev(Nn, psi_target, nrand):
    Qs, m1 = charges(Nn), m_first(Nn)
    secs = []
    for Q in range(Nn + 1):
        idx = np.where(Qs == Q)[0]
        a = 1.0 if Q % 2 == 0 else 0.0
        secs.append((idx, np.exp(1j * np.pi * a * m1[idx])))
    xi_pi = 2 / (3 * psi_target) + 0.0013
    psi = psi_of_xi(xi_pi)
    Rp = setup(q_of_xi(xi_pi))
    rng = np.random.default_rng(9)
    ts = np.concatenate([-np.linspace(0.002, 2.5, 600), np.linspace(0.002, 2.5, 300)])
    off = np.resize(np.array([0.05, -0.05, 0.02, -0.02]), Nn - 1)
    cfgs = [t * np.ones(Nn - 1) for t in ts] + [t * np.ones(Nn - 1) + off for t in ts[::3]] \
        + [rng.uniform(-3, 3, Nn - 1) for _ in range(nrand)]
    top = 0.0
    for ls in cfgs:
        T = transferN(Rp, 3, ls)
        for idx, Om in secs:
            d = np.abs(np.abs(np.linalg.eigvals(Om[:, None] * T[np.ix_(idx, idx)])) - 1).max()
            top = max(top, d)
    return psi, top


for Nn, lo, hi in ((4, 0.497, 0.5025), (5, 0.41, 0.4125)):
    p_lo, d_lo = max_dev(Nn, lo, 300 if Nn == 4 else 100)
    p_hi, d_hi = max_dev(Nn, hi, 300 if Nn == 4 else 100)
    check(f"N={Nn}: unbroken at psi/pi = {p_lo:.4f} (max dev {d_lo:.1e}), broken at psi/pi = {p_hi:.4f} "
          f"(max dev {d_hi:.1e})", d_lo < 1e-8 and d_hi > 1e-4)
sys.exit(0 if ok else 1)
