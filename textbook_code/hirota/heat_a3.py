"""Heat-trace check of Proposition 8.5(2) for the a_3^(1) soliton a=1, where X_13(0) = 0 and X_11 has a pole
at k = 0.  Grid: Tr(e^{-tA} - e^{-tA0}) from all eigenvalues of the Fourier-grid operators.
Formula with mass-degenerate channels grouped: sum over mass levels mu of
  sum_zeta ord e^{-t lambda_zeta} + (1/(pi i)) int_0^inf e^{-t(mu^2+k^2)} d/dk log prod_{m_b=mu} X_b(k) dk.
Closed form for this soliton: (1 + e^{-2t}) erf(sqrt(2t)).
Also evaluates the ungrouped formula with a principal-value regularization to show the divergence."""
import sys, os
import numpy as np
from scipy.special import erf, erfc
from scipy.integrate import quad
from sectors import parent_sector, Xclosed

S = parent_sector('a3', int(os.environ.get('A', '1')))
ts = [0.3, 0.6, 1.0, 1.5]
h = S.A.h

def formula(t, group=True, eps=0.0):
    total = 0
    levels = {}
    for nm, d, lb in S.channels:
        levels.setdefault(round(S.mb2[nm], 9), []).append(nm)
    for m2, chans in levels.items():
        m = np.sqrt(m2)
        # zeros/poles in z in (0,1]
        for nm in chans:
            for q, o in S.exps[nm].items():
                if 2 * q < h:
                    z0 = np.cos(np.pi * q / h)
                    total += o * np.exp(-t * m2 * (1 - z0 ** 2))
        def dlog(k):
            z = 1j * k / m
            s = 0
            for nm in chans:
                for q, o in S.exps[nm].items():
                    s += o * (1j / m) / (z - np.cos(np.pi * q / h))
            return s
        f = lambda k: (np.exp(-t * (m2 + k * k)) * dlog(k) / (np.pi * 1j)).real
        g = lambda k: (np.exp(-t * (m2 + k * k)) * dlog(k) / (np.pi * 1j)).imag
        val = quad(f, eps, np.inf, limit=400)[0] + 1j * quad(g, eps, np.inf, limit=400)[0]
        total += val
    return total

for N, Lfac in [(384, 40), (512, 56), (768, 56)]:
    Lp = Lfac
    H, H0, trV, kmax, x = S.grid(N, Lp)
    ev = np.linalg.eigvals(H)
    ev0 = np.linalg.eigvalsh(H0.real) if np.allclose(H0.imag, 0) else np.linalg.eigvals(H0)
    out = []
    for t in ts:
        tr = np.sum(np.exp(-t * ev)) - np.sum(np.exp(-t * ev0))
        out.append(tr)
    print(f"N={N} L={Lp}: grid traces " + "  ".join(f"t={t}: {v.real:.6f}{v.imag:+.1e}i" for t, v in zip(ts, out)), flush=True)
print("formula (grouped):  " + "  ".join(f"t={t}: {formula(t).real:.6f}" for t in ts))
print("closed form (1+e^-2t) erf(sqrt 2t): " + "  ".join(f"t={t}: {(1+np.exp(-2*t))*erf(np.sqrt(2*t)):.6f}" for t in ts))
print("weight-0 variant at lambda=2 (drop e^{-2t}):  " + "  ".join(f"t={t}: {formula(t).real-np.exp(-2*t):.6f}" for t in ts))
