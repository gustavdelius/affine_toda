"""Check that the Hirota soliton solves the static field equation phi'' = -i sum_j n_j alpha_j e^{i alpha_j.phi}
and that exact fluctuation solutions solve -u'' + M u = lambda u (finite differences, 4th order)."""
import sys
import numpy as np, mpmath as mp
from lie import Algebra
from fluctop import single_soliton, Channels, exact_solution

def d2(f, hgrid):
    return (-f[..., 4:] + 16 * f[..., 3:-1] - 30 * f[..., 2:-2] + 16 * f[..., 1:-3] - f[..., :-4]) / (12 * hgrid ** 2)

for nm in sys.argv[1:]:
    A = Algebra(nm[0], int(nm[1:]))
    ch = Channels(A)
    for a in sorted(A.mass):
        bg, t, mu = single_soliton(A, a)
        x = np.linspace(-15, 15, 6001) / float(mu)
        hg = x[1] - x[0]
        ph = bg.phi(x)
        lhs = d2(ph, hg)
        E = bg.expo(x)
        rhs = -1j * (A.alpha.T @ (A.n[:, None] * E))[:, 2:-2]
        # consistency of expo with phi: e^{i alpha_j.phi}
        E2 = np.exp(1j * (A.alpha @ ph))
        res1 = np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs))
        res0 = np.max(np.abs(E2 - E))
        worst = 0
        for b in sorted(A.mass):
            for z in [0.3 + 0.2j, -0.7 + 1.1j]:
                u, up, X, _ = exact_solution(A, t, mu, ch, b, z, bg.xi, x)
                lam = float(ch.lam[b]) * (1 - z * z)
                Mx = bg.M(x)
                Mu = np.einsum('xij,jx->ix', Mx, u)
                r = -d2(u, hg) + (Mu - lam * u)[:, 2:-2]
                worst = max(worst, np.max(np.abs(r)) / np.max(np.abs(Mu)))
        print(f"{nm} a={a}: xi={bg.xi:.3f} min|tau|={bg.min_tau:.3f}  e^(i a.phi) consistency {res0:.1e}  field eq {res1:.1e}  fluct eq (all b) {worst:.1e}")
