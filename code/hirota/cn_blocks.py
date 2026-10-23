"""The c_n^(1) threshold blocks: the Numerics after prp-cn-threshold-blocks (Chapter 18).

Run from this directory:  python3 cn_blocks.py     (about 2 min)

For every soliton a < n of c_2^(1), c_3^(1) and c_4^(1), realized as the sigma-invariant sector of a_{2n-1}^(1):
 (1) the soliton solves the static field equation and the channel-n solution at lambda_0 solves the fluctuation
     equation (4th-order finite differences);
 (2) X_n has a double zero at z_0 = sin(pi a/2n), that is at lambda_0 = m_{n-a}^2, the threshold of channel n-a;
 (3) psi_n(lambda_0) decays at both ends and is self-orthogonal: |int psi^T psi dx| << int |psi|^2 dx;
 (4) d psi_n / d lambda at lambda_0 tends to a constant vector as x -> +infinity, along channel n-a. The
     approach is like x e^{-m_a x}, slowest for a = 1: the tail is constant to 1e-6 from x = 30, not from x = 10.
"""
import numpy as np, mpmath as mp
from sectors import fold_sector

def d2(f, h):
    return (-f[..., 4:] + 16*f[..., 3:-1] - 30*f[..., 2:-2] + 16*f[..., 1:-3] - f[..., :-4]) / (12*h**2)

import sys
NS = [int(v) for v in sys.argv[1:]] or [2, 3, 4]
for n in NS:
    for a in range(1, n):
        S = fold_sector(f'c{n}(1)<a{2*n-1}', a, dps=80)
        A, bg = S.A, S.bg
        mp.mp.dps = 80
        z0m = mp.sin(mp.pi*a/(2*n)); z0 = float(z0m)
        mn2 = S.mb2[n]
        lam0 = mn2*(1 - z0**2)
        tag = f"c{n}(1) soliton {a}"

        # (2) double zero of X_n at z0, at the threshold of channel n-a
        X = lambda z: S.coeffs(n, z)[1]
        r1, r2 = abs(X(z0m + mp.mpf('1e-3')))/1e-6, abs(X(z0m + mp.mpf('1e-4')))/1e-8
        ok2 = abs(X(z0m)) < 1e-40 and abs(r1/r2 - 1) < 1e-2 and abs(lam0 - S.mb2[n - a]) < 1e-12
        print(f"{'PASS' if ok2 else 'FAIL'} {tag}: X_{n}(z0) = {abs(X(z0m)):.0e}, |X_{n}(z0+d)|/d^2 = {r1:.5f}, {r2:.5f} "
              f"(d = 1e-3, 1e-4); lambda_0 = {lam0:.6f} = m_{n-a}^2 = {S.mb2[n - a]:.6f}")

        # (1) residuals of the field and fluctuation equations
        s0 = S.coeffs(n, z0m)[0]
        res = []
        for npts in (3001, 6001):                               # halving h must cut a 4th-order residual by 16
            x = np.linspace(-15, 15, npts)/float(S.mu); h = x[1] - x[0]
            ph = bg.phi(x)
            rhs = -1j*(A.alpha.T @ (A.n[:, None]*bg.expo(x)))[:, 2:-2]
            u, _ = S.evaluate(s0, z0, n, x)
            Mu = np.einsum('xij,jx->ix', bg.M(x), u)
            res.append((np.max(np.abs(d2(ph, h) - rhs))/np.max(np.abs(rhs)),
                        np.max(np.abs(-d2(u, h) + Mu[:, 2:-2] - lam0*u[:, 2:-2]))/np.max(np.abs(lam0*u))))
        ratios = [res[0][i]/res[1][i] for i in (0, 1)]
        ok1 = all(12 < q < 20 or r < 1e-8 for q, r in zip(ratios, res[1]))      # 4th order, or at round-off
        print(f"{'PASS' if ok1 else 'FAIL'} {tag}: residuals of the field and fluctuation equations "
              f"{res[1][0]:.1e}, {res[1][1]:.1e}; ratios on halving h {ratios[0]:.1f}, {ratios[1]:.1f}")

        # (3) decay and self-orthogonality of psi_n(lambda_0)
        xl = np.linspace(-60, 60, 120001)/float(S.mu); hl = xl[1] - xl[0]
        u, _ = S.evaluate(s0, z0, n, xl)
        ends = max(np.max(np.abs(u[:, :10])), np.max(np.abs(u[:, -10:])))
        bil = abs(np.sum(np.einsum('ix,ix->x', u, u))*hl)
        norm = np.sum(np.abs(u)**2)*hl
        ok3 = ends < 1e-10*np.max(np.abs(u)) and bil < 1e-12*norm
        print(f"{'PASS' if ok3 else 'FAIL'} {tag}: |psi| at the ends {ends:.0e}; |<psi,psi>| = {bil:.1e}, "
              f"||psi||^2 = {norm:.2f}")

        # (4) tail of d psi / d lambda along channel n-a
        # d psi/dz = x m_n psi + e^{kappa x} (dS/dz)/T, with dS/dz by a central difference at 80 digits
        dz = mp.mpf('1e-30')
        sp, sm = S.coeffs(n, z0m + dz)[0], S.coeffs(n, z0m - dz)[0]
        ds = [[(p - q)/(2*dz) for p, q in zip(sj, tj)] for sj, tj in zip(sp, sm)]
        xt = np.array([10.0, 20.0, 30.0, 40.0])
        dpsi = xt*np.sqrt(mn2)*S.evaluate(s0, z0, n, xt)[0] + S.evaluate(ds, z0, n, xt)[0]
        dpsi *= -1/(2*mn2*z0)                                   # d/d lambda = (dz/d lambda) d/dz
        c = dpsi[:, -1]
        var = [np.linalg.norm(dpsi[:, i] - c)/np.linalg.norm(c) for i in range(3)]   # x = 10, 20, 30 against x = 40
        e = S.raw[n - a]/np.linalg.norm(S.raw[n - a])
        perp = np.linalg.norm(c - np.vdot(e, c)*e)/np.linalg.norm(c)
        ok4 = var[2] < 1e-6 and var[0] > var[1] > var[2] and perp < 1e-10 and np.linalg.norm(c) > 1e-3
        print(f"{'PASS' if ok4 else 'FAIL'} {tag}: tail of d psi/d lambda differs from its value at x = 40 by "
              f"{var[0]:.1e}, {var[1]:.1e}, {var[2]:.1e} at x = 10, 20, 30; component off channel {n-a} {perp:.1e}; |tail| = {np.linalg.norm(c):.3f}",
              flush=True)
