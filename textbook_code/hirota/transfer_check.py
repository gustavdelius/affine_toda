"""Independent check of (H)/(R) by direct ODE integration: full transfer matrix between all channels at real
lambda above all thresholds.  For each incoming channel c (pure e^{ik_c x} u_c at x=-L), decompose the solution at
x=+L into e^{+-ik_d x} u_d.  Checks: reflection R_cd = 0, transmission T_cd = X_c delta_cd, with X_c from the
exact Hirota solution.  Usage: python3 transfer_check.py d4 e6 ..."""
import sys
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
from lie import Algebra
from fluctop import single_soliton, Channels, exact_solution
from hirota import MpCtx, fluct_deg


def transfer(bg, A, lam, L, chans):
    r = A.r
    U = {c: A.u[c] for c in chans}            # channel vectors, u_c^T u_cbar = 1
    k = {c: np.sqrt(lam - A.mass[c] ** 2 + 0j) for c in chans}
    def rhs(x, y):
        u = y[:r]; up = y[r:]
        Mx = bg.M(np.array([x]))[0]
        return np.concatenate([up, (Mx - lam * np.eye(r)) @ u])
    T = np.zeros((len(chans), len(chans)), complex)
    R = np.zeros_like(T)
    for i, c in enumerate(chans):
        y0 = np.concatenate([np.exp(1j * k[c] * (-L)) * U[c], 1j * k[c] * np.exp(1j * k[c] * (-L)) * U[c]])
        sol = solve_ivp(rhs, (-L, L), y0, method='DOP853', rtol=1e-11, atol=1e-13)
        u, up = sol.y[:r, -1], sol.y[r:, -1]
        for jx, d in enumerate(chans):
            db = A.conj[d]
            pu, pup = U[db] @ u, U[db] @ up
            Ad = (1j * k[d] * pu + pup) / (2j * k[d]) * np.exp(-1j * k[d] * L)
            Bd = (1j * k[d] * pu - pup) / (2j * k[d]) * np.exp(1j * k[d] * L)
            T[i, jx], R[i, jx] = Ad, Bd
    return T, R, k


if __name__ == '__main__':
    for nm in sys.argv[1:]:
        A = Algebra(nm[0], int(nm[1:]))
        ch = Channels(A)
        chans = sorted(A.mass)
        mmax2 = max(A.mass[c] ** 2 for c in chans)
        for a in chans:
            bg, t, mu = single_soliton(A, a)
            L = 28 / float(mu)
            for lam in [mmax2 + 0.7, mmax2 + 3.1]:
                T, R, k = transfer(bg, A, lam, L, chans)
                Xh = []
                for c in chans:
                    z = 1j * k[c] / A.mass[c]
                    s, dr, res = fluct_deg(A, t, mu, ch.d[c], ch.lam[c], mp.mpc(z.real, z.imag), ctx=MpCtx(30))
                    D = [len(tj) - 1 for tj in t]
                    w = np.array([complex(s[j][D[j]] / t[j][D[j]]) for j in range(A.r + 1)])
                    d0 = np.array([complex(x) for x in ch.d[c]])
                    Xh.append(np.vdot(d0, w) / np.vdot(d0, d0))
                Xh = np.array(Xh)
                off = np.max(np.abs(T - np.diag(np.diag(T))))
                print(f"{nm} a={a} lam={lam:.3f}: max|R|={np.max(np.abs(R)):.1e} max|T_offdiag|={off:.1e} "
                      f"max|T_cc - X_c|={np.max(np.abs(np.diag(T) - Xh)):.1e}  |X_c|={np.round(np.abs(Xh),6).tolist()}", flush=True)
