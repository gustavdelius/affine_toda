"""Jordan-chain / Hypothesis (R) analysis at every zero of every X_b in the bound-state half-plane.

At +infinity the Hirota solution is psi_b = e^{kappa x} sum_{p>=0} sigma_p(z) E^{-p}, E = e^{mu x + xi}, with
sigma_0 = X_b(z) u_b.  For a zero z0 of order m, the chain chi_j = d^j psi/dz^j / j! (j < m) lies in L^2(R_+) iff
d^i sigma_p(z0) = 0 for all i <= j and all p >= 1 with kappa0 - p mu >= 0 (at -infinity it always decays).
Exponent kappa0 - p mu = 0 with nonzero coefficient means a bounded threshold tail (as in Proposition 8.7).
We report, for each zero, the first j for which chi_j fails to be in L^2 and the exponent responsible.
Usage: python3 jordan_check.py 'b3(1)<d4:3' 'g2(1)<d4:1' 'd4:1' ..."""
import sys
import numpy as np, mpmath as mp
from sectors import parent_sector, fold_sector
from hirota import MpCtx, fluct_deg

DPS = 40


def tail_coeffs(S, nm, z, P):
    """sigma_p (delta-phi vectors, p = 0..P) of the +inf expansion at z (mp complex)"""
    d, lb = [(dd, l) for n2, dd, l in S.channels if n2 == nm][0]
    s, dr, res = fluct_deg(S.A, S.t, S.mu, d, lb, z, ctx=MpCtx(DPS))
    out = []
    r1 = S.A.r + 1
    sig = []
    for j in range(r1):
        T = S.t[j][::-1]          # coefficients of E^{D}, E^{D-1}, ... -> series in E^{-1}
        Sx = s[j][::-1]
        Sx = Sx + [mp.mpc(0)] * (P + 1)
        T = T + [mp.mpc(0)] * (P + 1)
        q = []
        for p in range(P + 1):
            acc = Sx[p] - sum(q[i] * T[p - i] for i in range(p))
            q.append(acc / T[0])
        sig.append(q)
    al = mp.matrix(S.A.alpha.tolist())
    for p in range(P + 1):
        v = [sum(al[j, c] * sig[j][p] for j in range(r1)) for c in range(S.A.r)]
        out.append(v)
    return out


def derivs(S, nm, z0, P, m, h=mp.mpf('1e-7')):
    """d^i sigma_p / dz^i at z0, i = 0..m-1, via Cauchy integral on a circle"""
    K = 2 * m + 6
    pts = [z0 + h * mp.expj(2 * mp.pi * k / K) for k in range(K)]
    vals = [tail_coeffs(S, nm, z, P) for z in pts]
    res = []
    for i in range(m):
        Di = []
        for p in range(P + 1):
            v = []
            for c in range(S.A.r):
                acc = sum(vals[k][p][c] * mp.expj(-2 * mp.pi * k * i / K) for k in range(K)) / K
                v.append(acc * mp.factorial(i) / h ** i)
            Di.append(v)
        res.append(Di)
    return res   # res[i][p][c]


def analyse(spec, only_multiple=False):
    nm, a = spec.rsplit(':', 1)
    a = int(a)
    S = fold_sector(nm, a, dps=DPS) if '<' in nm else parent_sector(nm, a, dps=DPS)
    mp.mp.dps = DPS
    mu = S.mu
    h = S.A.h
    print(f"### {S.label}: mu={float(mu):.5f}; thresholds m^2 = { {k: round(v,5) for k,v in S.mb2.items()} }")
    rows = []
    for nmc, d, lb in S.channels:
        for q, m in S.exps[nmc].items():
            if m <= 0 or 2 * q >= h:
                continue
            if only_multiple and m < 2:
                continue
            z0 = mp.cos(mp.pi * q / h)
            mb = mp.sqrt(lb)
            kap0 = z0 * mb
            lam0 = float(lb * (1 - z0 ** 2))
            P = int(mp.floor(kap0 / mu + mp.mpf('1e-9')))
            D = derivs(S, nmc, mp.mpc(z0), P + 1, m)
            scale = max(abs(x) for p in range(P + 2) for x in D[0][p]) + max(abs(x) for x in D[min(1, m - 1)][1]) + 1e-30
            # first failing order
            fail = None
            for j in range(m):
                for i in range(j + 1):
                    for p in range(1, P + 1):
                        mag = max(abs(x) for x in D[i][p])
                        if mag > 1e-12 * (1 + float(scale)):
                            ex = float(kap0 - p * mu)
                            fail = (j, p, ex, float(mag))
                            break
                    if fail:
                        break
                if fail:
                    break
            sig0 = [float(max(abs(x) for x in D[i][0])) for i in range(m)]
            iso = lam0 < min(S.mb2.values()) - 1e-9
            thr = [k for k, v in S.mb2.items() if abs(v - lam0) < 1e-9]
            desc = f"  X_{nmc} zero q={q} order {m}: lambda0={lam0:.6f} ({'isolated' if iso else 'embedded'}{', AT THRESHOLD of '+str(thr) if thr else ''}); kappa0={float(kap0):.5f}, kappa0/mu={float(kap0/mu):.4f}"
            if fail is None:
                desc += f" -> chain of length {m} in L^2 (all tails decay)"
            else:
                j, p, ex, mag = fail
                kind = 'bounded (threshold tail)' if abs(ex) < 1e-9 else 'growing'
                desc += f" -> chi_{j} NOT in L^2: tail exponent kappa0-{p}mu={ex:+.5f} ({kind}), coeff {mag:.2e}"
                # direction of the offending tail coefficient (highest derivative i <= j that is nonzero)
                vecs = [np.array([complex(x) for x in D[i][p]]) for i in range(j + 1)]
                v = [vv for vv in vecs if np.max(np.abs(vv)) > 1e-12 * (1 + float(scale))][0]
                comps = {}
                for c2, d2, l2 in S.channels:
                    kc = np.sqrt(max(S.mb2[c2] - lam0, 0))
                    if abs(kc - abs(ex)) < 1e-7:
                        cb = [c3 for c3, _, _ in S.channels if abs(S.mb2[c3] - S.mb2[c2]) < 1e-9 and abs(S.raw[c3] @ S.raw[c2]) > 1e-8][0]
                        comps[c2] = complex(S.raw[cb] @ v / (S.raw[cb] @ S.raw[c2]))
                resid = v - sum(comps[c2] * S.raw[c2] for c2 in comps) if comps else v
                desc += f"; tail along channels {{{', '.join(f'{c2}: {abs(comps[c2]):.2e}' for c2 in comps)}}} (residual {np.max(np.abs(resid)):.1e})"
            print(desc, flush=True)
            rows.append((nmc, q, m, lam0, fail))
    return S, rows


if __name__ == '__main__':
    om = '--multiple' in sys.argv
    for spec in [s for s in sys.argv[1:] if not s.startswith('--')]:
        analyse(spec, only_multiple=om)
