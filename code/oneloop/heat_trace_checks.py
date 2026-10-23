"""Checks of the proofs of prp-net-count parts 2-3, lem-born and prp-one-loop-mass (Chapters 18-19).

Run from this directory:  python3 heat_trace_checks.py            (analytic checks and 19 sectors, no grids)
                          python3 heat_trace_checks.py --grid 768 (adds Fourier-grid heat traces, slow)

(a) one root: -1/(4 sqrt pi) int t^{-3/2}[e^{-t(1-w^2)} erf(w sqrt t) - 2w sqrt(t/pi) e^{-t}] dt = psi(arcsin w)
(b) one factor (z-w)/(z+w): residue + continuum term of part 2 = e^{-t mu^2(1-w^2)} erf(mu w sqrt t)   (part 3)
(c) Hollowood's a_n^(1) masses from the closed form
(d) remarks after prp-one-loop-mass: phase-shift form + sum_b c_b/4pi, and the unpaired bookkeeping
(e) per soliton sector: pairing antisymmetry of each mass level; Born constant from the roots versus
    c_mu = int tr(P_mu (M - M0)) dx computed from the soliton; closed form versus the Mellin integral
    (eq-one-loop-mellin) of part 3 plus the directly computed counterterm
(f) with --grid N: Fourier-grid Tr(e^{-tA} - e^{-tA0}) versus part 3 at t = 0.3, 0.6, 1, 1.5
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'hirota'))
import numpy as np
import mpmath as mp
from scipy.special import erf
from scipy.integrate import quad

pi = np.pi


def psi(w):
    th = np.arcsin(np.clip(w, -1, 1))
    return (th*np.cos(th) - np.sin(th))/pi


# ---------------------------------------------------------------- (a), (b), (c), (d)
def check_one_root():
    mp.mp.dps = 30
    def I_mellin(w):
        w = mp.mpf(w)
        f = lambda t: t**-1.5*(mp.exp(-t*(1-w*w))*mp.erf(w*mp.sqrt(t)) - 2*w*mp.sqrt(t/mp.pi)*mp.exp(-t))
        return -mp.quad(f, [0, 1, 10, 100, mp.inf])/(4*mp.sqrt(mp.pi))
    worst = 0
    for w in ['0.001', '0.1', '0.37', '0.5', '0.70710678118654752440', '0.9', '0.999', '1']:
        worst = max(worst, abs(I_mellin(w) - psi(float(w))))
    print(f"(a) one root, Mellin integral vs psi(arcsin w), 8 values of w in [1e-3, 1]: max diff {float(worst):.1e}")


def check_one_factor():
    worst = 0
    for mu in (0.7, 1.0, 1.9):
        for w in (0.2, 0.6, 1.0):
            for t in (0.05, 0.3, 1.0, 3.0):
                dlog = lambda k: (1j/mu)/(1j*k/mu - w) - (1j/mu)/(1j*k/mu + w)
                I = quad(lambda k: (np.exp(-t*(mu*mu+k*k))*dlog(k)/(pi*1j)).real, 0, np.inf,
                         limit=400, epsabs=1e-14)[0]
                lhs = np.exp(-t*mu*mu*(1-w*w)) + I
                rhs = np.exp(-t*mu*mu*(1-w*w))*erf(mu*w*np.sqrt(t))
                worst = max(worst, abs(lhs-rhs))
    print(f"(b) one factor, part 2 vs part 3, 36 values of (mu, w, t): max diff {worst:.1e}")


def an_channels(h, a):
    A = pi*a/h
    return [(2*np.sin(pi*b/h), [(np.cos(A-pi*b/h), 1), (np.cos(A+pi*b/h), -1)]) for b in range(1, h)]


def check_hollowood(hmax=30):
    worst = 0
    for h in range(2, hmax+1):
        for a in range(1, h):
            cf = sum(0.5*m*sum(o*psi(w) for w, o in r) for m, r in an_channels(h, a))
            ref = -2*np.sin(pi*a/h)*(h/(2*pi) - 0.25/np.tan(pi/h))
            worst = max(worst, abs(cf-ref))
    print(f"(c) Hollowood a_n^(1), all a, h <= {hmax}: max |closed form - (eq-hollowood-mass)| = {worst:.1e}")


def check_remarks(hmax=10):
    w1 = w2 = 0
    for h in range(3, hmax+1):
        for a in range(1, h):
            chans = an_channels(h, a)
            closed = sum(0.5*m*sum(o*psi(w) for w, o in r) for m, r in chans)
            lev = {}
            for m, r in chans:
                lev.setdefault(round(m, 12), []).extend(r)
            tot = csum = 0
            for m, r in lev.items():
                ordw = {}
                for w, o in r:
                    ordw[round(w, 13)] = ordw.get(round(w, 13), 0) + o
                r = [(w, o) for w, o in ordw.items() if o]
                c = -2*m*sum(o*w for w, o in r); csum += c
                bound = sum(o*0.5*m*np.sqrt(1-w*w) for w, o in r if w > 0)
                dlog = lambda k: sum(o*1j/(1j*k - m*w) for w, o in r)
                f = lambda k: ((np.sqrt(m*m+k*k)*dlog(k) - 1j*c/(2*np.sqrt(m*m+k*k)))/(2j*pi)).real
                tot += bound + quad(f, 0, np.inf, limit=500, epsabs=1e-13)[0]
            w1 = max(w1, abs(tot + csum/(4*pi) - closed))
            if all(abs(w) > 1e-12 for m, r in chans for w, o in r):
                unp = sum(0.5*m*sum(o*(psi(w) + 0.5*np.sqrt(1-w*w)) for w, o in r) for m, r in chans)
                w2 = max(w2, abs(unp - closed))
    print(f"(d) a_n^(1), h <= {hmax}: phase-shift form + sum c_b/4pi vs closed form {w1:.1e};"
          f" unpaired bookkeeping (no threshold roots) vs closed form {w2:.1e}")


# ---------------------------------------------------------------- (e), (f): soliton sectors
def levels_of(S):
    """mass levels: mu and the net orders ord_w X_mu, roots w = cos(pi q/h)"""
    h = S.A.h
    lev = {}
    for nm, d, lb in S.channels:
        L = lev.setdefault(round(S.mb2[nm], 9), {'mu': np.sqrt(S.mb2[nm]), 'ord': {}})
        for q, o in S.exps[nm].items():
            w = round(np.cos(pi*q/h), 12)
            L['ord'][w] = L['ord'].get(w, 0) + o
    for L in lev.values():
        L['ord'] = {w: o for w, o in L['ord'].items() if o}
    return lev


def c_direct(S, Lx=80.0, N=16001):
    """c_mu = int tr(P_mu (M - M0)) dx, P_mu the spectral projector of M0 (field subspace)"""
    x = np.linspace(-Lx, Lx, N)
    Mf = np.einsum('ia,xij,jb->xab', S.P, S.bg.M(x), S.P)
    M0 = S.P.T @ S.A.M0 @ S.P
    ev, U = np.linalg.eigh(M0)
    out = {}
    for key in sorted(set(np.round(ev, 9))):
        cols = U[:, np.abs(ev - key) < 1e-7]
        integrand = np.einsum('ab,xba->x', cols @ cols.T, Mf - M0[None])
        out[key] = np.real(np.sum((integrand[1:] + integrand[:-1])/2)*(x[1]-x[0]))
    return out


def H_PT(lev, t):
    """part 3 of prp-net-count"""
    return sum(o*np.exp(-t*L['mu']**2*(1-w*w))*erf(L['mu']*w*np.sqrt(t))
               for L in lev.values() for w, o in L['ord'].items() if w > 0)


def closed_form(S):
    h = S.A.h
    return sum(0.5*np.sqrt(S.mb2[nm])*sum(o*psi(np.cos(pi*q/h)) for q, o in S.exps[nm].items())
               for nm, d, lb in S.channels)


def mellin_mass(lev, cdir):
    def F(t):
        C = sum(cdir[k]*np.exp(-t*k) for k in cdir)*np.sqrt(t)/(2*np.sqrt(pi))
        return t**-1.5*(H_PT(lev, t) + C)
    v = sum(quad(F, a, b, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
            for a, b in [(0, 1e-4), (1e-4, 1), (1, 20), (20, np.inf)])
    return -v/(4*np.sqrt(pi))


def check_sector(S, grid=None):
    lev = levels_of(S)
    paired = all(abs(w) > 1e-12 and L['ord'].get(round(-w, 12), 0) == -o
                 for L in lev.values() for w, o in L['ord'].items())
    cd = c_direct(S)
    born = max(abs(cd[k] + 2*L['mu']*sum(o*w for w, o in L['ord'].items())) for k, L in lev.items())
    cf, mm = closed_form(S), mellin_mass(lev, cd)
    line = (f"   {S.label:24s} pairing {str(paired):5s}  Born {born:.1e}  "
            f"Delta M {cf:+.10f}  Mellin diff {abs(cf-mm):.1e}")
    if grid:
        Hm, H0, trV, kmax, x = S.grid(grid, 48)
        ev = np.linalg.eigvals(Hm); ev0 = np.linalg.eigvalsh(H0.real)
        g = max(abs(np.sum(np.exp(-t*ev)) - np.sum(np.exp(-t*ev0)) - H_PT(lev, t)) for t in (0.3, 0.6, 1.0, 1.5))
        line += f"  grid(N={grid}, L=48) {g:.1e}"
    print(line, flush=True)


SECTORS = ['a3:1', 'a4:1', 'a4:2', 'd4:1', 'd4:3', 'd4:2', 'd5:1', 'd5:4', 'e6:1', 'e6:2',
           'c3(1)<a5:1', 'c3(1)<a5:2', 'b3(1)<d4:3', 'g2(1)<d4:1', 'f4(1)<e6:1', 'f4(1)<e6:3',
           'a4(2)<a4:1', 'd4(3)<e6:2', 'e6(2)<e7:1']
GRID_SECTORS = ['a4:1', 'c3(1)<a5:1', 'd4:3', 'b3(1)<d4:3', 'g2(1)<d4:1', 'f4(1)<e6:1']

if __name__ == "__main__":
    import warnings
    warnings.filterwarnings('ignore')
    grid = int(sys.argv[sys.argv.index('--grid')+1]) if '--grid' in sys.argv else None
    check_one_root(); check_one_factor(); check_hollowood(); check_remarks()
    from sectors import parent_sector, fold_sector
    print("(e) soliton sectors: pairing; max |c_mu(direct) - c_mu(roots)|; closed form; |closed form - Mellin|"
          + ("; (f) grid heat trace vs part 3" if grid else ""))
    for c in (GRID_SECTORS if grid else SECTORS):
        nm, a = c.rsplit(':', 1)
        S = fold_sector(nm, int(a)) if '<' in nm else parent_sector(nm, int(a))
        check_sector(S, grid)
