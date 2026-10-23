"""The one-loop mass (eq-one-loop-mass) on a Dirichlet lattice (sec-resolution-semiclassical-mass-disagreements).

Run from this directory:  python3 dirichlet_lattice.py [case ...]     (about 30 min for all cases)

The fluctuation operator is discretized by second-order finite differences with spacing eps on [-L/2, L/2], with
u = 0 at both ends. On the lattice

    Delta M(eps) = 1/2 Tr(sqrt A - sqrt A0) - 1/4 sum_k tr[(M^2(x_k) - M^2) (A0^{-1/2})_kk],

where the second term is the normal-ordering counterterm built from the lattice equal-time propagator. The full
complex spectrum of A is used, and the eigenvalue of the zero mode, which the lattice moves away from 0 by O(eps^2),
is removed. The error is eps^2 (a + b log eps), with the logarithm coming from the log-divergent counterterm, so
Delta M is extrapolated from eps = 0.075, 0.05, 0.035 with that form. A pure eps^2 extrapolation leaves 1e-4 to 3e-4,
and for the heaviest c_4^(1) solitons eps = 0.1 is still outside the asymptotic regime.

Cases: the sine-Gordon kink (A = -d^2 + 1 - 2 sech^2 x, Delta M = -1/pi); the a_n^(1) solitons for n <= 5, against
Hollowood's masses (eq-hollowood-mass); the c_n^(1) solitons a < n for n <= 5, realized as sigma-invariant sectors of
a_{2n-1}^(1), against the closed form of prp-one-loop-mass, which counts the threshold blocks twice. Each case is
also compared with the closed form that counts every double zero once, which the lattice must miss.
L is chosen so that M^2(x) - M^2 is below 1e-6 at the ends; the result does not depend on L.
"""
import sys, os, time, warnings
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'hirota'))
from heat_trace_checks import closed_form, levels_of
import fold
from sectors import parent_sector, fold_sector

warnings.filterwarnings('ignore')
fold.FOLDS.setdefault('c5(1)<a9', ('a9', 'refl_a'))
EPS, TOL = (0.075, 0.05, 0.035), 2e-4


def lattice_dM(Mf, M0, L, eps):
    """Delta M on the Dirichlet lattice; Mf(x) -> (len(x), r, r), M0 (r, r) real symmetric."""
    Ns = int(round(L/eps)) - 1
    x = -L/2 + eps*np.arange(1, Ns + 1)
    Mx = Mf(x)
    r = M0.shape[0]
    D = (np.diag(np.full(Ns, 2.0)) - np.diag(np.ones(Ns - 1), 1) - np.diag(np.ones(Ns - 1), -1))/eps**2
    A = np.kron(D, np.eye(r)).astype(complex)
    for i in range(Ns):
        A[i*r:(i + 1)*r, i*r:(i + 1)*r] += Mx[i]
    lam = np.linalg.eigvals(A)
    iz = np.argmin(np.abs(lam)); zero = lam[iz]; lam = np.delete(lam, iz)
    mu, Sv = np.linalg.eigh(D)                                     # sine modes of the Dirichlet Laplacian
    m2, V = np.linalg.eigh(M0)
    lam0 = (mu[:, None] + m2[None, :]).ravel()
    g = (Sv**2) @ (1/np.sqrt(mu[:, None] + m2[None, :]))         # [A0^{-1/2}]_kk in the eigenbasis of M0
    G = np.einsum('aj,kj,bj->kab', V, g, V)
    ct = -0.25*np.einsum('kab,kba->', Mx - M0[None], G)
    return 0.5*(np.sum(np.sqrt(lam)) - np.sum(np.sqrt(lam0))) + ct, zero


def extrapolate(eps, d):
    e = np.array(eps)
    A = np.column_stack([np.ones_like(e), e**2, e**2*np.log(e)])
    return np.linalg.solve(A, np.array(d))[0]


def once_form(S):
    """each bound-state zero of order o >= 2 counted once: subtract (o-1) omega/2, omega = sqrt(lambda_zeta)"""
    lev = levels_of(S)
    extra = sum((o - 1)*0.5*L['mu']*np.sqrt(1 - w*w) for L in lev.values() for w, o in L['ord'].items() if 0 < w < 1 and o >= 2)
    return closed_form(S) - extra


def case(spec):
    if spec == 'sg':
        Mf = lambda x: (1 - 2/np.cosh(x)**2)[:, None, None].astype(complex)
        return 'sine-Gordon kink', Mf, np.eye(1), -1/np.pi, None, 30.0
    nm, a = spec.rsplit(':', 1)
    S = fold_sector(nm, int(a)) if '<' in nm else parent_sector(nm, int(a))
    Mf = lambda x: np.einsum('ia,xij,jb->xab', S.P, S.bg.M(x), S.P)
    M0 = (S.P.T @ S.A.M0 @ S.P).real
    once = once_form(S) if '<' in nm else None
    return S.label, Mf, M0, closed_form(S), once, max(24.0, 28/float(S.mu))


CASES = (['sg', 'a2:1', 'a3:1', 'a3:2', 'a4:1', 'a4:2', 'a5:1', 'a5:2', 'a5:3']
         + [f'c{n}(1)<a{2*n-1}:{a}' for n in (2, 3, 4, 5) for a in range(1, n)])

if __name__ == "__main__":
    for spec in sys.argv[1:] or CASES:
        t = time.time()
        label, Mf, M0, cf, once, L = case(spec)
        ds, zeros = zip(*[lattice_dM(Mf, M0, L, e) for e in EPS])
        d0 = extrapolate(EPS, [d.real for d in ds])
        ok = abs(d0 - cf) < TOL and max(abs(d.imag) for d in ds) < 1e-4 and max(abs(z) for z in zeros) < 0.05
        if once is not None:
            ok = ok and abs(d0 - once) > 10*TOL
        tail = "" if once is None else f"; counting each double zero once gives {once.real:.6f}"
        print(f"{'PASS' if ok else 'FAIL'} {label} (L = {L:.0f}): lattice "
              + ", ".join(f"{d.real:.6f}" for d in ds) + f" at eps = {EPS}; extrapolated {d0:.6f}, "
              f"closed form {cf:.6f}, diff {abs(d0 - cf):.1e}{tail} ({time.time() - t:.0f} s)", flush=True)
