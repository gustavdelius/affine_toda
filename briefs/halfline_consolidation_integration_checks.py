"""Additional re-derivation of referee B5 and the unchecked item O.4.

Run with OPENBLAS_NUM_THREADS=1 python3 -u
briefs/halfline_consolidation_integration_checks.py [krein|semiclassical].
Uses the referee's independently implemented matrices and Gamma products;
the finite-B, exact charge and d5 Hessian regressions are in code/.
"""
from pathlib import Path
from functools import lru_cache
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / 'halfline_consolidation_referee_scripts'))


def stable_rcheck(A, B, bc):
    """Solve only weight-allowed entries, equilibrating the nonzero rows.

    Dense solves leak into forbidden weights at large rapidity, producing
    apparent 1e-8 departures from unit modulus on uniform-wall grids.
    No Gram matrix is formed: SVD acts on the equilibrated equations.
    """
    d = A['dim']*B['dim']
    js = list(bc.nodes_of(A))
    win = np.array([np.diag(bc.cop(A, B, 'k', j)) for j in js]).T
    wout = np.array([np.diag(bc.cop(B, A, 'k', j)) for j in js]).T
    allowed = [(r, c) for r in range(d) for c in range(d)
               if np.max(np.abs(wout[r]-win[c])) < 1e-10]
    basis = np.zeros((len(allowed), d, d), complex)
    for i, (r, c) in enumerate(allowed):
        basis[i, r, c] = 1
    rows = []
    for j in js:
        for kind in 'ef':
            lhs, rhs = bc.cop(A, B, kind, j), bc.cop(B, A, kind, j)
            rows.append((basis@lhs-rhs@basis).reshape(len(allowed), -1).T)
    equations = np.vstack(rows)
    norms = np.linalg.norm(equations, axis=1)
    keep = norms > 0
    equations = equations[keep]/norms[keep, None]
    _, singular, vh = np.linalg.svd(equations, full_matrices=False)
    assert singular[-2] > 1e-5 and singular[-1] < 1e-12
    R = np.einsum('i,ijk->jk', vh[-1].conj(), basis)
    return R/R[B['hw']*A['dim']+A['hw'], A['hw']*B['dim']+B['hw']]


def krein():
    import b5_boundary_krein as k
    @lru_cache(None)
    def stable(mu, nu, om, t1, t2):
        return stable_rcheck(k.rep(mu, om, t1), k.rep(nu, om, t2), k.bc)
    for mu in ('V', 'W'):
        for nu in ('V', 'W'):
            assert np.linalg.norm(stable(mu, nu, .57, .91, -.27)-k.R(mu, nu, .57, .91, -.27)) < 1e-12
    k.R = stable
    worst = 0
    for mus in (('V', 'V'), ('V', 'W'), ('V', 'V', 'W')):
        ths = [.91, .27, -.17][:len(mus)]
        C = k.Cmat
        for _ in mus[1:]:
            C = np.kron(C, k.Cmat)
        for left, right in (('U', 'N'), ('N', 'U'), ('N', 'N')):
            matrices = [k.transfer(j, mus, ths, .57, left, right) for j in range(len(mus))]
            comm = max(np.linalg.norm(a@b-b@a) for a in matrices for b in matrices)
            residual = max(np.linalg.norm(a.conj().T@C@a-C) for a in matrices)
            symmetry = max(np.linalg.norm(a@k.Zop(mus)-k.Zop(mus)@a) for a in matrices)
            other = [k.transfer(j, mus, ths, .57, left, right, 'other') for j in range(len(mus))]
            wrong = max(np.linalg.norm(a@b-b@a) for a in other for b in other)
            assert wrong > .01 if left == 'N' else wrong < 1e-10
            worst = max(worst, comm, residual, symmetry)
    assert worst < 1e-10
    print(f'PASS mirror convention, cyclic sectors, and Krein form: {worst:.2e}', flush=True)
    for walls, bracket, target in (('UN', (.69, .71), .6982), ('NN', (.53, .54), .5362)):
        lo, hi = bracket
        for _ in range(15):
            mid = (lo+hi)/2
            if k.broken(('V', 'V'), mid, *walls, [12., .001]) > 1e-8:
                hi = mid
            else:
                lo = mid
        threshold = (lo+hi)/2
        assert abs(threshold-target) < .0002
        print(f'PASS VV {walls} threshold: {threshold:.6f}', flush=True)
    for mus in (('V', 'V'), ('V', 'W')):
        grid = [.0001, .05, .4, 1.3, 3., 7., 12.]
        for om in (.15, .45, .9, 1.3):
            dev, pair = max((k.broken(mus, om, 'U', 'U', [a, b]), (a, b))
                            for a in grid for b in grid if a != b)
            print(f'UU {mus} omega={om}: {dev:.3e} at {pair}', flush=True)
            assert dev < 1e-8, (mus, om, dev, pair)
        dev = max(k.broken(mus, .45, 'U', 'N', [a, b]) for a in grid for b in grid if a != b)
        assert dev < 1e-8
    print('PASS uniform-wall grids and UN below threshold', flush=True)


def semiclassical():
    import mpmath as mp
    import b1_dn_scalar as scalar
    mp.mp.dps = 30
    theta = mp.mpf('.83'); delta = mp.mpf('1e-7')
    values = []
    for scale in (40, 80, 160):
        om = mp.mpf(scale) + mp.mpf('.37'); T = 3*om
        mup = -1j*T*(theta+delta)/mp.pi
        mum = -1j*T*(theta-delta)/mp.pi
        row = []
        for sg in (1, -1):
            ratio = scalar.A_s(mup, T, om, sg)/scalar.A_s(mum, T, om, sg)
            row.append(mp.im(mp.log(ratio))/(2*delta*om))
        bulk = mp.im(scalar.logF(2*mup, T, om)-scalar.logF(2*mum, T, om))/(4*delta*om)
        row.append(bulk)
        values.append(row)
        print(f'omega={om}: phase slopes/omega = ' + ', '.join(mp.nstr(x, 10) for x in row), flush=True)
    extrapolated = [2*values[-1][i]-values[-2][i] for i in range(3)]
    errors = [abs(extrapolated[i]-extrapolated[2]) for i in (0, 1)]
    assert max(errors) < mp.mpf('.002')
    print('PASS equal leading semiclassical phase derivatives; extrapolated differences '
          + ', '.join(mp.nstr(x, 5) for x in errors), flush=True)


if __name__ == '__main__':
    {'krein': krein, 'semiclassical': semiclassical}[sys.argv[1]]()
