"""A2 library: soliton-preserving (SP) a_n^(1) reflection data in book conventions.

Built on code/boundary_common.py (an_vector, rcheck, cop, q_phys) only.
K_V, K_Vbar: the formulas of code/soliton_preserving.py, with omega as a parameter.
Book notation: h = N = n+1, T = h omega/2, mu = -i T theta/pi, (x)_I = sin(pi(mu+x)/2T)/sin(pi(mu-x)/2T).
"""
import sys
import numpy as np
sys.path.insert(0, '/home/gustav/Git/affine_toda/code')
from boundary_common import an_vector, rcheck, q_phys, cop, nodes_of, kron, report, prop_residual  # noqa


def K_V(n, om, l, eta, th):
    N, T = n + 1, (n + 1)*om/2
    return np.diag([np.exp(2*j*om*th)*np.sinh(eta + T*th) if j < l
                    else np.exp((2*j - N)*om*th)*np.sinh(eta - T*th) for j in range(N)])


def K_Vbar(n, om, l, etab, th):
    N, T = n + 1, (n + 1)*om/2
    return np.diag([np.exp(-2*j*om*th)*np.sinh(etab - T*th) if j < l
                    else np.exp(-(2*j - N)*om*th)*np.sinh(etab + T*th) for j in range(N)])


def etabar(n, om, l, eta):
    return eta + 1j*np.pi*om*(l - (n + 1)/2)


_cache = {}


def Rc(n, om, a, b, z):
    """R-check a(z) (x) b(1) -> b (x) a, a, b in {'V','W'} (W = Vbar), book rcheck (normalized on hw (x) hw)."""
    key = (n, om, a, b, complex(np.round(z, 13)))
    if key not in _cache:
        q = q_phys(om)
        A, B = an_vector(n, q, z, conj=(a == 'W')), an_vector(n, q, 1.0, conj=(b == 'W'))
        try:
            _cache[key] = rcheck(A, B)[0]
        except np.linalg.LinAlgError:          # rare LAPACK gesdd non-convergence: same equations, Gram + eigh
            _cache[key] = rcheck_gram(A, B)
    return _cache[key]


def rcheck_gram(A, B):
    """Same intertwining equations as boundary_common.rcheck; null vector by scipy's gesvd driver."""
    import scipy.linalg as sl
    dA, dB = A['dim'], B['dim']
    d = dA*dB
    M = np.vstack([np.kron(cop(A, B, kd, j).T, np.eye(d)) - np.kron(np.eye(d), cop(B, A, kd, j))
                   for j in nodes_of(A) for kd in 'efk'])
    _, s, Vh = sl.svd(M, full_matrices=False, lapack_driver='gesvd')
    assert s[-1] < 1e-9*s[0] and s[-2] > 1e-7*s[0]
    R = Vh[-1].conj().reshape(d, d).T
    i_in = A['hw']*dB + B['hw']
    i_out = B['hw']*dA + A['hw']
    return R/R[i_out, i_in]


def tensor_rep(reps):
    """Iterated Part VI coproduct on reps[0] (x) ... (x) reps[-1], built from boundary_common.cop."""
    A = reps[0]
    for B in reps[1:]:
        g = {}
        for j in nodes_of(A):
            for kd in 'efk':
                g[(kd, j)] = cop(A, B, kd, j)
        g['dim'] = A['dim']*B['dim']
        g['hw'] = None
        A = g
    return A
