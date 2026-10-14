"""A6 task 2 helpers: d_n vector K-matrices of prp-dn-k from the book's coideal (boundary_common), R-check from the
book's weight-filtered solver (reflection_matrices.rcheck_w, definition only)."""
import sys
import numpy as np
TC = '/home/gustav/Git/affine_toda/textbook_code/'
sys.path.insert(0, TC)
from boundary_common import (dn_vector, q_phys, k_intertwiner, nodes_of, cop, kron, coideal_gen, _mk, intertwiners,
                             report)

_src = open(TC + 'reflection_matrices.py').read()
_ns = {}
exec(_src.split('\nrng = np.random.default_rng(11)')[0], _ns)          # imports
_body = _src.split('def rcheck_w(A, B):')[1].split('\n# ----')[0]
exec('def rcheck_w(A, B):' + _body, _ns)
rcheck_w = _ns['rcheck_w']


def chevalley(V, y):
    """Vbar(y): V composed with e <-> f, k -> k^-1 (as reflection_matrices.chevalley); V given at y = 1."""
    g = {}
    for j in nodes_of(V):
        g[('e', j)] = y*V[('f', j)]; g[('f', j)] = V[('e', j)]/y; g[('k', j)] = np.linalg.inv(V[('k', j)])
    return _mk(g, 0)


def eps_star(om):
    q = q_phys(om)
    return 1/np.sqrt((1 - q)*(1 - 1/q))


def eps_vec(n, om, s):
    e = eps_star(om)*np.ones(n + 1)
    if s == -1:
        e[0] *= -1
    if s == 0:
        e[:] = 0
    return e


def K_raw(n, om, y, s):
    q = q_phys(om)
    sols, _ = k_intertwiner(dn_vector(n, q, y), dn_vector(n, q, 1/y), [1/q]*(n + 1), eps_vec(n, om, s))
    assert len(sols) == 1, len(sols)
    return sols[0]
