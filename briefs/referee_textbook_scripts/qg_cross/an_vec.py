"""sl3 vector: Part VI coproduct (rsolve.intertwiner); eigenvalues of R_VI(x) are 1 and (1-q^-2 x)/(x-q^-2) = rho_7(x;1/q)."""
import numpy as np
from rsolve import intertwiner
q = 0.83*np.exp(0.61j)
def Em(a, b):
    M = np.zeros((3, 3), complex); M[a, b] = 1; return M
def rep(x):
    e = {1: Em(0, 1), 2: Em(1, 2), 0: x*Em(2, 0)}; f = {1: Em(1, 0), 2: Em(2, 1), 0: Em(0, 2)/x}
    k = {1: np.diag([q, 1/q, 1]), 2: np.diag([1, q, 1/q]), 0: np.diag([1/q, 1, q])}
    return e, f, k, {0: q, 1: q, 2: q}
W = np.array([[1, 0], [-1, 1], [0, -1]], float)   # weights in a basis where alpha_1=(2,-1)... only equality matters
W = np.array([[2/3, 1/3], [-1/3, 1/3], [-1/3, -2/3]])
x = 1.7 + 0.4j
Rc, null, _, _ = intertwiner(rep(x), rep(1.0), [0, 1, 2], W, W); Rc = Rc/Rc[0, 0]
ev = np.linalg.eigvals(Rc); r = (1 - q**-2*x)/(x - q**-2)
print('null space dim', null, '; eigenvalues', np.round(np.sort_complex(ev), 6), '; rho_7(x;1/q) =', np.round(r, 6))
