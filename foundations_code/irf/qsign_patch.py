"""Physical sign of q for the a_2^(1) RSOS kinks.
At fixed spectral parameter x, q -> -q leaves [A-u]/[A], [1+u]/[1] and every quantum dimension
unchanged (integer A), and reverses [u]/[1]: it flips the sign of the off-diagonal face weight.
QSIGN=-1 applies this flip, i.e. uses q = -e^{i pi lam}, the value -e^{i pi lambda} of [tw1999]."""
import os
import irf_a2
QS = int(os.environ.get('QSIGN', '1'))
_W33 = irf_a2.A2RSOS.W33
def W33(self, a, b, c, d, u, sqrt_branch=None):
    w = _W33(self, a, b, c, d, u, sqrt_branch)
    if QS == -1 and d != b: w = -w
    return w
irf_a2.A2RSOS.W33 = W33
