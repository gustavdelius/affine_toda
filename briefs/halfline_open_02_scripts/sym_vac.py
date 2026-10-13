"""Uniform (C = +1) d_n vacua in the diagram-automorphism-invariant sector (d_4: Tzitzeica f'' = 2e^{2f}-2e^{-f};
d_5: sinh-Gordon f'' = 2e^f - 2e^{-f}), phi = f v; and the linearized reflection about them."""
import numpy as np
from scipy.integrate import quad, solve_ivp
from dn_classical import droots, masses, linear_K

def sym_data(n):
    E = np.eye(n)
    if n == 4:
        v = E[1] - E[2]
        V = lambda f: 4*(np.exp(-f) - 1) + 2*(np.exp(2*f) - 1)
        f0 = np.log(9/8)
    elif n == 5:
        v = E[1] - E[3]
        V = lambda f: 4*(np.exp(-f) - 1) + 4*(np.exp(f) - 1)
        f0 = np.log(4 - 2*np.sqrt(2))
    return v, V, f0

for n in (4, 5):
    al, nm = droots(n)
    v, V, f0 = sym_data(n)
    phi0 = f0*v
    # BC check
    fp0 = np.sqrt(V(f0))          # |v|^2 f'^2 / 2 = V  -> f' = sqrt(V) since |v|^2 = 2
    bc = fp0*v + (al.T*np.sqrt(nm)) @ np.exp(al @ phi0/2)
    W = quad(lambda f: 2*np.sqrt(V(f)), 0, f0)[0]       # int |phi'|^2 dx = int 2 f'^2 dx = int 2 f' df
    Bv = 2*np.sum(np.sqrt(nm)*(np.exp(al @ phi0/2) - 1))
    print(n, 'f0', f0, 'BC residual', np.abs(bc).max(), 'E =', W + Bv)


class Prof:
    """profile f(x) v on [-X, 0], integrated forward from the asymptotic data f = c e^{m x} (stable direction),
    c fixed by f(0) = f0."""
    def __init__(self, n, X=16.0):
        from scipy.optimize import brentq
        v, V, f0 = sym_data(n)
        dV = lambda f: (V(f + 1e-6) - V(f - 1e-6))/2e-6
        m = np.sqrt((V(1e-3) + V(-1e-3))/1e-6/2)
        def big(x, Y):
            return 5.0 - abs(Y[0])
        big.terminal = True
        def run(c, dense=False):
            return solve_ivp(lambda x, Y: [Y[1], dV(Y[0])/2], [-X, 0.0], [c*np.exp(-m*X), m*c*np.exp(-m*X)],
                             method='DOP853', rtol=1e-12, atol=1e-16, dense_output=dense, events=big)
        def g(c):
            r = run(c)
            return (r.y[0, -1] - f0) if r.status == 0 else 10.0
        cs = np.exp(np.linspace(np.log(1e-3), np.log(50), 60))
        vals = [g(c) for c in cs]
        k = next(i for i in range(len(cs) - 1) if vals[i] < 0 < vals[i + 1])
        c = brentq(g, cs[k], cs[k + 1], xtol=1e-14)
        self.s, self.v, self.n, self.m, self.c = run(c, True), v, n, m, c
        self.t = np.array([-X, 0.0])

    def sol(self, x):
        f = self.s.sol(x)[0]
        return np.concatenate([f*self.v, np.zeros(self.n)])


def bl(x, th, h):
    return np.sinh(th/2 + 1j*np.pi*x/(2*h))/np.sinh(th/2 - 1j*np.pi*x/(2*h))


def K1cl(th, n):
    h = 2*n - 2
    return -1/(bl(1, th, h)*bl(h - 1, th, h))*bl(h/2 - 1, th, h)*bl(h/2 + 1, th, h)/bl(h/2, th, h)**2


if __name__ == '__main__':
    np.set_printoptions(linewidth=180, precision=8, suppress=True)
    for n in (5, 4):
        m, U = masses(n)
        P = Prof(n)
        for th in (0.6, 0.3):
            Om = m[0]*np.cosh(th)
            print('  profile mass', P.m, 'c', P.c, 'f(0)', P.s.y[0, -1], flush=True)
            K, op = linear_K(n, np.ones(n + 1), P, Om)
            print(n, 'theta', th, 'open', op, '\n K =', K, '\n eig', np.linalg.eigvals(K), ' K1cl =', K1cl(th, n),
                  ' a_n-type -1/((1)(h-1)) =', -1/(bl(1, th, 2*n - 2)*bl(2*n - 3, th, 2*n - 2)))
