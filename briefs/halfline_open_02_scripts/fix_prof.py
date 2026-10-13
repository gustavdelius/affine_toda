s = open('sym_vac.py').read()
old = s[s.index('class Prof:'):s.index('def bl(')]
new = '''class Prof:
    """profile f(x) v on [-X, 0], integrated forward from the asymptotic data f = c e^{m x} (stable direction),
    c fixed by f(0) = f0."""
    def __init__(self, n, X=16.0):
        from scipy.optimize import brentq
        v, V, f0 = sym_data(n)
        dV = lambda f: (V(f + 1e-6) - V(f - 1e-6))/2e-6
        m = np.sqrt((V(1e-3) + V(-1e-3))/1e-6/2)
        def run(c, dense=False):
            return solve_ivp(lambda x, Y: [Y[1], dV(Y[0])/2], [-X, 0.0], [c*np.exp(-m*X), m*c*np.exp(-m*X)],
                             method='DOP853', rtol=1e-12, atol=1e-16, dense_output=dense)
        c = brentq(lambda c: run(c).y[0, -1] - f0, 1e-3, 50.0, xtol=1e-14)
        self.s, self.v, self.n, self.m, self.c = run(c, True), v, n, m, c
        self.t = np.array([-X, 0.0])

    def sol(self, x):
        f = self.s.sol(x)[0]
        return np.concatenate([f*self.v, np.zeros(self.n)])


'''
s = s.replace(old, new)
s = s.replace("            Om = m[0]*np.cosh(th)\n", "            Om = m[0]*np.cosh(th)\n            print('  profile mass', P.m, 'c', P.c, 'f(0)', P.s.y[0, -1], flush=True)\n")
open('sym_vac.py', 'w').write(s)
