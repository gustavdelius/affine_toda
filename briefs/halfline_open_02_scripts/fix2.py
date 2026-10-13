s = open('sym_vac.py').read()
s = s.replace("""        def run(c, dense=False):
            return solve_ivp(lambda x, Y: [Y[1], dV(Y[0])/2], [-X, 0.0], [c*np.exp(-m*X), m*c*np.exp(-m*X)],
                             method='DOP853', rtol=1e-12, atol=1e-16, dense_output=dense)
        c = brentq(lambda c: run(c).y[0, -1] - f0, 1e-3, 50.0, xtol=1e-14)""",
"""        def big(x, Y):
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
        c = brentq(g, cs[k], cs[k + 1], xtol=1e-14)""")
open('sym_vac.py', 'w').write(s)
