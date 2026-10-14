import sympy as sp
import integration_a1_dirichlet as m
# At lambda^{2h} = 1 the factor (lambda^{2h} - 1) vanishes; the obstruction should then disappear.
for name, rep, h in [('a2 V', m.sl_vector(2), 3), ('d4 V', m.so_vector(8), 6)]:
    w = m.generic_w(rep[6])
    for lam in (sp.Integer(1), sp.Integer(-1), sp.exp(sp.I*sp.pi/h), sp.Rational(3, 7)):
        m.lam = lam
        print(name, 'lambda =', lam, 'lambda^(2h) =', sp.simplify(lam**(2*h)),
              'dims:', [m.solve_linear(rep, p, w) for p in ('antipode', 'transpose')])
