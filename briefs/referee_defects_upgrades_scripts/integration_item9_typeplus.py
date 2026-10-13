"""Item 9 (ch. 43, prp-an-particle-transmission, type I_+): T^+_a = (e^{theta-theta_d} - u)/(e^{theta-theta_d} - 1/u) in
the type-I_- variable x = e^{-(theta-theta_d)} is F(x,u) = (1 - x u)/(1 - x/u).  Exact sympy checks of unitarity/crossing,
both fusings with the explicit cancelling factors, Hermitian analyticity and the pole.  No book or referee code."""
import sympy as sp
x, u, v, w = sp.symbols('x u v w', nonzero=True)
F = lambda X, U: (1 - X*U)/(1 - X/U)
ok = sp.simplify(((sp.exp(-sp.Symbol('D')) ** -1 - u)/(sp.exp(-sp.Symbol('D')) ** -1 - 1/u)).subs(sp.exp(-sp.Symbol('D')), x) - F(x, u)) == 0
ok &= sp.simplify(F(x, u)*F(-x, -1/u) - 1) == 0                                  # T_a(theta) T_abar(theta + i pi) = 1
ok &= sp.simplify(F(x/v, u)*F(x*u, v) - F(x, u*v)) == 0                           # a + b < h
ok &= sp.simplify((1 - x*u/v)/(1 - x/(u*v))*(1 - x*u*v)/(1 - x*u/v) - F(x, u*v)) == 0
ok &= sp.simplify(F(-x*v, u)*F(-x/u, v) - F(x, -u*v)) == 0                        # a + b > h
ok &= sp.simplify((1 + x*u*v)/(1 + x*v/u)*(1 + x*v/u)/(1 + x/(u*v)) - F(x, -u*v)) == 0
# Hermitian analyticity: Tbar^+_a(theta) = T^+_abar(i pi - theta): x -> -w (w = e^{theta+theta_d}), u -> -1/u;
# conj T^+_a(-theta) for real w, |u| = 1: (1 - w/u)/(1 - w u)
ok &= sp.simplify(F(-w, -1/u) - (1 - w/u)/(1 - w*u)) == 0
# pole at x = u, i.e. theta - theta_d = -i gamma_a
ok &= sp.simplify(sp.denom(sp.together(F(x, u))).subs(x, u)) == 0
print('PASS' if ok else 'FAIL', ': item 9 - type-I_+ factors in x: crossing, both fusings (cancelling 1 - xu/v and 1 + xv/u), Hermitian analyticity, pole at x = u')
