import sympy as sp
q, e, lam, x = sp.symbols('q epsilon lambda x')
cub = (x - q**4)*(x + q**6)*lam**3 + q**6*(2 + q**2)*(lam**2 + x**2*lam) + x*(q**2 - 1)*(1 - 3*q**4 + q**8)*(lam**2 + lam) - q**2*(1 + 2*q**2)*(x**2*lam**2 + lam) + (1 - q**4*x)*(1 + q**6*x)
P2 = sp.Poly(sp.expand(cub), x)        # coefficients in x: x^2 c2(lam) + x c1(lam) + c0(lam); large x -> c2(lam) dominant
c2, c1, c0 = P2.all_coeffs()
print("leading c2(lam) =", sp.factor(c2))
# lam = q^4(1+delta): c2 ~ A delta^2, c1 ~ B  => x A delta^2 + B = 0 => delta^2 = -B/(A x)
A = sp.cancel(sp.diff(c2, lam, 2).subs(lam, q**4)*q**8/2)
B = sp.cancel(c1.subs(lam, q**4))
print("c2'(q^4) =", sp.simplify(sp.diff(c2, lam).subs(lam, q**4)))
r = sp.cancel(-B/A); print("delta^2 * x =", sp.factor(sp.numer(r)), "/", sp.factor(sp.denom(r)))
