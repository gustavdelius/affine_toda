"""Which combination is conserved if d_- T = d_+ Theta (x^pm = x^0 +- x^1)?  Free massive field, standing wave."""
import sympy as sp
xp, xm, x0, x1 = sp.symbols('x_p x_m x0 x1', real=True)
m = 2
modes = [sp.Rational(2), sp.Rational(1, 2)]          # k=2 and k=1/2: spatial wavenumbers +-3/2, frequency 5/2
phi = sum(sp.cos(k*xp + sp.Rational(m**2, 4)/k*xm) for k in modes)
print('EOM 4 d+d- phi + m^2 phi =', sp.simplify(4*sp.diff(phi, xp, xm) + m**2*phi))
T = sp.diff(phi, xp)**2; Th = -sp.Rational(m**2, 4)*phi**2
print('d_- T - d_+ Theta =', sp.simplify(sp.diff(T, xm) - sp.diff(Th, xp)))
sub = {xp: x0 + x1, xm: x0 - x1}
L = 4*sp.pi/3                                            # spatial period
for sgn, nm in [(1, 'T + Theta (book)'), (-1, 'T - Theta')]:
    Q = sp.integrate(sp.expand(sp.expand_trig((T + sgn*Th).subs(sub))), (x1, 0, L))
    print(f'Q(t) = int (' + nm + ') dx1 =', sp.simplify(Q))
