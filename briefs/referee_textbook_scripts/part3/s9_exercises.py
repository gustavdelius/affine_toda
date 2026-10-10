"""Exercise spot checks: exr-lax (a_1, 2x2), exr-spin3 (Q_3 on the sG kink), exr-shg-cubic (a_2 C_111), exr-fusion."""
import sympy as sp, numpy as np
xp, xm, lam, m, b = sp.symbols('xp xm lambda m beta', positive=True)
f = sp.Function('f')(xp, xm)
E1 = sp.Matrix([[0, 1], [0, 0]]); E0 = sp.Matrix([[0, 0], [1, 0]]); H = sp.Matrix([[1, 0], [0, -1]])/sp.sqrt(2)
a1, a0 = sp.sqrt(2), -sp.sqrt(2)
Ap = b/2*sp.diff(f, xp)*H + m*lam*(sp.exp(b*a1*f/2)*E1 + sp.exp(b*a0*f/2)*E0)
Am = -b/2*sp.diff(f, xm)*H - m/lam*(sp.exp(b*a1*f/2)*E1.T + sp.exp(b*a0*f/2)*E0.T)
Fc = sp.diff(Am, xp) - sp.diff(Ap, xm) + Ap*Am - Am*Ap
eom = -(m**2/b)*(a1*sp.exp(b*a1*f) + a0*sp.exp(b*a0*f))       # = -(2 sqrt2 m^2/beta) sinh(sqrt2 beta f)
print('exr-lax: zero curvature with eom substituted:', sp.simplify(Fc.subs(sp.Derivative(f, xp, xm), eom)) == sp.zeros(2))
print('         eom =', sp.simplify(eom.rewrite(sp.sinh)))
# exr-spin3: sine-Gordon (beta^2 -> -beta^2): T4=(phi'')^2-(b^2/4)(phi')^4, Theta2=m^2 cos(b phi)(phi')^2 on kink phi=(4/b)atan(e^{mx})
x = sp.symbols('x', real=True)
ph = 4/b*sp.atan(sp.exp(m*x))
dens = sp.diff(ph, x, 2)**2 - b**2/4*sp.diff(ph, x)**4 + m**2*sp.cos(b*ph)*sp.diff(ph, x)**2
val = sp.Integral(sp.simplify(dens.subs({m: 1, b: 1})), (x, -sp.oo, sp.oo)).evalf(30)
print('exr-spin3: Q_3(static kink) at m=b=1:', val, ' (-16/3 =', sp.N(-sp.Rational(16, 3), 20), ')')
# exr-shg-cubic: a_2, Z_3 basis
w = np.exp(2j*np.pi/3); e1 = np.array([1, w, w**2])/np.sqrt(3)
al = [np.array(v, float) for v in ([1, -1, 0], [0, 1, -1], [-1, 0, 1])]
C111 = sum((a @ e1)**3 for a in al); C112 = sum((a @ e1)**2*(a @ e1.conj()) for a in al)
area = np.sqrt(3)/4*3
print('exr-shg-cubic: C_111 =', np.round(C111, 12), ' C_112 =', np.round(C112, 12), ' |C111|/((4/sqrt3) area) =', abs(C111)/(4/np.sqrt(3)*area))
# exr-fusion identity
ok = all(abs(np.sin(np.pi*a/h)*np.exp(1j*np.pi*bb/h) + np.sin(np.pi*bb/h)*np.exp(-1j*np.pi*a/h) - np.sin(np.pi*(a+bb)/h)) < 1e-14
         for h in range(2, 12) for a in range(1, h) for bb in range(1, h))
print('exr-fusion identity m_a e^{i pi b/h} + m_b e^{-i pi a/h} = m_{a+b}:', ok)
