import numpy as np, sympy as sp, cmath
th = sp.symbols('theta')
f = lambda s: (sp.sinh(th) + sp.I*s)/(sp.sinh(th) - sp.I*s)      # f_x with s = sin(pi x)
def check(S, name):
    pts = [0.31+0.2j, -0.77+1.1j, 1.3+2.4j]
    ev = lambda e, t: complex(e.subs(th, t).evalf(20))
    uni = max(abs(ev(S, t)*ev(S, -t) - 1) for t in pts)
    cr = max(abs(ev(S, 1j*np.pi - t) - ev(S, t)) for t in pts)
    ra = max(abs(np.conj(ev(S, t)) - ev(S, -np.conj(t))) for t in pts)
    print(f'{name}: unitarity {uni:.1e}, crossing {cr:.1e}, real analyticity {ra:.1e}')
# prp-cdd: solutions not of the listed form
y = sp.Rational(7, 10)
check(f(-sp.cosh(sp.pi*y)), 'staircase f_{-1/2+iy} (sin pi x = -cosh pi y <= -1)')
x = sp.Rational(3, 10) + sp.I*sp.Rational(1, 5)
check(f(sp.sin(sp.pi*x))*f(sp.sin(sp.pi*sp.conjugate(x))), 'pair f_x f_xbar, x = 0.3+0.2i')
check(f(sp.sin(sp.pi*x)), 'single f_x, x = 0.3+0.2i (should FAIL real analyticity)')
check(sp.exp(sp.I*sp.Rational(3, 2)*sp.sinh(th)), 'exp(1.5 i sinh theta) (no growth bound)')
# Lee-Yang residue: Res_{2 pi i/3} f_{2/3} = 2 i tan(2pi/3) = -2 sqrt3 i  => Gamma^2 = -2 sqrt 3
s = sp.sin(2*sp.pi/3)
print('Res f_{2/3} at 2 pi i/3 =', sp.nsimplify(sp.limit((th - 2*sp.pi*sp.I/3)*f(s), th, 2*sp.pi*sp.I/3)))
# exr-sg-expansion
b = sp.symbols('b', positive=True); m0 = sp.symbols('m0', positive=True)
xi = b**2/8/(1 - b**2/(8*sp.pi))
M = m0/(2*sp.sin(xi/2))
print('M series:', sp.series(M, b, 0, 2))
# Mostafazadeh: diag(i,i,-i) spectrum closed under conjugation as a set, but H^dagger not similar to H
H = np.diag([1j, 1j, -1j]); print('eig H:', np.linalg.eigvals(H), ' eig H^dagger:', np.linalg.eigvals(H.conj().T))
# exr-krein-collision counterexample: J = diag(1,-1), H = 0 has eigenvector (1,0) with [x,x] = 1 > 0
J = np.diag([1., -1.]); Hp = np.array([[0, 1.], [-1., 0]])
print('Hp J-self-adjoint:', np.allclose(J @ Hp.conj().T @ J, Hp), ' eig(0 + 0.01 Hp):', np.linalg.eigvals(0.01*Hp))
# free-field check of the conserved-charge sign: d_- T = d_+ Theta with T=(d_+ phi)^2, Theta=-m^2 phi^2/4
xp, xm, m = sp.symbols('x_p x_m m', positive=True)
k = sp.symbols('k', positive=True)
# plane-wave solution of 4 d_+ d_- phi + m^2 phi = 0: phi = cos(k x_p + m^2/(4k) x_m)
phi = sp.cos(k*xp + m**2/(4*k)*xm)
T = sp.diff(phi, xp)**2; Th = -m**2*phi**2/4
print('d_- T - d_+ Theta =', sp.simplify(sp.diff(T, xm) - sp.diff(Th, xp)))
x0, x1 = sp.symbols('x0 x1')
sub = {xp: x0 + x1, xm: x0 - x1}
for sgn, nm in [(1, 'T+Theta'), (-1, 'T-Theta')]:
    dens = (T + sgn*Th).subs(sub)
    # conserved iff d_0 dens is a total x1-derivative: integrate over a period in x1 numerically
    period = 2*sp.pi/(k - m**2/(4*k))
    val = sp.integrate(sp.diff(dens, x0).subs({k: 2, m: 1}), (x1, 0, period.subs({k: 2, m: 1})))
    print(f'd/dt int over a period of {nm}:', sp.simplify(val.subs(x0, sp.Rational(3, 10))).evalf())
