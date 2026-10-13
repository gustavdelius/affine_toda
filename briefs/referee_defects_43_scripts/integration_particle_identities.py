"""Integration check, ch. 43 item 2: prp-an-particle-transmission items 1-4 as exact identities for every h, both types.
Own sympy code from eq-dpt-particle.  x = e^{-(theta - theta_d)} (type I_-), y = e^{theta - theta_d} (type I_+),
u = e^{i gamma_a}, v = e^{i gamma_b}; a shift theta -> theta + i c multiplies x by e^{-i c}, y by e^{i c}."""
import sympy as sp
x, y = sp.symbols('x y', positive=True)
u, v = sp.symbols('u v', nonzero=True)
Tm = lambda X, U: (X + U)/(X + 1/U)          # type I_-: (e^{-D} + e^{i g})/(e^{-D} + e^{-i g})
Tp = lambda Y, U: (Y - U)/(Y - 1/U)          # type I_+: (e^{D} - e^{i g})/(e^{D} - e^{-i g})
res = {}
# bootstrap a + b < h: ub_a = gamma_b, ub_b = gamma_a; c = a + b, e^{i gamma_c} = u v
res['I- fuse a+b<h'] = sp.simplify(Tm(x/v, u)*Tm(x*u, v) - Tm(x, u*v))
res['I+ fuse a+b<h'] = sp.simplify(Tp(y*v, u)*Tp(y/u, v) - Tp(y, u*v))
# a + b > h: ub_a = pi - gamma_b, ub_b = pi - gamma_a; c = a + b - h, e^{i gamma_c} = -u v
# theta + i(pi - gamma_b): x -> x e^{-i pi} e^{i gamma_b} = -x v; theta - i(pi - gamma_a): x -> -x/u
res['I- fuse a+b>h'] = sp.simplify(Tm(-x*v, u)*Tm(-x/u, v) - Tm(x, -u*v))
res['I+ fuse a+b>h'] = sp.simplify(Tp(-y/v, u)*Tp(-y*u, v) - Tp(y, -u*v))
# unitarity + crossing: T_a(theta) T_abar(theta + i pi) = 1, e^{i gamma_abar} = -1/u, theta + i pi: x -> -x, y -> -y
res['I- unit/cross'] = sp.simplify(Tm(x, u)*Tm(-x, -1/u) - 1)
res['I+ unit/cross'] = sp.simplify(Tp(y, u)*Tp(-y, -1/u) - 1)
# Hermitian analyticity: Tbar_a(theta) = T_abar(i pi - theta) (a function of theta + theta_d), real theta, theta_d,
# |u| = 1.  With w = e^{theta + theta_d} > 0: T_abar(i pi - theta) has e^{-(i pi - theta - theta_d)} = -w;
# conj T_a(-theta) has x = w and conj(u) = 1/u.
w = sp.symbols('w', positive=True)
res['I- herm'] = sp.simplify(Tm(-w, -1/u) - Tm(w, 1/u))
res['I+ herm'] = sp.simplify(Tp(-1/w, -1/u) - Tp(1/w, 1/u))
# phase: for real x and |u| = 1 the denominator is the conjugate of the numerator
res['I- phase'] = sp.simplify((x + 1/u) - (x + u).subs(u, 1/u))
for k, val in res.items():
    print(f'{k}: {val}')
assert all(val == 0 for val in res.values())
# poles mod 2 pi i: I_-: x = -1/u -> D = -i(pi - gamma); I_+: y = 1/u -> D = -i gamma.  Physical strip 0 < Im < pi/2
# (and even 0 < Im < pi): the representatives -(pi - g), pi + g, -g, 2 pi - g lie outside for 0 < g < pi.
import numpy as np
ok = True
for h in range(2, 30):
    for a in range(1, h):
        g = np.pi*a/h
        for im in (-(np.pi - g), np.pi + g, -g, 2*np.pi - g):
            ok &= not (0 < im < np.pi)
print('no pole with 0 < Im theta < pi for real theta_d, h < 30:', ok)
print('ALL EXACT ZERO')
