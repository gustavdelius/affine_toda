"""Sine-Gordon checks for chapter 3.

1. The integral representation eq-sg-s0 of S_0 satisfies crossing,
   S_0(theta) = S_T(i pi - theta) with S_T from eq-sg-smatrix, tested at
   theta = i pi/2 +- x, which stays inside the strip of convergence for lambda < 2.
2. The breather eq-sg-breather solves the sine-Gordon equation.
"""
import random
from mpmath import mp, mpf, quad, sinh, cosh, sin, exp, pi, inf
import sympy as sp

mp.dps = 25


def S0(th, xi):
    f = lambda t: sinh(t*(pi - xi)/2)*sin(th*t)/(t*sinh(xi*t/2)*cosh(pi*t/2))
    return -exp(-1j*quad(f, [0, 1, 10, inf]))


worst = 0
for lam in (mpf('1.5'), mpf('1.2')):
    xi = pi/lam
    for xx in (mpf('0.3'), mpf('1.1')):
        th = 1j*pi/2 + xx
        lhs = S0(1j*pi/2 - xx, xi)/S0(th, xi)
        rhs = sinh(lam*th)/sinh(lam*(1j*pi - th))
        worst = max(worst, abs(lhs - rhs))
print(('PASS ' if worst < 1e-20 else 'FAIL ') + f'crossing of S_0, max error {float(worst):.1e}')

x, t, u, m, b = sp.symbols('x t u m b', positive=True)
phi = 4/b*sp.atan(sp.tan(u)*sp.sin(m*t*sp.cos(u))/sp.cosh(m*x*sp.sin(u)))
eq = sp.diff(phi, t, 2) - sp.diff(phi, x, 2) + m**2/b*sp.sin(b*phi)
random.seed(1)
err = max(abs(float(sp.N(eq.subs({x: random.random(), t: random.random(), u: random.random()*1.5,
                                    m: 1.3, b: 0.7})))) for _ in range(5))
print(('PASS ' if err < 1e-12 else 'FAIL ') + f'breather solves sine-Gordon equation, max residual {err:.1e}')
