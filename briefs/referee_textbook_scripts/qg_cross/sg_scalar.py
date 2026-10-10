"""Part VI scalar factor F = c f for n = 1 versus Zamolodchikov's S_0 (eq-sg-s0), at lambda = omega."""
import os; os.environ['NSPIN'] = '1'; os.environ['OMEGA'] = '1.5'; os.environ['JTRUNC'] = '20000'
import numpy as np
from mpmath import mp, quad, sinh, cosh, sin, exp, pi, inf, mpf
import spinorn as S
mp.dps = 20
lam = 1.5; xi = pi/lam
def S0(th):
    f = lambda t: sinh(t*(pi - xi)/2)*sin(th*t)/(t*sinh(xi*t/2)*cosh(pi*t/2))
    return complex(-exp(-1j*quad(f, [0, 1, 10, inf])))
for th in (0.3, 0.9, 0.4+0.5j, 1.2-0.8j):
    print(f"theta={th}:  F/S_0 = {np.round(S.F(th)/S0(th), 6)}")
