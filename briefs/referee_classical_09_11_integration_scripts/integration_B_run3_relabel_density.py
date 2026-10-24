import numpy as np, itertools
from integration_B_evalq import energy_momentum, fields, tau_terms
h = 3; M = 6*np.sqrt(3)
# two-soliton with one constituent at Re cosh < 0
th1, th2 = 0.3 + 1j*(np.pi - 0.4), -0.2 + 0.3j
sols = [(1, th1, 0.2+0.5j), (2, th2, -0.3+0.9j)]
E, P, eom, mt = energy_momentum(h, sols, t=0.1)
print('two-soliton, Re cosh th1 < 0: E', E, 'pred', -M*np.cosh(th1) + M*np.cosh(th2), 'P', P, 'pred', -M*np.sinh(th1)+M*np.sinh(th2), 'min|tau|', mt)
# energy density = 2 sum_j d_x^2 ln tau_j (m = beta = 1, alpha^2 = 2)
sols = [(1, 0.23+1j*np.pi/3, 0.1+0.2j), (2, 0.23-1j*np.pi/3, 0.1+0.8j)]
x = np.linspace(-6, 6, 13); t = 0.4
F = fields(h, sols, x, t)
sq = lambda v: np.sum(v*v, axis=0)
eps = -0.25*(sq(F['qp']) + sq(F['q'][1])) + np.sum(1 - F['Ej'], axis=0)
T0 = np.zeros((h, len(x)), complex); T1 = np.zeros_like(T0); T2 = np.zeros_like(T0)
for c, k, w in tau_terms(h, sols):
    ex = np.exp(k*x + w*t); T0 += np.outer(c, ex); T1 += np.outer(c, k*ex); T2 += np.outer(c, k*k*ex)
d2 = T2/T0 - (T1/T0)**2
print('density identity err', np.max(np.abs(eps - 2*np.sum(d2, axis=0))))
