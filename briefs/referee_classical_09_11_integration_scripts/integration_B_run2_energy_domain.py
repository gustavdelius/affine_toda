import numpy as np
from integration_B_evalq import energy_momentum, A_int
# SG antikink theta = i pi, xi = i pi/2 (h = 2, M = 8)
E, P, eom, mt = energy_momentum(2, [(1, 1j*np.pi, 1j*np.pi/2)], t=0.0)
print('SG theta=i pi: E', E, 'formula', 8*np.cosh(1j*np.pi), 'P', P, 'eom', eom, 'min|tau|', mt)
# a2 single soliton with Re cosh < 0
h = 3; M = 6*np.sqrt(3)
th = 0.3 + 1j*(np.pi - 0.4)
E, P, eom, mt = energy_momentum(h, [(1, th, 0.2+0.5j)], t=0.1)
print('a2 Re cosh<0: E', E, 'M cosh', M*np.cosh(th), 'P', P, 'M sinh', M*np.sinh(th), 'min|tau|', mt)
# generic non-conjugate complex pair
sols = [(1, 0.4+0.3j, 0.1+0.7j), (2, -0.5+0.2j, -0.3+0.4j)]
E, P, eom, mt = energy_momentum(h, sols, t=0.05)
print('generic pair: E', E, 'pred', sum(M*np.cosh(s[1]) for s in sols), 'P', P, 'pred', sum(M*np.sinh(s[1]) for s in sols), 'min|tau|', mt)
# singular-band breather u = 0.2 (A > 0): regular slices give the formula?
u = 0.2; print('A(2iu) =', A_int(3, 1, 2, 2j*u), '1 - sin^2/cos^2 =', 1 - np.sin(np.pi/3)**2/np.cos(u)**2)
for t in [0.0, 0.7, 1.9, 3.3, 5.1]:
    sols = [(1, 0.1+1j*u, 0.3+0.2j), (2, 0.1-1j*u, -0.1+1.1j)]
    E, P, eom, mt = energy_momentum(h, sols, t=t, N=80001)
    print(f'  t={t}: E {E:.12f} formula {2*M*np.cosh(0.1)*np.cos(u):.12f} min|tau| {mt:.3e}')
