import numpy as np
from integration_B_evalq import energy_momentum, charge
h = 3; M = 2*h*2*np.sin(np.pi/3)
for sols, lab in [([(1, 0.0, 0.3j)], 'sp1 rest'), ([(1, 0.4, 0.2+0.5j)], 'sp1 th=0.4'),
                  ([(1, 0.3+0.4j, 0.2+0.5j)], 'sp1 complex'),
                  ([(1, 0.23+1j*np.pi/3, 1j*np.pi/2), (2, 0.23-1j*np.pi/3, 1j*np.pi/2)], 'breather')]:
    E, P, eom, mt = energy_momentum(h, sols, t=0.17)
    pred = sum(M*np.cosh(s[1]) for s in sols), sum(M*np.sinh(s[1]) for s in sols)
    print(lab, 'E', E, 'pred', pred[0], 'P', P, 'pred', pred[1], 'eom', eom, 'min|tau|', mt)
for s in (2, 4):
    q1, _ = charge(s, [(1, 0.0, 0.3j)]); q2, _ = charge(s, [(2, 0.0, 0.3j)])
    qc, _ = charge(s, [(1, 0.3+0.4j, 0.2+0.5j)])
    print('spin', s, 'q1', q1, 'q2/q1', q2/q1, 'complex/(q1 e^{s th})', qc/(q1*np.exp(s*(0.3+0.4j))))
    for xis, t in [((1j*np.pi/2, 1j*np.pi/2), 0.17), ((0.2+1j*np.pi/2, -0.2+1j*np.pi/2), 0.43), ((0.1+0.2j, 0.1-0.2j+1j*np.pi), 0.0)]:
        th0, u = 0.23, np.pi/3
        Qb, mt = charge(s, [(1, th0+1j*u, xis[0]), (2, th0-1j*u, xis[1])], t=t)
        print('   breather Q/q1', Qb/q1, 'pred', np.exp(s*(th0+1j*u)) - np.exp(s*(th0-1j*u)), 'min|tau|', mt)
