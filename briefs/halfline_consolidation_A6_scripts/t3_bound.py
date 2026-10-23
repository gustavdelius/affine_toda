"""Bound states of the linearized problem about the C = +1 kinks (zeros of the BC determinant below threshold),
compared with the physical-strip poles of prod_j S_aj(theta, B = 1); and no negative modes (Om^2 < 0)."""
import numpy as np
import cmath
from scipy.optimize import brentq
from cl_lin import Kink, bstate_det

TC = '/home/gustav/Git/affine_toda/code/'
ns = {}
src = open(TC + 'smatrix_simply_laced.py').read().split('\nB = 0.37')[0]
src = src.replace("open('toda_classical.py')", "open('" + TC + "toda_classical.py')")
exec(src, ns)
setup, S = ns['setup'], ns['S']

for n in (4, 5):
    K = Kink(n)
    M0 = (K.A.T*K.nj) @ K.A
    ma = np.sqrt(np.linalg.eigvalsh(M0))
    mmin = ma.min()
    Oms = np.linspace(1e-3, mmin*(1 - 1e-6), 400)
    d = np.array([bstate_det(K, Om).real for Om in Oms])
    roots = []
    for i in range(len(Oms) - 1):
        if d[i]*d[i + 1] < 0:
            roots.append(brentq(lambda o: bstate_det(K, o).real, Oms[i], Oms[i + 1], xtol=1e-14))
    # minimum of |det| to detect even-order zeros
    absd = np.abs(d)
    loc = [Oms[i] for i in range(1, len(Oms) - 1) if absd[i] < absd[i - 1] and absd[i] < absd[i + 1]]
    # predicted: poles of prod_j S_aj(theta, 1) at theta = i u, 0 < u < pi/2, for each particle a, as frequencies m_a cos u
    h, blocks, mass, _ = setup('d', n)
    scale = mmin/mass.min()
    pred = set()
    for a in range(n):
        for x2 in range(1, h):          # candidate poles at theta = i pi x/h
            u = np.pi*x2/h
            if u >= np.pi/2:
                continue
            th = 1j*u
            val = abs(np.prod([S(blocks[(a, j)], h, 1.0, th*(1 + 1e-9)) for j in range(n)]))
            if val > 1e6:
                pred.add(round(mass[a]*scale*np.cos(u), 10))
    # negative modes: Om = i w
    neg = min(abs(bstate_det(K, 1j*w)) for w in np.linspace(0.01, 5, 60))
    print(f'd_{n}: bound-state frequencies (sign changes of det): {np.round(roots, 12)}; local minima of |det|: '
          f'{np.round(loc, 4)}; predicted from poles of prod_j S_aj(theta,1): {sorted(pred)}; '
          f'm_1 sin(pi/h) = {mmin*np.sin(np.pi/h):.12f}; min |det| for Om^2 < 0: {neg:.2e}')
