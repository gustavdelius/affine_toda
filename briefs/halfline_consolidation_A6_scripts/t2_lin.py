"""Linearized reflection about the d_4, d_5 C = +1 kinks vs prod_j S_aj(theta, B = 1) (book S-matrix) and vs the
B -> 0 limit of the book's eq-dn-uniform-amplitude."""
import numpy as np
import cmath
from cl_lin import Kink, kmatrix

TC = '/home/gustav/Git/affine_toda/textbook_code/'
ns = {}
src = open(TC + 'smatrix_simply_laced.py').read().split('\nB = 0.37')[0]
src = src.replace("open('toda_classical.py')", "open('" + TC + "toda_classical.py')")
exec(src, ns)
setup, S = ns['setup'], ns['S']


def bl(x, th, h):
    return cmath.sinh(th/2 + 1j*cmath.pi*x/(2*h))/cmath.sinh(th/2 - 1j*cmath.pi*x/(2*h))


def K1_book(th, h, B):
    return (bl(h - B/2, th, h)*bl(h/2 - 1, th, h)*bl(h/2 + 1 - B/2, th, h)
            / (bl(h - 1, th, h)*bl(1 - B/2, th, h)*bl(h/2, th, h)*bl(h/2 - B/2, th, h)))


np.set_printoptions(linewidth=200, precision=6)
for n in [int(a) for a in __import__("sys").argv[1:]] or (4, 5):
    K = Kink(n)
    h, blocks, mass, _ = setup('d', n)
    # absolute masses of the book's particles: m_a = sqrt(2) * pf... match by ratio with the classical masses
    M0 = (K.A.T*K.nj) @ K.A
    ev, U = np.linalg.eigh(M0)
    ma = np.sqrt(ev)
    print(f'd_{n}: classical masses {np.round(ma, 8)}; book PF ratios {np.round(np.sort(mass), 8)}; '
          f'ratio check {np.allclose(np.sort(ma)/ma.min(), np.sort(mass))}')
    # predicted classical factor for each node a of setup: prod_j S_aj(theta, 1)
    pred = lambda a, th: np.prod([S(blocks[(a, j)], h, 1.0, th) for j in range(n)])
    scale = ma.min()/mass.min()
    worst_diag, worst_off, worst_K1 = 0, 0, 0
    for Om in (1.1*ma.min(), 1.5*ma.min(), 0.99*ma.max(), 1.3*ma.max(), 2.7*ma.max()):
        mm, Uc, op, Km = kmatrix(K, Om)
        # off-diagonal between channels
        off = np.abs(Km - np.diag(np.diag(Km))).max()
        worst_off = max(worst_off, off)
        line = []
        for i, a in enumerate(op):
            th = np.arccosh(Om/mm[a])
            # book particles with this mass
            cands = [b for b in range(n) if abs(mass[b]*scale - mm[a]) < 1e-9]
            vals = [pred(b, th) for b in cands]
            assert max(abs(v - vals[0]) for v in vals) < 1e-12
            worst_diag = max(worst_diag, abs(Km[i, i]/vals[0] - 1))
            line.append(f'm={mm[a]:.4f} th={th:.3f} K={Km[i, i]:.10f} pred={vals[0]:.10f}')
            if abs(mm[a] - ma.min()) < 1e-9 and (n == 5 or True):
                if n == 5 or True:
                    worst_K1 = max(worst_K1, abs(Km[i, i]/K1_book(th, h, 1e-12) - 1)) if abs(mm[a] - mass[0]*scale) < 1e-9 else worst_K1
        print(f'  Om={Om:.4f} open={len(op)} max|offdiag|={off:.1e}')
        for l in line:
            print('     ' + l)
    print(f'  d_{n}: max |K_aa/prod_j S_aj(th,1) - 1| = {worst_diag:.1e}; max |offdiag| = {worst_off:.1e}; '
          f'max |K_1/K_1^book(B->0) - 1| = {worst_K1:.1e}')
