"""a_3 middle-soliton endpoint d -> 1 approached along non-solution paths (species-1 pair admixture)."""
import sys
import mpmath as mp
from a5lib import *
src = open('/home/gustav/Git/affine_toda/code/particle_reflection.py').read()
PR = {}; exec(src.split('\n# ---------------------------------------------------------------- 1. blocks')[0], PR)
mp.mp.dps = 60
n = 3; T = Theory(n); h = 4
th = mp.mpf('0.6'); Om = T.mass(1)*mp.cosh(th)
DG = (PR['K_uniform'](1, 3)*PR['M_sol'](1, 2, 3))(th, 0)
print('DG b=2 K_1 =', mp.nstr(DG, 10))
xs = [mp.mpf(0)] + [-mp.mpf(k)/40 for k in range(1, 400)]
def regular(bg): return min(bg.tau[j](x).real for j in range(h) for x in xs) > 0
for C in [(1, -1, 1, -1), (1, -1, -1, -1), (-1, -1, -1, -1)]:
    for etaf, chi in (('1', '0.3'), ('1', '2.0'), ('0.5', '1.1'), ('2', '4.0')):
        out = []
        for e in ('1e-2', '1e-4', '1e-6'):
            ee = mp.mpf(e); eta = ee**mp.mpf(etaf) if etaf != '0.5' else mp.sqrt(ee)
            bg = Background(T, [(2, 1 - ee), (1, eta*mp.exp(I*mp.mpf(chi))), (3, eta*mp.exp(-I*mp.mpf(chi)))])
            if not regular(bg):
                out.append('irregular'); continue
            opn, K, Kf, res = kmatrix(bg, C, Om)
            out.append(f'|K11-DG|={mp.nstr(abs(K[(1,1)]-DG),3)} |K31|={mp.nstr(abs(K[(3,1)]),3)}')
        print(''.join('+' if c > 0 else '-' for c in C), f'd=1-e, species-1 pair amplitude e^{etaf} chi={chi}:', ' | '.join(out), flush=True)
