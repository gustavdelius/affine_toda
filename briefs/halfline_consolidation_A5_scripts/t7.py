"""(i) a_5 middle-soliton endpoint with species-1 / species-2 pair admixture; (ii) a_2 pair endpoint chi = pi along
non-solution paths d != 1/(1+s)."""
import mpmath as mp
from a5lib import *
src = open('/home/gustav/Git/affine_toda/textbook_code/particle_reflection.py').read()
PR = {}; exec(src.split('\n# ---------------------------------------------------------------- 1. blocks')[0], PR)
mp.mp.dps = 60
th = mp.mpf('0.6')
xs = [mp.mpf(0)] + [-mp.mpf(k)/40 for k in range(1, 400)]
def regular(bg): h = bg.th.h; return min(bg.tau[j](x).real for j in range(h) for x in xs) > 0
n = 5; T = Theory(n); h = 6; Om = T.mass(1)*mp.cosh(th)
DG = (PR['K_uniform'](1, 5)*PR['M_sol'](1, 3, 5))(th, 0)
print('a_5 DG b=3 K_1 =', mp.nstr(DG, 10))
for C in [(1, -1, 1, -1, -1, -1), (1, -1, -1, -1, -1, -1)]:
    for sp, etaf, chi in ((1, 2, '0.3'), (1, 1, '0.3'), (2, 1, '0.3'), (2, 2, '1.0'), (1, 1, '3.5'), (2, 1, '3.5')):
        out = []
        for e in ('1e-2', '1e-4', '1e-6'):
            ee = mp.mpf(e); eta = ee**etaf
            bg = Background(T, [(3, 1 - ee), (sp, eta*mp.exp(I*mp.mpf(chi))), (h - sp, eta*mp.exp(-I*mp.mpf(chi)))])
            if not regular(bg):
                out.append('irregular'); continue
            opn, K, Kf, res = kmatrix(bg, C, Om)
            out.append(f'|K11-DG|={mp.nstr(abs(K[(1,1)]-DG),3)} |K51|={mp.nstr(abs(K[(5,1)]),3)}')
        print(''.join('+' if c > 0 else '-' for c in C), f'species-{sp} pair amp e^{etaf} chi={chi}:', ' | '.join(out), flush=True)
n = 2; T = Theory(n); h = 3; Om = T.mass(1)*mp.cosh(th); s = mp.sin(mp.pi/3); d0 = 1/(1+s)
for C in [(-1, 1, 1), (-1, -1, -1), (-1, 1, -1)]:
    for dfac, cfac in (('0', '1'), ('0', '-1'), ('-1', '0'), ('-1', '2')):
        out = []
        for e in ('1e-3', '1e-5', '1e-7'):
            ee = mp.mpf(e)
            bg = Background(T, [(1, d0*(1 + mp.mpf(dfac)*ee)*mp.exp(I*(mp.pi + mp.mpf(cfac)*ee))), (2, d0*(1 + mp.mpf(dfac)*ee)*mp.exp(-I*(mp.pi + mp.mpf(cfac)*ee)))])
            if not regular(bg):
                out.append('irregular'); continue
            opn, K, Kf, res = kmatrix(bg, C, Om)
            out.append(f'K11={mp.nstr(K[(1,1)],7)} |K21|={mp.nstr(abs(K[(2,1)]),6)}')
        print('a_2', ''.join('+' if c > 0 else '-' for c in C), f'd*(1+{dfac}e), chi=pi+{cfac}e:', ' | '.join(out), flush=True)
