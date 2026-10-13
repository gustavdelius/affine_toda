"""a_5 species-2 pair endpoint: other approach paths (non-solution configurations: d != 1/(1+s))."""
import sys
import mpmath as mp
from a5lib import *
mp.mp.dps = 70
n = 5; T = Theory(n); h = 6; s = mp.sin(mp.pi/3); d0 = 1/(1+s); chi_end = 5*mp.pi/3
Om = T.mass(1)*mp.cosh(mp.mpf('0.6'))
xs = [mp.mpf(0)] + [-mp.mpf(k)/40 for k in range(1, 500)]
def bg_of(d, chi): return Background(T, [(2, d*mp.exp(I*chi)), (4, d*mp.exp(-I*chi))])
def regular(bg): return min(bg.tau[j](x).real for j in range(h) for x in xs) > 0
S = lambda C: ''.join('+' if x > 0 else '-' for x in C)
CS = [(1, 1, -1, -1, -1, -1), (1, 1, -1, -1, 1, -1)] if len(sys.argv) < 2 else [(1, -1, -1, 1, -1, -1), (-1, 1, -1, -1, 1, -1), (-1, -1, -1, -1, -1, -1)]
for C in CS:
    for dfac, cfac in (('-1', '0'), ('1', '0'), ('-1', '1'), ('-1', '-1'), ('-1', '3'), ('-1', '-3'), ('0', '1'), ('0', '-1')):
        out = []
        for e in ('1e-3', '1e-5', '1e-7'):
            ee = mp.mpf(e)
            bg = bg_of(d0*(1 + mp.mpf(dfac)*ee), chi_end + mp.mpf(cfac)*ee)
            if not regular(bg):
                out.append('irregular'); continue
            opn, K, Kf, res = kmatrix(bg, C, Om)
            out.append(f'K11={mp.nstr(K[(1,1)],7)} |K51|={mp.nstr(abs(K[(5,1)]),6)} t2/t5={mp.nstr(bg.tau[2](0).real/bg.tau[5](0).real,3)}')
        print(S(C), f'd*(1+{dfac}e), chi_end+{cfac}e :', ' | '.join(out), flush=True)
