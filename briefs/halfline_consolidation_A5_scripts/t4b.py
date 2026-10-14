"""a_4 exceptional (BP 'overlapping charge sectors') background: locate it, then K with target Hessian along
several approach paths of regular static Hirota configurations."""
import sys
import mpmath as mp
from a5lib import *
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
n = 4; T = Theory(n); h = 5
P = (1, 1, -1, 1, -1)
def bg_of(d1, d2, chi1=9*mp.pi/5, chi2=18*mp.pi/5):
    c1 = d1*mp.exp(I*chi1); c2 = d2*mp.exp(I*chi2)
    return Background(T, [(1, c1), (4, mp.conj(c1)), (2, c2), (3, mp.conj(c2))])
def d2_of(d1, guess=mp.mpf('0.5637')):
    return mp.findroot(lambda d2: bg_of(d1, d2).tau[2](0).real, guess)
def Rfin(d1, eta):
    d2 = d2_of(d1); bg = bg_of(d1, d2); x = -eta
    t = [(bg.tau[j](x), bg.tau[j](x, 1)) for j in range(h)]
    eh = [mp.sqrt(v) for v in bg.psi(x)]
    dphi = [t[(k-1) % h][1]/t[(k-1) % h][0] - t[k][1]/t[k][0] for k in range(h)]
    return (dphi[0] + sum(P[j]*eh[j]*((0 == j) - (0 == (j+1) % h)) for j in range(h))).real
eta = mp.mpf(10)**(-mp.mp.dps//2)
d1s = mp.findroot(lambda d1: Rfin(d1, eta), mp.mpf('0.2579'))
d2s = d2_of(d1s)
print('exceptional: d1 =', mp.nstr(d1s, 15), 'd2 =', mp.nstr(d2s, 15))
bgs = bg_of(d1s, d2s)
print('tau_j(0):', [mp.nstr(bgs.tau[j](0).real, 5) for j in range(h)])
# regular? min tau_j on x <= 0 for the exact background (x<0)
xs = [-mp.mpf(k)/50 for k in range(1, 600)]
print('min tau_j(x<0):', mp.nstr(min(bgs.tau[j](x).real for j in range(h) for x in xs), 5))
th = mp.mpf('0.6'); Om = T.mass(1)*mp.cosh(th)
def regular(bg):
    return min(bg.tau[j](x).real for j in range(h) for x in [mp.mpf(0)] + xs) > 0
paths = {
  'd2 up':   lambda e: bg_of(d1s, d2s*(1+e)),
  'd2 down': lambda e: bg_of(d1s, d2s*(1-e)),
  'd1 up':   lambda e: bg_of(d1s*(1+e), d2s),
  'd1 down': lambda e: bg_of(d1s*(1-e), d2s),
  'chi1 +':  lambda e: bg_of(d1s, d2s, chi1=9*mp.pi/5 + e),
  'chi1 -':  lambda e: bg_of(d1s, d2s, chi1=9*mp.pi/5 - e),
  'chi2 +':  lambda e: bg_of(d1s, d2s, chi2=18*mp.pi/5 + e),
}
for name, f in paths.items():
    for e in ('1e-2', '1e-3', '1e-4', '1e-6'):
        bg = f(mp.mpf(e))
        if not regular(bg):
            print(f'{name:8s} e={e}: not regular on x<=0'); continue
        opn, K, Kf, res = kmatrix(bg, P, Om)
        print(f'{name:8s} e={e}: tau2(0)={mp.nstr(bg.tau[2](0).real,3)} tau4(0)={mp.nstr(bg.tau[4](0).real,3)} K11={mp.nstr(K[(1,1)],8)} K44={mp.nstr(K[(4,4)],8)} |K41|={mp.nstr(abs(K[(4,1)]),6)} |K14|={mp.nstr(abs(K[(1,4)]),6)} unit={mp.nstr(unitarity_defect(opn,Kf),2)}', flush=True)
