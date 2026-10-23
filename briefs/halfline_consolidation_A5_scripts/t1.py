import sys, io, contextlib
sys.path.insert(0, '/home/gustav/Git/affine_toda/code')
import mpmath as mp
from a5lib import *
src = open('/home/gustav/Git/affine_toda/code/particle_reflection.py').read()
PR = {}; exec(src.split('\n# ---------------------------------------------------------------- 1. blocks')[0], PR)
mp.mp.dps = 30
def Kdg(c, b, n, th): return (PR['K_uniform'](c, n)*PR['M_sol'](c, b, n))(th, 0)
def blk(x, th, h): return PR['blk'](x, th, h)
for n in (2,3):
    T = Theory(n); h = n+1
    bg = Background(T, [])
    Om = mp.mpf('2.1')
    for C in (1, -1):
        opn, K, Kf, res = kmatrix(bg, [C]*h, Om)
        print(n, C, opn, mp.nstr(res,3))
        for c in opn:
            th = mp.acosh(Om/T.mass(c))
            print('  ', c, mp.nstr(K[(c,c)],12), mp.nstr(-(blk(c,th,h)*blk(h-c,th,h))**(-C),12), mp.nstr(max(abs(K[(c2,c)]) for c2 in opn if c2!=c) if len(opn)>1 else 0,3))
# pair a_2
n=2; T=Theory(n); h=3; a=1; s=mp.sin(mp.pi*a/h); d=1/(1+s)
for chi in (mp.mpf('0.3'), mp.mpf('2.5'), mp.mpf('4.6')):
    bg = Background(T, [(a, d*mp.exp(I*chi)), (h-a, d*mp.exp(-I*chi))])
    # determine pattern from background BC
    t = bg.taus(mp.mpf(0)); eh=[mp.sqrt(v) for v in bg.psi(0)]
    # phi' = -sum alpha_j tau_j'/tau_j ; BC phi' = -sum C_j alpha_j e^{psi_j/2}
    dphi = [t[(k-1)%h][1]/t[(k-1)%h][0] - t[k][1]/t[k][0] for k in range(h)]
    import itertools
    good=[C for C in itertools.product([1,-1],repeat=h) if max(abs(dphi[k] + sum(C[j]*eh[j]*((1 if k==j else 0)-(1 if k==(j+1)%h else 0)) for j in range(h))) for k in range(h))<1e-20]
    print('chi',chi,'patterns',good, [mp.nstr(v[0].real,5) for v in t])
    for Om in (mp.cosh(mp.mpf('0.6'))*T.mass(1),):
        opn,K,Kf,res = kmatrix(bg, good[0], Om)
        th = mp.mpf('0.6')
        print(' res',mp.nstr(res,3),{k:mp.nstr(v,10) for k,v in K.items()}, 'DG', mp.nstr(Kdg(1,1,2,th),10), 'unit', mp.nstr(unitarity_defect(opn,Kf),3))
        print(' pde', mp.nstr(pde_residual(bg, bg.wave(1, I*mp.sqrt(Om**2-3)), mp.mpf('-0.4'), Om),3))
