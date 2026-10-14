import os
import numpy as np, mpmath as mp, os
TC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'textbook_code')
cwd = os.getcwd(); os.chdir(TC); ns = {}
exec(compile(open('smatrix_simply_laced.py').read().split('\nB = 0.37')[0], 'smatrix_simply_laced.py', 'exec'), ns); os.chdir(cwd)
mp.mp.dps = 30
blk = lambda x, h, th: mp.sinh(th/2 + 1j*mp.pi*x/(2*h))/mp.sinh(th/2 - 1j*mp.pi*x/(2*h))
n = 6
h, bl, mass, _ = ns['setup']('d', n)
scale = 2*np.sqrt(2)*np.sin(np.pi/h)/mass[0]
ms = mass*scale
m2, m4 = ms[1], ms[3]
print('d6 masses', np.round(ms, 6), 'm2', m2, 'm4', m4)
def prodS(a, th):
    r = mp.mpc(1)
    for b in range(n):
        for x, k in bl[(a, b)].items():
            r *= (blk(x-1, h, th)*blk(x+1, h, th)/(blk(x, h, th)**2))**k
    return r
for a in range(n):
    out = []
    for z in range(1, h//2 + 1):
        t0 = 1j*mp.pi*z/h
        v = [abs(prodS(a, t0 + e)) for e in (mp.mpf('1e-6'), mp.mpf('1e-8'))]
        o = int(round(float(mp.log(v[1]/v[0])/mp.log(100))))
        if o > 0:
            p = ms[a]*np.sin(np.pi*z/h)
            out.append(f'z={z} order {o} p/(i) = {p:.6f} (2p: {2*p:.6f})')
    print('particle', a+1, out)
