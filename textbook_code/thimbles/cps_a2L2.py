import numpy as np, time, sys
from lie import An
from model import Chain
from newton_batch import find_cps
g = An(2); lam = g.fundamental_weights_rep(1)[(1,)]
kappa = float(sys.argv[1]) if len(sys.argv)>1 else 1.0
L = int(sys.argv[2]) if len(sys.argv)>2 else 2
ch = Chain(g, L, kappa, lam)
t=time.time()
cps, cnt = find_cps(ch, nstart=40000, re_spread=4.0, im_spread=1.5)
print("found", len(cps), "in %.1fs"%(time.time()-t))
S = np.array([ch.action(u) for u in cps])
o = np.argsort(S.real); cps=[cps[i] for i in o]; cnt=[cnt[i] for i in o]; S=S[o]
np.save('cps_a2L%d_k%.2f.npy'%(L,kappa), np.array(cps))
print("Imax", ch.Imax, "Kmin", ch.Kmin)
for s,c in zip(S,cnt):
    flag = '*' if abs(s.imag) <= ch.Imax else ' '
    if s.real < 40: print("%s %10.5f %+10.5f  hits=%d"%(flag, s.real, s.imag, c))
