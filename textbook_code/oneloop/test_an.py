import numpy as np
from toda import *
pi=np.pi
r=4; h=5
sp, C, n, autos = species('a', r)
for m, v in sp:
    print(round(m,6), np.round(v/v[0],4))
# soliton of first species
z = 0.3+0.2j
for ia,(ma,va) in enumerate(sp):
    c = soliton(C, n, ma, va/va[0])
    print("species mass", round(ma,5), "degrees", degree(c))
    for ib,(mb,vb) in enumerate(sp):
        X, res, D, Dp, p = Xfactor(C, n, c, ma, mb, vb, z)
        # identify a, b: m = 2 sin(pi a/h)
        print("   channel", round(mb,5), "X =", np.round(X,8), "resid", f"{res:.1e}", "deg P", Dp)
