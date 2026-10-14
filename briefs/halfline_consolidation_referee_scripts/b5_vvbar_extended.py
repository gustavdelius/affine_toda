"""V Vbar: below the thresholds found by b5_vvbar_threshold.py (0.6667 UN, 0.5001 NN), is the sector zeta=1
still unbroken when both rapidities range up to 20?  And where does the breaking sit just above?"""
import numpy as np
from b5_boundary_krein import broken
A=[2.0,4.0,6.0,8.0,12.0,16.0,20.0]; Bs=np.linspace(0.5,20,40)
for lw,rw,below,above in (('U','N',0.660,0.672),('N','N',0.495,0.505)):
    for om in (below,above):
        m=(0,None)
        for a in A:
            for b in Bs:
                if abs(a-b)<1e-9: continue
                v=broken(('V','W'),om,lw,rw,[a,b])
                if v>m[0]: m=(v,(a,b))
        print(f"V Vbar {lw}{rw} omega={om}: max ||Lambda|-1| = {m[0]:.2e} at {m[1]}",flush=True)
