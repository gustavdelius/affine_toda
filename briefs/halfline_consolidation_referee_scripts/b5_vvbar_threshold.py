"""V Vbar two-soliton thresholds: minimal omega at which a level leaves the unit circle, searching
theta_1 in {8, 10, 12} (either particle fast) and the slow rapidity on a fine grid in [1e-4, 4]."""
import numpy as np
from b5_boundary_krein import broken
slow=np.concatenate([[1e-4,1e-2],np.linspace(0.2,4.0,39)])
def worst(om,lw,rw):
    m=(0,None)
    for t1 in (8.0,10.0,12.0):
        for t2 in slow:
            for p in ((t1,t2),(t2,t1)):
                b=broken(('V','W'),om,lw,rw,list(p))
                if b>m[0]: m=(b,p)
    return m
for lw,rw,lo,hi in (('U','N',0.60,0.70),('N','N',0.45,0.54)):
    assert worst(lo,lw,rw)[0]<1e-8
    for _ in range(14):
        mid=0.5*(lo+hi)
        if worst(mid,lw,rw)[0]>1e-8: hi=mid
        else: lo=mid
    print(f"V Vbar {lw}{rw}: threshold omega = {0.5*(lo+hi):.4f}; first breaking at (theta_1, theta_2) = {worst(hi,lw,rw)[1]}")
