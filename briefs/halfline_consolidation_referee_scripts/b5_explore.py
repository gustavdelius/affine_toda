import numpy as np
from b5_boundary_krein import broken
def thr(mus, lw, rw, pairs):
    f=lambda om: max(broken(mus, om, lw, rw, list(p)) for p in pairs)
    lo,hi=0.3,1.0
    if f(lo)>1e-8: return None
    for _ in range(25):
        mid=0.5*(lo+hi)
        if f(mid)>1e-8: hi=mid
        else: lo=mid
    return 0.5*(lo+hi)
big=(10.0,12.0); small=(1e-4,1e-3,1e-2)
for lw,rw in (('U','N'),('N','N')):
    a=thr(('V','W'),lw,rw,[(t1,t2) for t1 in big for t2 in small])
    b=thr(('V','W'),lw,rw,[(t2,t1) for t1 in big for t2 in small])
    c=thr(('W','V'),lw,rw,[(t1,t2) for t1 in big for t2 in small])
    print(f"V Vbar {lw}{rw}: soliton large/antisoliton small {a:.4f}; soliton small/antisoliton large {b:.4f}; (W,V) order large/small {c:.4f}")
