"""Find regular phases for the three-soliton images (helper for a2_half_coweights.py):
minimum over a space-time grid of min_j |tau_j| / max(1, |largest term|)."""
import numpy as np, itertools
from a2_half_coweights import taus, Acoef
def minrel(n,a,b,th,xi1,im3):
    h=n+1
    A_ab=Acoef(a,b,th,h)
    sols=[(a,th,xi1),(a,-th,-np.conj(xi1)-np.log(Acoef(a,a,2*th,h))-np.log(A_ab)),(b,0.0,-np.log(A_ab)+1j*im3)]
    m=np.inf
    for x in np.linspace(-6,0,61):
        for t in np.linspace(-4,4,41):
            T=taus(sols,x,t,h)
            # scale: sum of moduli of the terms ~ use |tau| relative to 1+|E| product bound
            sc=np.max(np.abs(taus([(s[0],s[1],np.real(s[2])) for s in sols],x,t,h)))
            m=min(m,np.min(np.abs(T))/sc)
    return m
for (n,a,b,th) in [(3,1,2,0.6),(4,2,1,0.9),(3,1,1,0.6),(2,1,1,0.7)]:
    best=None
    for im1 in np.linspace(-3,3,13):
        for im3 in np.linspace(-3,3,13):
            v=minrel(n,a,b,th,0.2+1j*im1,im3)
            if best is None or v>best[0]: best=(v,im1,im3)
    print(n,a,b,th,"best min|tau|/scale = %.3f at Im xi1=%.2f Im xi3=%.2f"%best)
