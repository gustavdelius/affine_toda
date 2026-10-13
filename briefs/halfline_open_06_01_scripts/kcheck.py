"""Linearized reflection factor about the regular solitonic vacua (continued stationary pair) vs DG (5.7) classical limit."""
from linrefl import *
from tauvac import nparams
from t5 import blk, Kdg, pattern_of
from bp_table import class_b
for n,a in [(2,1),(3,1),(4,1),(4,2),(5,1),(5,2)]:
    h=n+1; T=Tau(n); s=np.sin(np.pi*a/h); d=1/(1+s)
    sing=sorted(((np.pi-2*np.pi*j*a/h)%(2*np.pi)) for j in range(h))
    gaps=[(sing[k],sing[(k+1)%h]+(2*np.pi if k+1==h else 0)) for k in range(h)]
    seen={}
    for lo,hi in gaps:
        if hi-lo<1e-9: continue
        chi=(lo+hi)/2
        p=np.zeros(nparams(n)); p[0::2][:(h-1)//2]=-60.0; p[2*a-2]=np.log(d); p[2*a-1]=chi
        if h%2==0: p[-1]=0.0
        cs=cs_from_p(n,p)
        tf=lambda x: T.tau(cs,np.array([x]))[0][:,0].real
        r,C=pattern_of(T,cs)
        if C in seen: continue
        b=class_b(C); th=0.6; Om=2*np.sin(np.pi/h)*np.cosh(th)
        K,lab=refl_matrix(n,tf,C,Om,rtol=1e-12,atol=1e-14)
        kd=Kdg(1,b,h,th)
        seen[C]=1
        print('a_%d pair a=%d chi=%.3f pattern %s bc %.0e b=%d  |K11-Kdg|=%.1e |K22-Kdg|=%.1e max|offdiag|=%.1e'%(n,a,chi,''.join('+' if c>0 else '-' for c in C),r,b,abs(K[0,0]-kd),abs(K[1,1]-kd),max(abs(K[0,1]),abs(K[1,0]))),flush=True)
