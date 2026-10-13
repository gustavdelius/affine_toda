from redmin2 import *
from scipy.optimize import brentq
n=4; T=Tau(n); P=(1,1,-1,1,-1); s1=0; s2=np.pi; k=4
chi1=(s2-2*s1+2*np.pi*k)/5; chi2=2*chi1+s1
def cs_of(l1,l2): return cs_from_p(n,np.array([l1,chi1,l2,chi2]))
def t2(l1,l2): return T.tau(cs_of(l1,l2),np.array([0.0]))[0][2,0].real
for l1 in np.linspace(-3,-0.2,8):
    ls=np.linspace(-3,1,400); v=[t2(l1,l) for l in ls]
    roots_=[brentq(lambda l:t2(l1,l),ls[i],ls[i+1]) for i in range(len(ls)-1) if v[i]*v[i+1]<0]
    out=[]
    for r in roots_:
        for eps in (1e-5,):
            # step slightly inside: find l2 with tau2(0)=eps on the positive side
            try: l2=brentq(lambda l:t2(l1,l)-eps,r-0.01,r+0.01)
            except Exception: continue
            cs=cs_of(l1,l2); m=tau_min(T,cs).min()
            if m>0:
                F,t0,eh=Fval(T,cs,P); out.append((round(r,4),round(F,4),np.round(t0,6).tolist(),'minτ %.1e'%m))
    print(round(l1,3),out)
from scipy.optimize import minimize_scalar
def Fcurve(l1,eps=1e-7):
    ls=np.linspace(-1.2,0.2,300); v=[t2(l1,l)-eps for l in ls]
    r=[brentq(lambda l:t2(l1,l)-eps,ls[i],ls[i+1]) for i in range(len(ls)-1) if v[i]*v[i+1]<0]
    cs=cs_of(l1,r[0]); return Fval(T,cs,P)[0], r[0]
res=minimize_scalar(lambda l:Fcurve(l)[0],bounds=(-2.0,-1.0),method='bounded',options=dict(xatol=1e-9))
l1=res.x; F,l2=Fcurve(l1); print('min F on singular curve: %.6f at ln d1=%.5f ln d2=%.5f (d1=%.5f d2=%.5f)'%(F,l1,l2,np.exp(l1),np.exp(l2)),' BP -10.41616')
