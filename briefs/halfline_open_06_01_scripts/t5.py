from linrefl import *
from bp_table import class_b
import itertools
def blk(x,th,h): return np.sinh(th/2+1j*np.pi*x/(2*h))/np.sinh(th/2-1j*np.pi*x/(2*h))
def Kdg(a,b,h,th): return -1/(blk(a,th,h)*blk(h-a,th,h))*blk(a+b-h/2,th,h)*blk(-a-b-h/2,th,h)*blk(a-b+h/2,th,h)*blk(b-a+h/2,th,h)
def pattern_of(T,cs):
    n=T.n; h=n+1; al=roots(n)
    t,d,dd=T.tau(cs,np.array([0.0])); t0=t[:,0].real; d0=d[:,0].real
    phi_dot=-(al.T@(d0/t0)); eh=np.sqrt(np.roll(t0,1)*np.roll(t0,-1)+0j)/t0
    res=[]
    for C in itertools.product([1,-1],repeat=h):
        r=phi_dot+al.T@(np.array(C)*eh).real
        res.append((np.abs(r).max(),C))
    res.sort(); return res[0]
if __name__=="__main__":
    n=2; h=3; T=Tau(n)
    for a in (1,):
        s=np.sin(np.pi*a/h); d=1/(1+s)
        for chi in (0.3,1.0,1.5,2.0,3.0):
            p=np.array([np.log(d),chi]); cs=cs_from_p(n,p)
            tf=lambda x: T.tau(cs,np.array([x]))[0][:,0].real
            t=T.tau(cs,np.linspace(-20,0,4000))[0].real
            r,C=pattern_of(T,cs)
            th=0.6; Om=2*np.sin(np.pi/h)*np.cosh(th)
            K,lab=refl_matrix(n,tf,C,Om)
            b=class_b(C) if np.prod(C)==1 else None
            print('chi',chi,'minTau %.3f'%t.min(),'pattern',C,'bcres %.1e'%r,'K',np.round(K,5).tolist(),'DG',np.round(Kdg(1,b,h,th),5) if b else '')
