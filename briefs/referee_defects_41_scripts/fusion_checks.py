"""Independent bulk fusion and invariant-kernel restriction."""
import itertools
import numpy as np
from independent_ar import weights,rep,tensor,finite_hom,transitions,sector_T,rcheck,solve_sector


def fusion(h,w,t):
    q=-np.exp(-1j*np.pi*w);r=h-1
    ys=[np.exp(w*(t+1j*np.pi*(2*m-r+1)/h)) for m in range(r)]
    A=rep(h,q,ys[0]);ww=weights(h)
    for y in ys[1:]:
        A=tensor(A,rep(h,q,y));ww=np.array([u+v for u in ww for v in weights(h)])
    B=rep(h,q,np.exp(w*t),True)
    F,sv=finite_hom(B,A,weights(h,True),ww)
    phase=np.exp(-1j*np.pi*w/2);err=0
    states=list(itertools.product(range(h),repeat=r))
    for N in [np.zeros(h),np.arange(h)*.13-.13*(h-1)/2]:
      for i in range(h):
        d={(st,tuple(np.round(N,8))):F[k,i] for k,st in enumerate(states) if abs(F[k,i])>1e-12}
        for m in reversed(range(r)):
            out={}
            for (st,nn),val in d.items():
              for j,nn2,c in transitions(st[m],np.array(nn),q,phase*ys[m]):
                st2=list(st);st2[m]=j;key=(tuple(st2),tuple(np.round(nn2,8)))
                out[key]=out.get(key,0)+val*c
            d=out
        expected={}
        for j,nn,c in transitions(i,N,q,phase*np.exp(w*t),True):
          for k,st in enumerate(states):
            if abs(F[k,j])>1e-12: expected[st,tuple(np.round(nn,8))]=F[k,j]*c
        scale=max(map(abs,expected.values()))
        err=max(err,max(abs(d.get(k,0)-expected.get(k,0)) for k in d.keys()|expected.keys())/scale)
    return sv[-2:],err


def kernel_vector(Q,q,XP):
    M=sector_T(Q,q,XP);u=np.ones(len(Q),complex)
    for a in range(len(Q)-1): u[a+1]=-M[a,a]*u[a]/XP
    assert np.linalg.norm(M@u)<1e-10
    return u


def kernel_transmission(Q,q,XP,X,R):
    h=len(Q); ww=weights(h);u=kernel_vector(Q,q,XP)
    result=np.zeros((h,h),complex);err=0
    for i in range(h):
        raw=np.zeros((h,h),complex) # outgoing passing j, bound b
        for a,b,c in itertools.product(range(h),repeat=3):
          for j,nn,t in transitions(c,Q-ww[a],q,X):
            raw[j,b]+=u[a]*R[b,c,i,a]*t
        for j in range(h):
            unew=kernel_vector(Q+ww[i]-ww[j],q,XP)
            coef=np.vdot(unew,raw[j])/np.vdot(unew,unew)
            result[j,i]=coef
            err=max(err,np.linalg.norm(raw[j]-coef*unew))
    return result,err


def run():
    w=1.37;q=-np.exp(-1j*np.pi*w)
    for h in [3,4]: print('bulk fusion h,small singular values,error',h,fusion(h,w,.23))
    for h in [3,4]:
        ww=weights(h); t=.34
        tp=1j*np.pi/2-1j*(np.pi*(h-1)/h+np.pi/(h*w))
        yp=np.exp(w*tp); y=np.exp(w*t)
        phase=np.exp(-1j*np.pi*w/2);XP=phase*yp;X=phase*y
        R=rcheck(h,q,y,yp).reshape(h,h,h,h)
        for s in [np.zeros(h),ww[0]]:
            Q=s+ww[0];M,err=kernel_transmission(Q,q,XP,X,R)
            lam=M[0,0]/q**(-Q[0])
            diagerr=max(abs(M[i,i]/q**(-Q[i])-lam) for i in range(h))
            zeros=max(abs(M[j,i]) for j,i in itertools.product(range(h),repeat=2) if j not in [i,(i-1)%h])
            loop=1;Qwalk=Q.copy()
            for i in range(h):
                Mi,ei=kernel_transmission(Qwalk,q,XP,X,R)
                loop*=Mi[(i-1)%h,i]/(lam*X)
                Qwalk+=ww[i]-ww[(i-1)%h];err=max(err,ei)
            print('kernel h,s,invariance,diag,forbidden,cycle,q^2 error',h,s,err,diagerr,zeros,loop,abs(loop-q*q))
    # Exact Cartan conservation fails at roots of unity: an extra lattice shift commutes.
    h=3;qroot=np.exp(2j*np.pi/5);v=5*(np.eye(h)[0]-np.eye(h)[1])
    N=np.array([1.,-2.,1.]);mu=.71
    kerr=max(abs(qroot**(N[a]-N[(a+1)%h])-qroot**((N+v)[a]-(N+v)[(a+1)%h])) for a in range(h))
    ferr=max(abs(mu*qroot**N[(a+1)%h]-mu*qroot**(N+v)[(a+1)%h]) for a in range(h))
    print('root-of-unity commuting S_(5alpha1): K,F errors',kerr,ferr,'nonzero label displacement',v)
    # Minus-orientation linear matrix: display fitted diagonal and off-diagonal patterns.
    q=-np.exp(-1j*np.pi*w);mu=-np.exp(-1j*np.pi*w/2)/(q-1/q)
    Q=np.array([.2,-.05,-.15]);y=1.23
    M,sv=solve_sector(Q,q,y,mu,bar=True,orient=-1)
    print('orientation minus Vbar normalized matrix',M/M[0,0])


if __name__=='__main__':run()
