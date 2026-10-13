import numpy as np
from nk_common import k_book, fast_rcheck
from boundary_common import q_phys, an_vector, k_intertwiner
def prop(A,B):
    c=np.vdot(B.ravel(),A.ravel())/np.vdot(B.ravel(),B.ravel()); return np.linalg.norm(A-c*B)/np.linalg.norm(A)
def Phif(n,s,om,th):
    T=(n+1)*om/2; a=1j*np.pi*om/2
    return -np.sinh(T*th-a)/np.sinh(T*th+a) if s==1 else np.cosh(T*th-a)/np.cosh(T*th+a)
def Uf(n,s,om,th):
    T=(n+1)*om/2
    return (s*np.cosh(2*T*th)-np.cos(np.pi*T))/(2*np.sin(np.pi*om/2)**2)
def Phinum(n,s,om,th):
    q=q_phys(om); N=n+1
    km=k_book(n,om,1j*np.pi/2-th,s).T.reshape(N*N); kp=k_book(n,om,1j*np.pi/2+th,s).T.reshape(N*N)
    R2=fast_rcheck(an_vector(n,q,np.exp(2*om*th),True),an_vector(n,q,1.0,True))
    v=R2@kp; c=np.vdot(v,km)/np.vdot(v,v); return c,np.linalg.norm(km-c*v)/np.linalg.norm(km)
if __name__=='__main__':
  wi=wP=wU=0
  for n in (1,2,3,4,5):
    N=n+1
    for s in (1,-1):
      for om in (0.37,1.83,2.6,4.3):
        q=q_phys(om); w=np.exp(-0.5j*np.pi*om); es=1/(w+1/w)
        e=es*np.ones(N); e[0]*=s
        for th in (0.31+0.2j,-0.7+0.1j):
          y=np.exp(om*th)
          null,_=k_intertwiner(an_vector(n,q,y),an_vector(n,q,1/y,True),[1/q]*N,e)
          wi=max(wi,prop(null[0],k_book(n,om,th,s)) if len(null)==1 else 9)
          p,r=Phinum(n,s,om,th); wP=max(wP,abs(p/Phif(n,s,om,th)-1),r)
          P=k_book(n,om,-th,s).T@k_book(n,om,th,s); wU=max(wU,np.abs(P-Uf(n,s,om,th)*np.eye(N)).max()/abs(Uf(n,s,om,th)))
  print('intertwiner',wi,'Phi',wP,'U',wU)
