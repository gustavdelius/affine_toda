import numpy as np
from math import pi,sqrt
from scipy.linalg import lstsq
om=sqrt(2);q=-np.exp(-1j*pi*om);w=np.exp(-1j*pi*om/2);ep=1/(2*np.cos(pi*om/2));T=1.5*om
f=lambda x:-np.cos((om+1)*np.arccos(x))/np.sin(pi*om)
print('17 intercept',f(0),'epsilon*',ep,'slope finite difference',(f(1e-6)-f(-1e-6))/2e-6,'predicted',-(om+1)/(2*np.sin(pi*om/2)))

def rep(y,conj=False):
    out=[]
    for a,b in [(0,1),(1,2),(2,0)]:
        E=np.zeros((3,3),complex);E[a,b]=1
        k=np.eye(3,dtype=complex);k[a,a]=q;k[b,b]=1/q
        if conj:out.append((y*E.T,E/y,np.linalg.inv(k)))
        else:out.append((y*E,E.T/y,k))
    return out

def cop(A,B):
    return [x for (e,f,k),(ee,ff,kk) in zip(A,B) for x in [np.kron(e,np.eye(3))+np.kron(k,ee),np.kron(f,np.linalg.inv(kk))+np.kron(np.eye(3),ff),np.kron(k,kk)]]
def braid(z):
    dom=cop(rep(z),rep(1));tar=cop(rep(1),rep(z))
    basis=[]
    for a in range(3):
        for b in range(3):
            for c,d in set([(a,b),(b,a)]):
                M=np.zeros((9,9),complex);M[3*c+d,3*a+b]=1;basis.append(M)
    mat=np.vstack([np.column_stack([(M@D-C@M).ravel() for M in basis]) for D,C in zip(dom,tar)])
    normal=np.array([M[0,0] for M in basis]);L=np.vstack([mat,normal]);rhs=np.zeros(L.shape[0],complex);rhs[-1]=1
    sol=lstsq(L,rhs)[0];R=sum(x*M for x,M in zip(sol,basis))
    assert np.linalg.norm(L@sol-rhs)<1e-12
    return R

def K(theta):
    logX=om*(theta-1j*pi/2);p=lambda x:np.exp(x*logX)
    k=(w*p(-1.5)-p(1.5)/w)/(ep*(1/q-q))
    return np.array([[k,p(-.5),p(.5)],[p(.5),k,p(-.5)],[p(-.5),p(.5),k]])
mu=.29-.08j; mua=mu-T/2+.5;mub=mu+T/2-.5
ta=1j*pi*mua/T;tb=1j*pi*mub/T;ya=np.exp(om*ta);yb=np.exp(om*tb)
r=-q*yb/ya;v=np.zeros(9,complex);v[[0,4,8]]=[1,r,r*r]
ops=cop(rep(ya),rep(yb,True))
singerr=max(np.linalg.norm(O@v-(v if j%3==2 else 0)) for j,O in enumerate(ops))
out=np.kron(np.eye(3),K(ta))@braid(ya*yb)@np.kron(np.eye(3),K(tb).T)@v
c=-np.cos(pi*mu)*np.cos(pi*(mu+om/2))*np.cos(pi*(mu-T))/(np.sin(pi*om/2)**2*np.cos(pi*(mu-om/2)))
print('18 affine singlet residual',singerr,'matrix fusion residual',np.linalg.norm(out-c*v)/np.linalg.norm(v))
assert singerr<1e-12 and np.linalg.norm(out-c*v)<1e-11
print('PASS 17-18: SG Taylor coefficients and independently solved a2 matrix fusion')
