import numpy as np
from cmath import sinh,sin
from math import pi,sqrt
# 7: continuation, with modulo 2h equivalence retained.
t=.63+.11j; B=.47;h=3;n=h-1;om=-2/B
blk=lambda x:sinh(t/2+1j*pi*x/(2*h))/sinh(t/2-1j*pi*x/(2*h))
labels=[2*(-n*om/2)/om,2*(.5-h*om/2)/om,2*(-om/2-.5)/om]
a=np.prod([blk(x) for x in labels]);b=blk(h-B/2)/(blk(h-1)*blk(1-B/2))
print('7 labels',labels,'continuation error',abs(a-b));assert abs(a-b)<1e-12
print('PASS 7: continuation correct using block period 2h')
# 8: independent Uq(a2) reps and coproduct; no author helpers.
q=-np.exp(-1j*pi*sqrt(2));d=3
E=lambda a,b:np.eye(3)[:,a:a+1]@np.eye(3)[b:b+1,:]
base=[]
for a,b in [(0,1),(1,2)]:
    k=np.eye(3,dtype=complex);k[a,a]=q;k[b,b]=1/q
    base.append((E(a,b),E(b,a),k))
reps=base
for size in [3,9,27]:
    bi=[];cs=[.3+.2j,.8-.4j];eps=[.47,.91]
    for (e,f,k),c,z in zip(reps,cs,eps):bi.append(e+c*f@k+z*k)
    X,Y=bi;ki=reps[0][2];kj=reps[1][2]
    lhs=X@X@Y-(q+1/q)*X@Y@X+Y@X@X
    rhs=q*cs[0]*Y+eps[1]*((1-q)*(1-1/q)*eps[0]**2-q*cs[0])*ki@ki@kj
    err=np.max(abs(lhs-rhs));print('8 dim',size,'Serre residual',err);assert err<1e-10
    if size<27:
        reps=[(np.kron(e,np.eye(3))+np.kron(k,ee),np.kron(f,np.linalg.inv(kk))+np.kron(np.eye(len(k)),ff),np.kron(k,kk)) for (e,f,k),(ee,ff,kk) in zip(reps,base)]
print('PASS 8: deformed Serre identity on fundamental and independent tensor modules')
# 9 fractional powers defined with the continuous log omega(theta-i*pi/2).
om=sqrt(2);w=np.exp(-1j*pi*om/2);eps=1/(2*np.cos(pi*om/2));th=.43
for s in [1,-1]:
    def km(th):
        z=om*(th-1j*pi/2); X=lambda p:np.exp(p*z)
        k=(w*X(-1.5)-s/w*X(1.5))/(eps*(1/q-q))
        return np.array([[k,X(-.5),X(.5)],[s*X(.5),k,X(-.5)],[s*X(-.5),s*X(.5),k]])
    U=(s*np.cosh(3*om*th)-np.cos(3*pi*om/2))/(2*np.sin(pi*om/2)**2)
    print('9 s',s,'unitarity error',np.max(abs(km(-th).T@km(th)-U*np.eye(3))))
    Z=np.array([[0,1,0],[0,0,1],[s,0,0]])
    print('9 Z^T-Z^-1',np.max(abs(Z.T-np.linalg.inv(Z))),'printed Z^T-sZ^-1',np.max(abs(Z.T-s*np.linalg.inv(Z))))
    assert np.max(abs(km(-th).T@km(th)-U*np.eye(3)))<1e-11
print('PASS 9: matrix unitarity; printed transpose identity false for s=-1')
