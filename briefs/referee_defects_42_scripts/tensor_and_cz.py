"""Independent finite-dimensional tensor K check and final literal CZ correction check."""
import numpy as np
import scipy.linalg as la
from independent_checks import q, LAM, rep

def checkR(y,z):
    rows=[]
    for j in (0,1):
        e,f,k=rep(y,j);E,F,K=rep(z,j)
        for A,B in ((np.kron(e,np.eye(2))+np.kron(k,E),np.kron(E,np.eye(2))+np.kron(K,e)),
                    (np.kron(f,la.inv(K))+np.kron(np.eye(2),F),np.kron(F,la.inv(k))+np.kron(np.eye(2),f)),
                    (np.kron(k,K),np.kron(K,k))):
            rows.append(np.kron(A.T,np.eye(4))-np.kron(np.eye(4),B))
    ns=la.null_space(np.vstack(rows));assert ns.shape[1]==1
    R=ns[:,0].reshape(4,4,order='F');return R/R[0,0]

P=np.eye(4)[[0,2,1,3]];eps={1:.34+.12j,0:-.47+.08j};d=np.exp(.29*LAM)
def K(y):
    return np.array([[(q-q**-1)*(eps[1]*y+eps[0]/y),y*y-y**-2],
                     [y*y-y**-2,(q-q**-1)*(eps[0]*y+eps[1]/y)]])
def X(y,j):
    e,f,k=rep(y,j);E,F,kk=rep(d,j)
    return np.kron(e+q**-1*f@k,np.eye(2))+np.kron(k,E+q**-1*F@kk+eps[j]*kk)
def fused(y):
    L=P@checkR(y,d);Lm=P@checkR(1/y,d)
    return la.solve(Lm,np.kron(K(y),np.eye(2))@L)
for th in (.37,-.28+.13j):
    y=np.exp(LAM*th);F=fused(y);rows=[];err=[]
    for j in (0,1):
        A=X(y,j);B=X(1/y,j)
        rows.append(np.kron(A.T,np.eye(4))-np.kron(np.eye(4),B))
        err.append(la.norm(F@A-B@F)/la.norm(F@A))
    print('tensor K theta/residual/intertwiner dimension:',th,max(err),la.null_space(np.vstack(rows)).shape[1])

def RE(t1,t2):
    y1,y2=np.exp(LAM*t1),np.exp(LAM*t2)
    A=np.kron(checkR(1/y2,1/y1),np.eye(2))@np.kron(np.eye(2),fused(y1))@np.kron(checkR(y1,1/y2),np.eye(2))@np.kron(np.eye(2),fused(y2))
    B=np.kron(np.eye(2),fused(y2))@np.kron(checkR(y2,1/y1),np.eye(2))@np.kron(np.eye(2),fused(y1))@np.kron(checkR(y1,y2),np.eye(2))
    return la.norm(A-B)/la.norm(A)
print('tensor K reflection equation:',RE(.37,-.28))

# Read source formulas from the independently transcribed module; it also prints diagnostic failures.
import source_and_limits as cz
for a in (0,2,4):
    F=la.solve(cz.T(1/cz.x,a),cz.kd@cz.T(cz.x,a))
    Pc=cz.printed(a,cz.hcorrect);Pc[0,1]*=-1
    scale=np.vdot(Pc,F)/np.vdot(Pc,Pc)
    print('CZ 2.19 corrected h and -s-minus, label/residual:',a,la.norm(F-scale*Pc)/la.norm(F))
