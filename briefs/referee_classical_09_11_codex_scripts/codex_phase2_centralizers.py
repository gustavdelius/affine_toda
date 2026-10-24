"""Independent homogeneous linear solves for principal loop centralizers.
Uses defining matrices, not exponent tables. Triality is built from Chevalley brackets.
"""
import numpy as np
from itertools import product

def E(n,i,j):
 a=np.zeros((n,n));a[i,j]=1;return a
def null(A):
 U,s,Vh=np.linalg.svd(A,full_matrices=True)
 rank=np.sum(s>1e-8)
 return Vh.conj().T[:,rank:],s

def loop_dims(B,heights,sigma,twist,L,R,Z,period,name):
 # X=sum c_b z^k B_b has principal degree height(B_b)+L*k.
 out=[];maxres=0.;gap=1.
 for degree in range(period):
  choices=[(i,(degree-h)//L) for i,h in enumerate(heights) if (degree-h)%L==0]
  ks=sorted(set(k for i,k in choices)|set(k+1 for i,k in choices))
  cols=[]
  for i,k in choices:
   blocks=[]
   for kk in ks:
    v=np.zeros_like(B[0],dtype=complex)
    if kk==k:v+=R@B[i]-B[i]@R
    if kk==k+1:v+=Z@B[i]-B[i]@Z
    blocks.append(v.ravel())
   # Each eigenspace condition lies in its own homogeneous loop coefficient.
   for kk in sorted(set(k0 for _,k0 in choices)):
    blocks.append((sigma(B[i])-twist**k*B[i]).ravel() if kk==k else np.zeros(B[i].size))
   cols.append(np.concatenate(blocks))
  A=np.array(cols).T if cols else np.zeros((1,0))
  N,sv=null(A);out.append(N.shape[1])
  if N.size:maxres=max(maxres,float(np.max(abs(A@N))))
  nz=sv[sv>1e-8]
  if len(nz):gap=min(gap,float(min(nz)))
 print(f'PASS {name}: degrees 0..{period-1} multiplicities={out}; residual={maxres:.2g}; min nonzero singular value={gap:.2g}',flush=True)
 return out

def typeA(N,twisted=False,symplectic=False,L=None,specialZ=False,name=''):
 I=np.eye(N);J=np.fliplr(np.diag([(-1.)**i for i in range(N)]))
 sig=lambda X:-J@X.T@np.linalg.inv(J)
 B=[];h=[]
 for i in range(N):
  for j in range(N):
   if i!=j:B.append(E(N,i,j));h.append(j-i)
 for i in range(N-1):B.append(E(N,i,i)-E(N,i+1,i+1));h.append(0)
 R=sum((E(N,i,i+1) for i in range(N-1)),np.zeros((N,N)))
 Z=E(N,N-1,0)
 if specialZ:Z=E(N,2,0)+E(N,3,1)
 if not twisted and not symplectic:sig=lambda X:X
 L=N if L is None else L
 return loop_dims(B,h,sig,-1 if twisted else 1,L,R,Z,2*L if twisted else L,name)
for N in [3,4,5]:
 out=typeA(N,name=f'A{N-1} untwisted')
 assert out==[0]+[1]*(N-1)
out=typeA(5,twisted=True,name='A4^(2)');assert out==[0,1,0,1,0,0,0,1,0,1]
out=typeA(4,symplectic=True,name='C2^(1)');assert out==[0,1,0,1]
out=typeA(4,twisted=True,L=3,specialZ=True,name='D3^(2) = A3^(2)');assert out==[0,1,0,1,0,1]
# D4 defining representation and triality.
N=8
es=[E(N,0,1)-E(N,5,4),E(N,1,2)-E(N,6,5),E(N,2,3)-E(N,7,6),E(N,2,7)-E(N,3,6)]
perm=[2,1,3,0] # 1 -> 3 -> 4 -> 1; 2 fixed
pos={}
for i,e in enumerate(es):
 key=tuple(int(i==j) for j in range(4));pos[key]=(e,es[perm[i]])
changed=True
while changed:
 changed=False
 for key,(v,img) in list(pos.items()):
  for i,e in enumerate(es):
   newkey=list(key);newkey[i]+=1;newkey=tuple(newkey)
   if newkey in pos:continue
   vv=e@v-v@e
   if np.linalg.norm(vv)>1e-9:
    pos[newkey]=(vv,es[perm[i]]@img-img@es[perm[i]]);changed=True
assert len(pos)==12
B=[];imgs=[];heights=[]
for i,e in enumerate(es):
 B.append(e@e.T-e.T@e);ee=es[perm[i]];imgs.append(ee@ee.T-ee.T@ee);heights.append(0)
for key,(v,img) in sorted(pos.items()):
 B.extend([v,v.T]);imgs.extend([img,img.T]);heights.extend([sum(key),-sum(key)])
flat=np.array([v.ravel() for v in B]).T
pinv=np.linalg.pinv(flat)
imgflat=np.array([v.ravel() for v in imgs]).T
sig=lambda X:(imgflat@pinv@X.ravel()).reshape(N,N)
assert np.max(abs(imgflat@pinv@imgflat@pinv@imgflat-flat))<1e-12
R=sum(es);theta=(1,2,1,1);Z=pos[theta][0].T
assert np.linalg.norm(sig(Z)-Z)<1e-12
out=loop_dims(B,heights,sig,1,6,R,Z,6,'G2^(1) as triality fixed D4');assert out==[0,1,0,0,0,1]
zeta=np.exp(2j*np.pi/3)
v=pos[(1,1,1,0)][0].T
Z=v+sig(v)/zeta+sig(sig(v))/zeta**2
assert np.linalg.norm(sig(Z)-zeta*Z)<1e-12
out=loop_dims(B,heights,sig,zeta,4,R,Z,12,'D4^(3)');assert out==[0,1,0,0,0,1,0,1,0,0,0,1]
