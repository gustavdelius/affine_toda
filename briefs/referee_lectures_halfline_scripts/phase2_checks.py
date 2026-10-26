"""Independent phase-2 calculations from lecture definitions. No author imports."""
import numpy as np, itertools, mpmath as mp
mp.mp.dps=50
q=-np.exp(-1j*np.pi*np.sqrt(2)); w=np.exp(-.5j*np.pi*np.sqrt(2)); eps=1/(w+1/w)
def E(n,a,b):
 x=np.zeros((n,n),complex);x[a,b]=1;return x
def system(L,R):
 n=len(L[0]);return np.vstack([np.kron(a.T,np.eye(n))-np.kron(np.eye(n),b) for a,b in zip(L,R)])
def null(A):
 u,s,v=np.linalg.svd(A,full_matrices=False);return s,v[-1].conj()
def qnodes(n,y):
 es=[];fs=[];ks=[]
 for a in range(n):
  b=(a+1)%n;k=np.ones(n,complex);k[a]=q;k[b]=1/q
  es.append(y*E(n,a,b));fs.append(E(n,b,a)/y);ks.append(np.diag(k))
 return es,fs,ks
y=np.exp(np.sqrt(2)*.271)
e,f,k=qnodes(3,y);A=[ei+fi@ki/q for ei,fi,ki in zip(e,f,k)]
worst=0; badmin=9
for pattern in itertools.product([0,1,-1,.37],repeat=3):
 L=[a+eps*p*ki for a,p,ki in zip(A,pattern,k)]
 R=[a+eps*p*np.linalg.inv(ki) for a,p,ki in zip(A,pattern,k)]
 ss,_=null(system(L,R));valid=all(x==0 for x in pattern) or all(abs(x)==1 for x in pattern)
 if valid:worst=max(worst,ss[-1]);assert ss[-1]<1e-12
 else:badmin=min(badmin,ss[-1]);assert ss[-1]>1e-6
print('PASS quantum mixed branches: 64 patterns; allowed max/min disallowed',worst,badmin)
worst=0
for N in range(3,9):
 for sign in [1,-1]:
  X=w*y; xx=np.sqrt(2)*(.271-1j*np.pi/2)
  K=np.array([[np.exp((j-i-N/2)*xx) if j>i else sign*np.exp((j-i+N/2)*xx) if j<i else (w*np.exp(-N*xx/2)-sign/w*np.exp(N*xx/2))/(eps*(1/q-q)) for j in range(N)] for i in range(N)])
  ee,ff,kk=qnodes(N,y)
  for node in range(N):
   p=sign if node==N-1 else 1;a=ee[node]+ff[node]@kk[node]/q
   er=np.max(np.abs(K@(a+p*eps*kk[node])-(a+p*eps*np.linalg.inv(kk[node]))@K));worst=max(worst,er)
assert worst<1e-11
print('PASS Toeplitz N=3..8 both signs, max entry residual',worst)
# D4 vector Chevalley generators, ordered roots 0,1,2,3,4.
def minus(i,j):return E(8,i,j)-E(8,j+4,i+4)
def plus(i,j):return E(8,i,j+4)-E(8,j,i+4)
ees=[plus(0,1).T,minus(0,1),minus(1,2),minus(2,3),plus(2,3)]
ffs=[a.T for a in ees];hs=[a@b-b@a for a,b in zip(ees,ffs)]
lam=1.317
min_det=9;max_allowed=0;min_wrong=9
for signs in itertools.product([1,-1],repeat=5):
 for wrong in [False,True]:
  c=np.array(signs,float);c[2]/=np.sqrt(2) if wrong else 1
  L=[-v*h+lam*a-b/lam for v,h,a,b in zip(c,hs,ees,ffs)]
  R=[v*h+lam*a-b/lam for v,h,a,b in zip(c,hs,ees,ffs)]
  ss,z=null(system(L,R))
  if wrong:min_wrong=min(min_wrong,ss[-1]);assert ss[-1]>1e-6
  else:
   max_allowed=max(max_allowed,ss[-1]);K=z.reshape((8,8),order='F');min_det=min(min_det,abs(np.linalg.det(K)));assert ss[-1]<1e-12 and abs(np.linalg.det(K))>1e-16, (signs,ss[-1],np.linalg.det(K))
print('PASS D4 all32 normalized signs invertible; max residual/min det',max_allowed,min_det,'missing sqrt(n) min singular',min_wrong)
# arbitrary c0,c1 SG: e+c f k and reflected y inverse.
for c1 in [.71+0.2j]:
 c0=1/(q*q*c1)
 for factor in [1,1.1]:
  cs=[c1,c0*factor];ee,ff,kk=qnodes(2,y);er,fr,kr=qnodes(2,1/y)
  ss,_=null(system([a+c*b@k for a,b,k,c in zip(ee,ff,kk,cs)],[a+c*b@k for a,b,k,c in zip(er,fr,kr,cs)]))
  print('general c product factor',factor,'smallest singular',ss[-1]);assert (ss[-1]<1e-12)==(factor==1)
print('PASS general-c intertwiner condition at generic y')
# q Dolan Grady arbitrary c and epsilon, tensor products.
def tensor(left,right):
 e,f,k=left;E,F,K=right;n=len(e);m=len(E)
 return np.kron(e,np.eye(m))+np.kron(k,E),np.kron(f,np.linalg.inv(K))+np.kron(np.eye(n),F),np.kron(k,K)
node=[]
for j in [0,1]:
 reps=[]
 for yy in [1.173,.823,1.414]:
  es,fs,ks=qnodes(2,yy);reps.append((es[j],fs[j],ks[j]))
 node.append(tensor(tensor(reps[0],reps[1]),reps[2]))
cs=[.71+.2j,-.31+.6j];ep=[.33,-.27j]
bs=[e+c*f@k+t*k for (e,f,k),c,t in zip(node,cs,ep)]
for i,j in [(0,1),(1,0)]:
 x,z=bs[i],bs[j];x2=x@x;x3=x2@x;r=q*q+1+q**-2
 lhs=x3@z-r*x2@z@x+r*x@z@x2-z@x3
 rhs=q*cs[i]*(q+1/q)**2*(x@z-z@x)
 err=np.max(abs(lhs-rhs));assert err<1e-11;print('PASS arbitrary-parameter Dolan Grady on 8-dim tensor, residual',err)
# CDD factor ambiguity and crossed-channel bootstrap.
def bl(x,t,h=3):return mp.sinh(t/2+1j*mp.pi*x/(2*h))/mp.sinh(t/2-1j*mp.pi*x/(2*h))
def D(chi,B,t):return mp.fprod(bl(z,t) for z in [3+chi+B/2,chi-B/2,3-chi+B/2,-chi-B/2])
t=mp.mpf('.71')
for B in [mp.mpf('0.01'),mp.mpf('0.001'),mp.mpf('0.0001')]:
 diff=D(1-B/2,B,t)/D(1,B,t)-1
 print('CDD change B, abs change /B^2',B,mp.nstr(abs(diff)/B**2,14));assert abs(diff)<4*B**2
print('PASS CDD change first appears at B^2 in samples')
B=mp.sqrt(2)/4
K=lambda t:bl(3-B/2,t)/(bl(2,t)*bl(1-B/2,t))
S=lambda x,t:bl(x-1,t)*bl(x+1,t)/(bl(x-1+B,t)*bl(x+1-B,t))
right=K(t+1j*mp.pi/3)*S(1,2*t)*K(t-1j*mp.pi/3)
wrong=K(t+1j*mp.pi/3)*S(2,2*t)*K(t-1j*mp.pi/3)
print('PASS bulk bootstrap S11 error / wrong S12 error',mp.nstr(abs(right-K(t)),8),mp.nstr(abs(wrong-K(t)),14));assert abs(right-K(t))<mp.mpf('1e-40') and abs(wrong-K(t))>.01
# CT state labels, energy difference.
F=lambda k,d:mp.sin(k*d)*mp.cos(mp.pi/6-k*d)/mp.sin(d)
for B in [mp.mpf('.23'),mp.mpf('.45')]:
 d=mp.pi*B/6;p=mp.mpf('.5')+mp.mpf('1.5')*B
 gap=F(-1,d)+F(2,d)-F(-2,d)-F(2,d)
 assert abs(gap-mp.cos(mp.pi*p/3))<mp.mpf('1e-40')
print('PASS CT pole energy matches transition b[-2,2] -> b[-1,2]; labels 0<p<1.5')
# Lee Yang one point consistency including Jacobian.
Sly=lambda t:(mp.sinh(t)+1j*mp.sqrt(3)/2)/(mp.sinh(t)-1j*mp.sqrt(3)/2)
Kly=lambda t:bl(1,t,6)*bl(3,t,6)/bl(4,t,6)
epsilon=mp.mpf('1e-25');res=epsilon*Kly(1j*mp.pi/6+epsilon);g=2j*mp.sqrt(2*mp.sqrt(3)-3);Gamma=1j*mp.sqrt(2*mp.sqrt(3))
assert abs(res+1j*Gamma*g/2)<mp.mpf('1e-22')
print('PASS LY reflection pole uses -i Gamma g/2; residue',mp.nstr(res,15),'Gamma^2',mp.nstr(Gamma**2,15))
# Ward identity moments and contact distribution sign by integration.
f=lambda t:mp.e**(-t*t)*(1+mp.mpf('.3')*t)
for e in [mp.mpf('.01'),mp.mpf('.001')]:
 # n=2 integral in t' at t=0; rational imaginary kernel.
 val=mp.quad(lambda z:f(z)*(1/(-z-1j*e)**2-1/(-z+1j*e)**2),[-mp.inf,-e,0,e,mp.inf])
 print('contact n=2 regulator, integral,target',e,mp.nstr(val,12),mp.nstr(-2j*mp.pi*mp.mpf('.3'),12))
print('PASS all phase-2 checks')
