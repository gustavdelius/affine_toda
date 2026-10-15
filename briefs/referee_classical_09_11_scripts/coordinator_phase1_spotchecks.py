import numpy as np, itertools, mpmath as mp
# 1. d4^(1): couplings in a real charge basis (Z2xZ2 acting by double transpositions of outer nodes)
e=np.eye(4)
al={1:e[0]-e[1],2:e[1]-e[2],3:e[2]-e[3],4:e[2]+e[3]}
al[0]=-(e[0]+e[1])           # minus highest root e1+e2
n={0:1,1:1,2:2,3:1,4:1}
M2=sum(n[j]*np.outer(al[j],al[j]) for j in n)
w,V=np.linalg.eigh(M2); print("d4 mass^2:",np.round(w,6))
# orthogonal maps realising the double transpositions of outer nodes {0,1,3,4}
def ortho(perm):  # find O with O al[j] = al[perm[j]] for simple j=1..4
    A=np.array([al[j] for j in (1,2,3,4)]).T; B=np.array([al[perm[j]] for j in (1,2,3,4)]).T
    return B@np.linalg.inv(A)
# perms fixing node 2, given on all nodes; check O al0 = al[perm0]
gens=[{0:1,1:0,3:4,4:3,2:2},{0:3,3:0,1:4,4:1,2:2}]
Os=[]
for p in gens:
    O=ortho(p); assert np.allclose(O@al[0],al[p[0]]) and np.allclose(O.T@O,np.eye(4)); Os.append(O)
light=V[:,w<3]; heavy=V[:,w>3]
# simultaneous eigenvectors in light space
P=light.T@(Os[0]+2.5*Os[1])@light; ev,U=np.linalg.eigh((P+P.T)/2)
basis=[heavy[:,0]]+[light@U[:,i] for i in range(3)]
chars=[(round(float(b@Os[0]@b),3),round(float(b@Os[1]@b),3)) for b in basis]
print("characters (heavy first):",chars)
C=lambda a,b,c: sum(n[j]*(al[j]@a)*(al[j]@b)*(al[j]@c) for j in n)
names=['H','L1','L2','L3']
for i,j,k in itertools.combinations_with_replacement(range(4),3):
    v=C(basis[i],basis[j],basis[k])
    if abs(v)>1e-9: print("nonzero C",names[i],names[j],names[k],round(abs(v),6))
# 2. A_aa(2iu): book formula vs Gemini's -tan^2(pi a/h)/tan^2 u
A=lambda th,a,b,h: (mp.cosh(th)-mp.cos(mp.pi*(a-b)/h))/(mp.cosh(th)-mp.cos(mp.pi*(a+b)/h))
for (a,h,u) in [(1,3,mp.pi/4),(1,3,0.2),(1,4,1.2)]:
    print("A_aa(2iu) a,h,u=",a,h,float(u)," book:",mp.nstr(mp.re(A(2j*u,a,a,h)),8)," Gemini:",mp.nstr(-mp.tan(mp.pi*a/h)**2/mp.tan(u)**2,8))
# 3. type-B band: A_{a,h-a}(2iu)=1-sin^2(pi a/h)/cos^2 u
for (a,h,u) in [(1,3,0.2),(1,3,mp.pi/3)]:
    print("A_{a,h-a}(2iu)",a,h,float(u),mp.nstr(mp.re(A(2j*u,a,h-a,h)),12),mp.nstr(1-mp.sin(mp.pi*a/h)**2/mp.cos(u)**2,12))
# 4. singular sine-Gordon static: energy density 2*d2 log(1-e^{4x}) ; integral with (-eps,eps) removed
mp.mp.dps=30
dens=lambda x: 2*mp.diff(lambda y: mp.log(1-mp.exp(4*y)),x,2)
print("density at 0.3 vs -8csch^2(0.6):",mp.nstr(dens(mp.mpf('0.3')),12),mp.nstr(-8/mp.sinh(0.6)**2,12))
for eps in ['0.01','0.001']:
    eps=mp.mpf(eps); Ee=2*mp.quad(lambda x:-8/mp.sinh(2*x)**2,[eps,1,mp.inf])
    print("E_eps",eps,mp.nstr(Ee,15),"  8(1-coth2eps)=",mp.nstr(8*(1-mp.coth(2*eps)),15))
# 5. Codex counterexample: SG theta=i pi, xi=i pi/2 -> static antikink, energy by quadrature
t0=lambda x:1+1j*mp.exp(-2*x); t1=lambda x:1-1j*mp.exp(-2*x)
phi=lambda x:(1j/1)*mp.sqrt(2)*(mp.log(t1(x))-mp.log(t0(x)))/2*2/ mp.sqrt(2)  # alpha1=-alpha0=sqrt2, coroot=alpha
# phi = (i/beta) sum_j alpha_j^vee ln tau_j with alpha_1^vee=sqrt2, alpha_0^vee=-sqrt2 (tau_0 at node 0)
phi=lambda x: 1j*mp.sqrt(2)*(mp.log(t1(x))-mp.log(t0(x)))
print("phi(0.2) =",mp.nstr(phi(mp.mpf('0.2')),12))
V=lambda f: -(1)*(mp.exp(1j*mp.sqrt(2)*f)-1)-(mp.exp(-1j*mp.sqrt(2)*f)-1)  # m=beta=1, n0=n1=1
Ed=lambda x: 0.5*mp.diff(phi,x)**2+V(phi(x))
import mpmath as mp
mp.mp.dps=30
# Codex counterexample: SG (a_1^(1), alpha_1=-alpha_0=sqrt2, m=beta=1), theta=i*pi, xi=i*pi/2
# tau_j = 1 + omega^j E, omega=-1, E=exp(m_1(x cosh th - t sinh th)+xi), m_1=2
E=lambda x: mp.exp(2*x*mp.cosh(1j*mp.pi)+1j*mp.pi/2)
t0=lambda x:1+E(x); t1=lambda x:1-E(x)
# phi = (i/beta) sum_j alpha_j^vee ln tau_j ; alpha^vee=alpha for |alpha|^2=2
phi=lambda x: 1j*mp.sqrt(2)*(mp.log(t1(x))-mp.log(t0(x)))
dphi=lambda x: mp.diff(phi,x)
V=lambda f: -((mp.exp(1j*mp.sqrt(2)*f)-1)+(mp.exp(-1j*mp.sqrt(2)*f)-1))
Ed=lambda x: dphi(x)**2/2+V(phi(x))
print("phi(-5),phi(5):",mp.nstr(phi(-5),8),mp.nstr(phi(5),8))
print("energy by quadrature:",mp.nstr(mp.quad(Ed,[-15,0,15]),12),"; formula M cosh(theta) =",mp.nstr(8*mp.cosh(1j*mp.pi),6))
