"""Small independent cross-checks, using only referee-owned root construction."""
import numpy as np
import sympy as S
from check_algebra import finite,cartan

b=S.symbols('b',positive=True);B=lambda y:2*y/(y+4*S.pi)
Hb=4-B(b)/2;Hc=4+B(b/2)
assert S.simplify(1/Hb+1/Hc-S.Rational(1,2))==0
print('PASS B2=C2: 1/H_B+1/H_C=1/2 exactly, so the two mass-ratio parametrizations agree.')

for N in range(3,9):
    w=np.exp(2j*np.pi/N);R=np.roll(np.eye(N),1,axis=0);F=np.eye(N)[[(-j)%N for j in range(N)]]
    for a in range(1,N):
        v=np.array([w**(a*j) for j in range(N)])
        assert np.allclose(R@v,w**(-a)*v) and np.allclose(F@v,v.conj())
        weight=np.array([1-a/N if j<a else -a/N for j in range(N)])
        lam1=np.array([1-1/N]+[-1/N]*(N-1))
        assert np.allclose(weight-a*lam1,np.rint(weight-a*lam1)) and abs(sum(weight-a*lam1))<1e-10
print('PASS A_n rotation acts by phases; reflection exchanges a and N-a; [lambda_a]=a[lambda_1].')

for typ,n in [('A',4),('D',5),('E',6),('E',7),('E',8)]:
    G,roots,Ga,marks=finite(typ,n);R=np.linalg.cholesky(np.array(G,float))
    Rp=np.vstack([-np.array(marks[1:],float).reshape(n)@R,R]);M=Rp.T@np.diag(np.array(marks,float).ravel())@Rp
    masses=np.sort(np.sqrt(np.linalg.eigvalsh(M)))
    ev,vec=np.linalg.eigh(np.array(cartan(G),float));pf=np.sort(abs(vec[:,0]));pf=pf/min(pf)
    assert np.allclose(masses/min(masses),pf)
    # Principal sl2 decomposition: multiplicity of weight k minus weight k+1.
    heights={k:sum(sum(r)==k for r in roots) for k in range(1,max(map(sum,roots))+2)}
    sl2=[k for k in range(1,max(heights)) for _ in range(heights[k]-heights[k+1])]
    print(f'PASS {typ}{n}: Hessian masses/PF; principal sl2 spins from root heights={sl2}')

h=30
for theta in [np.sqrt(2),np.e/2,np.sqrt(3)+.2j]:
    block=lambda x:np.sinh(theta/2+1j*np.pi*x/(2*h))/np.sinh(theta/2-1j*np.pi*x/(2*h))
    minimal=np.prod([block(x-1)*block(x+1) for x in [1,11,19,29]])
    ising=np.prod([(np.sinh(theta)+1j*np.sin(np.pi*a))/(np.sinh(theta)-1j*np.sin(np.pi*a)) for a in [2/3,2/5,1/15]])
    assert abs(minimal-ising)<1e-12
print('PASS E8 minimal S11: product blocks at {1,11,19,29} equals f_(2/3) f_(2/5) f_(1/15).')

G,roots,Ga,marks=finite('D',4)
ev,vec=np.linalg.eigh(np.array(G,float));v=abs(vec[:,0]);v=v/min(v)
assert np.allclose(sorted(v),[1,1,1,np.sqrt(3)])
print('PASS D4 exercise PF entries:',v)
v=S.Matrix([S.cos(5*S.pi/12),S.sin(5*S.pi/12)])
a=S.Matrix([S.sqrt(S.Rational(2,3)),0]);c=S.Matrix([-S.sqrt(S.Rational(3,2)),1/S.sqrt(2)])
assert float(a.dot(v))>0 and float(c.dot(v))>0
print('PASS G2 chamber: 60 degrees < angle < 90 degrees; interior 75 degrees has both simple-root pairings positive.')
