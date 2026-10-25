"""Independent PR60 checks. No imports from the book or original referees."""
import itertools
from collections import defaultdict
import numpy as np
import sympy as s

q,y,e0,e1=s.symbols('q y e0 e1', nonzero=True)
K=s.Matrix([[(q-1/q)*(e1*y+e0/y),y*y-y**-2],
            [y*y-y**-2,(q-1/q)*(e0*y+e1/y)]])
for j,e in [(0,e0),(1,e1)]:
    raising=s.Matrix([[0,1],[0,0]]) if j==1 else s.Matrix([[0,0],[1,0]])
    k=s.diag(q,1/q) if j==1 else s.diag(1/q,q)
    b=y*raising+raising.T*k/(q*y)+e*k
    assert (K*b-b.subs(y,1/y)*K).applyfunc(s.simplify)==s.zeros(2)
print('PASS O1: symbolic coideal intertwining, both nodes')

# Independent Bourbaki E8 roots; Kac labels from solving highest-root relation.
E=np.eye(8)
roots=[(E[0]+E[7]-sum(E[1:7]))/2,E[0]+E[1],E[1]-E[0]]
roots += [E[j]-E[j-1] for j in range(2,7)]
marks=np.array([2,3,4,6,5,4,3,2])
affine=-marks@np.array(roots)
assert abs(affine@affine-2)<1e-14
A=np.vstack([roots,affine]); n=np.r_[marks,1]
H=A.T@np.diag(n)@A
m=np.sqrt(np.linalg.eigvalsh(H)); m/=m[0]
c=lambda x: np.cos(np.pi*x)
m2=2*c(1/5)
expected=np.array([1,m2,2*c(1/30),2*m2*c(7/30),2*m2*c(2/15),
                   2*m2*c(1/30),4*m2*c(1/5)*c(7/30),4*m2*c(1/5)*c(2/15)])
err=max(abs(m-expected)); assert err<1e-12
print('PASS O7a: independent E8 Hessian masses',m,'max error',err)

# Exact Laurent quantum torus U V = q^2 V U, normal form U^a V^b.
# Rational q=2 is generic (not a root of unity), and exact evaluation suffices
# here because the symbolic q identity is computed as well.
def add(*polys):
    d=defaultdict(lambda:s.S(0))
    for p in polys:
        for k,v in p.items(): d[k]+=v
    return {k:s.factor(v) for k,v in d.items() if s.factor(v)!=0}
def scale(p,c): return {k:s.factor(c*v) for k,v in p.items()}
def mul(a,b):
    d=defaultdict(lambda:s.S(0))
    for (u,v),coef in a.items():
        for (w,z),other in b.items(): d[u+w,v+z]+=coef*other*q**(-2*v*w)
    return add(d)
def comm(a,b): return add(mul(a,b),scale(mul(b,a),-1))
def qc(a,b,t): return add(scale(mul(a,b),t),scale(mul(b,a),-1/t))
c,alpha,beta=s.symbols('c alpha beta',nonzero=True)
fac=s.I*c/(q-1/q)
P={(1,0):fac*alpha,(-1,0):fac/alpha}
Q={(0,1):fac*beta,(0,-1):fac/beta}
rho=c*c*(q+1/q)**2
bk=add(comm(mul(P,Q),mul(Q,P)),scale(add(mul(P,P),scale(mul(Q,Q),-1)),-rho))
assert not bk
print('PASS O7b: exact symbolic BK (19), arbitrary q,c,alpha,beta')
for X,Y in [(P,Q),(Q,P)]:
    dg=add(comm(X,qc(X,qc(X,Y,q),1/q)),scale(comm(X,Y),-rho))
    assert not dg
print('PASS O7b: exact symbolic q-Dolan-Grady, both signs')

# A direct counterexample to the OLD REPORT's exponent inequality, not the book.
mm,aa=s.Integer(2),s.Integer(1)
gamma=(2*mm-aa)/(mm+1)
assert 1-gamma < aa/(2*mm)
print('PASS O6: report inequality is false: m=2 alpha=1, 1-gamma=0 < alpha/(2m)=1/4')

# Enumerate compatible dendrograms for every labelled tree n <= 5.
def tree_from_prufer(seq,n):
    deg=[1]*n
    for v in seq: deg[v]+=1
    edges=[]
    for v in seq:
        u=next(i for i in range(n) if deg[i]==1)
        edges.append((u,v)); deg[u]-=1; deg[v]-=1
    edges.append(tuple(i for i in range(n) if deg[i]==1))
    return edges
count=0
for nn in range(2,6):
    for seq in itertools.product(range(nn),repeat=nn-2):
        edges=tree_from_prufer(seq,nn); topologies={}
        for order in itertools.permutations(range(nn-1)):
            comp=[frozenset([i]) for i in range(nn)]; clusters=[]
            for ix in order:
                u,v=edges[ix]; C=comp[u]|comp[v]; clusters.append(C)
                for k in C: comp[k]=C
            key=tuple(sorted(tuple(sorted(C)) for C in clusters))
            topologies[key]=topologies.get(key,0)+1
        weights=0
        for key,num in topologies.items():
            hs=[len(C)-1 for C in key]; prod=s.prod(hs)
            assert num==s.factorial(nn-1)/prod
            weights+=s.Rational(1,prod)
            for mm in range(1,nn):
                gamma=s.Rational(2*mm,mm+1)*s.Rational(9,10)
                aa=2*mm-(mm+1)*gamma
                for h in hs:
                    C=h+1 if h>=mm else 0
                    exponent=2*h-gamma*C
                    assert exponent>0 and 2*h/exponent<=2*mm/aa
        assert weights==1
        count+=1
print('PASS O6: hook counts, sum=1, corrected exponent bound for',count,'labelled trees')

# Nested caterpillar integral, m=2, gamma=6/5: c=(0,3,1), E=(2,2/5,6/5).
u,v,w=s.symbols('u v w', positive=True)
integral=s.integrate(s.integrate(s.integrate(u*v**s.Rational(-13,5)*w**s.Rational(-1,5),
                  (u,0,v)),(v,0,w)),(w,0,1))
assert integral==s.Rational(25,24)
print('PASS O6: exact nested integral',integral)
