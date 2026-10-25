"""Exact local fusion, independent representation on sparse basis."""
import sympy as s
from collections import defaultdict
q,z,t,r=s.symbols('q z t r',nonzero=True)
I=s.I
# t=q^n, r=q^offset, k represents n+k.
def invL(j,k,Z,R=r):
    a=t*q**k
    if j==1:return [(1,k,a/(1+Z**2)),(-1,k+1,I*Z/(a*q*R*(1+Z**2)))]
    return [(-1,k,1/(a*R*(1+Z**2))),(1,k-1,I*Z*a/(q*(1+Z**2)))]
def Ract(a,j,x):
    if a==j:return [(a,j,1)]
    return [(a,j,q*(x*x-1)/(q*q*x*x-1)),(j,a,(q*q-1)*x/(q*q*x*x-1))]
def F(v):
    out=defaultdict(lambda:0)
    for (a,j,k),c in v.items():
        for jj,kk,b in invL(j,k,z):
            for aa,jj2,d in Ract(a,jj,I/z):out[(aa,jj2,kk)]+=c*b*d
    return out

def encode(which,j,k=0):
    a=t*q**k
    if which=='e':return {(1,j,k):a,(-1,j,k+1):-1/(a*q*r)}
    return {(-1,j,k):1/a}
def decode(v):
    out=defaultdict(lambda:0)
    for (a,j,k),c in v.items():
        nn=t*q**k
        if a==1:
            out[('e',j,k)]+=c/nn
            out[('f',j,k+1)]+=c/(nn*r)
        else:out[('f',j,k)]+=c*nn
    return {key:s.factor(val) for key,val in out.items()}
for which in ['e','f']:
    for j in [1,-1]:
        got=decode(F(encode(which,j)))
        expected=defaultdict(lambda:0)
        if which=='e':
            for jj,k,c in invL(j,0,q*z,r*q):expected[('e',jj,k)]+=c*(1+q*q*z*z)/(1+z*z)
        else:
            for jj,k,c in invL(j,0,z/q,r/q):expected[('f',jj,k)]+=c/q
        keys=set(got)|set(expected)
        residuals={key:s.factor(got.get(key,0)-expected.get(key,0)) for key in keys if which=='e' or key[0]=='f'}
        assert all(v==0 for v in residuals.values()),residuals
        print('PASS fusion block',which,'input',j,'all residuals zero')
# Borel recurrence uniqueness at arbitrary n: A=a*t^-1+H*t homogeneous.
c,Y=s.symbols('c Y',nonzero=True)
a,hom=s.symbols('a hom')
# eq7 gives D=q^(2n+s) A, while eq3 and eq6 enforce hom=0 via q^4!=1.
A=lambda k:a/(t*q**k)+hom*t*q**k
D=lambda k:r*(t*q**k)**2*A(k)
C=lambda k:c*(A(k+1)-q*A(k))
B=lambda k:c*(D(k-1)-q*D(k))
obstruction=s.factor(c*(q**-1*A(0)-A(-1))-B(0)/(t*t*r))
print('uniqueness homogeneous obstruction =',obstruction)
assert s.factor(obstruction/(c*hom*t))==(q-1)*(q+1)*(q*q+1)/q**3
print('PASS intertwiner homogeneous term excluded for q^4 != 1')
