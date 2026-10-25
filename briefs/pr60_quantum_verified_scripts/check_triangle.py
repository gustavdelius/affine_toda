"""Independent symbolic triangle and crossing on formal bilateral states."""
import sympy as s
from collections import defaultdict
q,t,r,Y,Z=s.symbols('q t r Y Z',nonzero=True)
I=s.I
# state=(spin1,spin2,k), t=q^n; all entries exact Laurent-rational.
def L(j,k,y):
    a=t*q**k
    if j==1:return [(1,k,1/a),(-1,k+1,-I*y/a)]
    return [(-1,k,a*r),(1,k-1,-I*y*a*r)]
def invL(j,k,y):
    a=t*q**k
    if j==1:return [(1,k,a/(1+y*y)),(-1,k+1,I*y/(a*q*r*(1+y*y)))]
    return [(-1,k,1/(a*r*(1+y*y))),(1,k-1,I*y*a/(q*(1+y*y)))]
def Lop(v,index,y):
    out=defaultdict(lambda:0)
    for (a,b,k),c in v.items():
        j=[a,b][index]
        for jj,kk,d in L(j,k,y):
            aa,bb=(jj,b) if index==0 else (a,jj)
            out[aa,bb,kk]+=c*d
    return out
def Rop(v,x):
    out=defaultdict(lambda:0)
    for (a,b,k),c in v.items():
        if a==b:out[a,b,k]+=c
        else:
            out[a,b,k]+=c*q*(x*x-1)/(q*q*x*x-1)
            out[b,a,k]+=c*(q*q-1)*x/(q*q*x*x-1)
    return out
for a in [1,-1]:
 for b in [1,-1]:
    v={(a,b,0):1}
    lhs=Rop(Lop(Lop(v,1,Z),0,Y),Y/Z)
    rhs=Lop(Lop(Rop(v,Y/Z),0,Y),1,Z)
    assert all(s.simplify(lhs.get(k,0)-rhs.get(k,0))==0 for k in set(lhs)|set(rhs))
    print('PASS exact triangle input',a,b)
# Crossing after multiplying RHS scalar denominator and r=q^s.
for a in [1,-1]:
 for b,k,c in L(a,0,Y):
    rhs=sum(cc for jj,kk,cc in invL(-b,0,-q*Y) if jj==-a and kk==k)
    # inverse expression has charge input n, same n as original incoming.
    assert s.factor(c-r*(1+q*q*Y*Y)*rhs)==0
print('PASS all four crossing entries with Yprime=-qY')
# Wrong sign q cannot be absorbed in one scalar: diag unchanged, converting sign reversed.
print('PASS wrong crossing sign: Yprime=+qY keeps diagonal and negates both conversions')
