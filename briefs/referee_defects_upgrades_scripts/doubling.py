"""Exact difference-operator arithmetic for chapter 42, symbolic q, q^n, q^s.
No lattice truncation and no imports from the book.
"""
import sympy as S
q,z,s,y,d=S.symbols('q z s y yd',nonzero=True)
def clean(a):return {k:S.cancel(v) for k,v in a.items() if S.cancel(v)!=0}
def mul(A,B):
    C={}
    for (i,j,k),a in A.items():
        for (jj,l,kk),b in B.items():
            if j==jj:
                key=(i,l,k+kk);C[key]=C.get(key,0)+a.subs(z,z*q**kk)*b
    return clean(C)
def add(A,B,c=1):
    C=dict(A)
    for k,b in B.items():C[k]=C.get(k,0)+c*b
    return clean(C)
def scale(A,c):return {k:c*v for k,v in A.items()}
def scalar(a):return {(0,0,0):a,(1,1,0):a}
def diag(a,b):return {(0,0,0):a,(1,1,0):b}
def pmat(i,j,a=1):return {(i,j,0):S.sympify(a)}
def tensor(i,j,k,a):return {(i,j,k):S.sympify(a)}
def wop(k,a):return {(0,0,k):a,(1,1,k):a}
L={(0,0,0):1/z,(1,1,0):z*s,(0,1,-1):-S.I*y*z*s/d,(1,0,1):-S.I*y/(d*z)}
one=scalar(1)
Y1=wop(1,S.I*d/(s*(q*q-1)*z*z))
Y0=wop(-1,S.I*d*s*z*z/(q*q-1))
K1=wop(2,1/(q*q*s*z*z))
K0=wop(-2,s*z*z/(q*q))
for j,Y,K,Kinv in [(1,Y1,K1,K0),(0,Y0,K0,K1)]:
    i0,i1=(0,1) if j==1 else (1,0)
    kp=diag(q,1/q) if j==1 else diag(1/q,q)
    kout=diag(1/q,q) if j==1 else diag(q,1/q)
    pin=pmat(i0,i1,y); pout=pmat(i1,i0,y)
    incoming=add(pin,mul(kp,Y))
    outgoing=add(Y,mul(pout,K))
    assert not add(mul(L,incoming),mul(outgoing,L),-1)
    assert not add(mul(L,mul(kp,K)),mul(mul(kout,K),L),-1)
    assert not add(mul(K,Kinv),one,-1)
    assert not add(mul(mul(K,Y),Kinv),scale(Y,q*q),-1)
    print(f'node {j}: e equation, k equation, inverse and own Cartan relation: exact zero')
for K,Kinv,Y in [(K1,K0,Y0),(K0,K1,Y1)]:
    assert not add(mul(mul(K,Y),Kinv),scale(Y,q**-2),-1)
assert not add(mul(Y1,Y0),mul(Y0,Y1),-1)
print('cross Cartan equations: exact zero; e1 e0=e0 e1=-yd^2 q^2/(q^2-1)^2; cubic Serre follows')

# Fully arbitrary noncommuting unknowns for the obstruction; no charge ansatz.
A,B,C,D,Y,K,t=S.symbols('A B C D Y K t',commutative=False)
Op=S.Matrix([[0,B],[C,0]]); Di=S.diag(A,D)
Ep=S.Matrix([[0,1],[0,0]]); Em=Ep.T
print('minus untwisted e, y^2:',Op*Ep-t*Ep*S.diag(K,K)*Op)
print('minus twisted f, y^-1:',Di*Em*S.diag(K,K)-t**-1*Ep*Di)
print('plus untwisted f, y^-2:',Op*Em*S.diag(K,K)-t**-1*Em*Op)
print('plus twisted e, y:',Di*Ep-t*Em*S.diag(K,K)*Di)
print('plus coideal isolated 1 tensor f0, y^-2:',-Ep*Op)
print('PASS. The report derives the forced inverse-Cartan contradictions in both limits.')
