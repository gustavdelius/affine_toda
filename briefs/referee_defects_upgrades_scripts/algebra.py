"""Independent exact lattice propagation. No imports from the book or old checks.
Node j connects coordinate j to j+1 (zero-based, cyclic).
"""
import sympy as S
from collections import defaultdict

def shift(n,a,b):
    m=list(n); m[a]+=1; m[b]-=1
    return tuple(m)
def linear(op,v):
    out=defaultdict(lambda:S.S(0))
    for s,c in v.items():
        for t,d in op(s).items(): out[t]+=c*d
    return {s:S.cancel(c) for s,c in out.items() if S.cancel(c)!=0}
def difference(v,w,c=1):
    out=dict(v)
    for s,d in w.items(): out[s]=out.get(s,0)-c*d
    return {s:S.cancel(d) for s,d in out.items() if S.cancel(d)!=0}
def add(v,w): return difference(v,w,-1)

for h in range(3,8):
    q=S.Integer(2)**h; y=S.Rational(3,2); yd=S.Rational(5,3); mu=S.Rational(7,4)
    X=-(q-1/q)*mu*y/yd
    Xp=(q-1/q)*mu*yd/y
    n0=tuple([S.Rational(i+1,h) for i in range(h-1)]+[-S.Rational(h-1,2)])
    assert sum(n0)==0
    def particle(kind,j,s,bar):
        i,n=s; k=(j+1)%h
        a,b=(j,k) if (kind=='e') != bar else (k,j)
        if i!=b:return {}
        return {(a,n):y if kind=='e' else 1/y}
    def pk(j,i,bar): return q**(((-1)**bar)*(int(i==j)-int(i==(j+1)%h)))
    def module(kind,j,s,orient=1):
        i,n=s; k=(j+1)%h
        nn=shift(n,j,k) if kind=='e' else shift(n,k,j)
        return {(i,nn):mu*(yd if kind=='e' else 1/yd)*q**(orient*n[k])}
    def copro(kind,j,s,bar,opposite=False,orient=1):
        i,n=s; k=(j+1)%h
        p=particle(kind,j,s,bar); w=module(kind,j,s,orient)
        if kind=='e':
            if opposite:p={t:c*q**(n[j]-n[k]) for t,c in p.items()}
            else:w={t:c*pk(j,i,bar) for t,c in w.items()}
        else:
            if opposite:w={t:c/pk(j,i,bar) for t,c in w.items()}
            else:p={t:c*q**(-n[j]+n[k]) for t,c in p.items()}
        return add(p,w)
    def matrix(s,kind):
        i,n=s
        if kind=='V':return {(i,n):q**(-n[i]),((i-1)%h,shift(n,i,(i-1)%h)):X}
        if kind=='bar':
            return {((i+d)%h,shift(n,(i+d)%h,i)):X**d*q**sum(n[(i+k)%h] for k in range(d+1)) for d in range(h)}
        if kind=='rev':return {(i,n):q**n[i],((i+1)%h,shift(n,(i+1)%h,i)):X*q**(n[i]-n[(i+1)%h])}
        if kind=='plusbar':return {(i,n):q**(-n[i]),((i-1)%h,shift(n,(i-1)%h,i)):Xp*q**(n[i]-n[(i-1)%h]-1)}
        raise ValueError(kind)
    for bar in (False,True):
        for kind in ('f','e'):
            fac=1/q if kind=='f' else q
            for j in range(h):
                mat=S.zeros(h)
                for i in range(h):
                    n=shift(n0,0,i) if not bar else shift(n0,i,0)
                    st=(i,n)
                    op1=lambda s:copro(kind,j,s,bar)
                    op2=lambda s:copro(kind,(j+1)%h,s,bar)
                    o=difference(linear(op1,op2(st)),linear(op2,op1(st)),fac)
                    for (ii,nn),c in o.items():mat[ii,i]+=c
                line_image=(kind=='f' and not bar) or (kind=='e' and bar)
                assert mat.rank()==1
                if line_image: assert all(mat[k,:]==S.zeros(1,h) for k in range(h) if k!=(j+2)%h)
                else: assert all(mat[:,k]==S.zeros(h,1) for k in range(h) if k!=(j+2)%h)
    for kind,bar,gen,orient in [('V',False,'f',1),('bar',True,'f',1),('rev',True,'f',-1),('plusbar',True,'e',1)]:
        for j in range(h):
            for i in range(h):
                st=(i,n0)
                lhs=linear(lambda s:matrix(s,kind),copro(gen,j,st,bar,orient=orient))
                rhs=linear(lambda s:copro(gen,j,s,bar,opposite=True,orient=orient),matrix(st,kind))
                assert not difference(lhs,rhs),(h,kind,j,i)
    # The opposite Borel on the second module, checked on every exact reachable label.
    t=y/X
    def Y(j,s):
        i,n=s; return {(i,shift(n,j,(j+1)%h)):t/(q-1/q)*q**(-n[j])}
    def K(j,s):
        i,n=s; return {(i,shift(shift(n,j,(j+1)%h),j,(j-1)%h)):S.S(1)}
    for bar,kind in [(False,'V'),(True,'bar')]:
        for j in range(h):
            for i in range(h):
                st=(i,n0)
                left_input=add(particle('e',j,st,bar),{s:c*pk(j,i,bar) for s,c in Y(j,st).items()})
                lhs=linear(lambda s:matrix(s,kind),left_input)
                def out(s): return add(Y(j,s),linear(lambda ss:particle('e',(j-1)%h,ss,bar),K(j,s)))
                rhs=linear(out,matrix(st,kind))
                assert not difference(lhs,rhs),(h,kind,'second module',j,i)
    print(f'h={h}: Omega image/kernel; four displayed intertwiners; second module on V and Vbar: exact PASS',flush=True)

# Arbitrary diagonal variables (product one): sector determinants and residue factorization.
X,q=S.symbols('X q',nonzero=True)
for h in range(3,7):
    z=list(S.symbols('z0:'+str(h-1),nonzero=True)); z.append(1/S.prod(z))
    # Use unprojected Q with sum +1/-1: diag factors a_i have product q^(h-1).
    a=[q*zz for zz in z]
    a[-1]/=q
    MV=S.diag(*a)
    for i in range(h):MV[(i-1)%h,i]=X
    assert S.factor(MV.det()-q**(h-1)-(-1)**(h-1)*X**h)==0
    # General Vbar sector: total charge projected to sum zero; choose q^Q_i=z_i.
    # Remove the common q^(-1/h); avoiding fractional powers, use b_i product=q^(h-1).
    b=[q*zz for zz in z]; b[-1]/=q
    U=S.zeros(h)
    for i in range(h):U[(i+1)%h,i]=b[(i+1)%h]/q
    A=S.diag(*b)
    M=sum(((X*U)**j*A for j in range(h)),S.zeros(h))
    assert S.simplify(U**h-S.eye(h)/q)==S.zeros(h)
    # det via (1-XU) M=(1-X^h/q) A; det(1-XU)=1-X^h/q.
    assert S.simplify((S.eye(h)-X*U)*M-(1-X**h/q)*A)==S.zeros(h)
    print(f'h={h}: symbolic cyclic inverse identity; det Mbar=(q-X^h)^(h-1)',flush=True)
print('PASS. Generic-r proofs and first-order universal-R calculations are written in the report.')
