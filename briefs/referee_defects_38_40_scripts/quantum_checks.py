"""Independent algebraic and finite-chain checks for chapter 40.
No imports from code. Infinite auxiliary paths are evaluated exactly
at each site; only the absolutely convergent Fock trace is cut off.
"""
import itertools
import numpy as np
import sympy as s

def kt_check():
    q,x=s.symbols('q x',nonzero=True)
    E=s.Matrix([[0,1],[0,0]]); F=E.T; K=s.diag(q,1/q); I=s.eye(2)
    # Solve the finite-dimensional intertwining equations independently.
    a,b,c,d=s.symbols('a b c d')
    R=s.Matrix([[1,0,0,0],[0,a,b,0],[0,c,d,0],[0,0,0,1]])
    equations=[]
    for e1,e2,k1,k2 in ((E/x,E,K,K),(F/x,F,K.inv(),K.inv())):
        delta=s.kronecker_product(e1,I)+s.kronecker_product(k1,e2)
        op=s.kronecker_product(I,e2)+s.kronecker_product(e1,k2)
        equations+=list(R*delta-op*R)
    solution=s.solve(equations,(a,b,c,d)); R=R.subs(solution)
    print('spin-half R from coproduct =',R,flush=True)
    ed=E*F-q**2*F*E; fd=E*F-q**-2*F*E
    ea=(E*ed-ed*E)/(q+1/q); e0=(ed*F-F*ed)/(q+1/q)
    fa=(fd*F-F*fd)/(q+1/q); f0=(E*fd-fd*E)/(q+1/q)
    kap=q-1/q; ident=s.eye(4); kron=s.kronecker_product
    real0=ident+kap*x*kron(E,F)
    real3=ident+kap*x**3*kron(f0,e0)
    imag=ident+kap*x*x*kron(fd,ed)/(q+1/q)
    reala3=ident+kap*x**3*kron(fa,ea)
    reala0=ident+kap*x*kron(F,E)
    # Divide Cartan factor by sqrt(q); ratios are all rational in q.
    cartan=s.diag(1,1/q,1/q,1)
    def trunc(m):return m.applyfunc(lambda e:s.series(e,x,0,4).removeO().expand())
    raw=trunc(cartan*real0*real3*imag*reala3*reala0)
    normal=trunc(raw/raw[0,0])
    residual=trunc(normal-R).applyfunc(s.simplify)
    print('KT through x^3, all coefficients zero:',residual==s.zeros(4),flush=True)
    assert residual==s.zeros(4)
    for label,m in [('no imaginary',cartan*real0*real3*reala3*reala0),('no height 3',cartan*real0*imag*reala0)]:
        err=trunc(m/m[0,0]-R).subs(q,s.Rational(7,5))
        powers=sorted({term[0] for e in err for term,coef in s.Poly(e,x).terms() if coef!=0})
        print(label,'first nonzero order',min(powers),flush=True)
    return R

def r_matrix(q,z):
    a=q*(z*z-1)/(q*q*z*z-1)
    b=(q*q-1)*z/(q*q*z*z-1)
    return np.array([[1,0,0,0],[0,a,b,0],[0,b,a,0],[0,0,0,1]],complex)

def local(q,z,offset,a,n,module='weyl',inverse=True):
    # a=+1/-1; z=y_particle/y_defect. Return b,n',amplitude.
    if module=='weyl':
        if inverse:
            den=1+z*z
            return [(a,n,q**(a*(n+offset/2))/den),(-a,n+a,1j*z*q**(-a*(n+a+offset/2))/den)]
        return [(a,n,q**(-a*(n+offset/2))),(-a,n+a,-1j*z*q**(-a*(n+offset/2)))]
    x=1/z; den=1+x*x
    if a==1:
        return [(1,n,q**(-n-offset/2)*(1+x*x*q**(2*n+2))/den),(-1,n+1,-1j*x*q**(-n-offset/2)/den)]
    out=[(-1,n,q**(n+offset/2)/den)]
    if n>0:out.append((1,n-1,-1j*x*q**(n+offset/2)*(1-q**(2*n))/den))
    return out

def d_operator(q,ys,yd,offset=.3,n0=0,module='weyl',inverse=True):
    N=len(ys); dim=2**N; out=np.zeros((dim,dim),complex)
    for col in range(dim):
        spins=[1-2*((col>>(N-1-j))&1) for j in range(N)]
        paths={(0,n0):1.+0j}
        for j,a in enumerate(spins):
            nxt={}
            for (bits,n),ampl in paths.items():
                for b,k,c in local(q,ys[j]/yd,offset,a,n,module,inverse):
                    key=(2*bits+(1-b)//2,k)
                    nxt[key]=nxt.get(key,0)+ampl*c
            paths=nxt
        for (row,n),ampl in paths.items():
            if n==n0:out[row,col]+=ampl
    return out

def transfer(q,ys,y,omega):
    N=len(ys); dim=2**N; out=np.zeros((dim,dim),complex)
    for col in range(dim):
        bits=[(col>>(N-1-j))&1 for j in range(N)]
        for aux0 in (0,1):
            paths={(0,aux0):1.+0j}
            for j,a in enumerate(bits):
                R=r_matrix(q,y/ys[j]); nxt={}
                for (row,aux),ampl in paths.items():
                    for b,ax in itertools.product((0,1),repeat=2):
                        c=R[2*ax+b,2*aux+a]
                        if abs(c):
                            key=(2*row+b,ax); nxt[key]=nxt.get(key,0)+c*ampl
                paths=nxt
            for (row,aux),ampl in paths.items():
                if aux==aux0:out[row,col]+=omega**(1-2*aux)*ampl
    return out

def relative(a,b):return np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b))

def chain_checks():
    for angle in (.71,-np.pi*.37):
        q=np.exp(1j*angle); offset=.3
        for N in (2,3,4,5,6):
            ys=np.exp(np.linspace(-.43,.51,N)+.09*np.sin(np.arange(N)))
            yd=.79
            D=d_operator(q,ys,yd,offset); Dp=d_operator(q,ys,q*yd,offset); Dm=d_operator(q,ys,yd/q,offset)
            other=d_operator(q,ys,1.21,offset)
            charges=np.array([N-2*i.bit_count() for i in range(2**N)])
            fac=np.diag(q**(2*charges))@D
            factor=relative(d_operator(q,ys,yd,offset,n0=2),fac)
            tq=comm=badcomm=0.
            coef=np.prod(q**.5*(q*(ys/yd)**2+1/q)/((ys/yd)**2+1))
            for m in np.unique(charges):
                ids=np.flatnonzero(charges==m); sl=np.ix_(ids,ids); om=q**(-m/2)
                t=transfer(q,ys,1j*yd,om)
                tq=max(tq,relative((t@D)[sl],(q**(-N/2)*Dp+coef*Dm)[sl]))
                t=transfer(q,ys,1.37,om)
                comm=max(comm,relative((t@D)[sl],(D@t)[sl]))
                if m:
                    tb=transfer(q,ys,1.37,1)
                    badcomm=max(badcomm,relative((tb@D)[sl],(D@tb)[sl]))
            dd=relative(D@other,other@D)
            direct=d_operator(q,ys,yd,offset,inverse=False)
            ids=np.flatnonzero(charges==(N%2)); sl=np.ix_(ids,ids)
            t=transfer(q,ys,1.37,q**(-(N%2)/2))
            badinverse=relative((t@direct)[sl],(direct@t)[sl])
            print(f'angle={angle:.4f} N={N}: factor={factor:.1e} TQ={tq:.1e} [t,D]={comm:.1e} [D,D]={dd:.1e}; wrong twist={badcomm:.3g}, direct L={badinverse:.3g}',flush=True)
            assert max(factor,tq,comm,dd)<1e-10
        # convergent Fock trace and its TQ identity
        ys=np.array([.81,1.07,1.32]); yd=.79; omega=.63
        def fock(yd,off,M):
            return sum(omega**(2*n+off)*d_operator(q,ys,yd,off,n,module='fock') for n in range(M+1))
        Q=fock(yd,offset,45); Q2=fock(yd,offset,65)
        T=transfer(q,ys,1.29,omega)
        A=q**(len(ys)/2)
        C=np.prod((q*(ys/yd)**2+1/q)/(q**.5*((ys/yd)**2+1)))
        tq=relative(transfer(q,ys,1j*yd,omega)@Q,A*fock(q*yd,offset+1,65)+C*fock(yd/q,offset-1,65))
        print('Fock: cutoff=',relative(Q,Q2),'comm=',relative(T@Q,Q@T),'TQ=',tq,flush=True)
        assert max(relative(Q,Q2),relative(T@Q,Q@T),tq)<1e-10
    # N=1 and extremal sectors are scalars and commute at every twist.
    q=np.exp(.71j); D=d_operator(q,[.9],.8); t=transfer(q,[.9],1.2,.37)
    print('N=1 counterexample to ONLY sector twist:',np.linalg.norm(t@D-D@t),flush=True)

def root_unity_counterexample():
    # The theorem assumes only q^4!=1, allowing q a primitive sixth root.
    # S^3 commutes with k and with both constant Weyl shifts, producing an
    # independent intertwiner (S^3 tensor 1) T.
    q=np.exp(1j*np.pi/3)
    print('q sixth root: |q^4-1|=',abs(q**4-1),'|q^6-1|=',abs(q**6-1),flush=True)
    print('Weyl shift S^3 commutator with k:',max(abs(q**(2*(n+3))-q**(2*n)) for n in range(-5,6)),flush=True)
    print('S^3 commutes with f0,f1 exactly (constant shifts); T and (S^3 x 1)T are independent.',flush=True)

if __name__=='__main__':
    kt_check()
    root_unity_counterexample()
    chain_checks()
