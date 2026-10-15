"""Independent checks from Euclidean simple roots and reflection closure.
No imports of book code. Run from any directory; output is written by shell.
"""
import itertools, math
from functools import reduce, lru_cache
import numpy as np
import sympy as S

def cartan(G):
    return S.Matrix(G.rows,G.cols,lambda i,j:2*G[i,j]/G[i,i])

@lru_cache(None)
def finite(t,n):
    if t in 'ABCD':
        dim=n+1 if t=='A' else n
        e=S.eye(dim)
        rows=[e[i,:]-e[i+1,:] for i in range(n-1)]
        if t=='A': rows.append(e[n-1,:]-e[n,:])
        if t=='B': rows.append(e[n-1,:])
        if t=='C': rows.append(2*e[n-1,:])
        if t=='D': rows.append(e[n-2,:]+e[n-1,:])
        R=S.Matrix.vstack(*rows); G=R*R.T
        if t=='C': G=G/2
    elif t=='F':
        R=S.Matrix([[0,1,-1,0],[0,0,1,-1],[0,0,0,1],[S.Rational(1,2),-S.Rational(1,2),-S.Rational(1,2),-S.Rational(1,2)]])
        G=R*R.T
    elif t=='G': G=S.Matrix([[S.Rational(2,3),-1],[-1,2]])
    elif t=='E':
        G=2*S.eye(n); idx=1
        for length in [n-4,2,1]:
            prev=0
            for _ in range(length):
                G[prev,idx]=G[idx,prev]=-1;prev=idx;idx+=1
    A=np.array(cartan(G),dtype=int)
    roots={tuple(int(i==j) for i in range(n)) for j in range(n)}
    todo=list(roots)
    while todo:
        b=todo.pop()
        for i in range(n):
            c=list(b);c[i]-=int(A[i]@b);c=tuple(c)
            if c not in roots: roots.add(c);todo.append(c)
    high=max(roots,key=sum)
    R=S.Matrix.vstack(-S.Matrix([high]),S.eye(n))
    Ga=R*G*R.T
    return G,sorted(roots),Ga,S.Matrix([1,*high])

def isomorphisms(A,B):
    if A.shape!=B.shape: return
    n=A.rows
    sig=lambda M,i:(M[i,i],tuple(sorted(M[i,j] for j in range(n))),tuple(sorted(M[j,i] for j in range(n))))
    opts={i:[j for j in range(n) if sig(A,i)==sig(B,j)] for i in range(n)}
    order=sorted(range(n),key=lambda i:len(opts[i])); p={}; used=set()
    def rec(k):
        if k==n: yield tuple(p[i] for i in range(n));return
        i=order[k]
        for j in opts[i]:
            if j in used or any(A[i,a]!=B[j,b] or A[a,i]!=B[b,j] for a,b in p.items()):continue
            p[i]=j;used.add(j);yield from rec(k+1);used.remove(j);del p[i]
    yield from rec(0)

def orbits(p):
    unseen=set(range(len(p)));out=[]
    while unseen:
        x=min(unseen);o=[]
        while x not in o: o.append(x);unseen.remove(x);x=p[x]
        out.append(o)
    return out

def primitive(v):
    v=v*S.ilcm(*[x.q for x in v]);return v/reduce(math.gcd,[int(x) for x in v])

def folded(Ga,marks,p):
    oo=orbits(p);P=S.zeros(len(oo),Ga.rows)
    for i,o in enumerate(oo):
        for j in o:P[i,j]=S.Rational(1,len(o))
    G=P*Ga*P.T;raw=S.Matrix([sum(marks[j] for j in o) for o in oo]);n=primitive(raw)
    return oo,G,raw,n

def target(t,n):
    if t=='Ae':
        R=-2*S.eye(n)[0,:]
        if n>1:R=R.col_join(S.Matrix.vstack(*[S.eye(n)[i,:]-S.eye(n)[i+1,:] for i in range(n-1)]))
        R=R.col_join(S.eye(n)[n-1,:]);return cartan(R*R.T/2)
    dual={'Ao':'B','Dt':'C','Et':'F','Gt':'G'}
    return cartan(finite(dual.get(t,t),n)[2]).T if t in dual else cartan(finite(t,n)[2])

def finite_checks():
    families=[('A',n) for n in range(1,7)]+[(t,n) for t in 'BC' for n in range(2,7)]+[('D',n) for n in range(4,9)]+[('E',n) for n in [6,7,8]]+[('F',4),('G',2)]
    for t,n in families:
        G,roots,Ga,marks=finite(t,n);A=cartan(G);h=int(sum(marks));dual=S.Matrix([marks[i]*Ga[i,i]/2 for i in range(n+1)])
        # Euclidean representation, independent of author numbering.
        R=np.linalg.cholesky(np.array(G,float));W=np.eye(n)
        for a in R: W=W@(np.eye(n)-2*np.outer(a,a)/(a@a))
        exps=sorted(np.rint(np.mod(np.angle(np.linalg.eigvals(W)),2*np.pi)*h/(2*np.pi)).astype(int))
        assert np.allclose(np.sort(np.linalg.eigvals(np.array(A,float))),np.sort([4*np.sin(np.pi*s/(2*h))**2 for s in exps]))
        order=next(k for k in range(1,h+1) if np.allclose(np.linalg.matrix_power(W,k),np.eye(n)))
        rr=np.array(roots)@R;unused=set(range(len(rr)));lens=[]
        while unused:
            i=next(iter(unused));x=rr[i];o=set()
            while i not in o:
                o.add(i);x=W@x;i=int(np.argmin(np.linalg.norm(rr-x,axis=1)))
            unused-=o;lens.append(len(o))
        weights=np.linalg.inv(2*R/np.diag(np.array(G,float))[:,None]).T
        rho=weights.sum(axis=0);half=sum((np.array(b)@R for b in roots if min(b)>=0),np.zeros(n))/2
        theta=-np.array(finite(t,n)[3][1:],float).reshape(n)@R
        cas=theta@(theta-2*rho) # theta here negative highest root
        assert len(roots)==n*h and order==h and lens==[h]*n and np.allclose(rho,half) and np.isclose(cas,2*float(sum(dual)))
        assert cartan(Ga)*marks==S.zeros(n+1,1) and cartan(Ga).T*dual==S.zeros(n+1,1)
        print(f'PASS finite {t}{n}: dim={len(roots)+n} marks={list(marks)} h={h} hdual={sum(dual)} exp={exps} Coxeter_orbits={n}x{h}')
        Aa=cartan(Ga);aa=list(isomorphisms(Aa,Aa));stab=[p for p in aa if p[0]==0]
        special=[i for i,v in enumerate(marks) if v==1]
        assert len(aa)==len(stab)*len(special) and {p[0] for p in aa}==set(special)
        # Find a normal abelian regular subgroup with the expected centre structure.
        ident=tuple(range(n+1));compose=lambda p,q:tuple(p[q[i]] for i in range(n+1))
        def subgroup(gs):
            out={ident};todo=[ident]
            while todo:
                p=todo.pop()
                for q in gs:
                    v=compose(p,q)
                    if v not in out:out.add(v);todo.append(v)
            return out
        reg=None
        for count in [0,1,2]:
            for gs in itertools.combinations(aa,count):
                H=subgroup(gs)
                if len(H)!=len(special) or len({p[0] for p in H})!=len(H):continue
                if any(compose(p,q)!=compose(q,p) for p in H for q in H):continue
                if t=='D' and n%2==0 and any(compose(p,p)!=ident for p in H):continue
                if any({compose(p,q) for q in H}!={compose(q,p) for q in H} for p in aa):continue
                reg=H;break
            if reg is not None:break
        assert reg is not None
        orders=sorted(len(subgroup([p])) for p in reg)
        print(f'PASS symmetry {t}{n}: Aut={len(aa)} finite_stabilizer={len(stab)} centre={len(reg)} element_orders={orders}')
        if n>1:
            ev,vv=np.linalg.eig(W.T);v=vv[:,np.argmin(abs(ev-np.exp(2j*np.pi/h)))]
            radii={}
            for root in rr:
                length=round(float(root@root),8)
                radii.setdefault(length,set()).add(round(float(abs(root@v)),8))
            if t in ['B','C','F','G']:print(f'INFO Coxeter radii {t}{n} by root length2: {radii}')

def fold_checks():
    specs=[]
    for n in range(2,5):specs += [('B',n,'D',n+1,1),('C',n,'A',2*n-1,1),('Ao',n,'D',2*n,2),('Dt',n,'D',n+2,2)]
    specs += [('F',4,'E',6,1),('G',2,'D',4,1),('Et',4,'E',7,2),('Gt',2,'E',6,3)]
    for n in [1,2,3]:specs += [('Ae',n,'A',2*n,2),('Ae',n,'D',2*n+2,2)]
    for t,n,pt,pn,k in specs:
        _,_,Ga,marks=finite(pt,pn);T=target(t,n);matches=[]
        for p in isomorphisms(Ga,Ga):
            if len(orbits(p))!=n+1:continue
            oo,G,raw,nn=folded(Ga,marks,p);A=cartan(G)
            if next(isomorphisms(A,T),None) is not None:matches.append((p,oo,G,raw,nn))
        assert matches,(t,n,pt,pn)
        p,oo,G,raw,nn=matches[0];A=cartan(G)
        Gfinite=finite(pt,pn)[0];Rp=S.Matrix.vstack(-S.Matrix([list(marks)[1:]]),S.eye(pn))
        Rp=np.array(Rp,float)@np.linalg.cholesky(np.array(Gfinite,float))
        Rf=np.array([Rp[o].mean(axis=0) for o in oo]);U,sv,_=np.linalg.svd(Rf.T,full_matrices=False);U=U[:,sv>1e-9]
        Mp=Rp.T@np.diag(np.array(marks,float).ravel())@Rp
        Mf=Rf.T@np.diag(np.array(raw,float).ravel())@Rf
        assert np.allclose(U.T@Mp@U,U.T@Mf@U)
        dn=S.Matrix([k*nn[i]*G[i,i]/2 for i in range(n+1)])
        assert dn==primitive(A.T.nullspace()[0])
        trans=next(isomorphisms(A,T.T),None) is not None
        print(f'PASS fold {pt}{pn}->{t}{n}: orbits={oo} lengths={list(G.diagonal())} raw={list(raw)} marks={list(nn)} h={sum(nn)} dual={list(dn)} hdual={sum(dn)} target=True transpose_isomorphic={trans}')
    for n in [1,2,3]:
        for pt,pn in [('A',2*n),('D',2*n+2)]:
            _,_,Ga,marks=finite(pt,pn)
            if pt=='A':p=tuple((-j)%(pn+1) for j in range(pn+1))
            else:
                p=list(range(pn+1));p[0]=pn-1;p[pn-1]=1;p[1]=pn;p[pn]=0
                for j in range(2,pn-1):p[j]=pn-j
                p=tuple(p)
            assert all(Ga[i,j]==Ga[p[i],p[j]] for i in range(pn+1) for j in range(pn+1))
            oo,G,raw,nn=folded(Ga,marks,p);fac=int(raw[0]/nn[0])
            # Nonzero eigenvalues of R^T diag(n) R = diag(sqrt(n)) G diag(sqrt(n)).
            mass=np.linalg.eigvalsh(np.diag(np.sqrt(np.array(nn,float).ravel()))@np.array(G,float)@np.diag(np.sqrt(np.array(nn,float).ravel())))[1:]
            print(f'PASS Aeven n={n} parent={pt}{pn}: permutation={p} orbits={oo} raw={list(raw)} gcd={fac} standard_mass2={np.round(mass,10).tolist()} parent_restricted_mass2={np.round(fac*mass,10).tolist()}')

def charges():
    for t,n in [('A',4),('D',5),('E',6),('D',4),('D',6),('E',8)]:
        G,roots,_,marks=finite(t,n);A=np.array(cartan(G),float);h=int(sum(marks));R=np.linalg.cholesky(np.array(G,float))
        colours={0:0};todo=[0]
        while todo:
            i=todo.pop()
            for j in range(n):
                if A[i,j]<0 and j not in colours:colours[j]=1-colours[i];todo.append(j)
        W=np.eye(n)
        for c in [0,1]:
            for i in range(n):
                if colours[i]==c:W=W@(np.eye(n)-np.outer(R[i],R[i]))
        eig,v=np.linalg.eig(W.T);v=v[:,np.argmin(abs(eig-np.exp(2j*np.pi/h)))]
        gamma=np.array([(-1)**colours[i]*R[i] for i in range(n)])
        base=gamma@v;rr=[];labels=[];phase=[]
        for i in range(n):
            x=gamma[i].copy()
            for j in range(h):
                rr.append(x);labels.append(i);phase.append(np.angle(x@v));x=W@x
        rr=np.array(rr);triples=[]
        for i in range(len(rr)):
            for j in range(i,len(rr)):
                ds=np.linalg.norm(rr+rr[i]+rr[j],axis=1);ks=np.where(ds<1e-7)[0]
                triples.extend((i,j,int(k)) for k in ks if k>=j)
        nullspaces=[]
        for spin in ([1,h+1] if t!='D' or n%2 else [n-1,n-1+h]):
            rows=[]
            for tri in triples:
                row=np.zeros(n,dtype=complex)
                for j in tri:row[labels[j]]+=np.exp(1j*spin*phase[j])
                rows.append(row)
            M=np.array(rows);U,sv,Vh=np.linalg.svd(np.vstack([M.real,M.imag]),full_matrices=False)
            null=Vh[sv<1e-8].T;lam=4*np.sin(np.pi*spin/(2*h))**2
            residual=np.linalg.norm(A@null-lam*null)
            assert null.shape[1] and residual<1e-8
            nullspaces.append(null)
            print(f'PASS fusion {t}{n} s={spin} h={h}: root_triples={len(triples)} conservation_nullity={null.shape[1]} Cartan_residual={residual:.2e}')
        J=np.diag([(-1)**colours[i] for i in range(n)])
        err=np.linalg.norm(nullspaces[1]@nullspaces[1].T-J@nullspaces[0]@nullspaces[0].T@J)
        assert err<1e-8
        print(f'PASS fusion colour flip {t}{n}: projector_error={err:.2e}')
        # E8 orbit radii are precisely PF entries in the same species ordering.
        vals,vec=np.linalg.eigh(A);pf=abs(vec[:,0]);radii=abs(base)
        assert np.ptp(radii/pf)<1e-10
        if t=='E' and n==8:
            ratios=np.sort(pf)/min(pf)
            print('PASS E8 radii/PF and masses:',np.round(ratios,10).tolist(),'golden_error=',abs(ratios[1]-2*np.cos(np.pi/5)))

def duality():
    b=S.symbols('b',positive=True);B=lambda x:2*x/(x+4*S.pi)
    for name,k,a,ad,end,delta in [('C2/D3t',2,2,1,6,1),('B2/A3t',2,1,2,3,-S.Rational(1,2)),('G2/D4t',3,3,1,12,3)]:
        bd=16*S.pi**2*k/b;res=S.simplify(2-B(b/a)-B(bd/ad))
        hu=end-2*delta+delta*B(b/a);hd=end-delta*B(bd/ad)
        assert res==0 and S.simplify(hu-hd)==0
        print(f'PASS duality {name}: B residual={res} H residual={S.simplify(hu-hd)}; beta2=sqrt(2): B={float(B(S.sqrt(2)/a)):.10f} Bdual={float(B(bd.subs(b,S.sqrt(2))/ad)):.10f}')
    # C2 and B2 are isomorphic; their H parameters do not use identical mass parametrizations.
    print('NOTE duality checks verify the stated B conjecture algebraically, not its all-orders physical validity.')

if __name__=='__main__':
    finite_checks();fold_checks();charges();duality()
