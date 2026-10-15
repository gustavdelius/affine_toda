"""Twisted-loop centralizers from explicit matrix commutators; no book code."""
import numpy as np
import sympy as S

def aeven(N):
    z=S.symbols('z',nonzero=True);E=S.zeros(N)
    for i in range(N-1):E[i,i+1]=1
    E[N-1,0]=z
    J=S.zeros(N)
    for i in range(N):J[i,N-1-i]=(-1)**i
    sigma=lambda X:-J*X.T*J.inv()
    assert sigma(E)==E.subs(z,-z)
    assert E**N==z*S.eye(N)
    degrees=[]
    for r in range(1,N):
        m=(r+1)%2;X=z**m*E**r
        assert sigma(X)==X.subs(z,-z) and E*X-X*E==S.zeros(N) and S.trace(X)==0
        degrees.append((m*N+r)%(2*N))
    print(f'PASS A{N-1}^(2) in sl({N}): E^{N}=z I, centralizer basis z^m E^r; residues={sorted(degrees)} mod {2*N}')

def d4triality():
    N=8
    def unit(i,j):
        a=np.zeros((N,N));a[i,j]=1;return a
    es=[unit(i,i+1)-unit(i+5,i+4) for i in range(3)]
    es.append(unit(2,7)-unit(3,6))
    fs=[x.T for x in es];hs=[e@f-f@e for e,f in zip(es,fs)]
    perm=[2,1,3,0];gen=es+fs+hs;img=[es[p] for p in perm]+[fs[p] for p in perm]+[hs[p] for p in perm]
    basis=[];images=[];Q=np.zeros((64,0));todo=[]
    def add(x,y):
        nonlocal Q
        v=x.ravel();v=v-Q@(Q.T@v)
        if np.linalg.norm(v)<1e-9:return
        Q=np.column_stack([Q,v/np.linalg.norm(v)]);basis.append(x);images.append(y);todo.append(len(basis)-1)
    for x,y in zip(gen,img):add(x,y)
    while todo:
        i=todo.pop(0)
        for x,y in zip(gen,img):add(x@basis[i]-basis[i]@x,y@images[i]-images[i]@y)
    assert len(basis)==28
    B=np.column_stack([x.ravel() for x in basis]);inv=np.linalg.pinv(B)
    sig=inv@np.column_stack([x.ravel() for x in images]);assert np.allclose(np.linalg.matrix_power(sig,3),np.eye(28))
    D=np.diag([3,2,1,0,-3,-2,-1,0])
    grade=inv@np.column_stack([(D@x-x@D).ravel() for x in basis])
    assert np.allclose(grade,np.diag(np.diag(grade)))
    grades=np.rint(np.diag(grade)).astype(int);omega=np.exp(2j*np.pi/3)
    def null(M):
        u,s,vh=np.linalg.svd(M,full_matrices=True);rank=sum(s>1e-8);return vh[rank:].conj().T
    ids=np.where(grades==-3)[0];V=null(sig[:,ids]-omega*np.eye(28)[:,ids]);assert V.shape[1]==1
    e0=(B[:,ids]@V[:,0]).reshape(N,N);ef=sum(es)
    assert np.linalg.norm(ef@e0-e0@ef)>1e-5
    dims=[]
    for d in range(12):
        candidates=[]
        for m in range(-2,5):
            ids=np.where(grades==d-4*m)[0]
            if not len(ids):continue
            V=null(sig[:,ids]-omega**m*np.eye(28)[:,ids])
            for v in V.T:candidates.append((m,(B[:,ids]@v).reshape(N,N)))
        if not candidates:dims.append(0);continue
        mmin=min(m for m,x in candidates);mmax=max(m for m,x in candidates)+1
        cols=[]
        for m,x in candidates:
            comm={m:ef@x-x@ef,m+1:e0@x-x@e0}
            cols.append(np.concatenate([comm.get(j,np.zeros((N,N))).ravel() for j in range(mmin,mmax+1)]))
        dims.append(null(np.column_stack(cols)).shape[1])
    assert dims==[0,1,0,0,0,1,0,1,0,0,0,1]
    print('PASS D4^(3) in twisted so(8) loop: sigma^3=1; degree(lambda)=4; centralizer dimensions degrees 0..11 =',dims)

aeven(3);aeven(5);d4triality()
