"""Solve ch43 doubled equations with both Y and kappa unknown.

No truncated transmission matrices: each input is propagated exactly and every
reached output equation retained. Unknown coefficients at every reached label
are shared between rapidities and between V and Vbar.
"""
import itertools
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr
from independent_breathers import transmit

def add(N,v): return tuple(np.array(N)+v)

def assemble(h,kind,rotation,node,extent=1):
    w=np.sqrt(2); q=-np.exp(-1j*np.pi*w); td=.27
    offset=np.linspace(-.231,.231,h)
    alpha=np.eye(h,dtype=int)[(node-1)%h]-np.eye(h,dtype=int)[node]
    jp=(node+rotation)%h
    alphap=np.eye(h,dtype=int)[(jp-1)%h]-np.eye(h,dtype=int)[jp]
    wy=alpha if kind=='e' else -alphap
    wk=alpha-alphap
    variables={}; rows=[];cols=[];vals=[];rhs=[]
    labels=[tuple(a)+(-sum(a),) for a in itertools.product(range(-extent,extent+1),repeat=h-1)]
    for conj in (False,True):
        sign=-1 if conj else 1
        def generator(j,i,t,typ):
            a,b=(j-1)%h,j
            if typ=='k': return i,q**(sign*((i==a)-(i==b)))
            src,dst=(b,a) if typ=='e' else (a,b)
            if conj: src,dst=dst,src
            return (dst,np.exp((1 if typ=='e' else -1)*w*t)) if i==src else (i,0)
        for t in (-.32,.41,.82+.11j):
            for N in labels:
                for i in range(h):
                    eq=defaultdict(lambda:defaultdict(complex))
                    def put(j,n,key,c):
                        if c: eq[(j,n)][key]+=c
                    if kind=='e':
                        ii,c=generator(node,i,t,'e')
                        for j,n,b in transmit(ii,N,offset,t,h,w,conj):put(j,n,None,c*b)
                        _,k=generator(node,i,t,'k')
                        for j,n,b in transmit(i,add(N,wy),offset,t,h,w,conj):put(j,n,('Y',N),k*b)
                        for j,n,b in transmit(i,N,offset,t,h,w,conj):
                            put(j,add(n,wy),('Y',n),-b)
                            jj,c=generator(jp,j,t,'e');put(jj,add(n,wk),('K',n),-b*c)
                    else:
                        ii,c=generator(node,i,t,'f')
                        for j,n,b in transmit(ii,add(N,wk),offset,t,h,w,conj):put(j,n,('K',N),c*b)
                        for j,n,b in transmit(i,add(N,wy),offset,t,h,w,conj):put(j,n,('Y',N),b)
                        for j,n,b in transmit(i,N,offset,t,h,w,conj):
                            _,k=generator(jp,j,t,'k');put(j,add(n,wy),('Y',n),-b/k)
                            jj,c=generator(jp,j,t,'f');put(jj,n,None,-b*c)
                    for terms in eq.values():
                        row=len(rhs);rhs.append(-terms.pop(None,0))
                        for key,c in terms.items():
                            if key not in variables:variables[key]=len(variables)
                            rows.append(row);cols.append(variables[key]);vals.append(c)
    A=coo_matrix((vals,(rows,cols)),shape=(len(rhs),len(variables))).tocsr()
    b=np.array(rhs)
    fit=lsmr(A,b,atol=3e-14,btol=3e-14,maxiter=3000)
    residual=np.linalg.norm(A@fit[0]-b)/np.linalg.norm(b)
    closed=None
    if (kind,rotation)==('e',-1):
        C=np.exp(w*(td+1j*np.pi/2))/(q-1/q)
        solution=np.array([1 if typ=='K' else C*q**(-N[(node-1)%h]-offset[(node-1)%h]) for typ,N in variables])
        closed=np.linalg.norm(A@solution-b)/np.linalg.norm(b)
    elif (kind,rotation)==('f',0):
        mu=-np.exp(-1j*np.pi*w/2)/(q-1/q)
        solution=np.array([q**(-np.dot(alpha,np.array(N)+offset)) if typ=='K' else mu*np.exp(-w*td)*q**(N[node]+offset[node]) for typ,N in variables])
        closed=np.linalg.norm(A@solution-b)/np.linalg.norm(b)
    return residual,closed,A.shape,fit[1]

def module_relations(h):
    q=-np.exp(-1j*np.pi*np.sqrt(2));offset=np.linspace(-.231,.231,h)
    alpha=[np.eye(h,dtype=int)[(j-1)%h]-np.eye(h,dtype=int)[j] for j in range(h)]
    kshift=[alpha[j]-alpha[(j-1)%h] for j in range(h)]
    def act(word,N):
        N=np.array(N);c=1
        for kind,j in reversed(word):
            if kind=='E':c*=q**(-N[(j-1)%h]-offset[(j-1)%h]);N=N+alpha[j]
            else:N=N+kshift[j]
        return tuple(N),c
    worst=0
    for N in (np.zeros(h,dtype=int),np.array([2,-1,-1]+[0]*(h-3))):
        for i,j in itertools.product(range(h),repeat=2):
            cartan=int(alpha[i]@alpha[j])
            x,c=act([('K',i),('E',j)],N);y,d=act([('E',j),('K',i)],N)
            assert x==y
            worst=max(worst,abs(c-q**cartan*d))
            if i!=j:
                n=1-cartan
                # q binomial for Serre: n=1 (nonadjacent), 2 (adjacent).
                coef=[1,1] if n==1 else [1,q+1/q,1]
                s=0
                for k in range(n+1):
                    _,c=act([('E',i)]*(n-k)+[('E',j)]+[('E',i)]*k,N)
                    s+=(-1)**k*coef[k]*c
                worst=max(worst,abs(s))
        assert np.all(sum(kshift)==0)
    return worst

if __name__=='__main__':
    for h in (3,4):
        for kind,rot in (('f',0),('e',-1),('e',0),('e',1),('f',-1),('f',1)):
            values=[assemble(h,kind,rot,j) for j in range(h)]
            print(f'h={h} {kind} rotation={rot}: LS residuals='+','.join(f'{v[0]:.3e}' for v in values)+
                  ' closed='+','.join('n/a' if v[1] is None else f'{v[1]:.3e}' for v in values)+f' shape={values[0][2]}',flush=True)
        print(f'h={h} second-module Cartan and Serre max residual={module_relations(h):.3e}',flush=True)
    for h in (5,6):print(f'h={h} second-module Cartan and Serre max residual={module_relations(h):.3e}',flush=True)
