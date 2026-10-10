# (1) one-loop particle mass shift (bubble only, normal-ordered scheme) for a_n^(1), d_4^(1)
# (2) Hollowood's bound-state sum (3.2) vs MacKay-Watts (4.29)
import numpy as np
from scipy.integrate import quad

def affine_simple_roots(alg, n):
    if alg=='a':   # a_n^(1), h=n+1, in R^{n+1} projected
        N=n+1; E=np.eye(N)
        roots=[E[i]-E[i+1] for i in range(n)]+[E[n]-E[0]]
        marks=[1]*(n+1)
        # orthonormal basis of the hyperplane sum=0
        Q,_=np.linalg.qr(np.vstack([np.ones(N)]+[E[i]-E[i+1] for i in range(n)]).T)
        B=Q[:,1:]
        roots=[B.T@r for r in roots]
        return roots,marks
    if alg=='d':   # d_n^(1)
        E=np.eye(n)
        roots=[E[i]-E[i+1] for i in range(n-1)]+[E[n-2]+E[n-1], -(E[0]+E[1])]
        marks=[1]+[2]*(n-3)+[1,1,1]
        return roots,marks

def I_bubble(p2, mb, mc):
    # Euclidean 2D: int d^2k/(2pi)^2 1/((k^2+mb^2)((k+p)^2+mc^2)), at p^2 = p2 (Euclidean)
    f=lambda x: 1.0/(x*mb**2+(1-x)*mc**2+x*(1-x)*p2)
    return quad(f,0,1,limit=200)[0]/(4*np.pi)

def particle_shift(alg,n,m=1.0,beta=1.0):
    roots,marks=affine_simple_roots(alg,n)
    M2=m**2*sum(nj*np.outer(a,a) for nj,a in zip(marks,roots))
    w,V=np.linalg.eigh(M2); r=len(w)
    masses=np.sqrt(w)
    # cubic couplings C_abc = m^2 beta sum_j n_j (a_j.e_a)(a_j.e_b)(a_j.e_c)
    P=np.array([[a@V[:,k] for k in range(r)] for a in roots])  # j x a
    C=m**2*beta*np.einsum('j,ja,jb,jc->abc',np.array(marks,float),P,P,P)
    # self-energy matrix at p^2=-m_a^2 within degenerate block: compute Sigma_ab(p^2)
    out=[]
    for a in range(r):
        p2=-masses[a]**2
        Sig=np.zeros((r,r))
        for c in range(r):
            for d in range(r):
                Sig+= -0.5*np.outer(C[:,c,d],C[:,c,d])*I_bubble(p2,masses[c],masses[d])
        # degenerate subspace of a
        deg=[b for b in range(r) if abs(masses[b]-masses[a])<1e-9]
        sub=Sig[np.ix_(deg,deg)]
        ev=np.linalg.eigvalsh(sub)
        out.append((masses[a], ev/masses[a]**2))
    return out

for alg,n,h in [('a',2,3),('a',3,4),('a',4,5),('a',5,6),('d',4,6),('d',5,8)]:
    res=particle_shift(alg,n)
    print(f"{alg}_{n}^(1) h={h}: predicted -beta^2/(4h) cot(pi/h) = {-1/(4*h)/np.tan(np.pi/h):.6f}")
    for ma,ev in res: print(f"   m_a={ma:.5f}  delta m^2/m^2 = {np.round(ev,6)}")

print("\nHollowood bound-state sum (3.2), in units of m, vs m_a cot and m_a cot/2 (n even)")
for n in [4,6,8,10]:
    for a in range(1,n//2+1):
        s=0.0
        for b in range(1,n):
            # net bound states: b<a, or n/2<=b<n/2+a (b-a<n/2), excluding threshold b=a+n/2
            if b<a or (n/2<=b<n/2+a):
                s+= np.sin(np.pi*b/n)*abs(np.sin(np.pi*(a-b)/n))
        ma=2*np.sin(np.pi*a/n)
        print(f" n={n} a={a}: sum_bound (1/2)omega = {s:.6f}, m_a cot = {ma/np.tan(np.pi/n):.6f}, m_a cot/2 = {ma/np.tan(np.pi/n)/2:.6f}")
