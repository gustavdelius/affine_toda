# One-loop (bubble, normal-ordered scheme) mass shifts for non-simply-laced theories in the book's
# normalization (long roots length^2 2, V = m^2/beta^2 sum n_j exp(beta alpha_j.phi)), and the implied
# floating Coxeter number H = h + dH*beta^2. Also classical masses of the twisted duals.
import numpy as np
from scipy.integrate import quad
def I_bubble(p2, mb, mc):
    f=lambda x: 1.0/(x*mb**2+(1-x)*mc**2+x*(1-x)*p2)
    return quad(f,0,1,limit=400,points=None)[0]/(4*np.pi)
def spectrum(roots,marks):
    M2=sum(nj*np.outer(a,a) for nj,a in zip(marks,roots))
    w,V=np.linalg.eigh(M2); return np.sqrt(np.abs(w)),V
def shifts(roots,marks,beta=1.0):
    ms,V=spectrum(roots,marks); r=len(ms)
    P=np.array([[a@V[:,k] for k in range(r)] for a in roots])
    C=beta*np.einsum('j,ja,jb,jc->abc',np.array(marks,float),P,P,P)
    d=[]
    for a in range(r):
        s=0.0
        for c in range(r):
            for e in range(r):
                if abs(C[a,c,e])>1e-12:
                    s+= -0.5*C[a,c,e]**2*I_bubble(-ms[a]**2,ms[c],ms[e])
        d.append(s)   # assumes nondegenerate masses (true below)
    return ms,np.array(d)   # delta m_a^2 per beta^2
E=lambda n: np.eye(n)
def bn(n):
    e=E(n); R=[e[i]-e[i+1] for i in range(n-1)]+[e[n-1], -(e[0]+e[1])]
    return R,[1]+[2]*(n-2)+[2,1]   # alpha_1 mark 1, alpha_2..alpha_{n-1} mark 2, alpha_n mark 2, alpha_0 mark 1
def cn(n):
    e=E(n); R=[(e[i]-e[i+1])/np.sqrt(2) for i in range(n-1)]+[np.sqrt(2)*e[n-1], -np.sqrt(2)*e[0]]
    return R,[2]*(n-1)+[1,1]
def g2():
    R=[np.array([1/np.sqrt(6),-1/np.sqrt(2)]), np.array([0,np.sqrt(2)]), np.array([-np.sqrt(1.5),-1/np.sqrt(2)])]
    return R,[3,2,1]
def check(name,R,marks,ratio_fn,h):
    for a,nj in zip(R,marks): pass
    assert np.allclose(sum(nj*a for nj,a in zip(marks,R)),0)
    ms,d=shifts(R,marks)
    order=np.argsort(ms); ms=ms[order]; d=d[order]
    print(name, "lengths^2:", [round(a@a,3) for a in R], "classical masses:", np.round(ms,5))
    # extract dH from each ratio via finite difference in H
    for (i,j,f) in ratio_fn:
        r2=ms[i]**2/ms[j]**2; dr2=r2*(d[i]/ms[i]**2-d[j]/ms[j]**2)
        eps=1e-6; dfdH=(f(h+eps)**2-f(h-eps)**2)/(2*eps)
        print(f"   ratio m{i+1}/m{j+1}: classical {np.sqrt(r2):.6f} vs formula {f(h):.6f};  dH/beta^2 = {dr2/dfdH:.6f}  (x 2pi = {2*np.pi*dr2/dfdH:.5f})")
for n in [3,4]:
    R,mk=bn(n)
    # b_n: m_a = 2 sin(pi a/H) * m_spinor ; spinor is lightest? order masses
    ms,_=spectrum(R,mk); print()
    srt=np.sort(ms); 
    # identify: spinor mass = 1*s, vectors 2 s sin(pi a/2n)
    fns=[]
    sp=np.argmin(np.abs(srt-np.sqrt(2)))  # spinor sqrt2 in book units
    vec=[i for i in range(n) if i!=sp]
    for k,i in enumerate(vec): fns.append((i,sp,(lambda a: (lambda H: 2*np.sin(np.pi*a/H)))(k+1)))
    check(f"b_{n}^(1)",R,mk,fns,2*n)
for n in [3,4]:
    R,mk=cn(n); print()
    fns=[(i,n-1,(lambda a:(lambda H: np.sin(np.pi*a/H)/np.sin(np.pi*n/H)))(i+1)) for i in range(n-1)]
    check(f"c_{n}^(1)",R,mk,fns,2*n)
print(); R,mk=g2(); check("g_2^(1)",R,mk,[(1,0,lambda H: 2*np.cos(np.pi/H))],6)
# classical masses of twisted duals via coroots
def dual(R,mk):
    Rv=[2*a/(a@a) for a in R]; mkv=[nj*(a@a)/2 for nj,a in zip(mk,R)]
    return Rv,mkv
print("\nclassical mass ratios of duals:")
for n in [3,4]:
    Rv,mv=dual(*bn(n)); ms,_=spectrum(Rv,mv); ms=np.sort(ms)
    print(f" a_{2*n-1}^(2): marks {mv}, masses/min {np.round(ms/ms[0],5)}; 2sin(pi a/(2n-1))*m_n:",np.round([2*np.sin(np.pi*a/(2*n-1)) for a in range(1,n)],5))
    Rv,mv=dual(*cn(n)); ms,_=spectrum(Rv,mv); ms=np.sort(ms)
    print(f" d_{n+1}^(2): masses/max {np.round(ms/ms[-1],5)}; sin(pi a/(2n+2))/sin(pi n/(2n+2)):",np.round([np.sin(np.pi*a/(2*n+2))/np.sin(np.pi*n/(2*n+2)) for a in range(1,n+1)],5))
Rv,mv=dual(*g2()); ms,_=spectrum(Rv,mv); ms=np.sort(ms)
print(" d_4^(3): m2/m1 =",ms[1]/ms[0]," 2cos(pi/12)=",2*np.cos(np.pi/12), " 2cos(pi/8)=",2*np.cos(np.pi/8))
