"""Independent small algebraic checks, exercises, and scalar asymptotics."""
import numpy as np
import sympy as s
from scipy.integrate import quad

def cft():
    rt=s.sqrt(2)
    S=s.Matrix([[1,1,rt],[1,1,-rt],[rt,-rt,0]])/2
    de=[s.simplify(S[1,j]/S[0,j]) for j in range(3)]
    ds=[s.simplify(S[2,j]/S[0,j]) for j in range(3)]
    print('Ising epsilon eigenvalues=',de,'sigma=',ds,flush=True)
    fusion=lambda a,i,j:s.simplify(sum(S[a,k]*S[i,k]*S[j,k]/S[0,k] for k in range(3)))
    labels=['1','epsilon','sigma']
    for a in (1,2):
        sectors=[(labels[i],labels[j]) for i in range(3) for j in range(3) if fusion(a,i,j)]
        print(labels[a]+'-twisted Hilbert sectors:',sectors,flush=True)
    # M(2,5) modular S in the order (vacuum, h=-1/5).
    SY=s.Matrix([[-s.sin(2*s.pi/5),s.sin(s.pi/5)],[s.sin(s.pi/5),s.sin(2*s.pi/5)]])*2/s.sqrt(5)
    ly=[s.simplify(SY[1,j]/SY[0,j]) for j in range(2)]
    print('Lee-Yang eigenvalues=',[x.evalf(12) for x in ly],'; each d^2-d-1=',[s.simplify(x*x-x-1) for x in ly],flush=True)
    # Local massless gluing equations give psiL'=phiL', psiR'=-phiR'.
    a,b,c,d=s.symbols('phiL phiR psiL psiR')
    print('massless gluing:',s.solve([a-b-c-d,c-d-a-b],(c,d)),flush=True)
    print('compact-boson same-radius lattice: (n,w)->(R^2 w,n/R^2); R^2=2 sends (1,0) to (0,1/2)',flush=True)

def classical_limit():
    ph,ps,b=s.symbols('phi psi b')
    u=s.exp(-b*(ph-ps)/(2*s.sqrt(2)))
    for sign in (1,-1):
        w=s.exp(sign*b*ph/s.sqrt(2))
        pb=-s.diff(u,ph)*s.diff(w,ps)+s.diff(u,ps)*s.diff(w,ph)
        print('Poisson ratio sign',sign,'=',s.simplify(pb/(u*w)),flush=True)
    q=np.exp(.47j);offset=.3;n=2;z=.73
    # Acting on |n>, W+ has coefficient q^(-2n-1-s); W- the inverse.
    wp=lambda n:q**(-2*n-1-offset);wm=lambda n:q**(2*n-1+offset)
    print('W+W- and W-W+ =',wp(n-1)*wm(n),wm(n+1)*wp(n),flush=True)
    # Diagonal(U^-1,U)[1-ix(E21 W+ + E12 W-)]
    print('L factorization offdiag residuals:',abs(q**(n+1+offset/2)*wp(n)-q**(-n-offset/2)),abs(q**(-(n-1)-offset/2)*wm(n)-q**(n+offset/2)),flush=True)

def tau_check():
    q=s.Rational(7,5);E=s.Matrix([[0,1],[0,0]]);I=s.eye(2);K=s.diag(q,1/q)
    # Delta(k^-1 e) = k^-1 e tensor k^-1 + 1 tensor k^-1 e
    lhs=s.kronecker_product(K.inv(),K.inv())*(s.kronecker_product(E,I)+s.kronecker_product(K,E))
    rhs=s.kronecker_product(K.inv()*E,K.inv())+s.kronecker_product(I,K.inv()*E)
    print('tau coproduct residual=',lhs-rhs,flush=True)
    # Complex s cannot give unitarity; crossing is an algebraic condition.
    q=np.exp(.71j);z=.8;n=0
    for offset in (.3,.3+.2j):
        L=np.array([[1,-1j*z],[-1j*z,1]])@np.diag([q**(-n-offset/2),q**(n+1+offset/2)])/np.sqrt(1+z*z)
        print('unitarity s=',offset,'residual=',np.linalg.norm(L.conj().T@L-np.eye(2)),flush=True)

def phi(delta,lam):
    # Derive the phase by summing the integral representation of log Gamma.
    # The linear gamma terms cancel between the two numerator/denominator pairs.
    t=lam*delta/(2*np.pi)
    def integrand(v):
        if v<1e-7:return (lam+1)*t/4
        H=(np.exp(-v/4)-np.exp(-(lam/2+.75)*v))/((1+np.exp(-lam*v/2))*(-np.expm1(-v)))
        return np.sin(t*v)*H/v
    phase=np.pi*(1+lam)/4+2*quad(integrand,0,np.inf,epsabs=3e-11,limit=400)[0]
    return np.exp(1j*phase)/(1+1j*np.exp(lam*delta))

def asymptotes():
    for lam in (.73,1.37,2.6):
        print('lambda',lam,'Phi(-6),Phi(-10) errors=',[abs(phi(d,lam)-1) for d in (-6,-10)],flush=True)
        print('lambda',lam,'positive Phi*Y / exp(i*pi*lambda/2) error=',abs(phi(10,lam)*np.exp(10*lam)/np.exp(1j*np.pi*lam/2)-1),flush=True)
    # At negative rapidity, use Phi from the integral and compare full block.
    lam=1.37;q=np.exp(1j*np.pi*(1-lam));offset=.3;n=2
    target=np.diag([q**(-n),q**(n+1+offset)])
    for d in (-6,-10):
        z=np.exp(lam*d)
        L=phi(d,lam)*np.array([[1,-1j*z],[-1j*z,1]])@target
        print('twist limit delta=',d,'full charge-block residual=',np.linalg.norm(L-target),flush=True)

if __name__=='__main__':cft();classical_limit();tau_check();asymptotes()
