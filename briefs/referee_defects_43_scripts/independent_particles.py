"""Independent scalar, classical sewing, spectrum, CDD and exercise checks."""
import itertools
import numpy as np
import sympy as sp
from independent_breathers import transmit,scalar

def T(z,a,h,typ=-1):
    A=np.exp(1j*np.pi*a/h)
    x=np.exp(-z if typ==-1 else z)
    return (x+A)/(x+1/A) if typ==-1 else (x-A)/(x-1/A)

def block(z,x,h):
    return np.sinh((z+1j*np.pi*x/h)/2)/np.sinh((z-1j*np.pi*x/h)/2)

def S(z,a,b,h,B=.73):
    out=1
    for p in range(abs(a-b)+1,a+b,2):
        out*=block(z,p-1,h)*block(z,p+1,h)/(block(z,p-1+B,h)*block(z,p+1-B,h))
    return out

def scalar_checks():
    errors=dict(unitarity=0,HA=0,phase=0,bootstrap=0,block=0,bulk_CDD=0)
    count=0
    for h in range(2,9):
        for typ in (-1,1):
            for a in range(1,h):
                for z in (.31-.27,-.42+.13j):
                    errors['unitarity']=max(errors['unitarity'],abs(T(z,a,h,typ)*T(z+1j*np.pi,h-a,h,typ)-1))
                    tb=(-np.exp(1j*np.pi*a/h)/block(z,h-a,h)) if typ==-1 else np.exp(1j*np.pi*a/h)/block(z,a,h)
                    errors['block']=max(errors['block'],abs(T(z,a,h,typ)-tb))
                theta=.53;td=.27
                errors['HA']=max(errors['HA'],abs(T(1j*np.pi-theta-td,h-a,h,typ)-np.conj(T(-theta-td,a,h,typ))))
                errors['phase']=max(errors['phase'],abs(abs(T(theta-td,a,h,typ))-1))
                for b in range(1,h):
                    if a+b==h:continue
                    c=(a+b)%h
                    ua,ub=(np.pi*b/h,np.pi*a/h) if a+b<h else (np.pi*(h-b)/h,np.pi*(h-a)/h)
                    z=.317+.103j
                    errors['bootstrap']=max(errors['bootstrap'],abs(T(z+1j*ua,a,h,typ)*T(z-1j*ub,b,h,typ)-T(z,c,h,typ)))
                    for e in range(1,h):
                        errors['bulk_CDD']=max(errors['bulk_CDD'],abs(S(z+1j*ua,a,e,h)*S(z-1j*ub,b,e,h)-S(z,c,e,h)))
                    count+=1
    print('scalar checks h=2..8, both types, fusion cases=',count)
    for key,val in errors.items():print(f'{key}: {val:.3e}')

def classical_sewing():
    worst=reflect=0
    for h in range(2,9):
        roots=np.eye(h)-np.roll(np.eye(h),-1,axis=1)
        # pi e_j=e_(j-1); inverse on the sum-zero root space.
        pi=np.roll(np.eye(h),-1,axis=0)
        D=2*np.linalg.pinv(np.eye(h)-pi)
        P=np.eye(h)-np.ones((h,h))/h
        E=P-D
        for typ in (-1,1):
            for td in (-.42,.27):
                sig=np.exp(-td);sigma=sig if typ==-1 else 1/sig
                f=np.concatenate([roots@D.T,roots@D],axis=1)/2
                g=np.concatenate([roots@D,-roots@D],axis=1)/2
                H=sigma*f.T@f+g.T@g/sigma
                Dkin,Ekin=(D,E) if typ==-1 else (-D,-E)
                for a in range(1,h):
                    gamma=np.pi*a/h;v=np.exp(2j*np.pi*a*np.arange(h)/h)
                    assert np.linalg.norm(pi@v-np.exp(2j*gamma)*v)<1e-12
                    mass=2*np.sin(gamma);energy=mass*np.cosh(.63);mom=mass*np.sinh(.63)
                    def residual(R,U):
                        phi=(1+R)*v;psi=U*v;grad=H@np.r_[phi,psi]
                        left=1j*mom*(1-R)*v+1j*energy*(Ekin@phi+Dkin@psi)+grad[:h]
                        right=1j*mom*psi+1j*energy*(Dkin.T@phi-Ekin@psi)-grad[h:]
                        return np.r_[left,right]
                    zero=residual(0,0)
                    mat=np.column_stack([residual(1,0)-zero,residual(0,1)-zero])
                    R,U=np.linalg.lstsq(mat,-zero,rcond=None)[0]
                    reflect=max(reflect,abs(R));worst=max(worst,abs(U-T(.63-td,a,h,typ)))
    print(f'independent potential Hessian and vector sewing h=2..8: reflection={reflect:.3e} factor_error={worst:.3e}')

def loop_and_exercises():
    x=sp.symbols('x',positive=True)
    tau=(1+sp.I*x)/(1-sp.I*x)
    derivative=sp.factor(x*sp.diff(tau,x)+1/(1-sp.I*(x-1/x)/2))
    print('BS derivative identity after x=exp(theta+mu):',derivative)
    z=.31+.17j
    cz1=np.exp(-1j*np.pi/3)*np.sinh(z/2-1j*np.pi/6)/np.sinh(z/2+1j*np.pi/6)
    cz2=np.exp(1j*np.pi/3)*np.sinh(z/2+2j*np.pi/3)/np.sinh(z/2+1j*np.pi/3)
    print('CZ(7.5) label swap max error:',max(abs(cz1-T(z,2,3)),abs(cz2-T(z,1,3))))
    v=np.pi;h=3;a=1;u=v-np.pi+np.pi*a/h
    print(f'exercise imaginary counterexample: h={h},a={a},v=pi,u/pi={u/np.pi:.6f},|T(.4)|={abs(T(.4-1j*v,a,h)):.15f},sigma=-1 real')
    print('higher SG breather B2 at theta=.4,td=.27:')
    for w in (2.6,3.1):print(w,T(.13+1j*np.pi/(2*w),1,2)*T(.13-1j*np.pi/(2*w),1,2))

def residue_checks():
    hlist=(3,4,5);w=np.sqrt(5);offsets={h:np.linspace(-.13,.13,h) for h in hlist}
    for h in hlist:
        for conj,a in ((False,1),(True,h-1)):
            theta=.27+1j*(np.pi*a/h-np.pi/2-np.pi/(h*w))
            Q=np.zeros(h,dtype=int);Q[0]=-1 if conj else 1
            # Incoming N=Q-weight(i), so every label is in the root lattice.
            matrix=np.zeros((h,h),complex)
            for i in range(h):
                N=Q.copy();N[i]+=1 if conj else -1
                for j,NN,coef in transmit(i,tuple(N),offsets[h],theta,h,w,conj):matrix[j,i]+=coef
            singular=np.linalg.svd(matrix,compute_uv=False)
            rank=int(sum(singular>1e-9*singular[0]))
            print(f'species pole h={h} a={a}: rank(M)={rank}, normalized singular values={singular/singular[0]}')
    print('a2 resonance windows from -1/2 < a/3-1/2-1/(3w) < 0: a1 w>1; a2 1/2<w<2')
    for w in (1.1,2.3,10):
        print(f'a2 omega={w} Imtheta/pi for a1,a2=',[a/3-.5-1/(3*w) for a in (1,2)])

def cdd_and_tba():
    # Entire, unitary, nonminimal scalar with no poles. It also satisfies
    # the full fused soliton crossing equation, not only f(D)f(D+i*pi)=1.
    eta=.19;D=.41;h=4;a=2
    f=lambda z:np.exp(1j*eta*np.sinh(z))
    for w in (np.sqrt(2),np.sqrt(5)):
        u=np.pi/2-np.pi/(h*w)
        shifts=[-u+np.pi*(2*m-(h-a)+1)/h for m in range(h-a)]
        shifts += [u+np.pi*(2*m-a+1)/h for m in range(a)]
        product=np.prod([f(D+1j*s) for s in shifts])
        predicted=np.exp(2j*eta*np.sin(np.pi*a/h)/np.sin(np.pi/h)*np.sin(np.pi/(h*w))*np.sinh(D))
        cross=f(D-1j*np.pi)*np.prod([f(D+1j*np.pi*(2*m-h+2)/h) for m in range(h-1)])
        print(f'admissible exponential CDD omega={w:.12g}: breather multiplier={product}, formula_error={abs(product-predicted):.3e}, full_crossing={abs(cross-1):.3e}')
    # The single-species condition printed in ch43 is insufficient for h>2.
    bad=lambda z:(np.sinh(z)+.3j)/(np.sinh(z)-.3j)
    short=bad(D)*bad(D+1j*np.pi)
    hb=3
    full=bad(D-1j*np.pi)*np.prod([bad(D+1j*np.pi*(2*m-hb+2)/hb) for m in range(hb-1)])
    print(f'ordinary sinh CDD at h=3: printed crossing residual={abs(short-1):.3e}, actual fused crossing residual={abs(full-1):.6f}')
    # A bounded, meromorphic CDD family obtained from a second particle factor.
    for w in (np.sqrt(2),np.sqrt(5)):
        u=np.pi/2-np.pi/(h*w)
        shifts=[-u+np.pi*(2*m-(h-a)+1)/h for m in range(h-a)]+[u+np.pi*(2*m-a+1)/h for m in range(a)]
        factor=np.prod([T(D+1j*s,1,h) for s in shifts])
        crossing=T(D-1j*np.pi,1,h)*np.prod([T(D+1j*np.pi*(2*m-h+2)/h,1,h) for m in range(h-1)])
        print(f'meromorphic CDD omega={w:.12g}: multiplier={factor}, full_crossing={abs(crossing-1):.3e}')
    worst=0
    for h in range(2,9):
        for a in range(1,h):
            F=T(1j*np.pi/2+.73-.27,h-a,h)
            Fbar=T(1j*np.pi/2+.73-.27,a,h)
            worst=max(worst,abs(F-np.conj(Fbar)))
    print(f'TBA crossed factor conjugacy={worst:.3e}; self-conjugate F at theta=.73,td=.27: {T(1j*np.pi/2+.73-.27,1,2)}')

if __name__=='__main__':
    scalar_checks();classical_sewing();loop_and_exercises();residue_checks();cdd_and_tba()
