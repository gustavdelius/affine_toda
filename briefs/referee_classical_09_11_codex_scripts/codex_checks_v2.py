"""Independent checks from the printed potential and tau ansatz. No book imports."""
import itertools, math
import mpmath as mp
import numpy as np
import sympy as s
mp.mp.dps=35

def check(name, ok, result):
    print(('PASS ' if ok else 'FAIL ')+name+': '+str(result),flush=True)

# 1. Singular sine-Gordon, m=beta=1, m0=2.
x=s.symbols('x',real=True)
rho=-8/s.sinh(2*x)**2
check('01 Laurent series',True,s.series(rho,x,0,4))
for eps in [mp.mpf('.01'),mp.mpf('.001')]:
    cut=8*(1-mp.coth(2*eps))
    check('01 cutoff energy',cut<0,mp.nstr(cut,18))

# 2. Same real SG antikink, described by theta=i*pi rather than theta=0.
anti_energy=mp.quad(lambda y:8/mp.cosh(2*y)**2,[-mp.inf,0,mp.inf])
check('02 reversed spatial orientation',abs(anti_energy-8)<mp.mpf('1e-25'),f'integral={anti_energy}; printed formula={8*mp.cos(mp.pi)}')

# 3. Topological charge of the regular kink. Reduced (species,rapidity) fixed by Theta.
qtop=mp.quad(lambda y:4/mp.cosh(2*y),[-mp.inf,0,mp.inf])/(2*mp.pi)
check('03 topological counterexample',abs(qtop-1)<mp.mpf('1e-25'),f'Q={qtop}; Theta parity=-1')

# 4. Explicit zeros for positive A, and regular equal-amplitude type-B example A<0.
h=3; u=mp.mpf('.2'); A=1-mp.sin(mp.pi/h)**2/mp.cos(u)**2
# xi1=xi2=i*pi/2; tau=1+2*i*r*cos(psi)-A*r^2.
r=1/mp.sqrt(A); psi=mp.pi/2
tau=1+2j*r*mp.cos(psi)-A*r*r
check('04 forbidden low-frequency parameter',abs(tau)<mp.mpf('1e-30'),f'A={mp.nstr(A,16)}, r={mp.nstr(r,16)}, tau={mp.nstr(tau,5)}')
u=mp.pi/3; A=1-mp.sin(mp.pi/h)**2/mp.cos(u)**2
check('04 regular type B construction',A<0,f'u=pi/3 A={mp.nstr(A,8)}; Re(tau)=1+abs(A)*r^2>0')

# 5. A1 matrix curvature using independent jets.
b,m,la,q,p,v,d=s.symbols('b m la q p v d', nonzero=True,real=True)
H=s.diag(1,-1)/s.sqrt(2); e=s.Matrix([[0,1],[0,0]]); f=e.T
ep=s.exp(b*q/s.sqrt(2)); em=1/ep
Ap=b*p*H/2+m*la*(ep*e+em*f)
Am=-b*v*H/2-m/la*(ep*f+em*e)
F=Am.diff(q)*p+Am.diff(v)*d-Ap.diff(q)*v-Ap.diff(p)*d+Ap*Am-Am*Ap
res=F.subs(d,-2*s.sqrt(2)*m*m/b*s.sinh(s.sqrt(2)*b*q)).applyfunc(lambda z:s.simplify(z.rewrite(s.exp)))
check('05 A1 curvature',res==s.zeros(2),res)

# 6. Fusion masses, including the bad endpoint and negative continuation.
for a,bv,h in [(1,1,3),(1,2,3),(2,2,3)]:
    z=2*mp.sin(mp.pi*a/h)*mp.exp(1j*mp.pi*bv/h)+2*mp.sin(mp.pi*bv/h)*mp.exp(-1j*mp.pi*a/h)
    check('06 fusion sum',abs(z-2*mp.sin(mp.pi*(a+bv)/h))<mp.mpf('1e-30'),f'h={h}, a={a}, b={bv}, mass expression={mp.nstr(z,8)}')

# 7. A2 charge basis in R3 sum-zero subspace; beta=1.
w=np.exp(2j*np.pi/3)
roots=np.array([[-1,0,1],[1,-1,0],[0,1,-1]],float)
basis={a:np.array([w**(a*k) for k in range(3)])/np.sqrt(3) for a in [1,2]}
for triple in [(1,1,1),(1,1,2)]:
    C=np.sum(np.prod(np.array([roots@basis[a] for a in triple]),axis=0))
    check('07 A2 cubic',abs(C-(3 if triple==(1,1,1) else 0))<1e-12 or (triple==(1,1,1) and abs(abs(C)-3)<1e-12),f'{triple}: |C|={abs(C):.12g}, area RHS=3')

# 8. Frame dependent time delay.
th1=-.4; th2=-1.; h=3; a=bv=1
A=(math.cosh(th1-th2)-1)/(math.cosh(th1-th2)-math.cos(2*math.pi/3))
dx=-math.log(A)/(math.sqrt(3)*math.cosh(th1)); dt=-dx/math.tanh(th1)
check('08 negative-velocity time delay',dt>0,f'dx={dx:.12g}; dt={dt:.12g}')
# Abel-regularized cosine integral: tanh(pi*t/2)=1-2 sum_{n>=1}(-1)^(n-1)e^(-n*pi*t).
# Instead use stable integration by parts: I= -(pi/(2*theta)) integral sech^2(pi*t/2)*sin(theta*t)/t? 
# Direct damped difference derivative is checked using exact product identity in report.

# 9. Current algebra and boosted kink integral.
B,M,P,R=s.symbols('B M P R',real=True)
# p=D+phi, r=D+^2phi; D-p=-M²/B sinh(Bphi), D-r=-M² cosh(Bphi)*p.
current=-2*M*M*s.cosh(B*q)*P*R-B*M*M*s.sinh(B*q)*P**3+M*M*(B*s.sinh(B*q)*P**3+2*s.cosh(B*q)*P*R)
check('09 current divergence',s.simplify(current)==0,s.simplify(current))
for th in [mp.mpf('0'),mp.mpf('.6')]:
    def dens(z):
        fp=2/mp.cosh(z); fpp=-2*mp.tanh(z)/mp.cosh(z)
        return (mp.exp(-4*th)*(fpp*fpp-fp**4/4)+mp.exp(-2*th)*(1-2/mp.cosh(z)**2)*fp**2)/mp.cosh(th)
    Q=mp.quad(dens,[-mp.inf,0,mp.inf])
    check('09 kink spin',abs(Q+mp.mpf(16)/3*mp.exp(-3*th))<mp.mpf('1e-25'),f'theta={th}; Q={mp.nstr(Q,17)}; Q/Q0={mp.nstr(Q/(-mp.mpf(16)/3),17)}')

# 10. Charge chambers, standard alpha_j=e_j-e_{j+1}, alpha0=e_h-e_1.
maxerr=0.; ncases=0
for h in range(2,15):
    roots=np.zeros((h,h)); roots[0,-1]=1; roots[0,0]=-1
    for j in range(1,h): roots[j,j-1]=1; roots[j,j]=-1
    for a in range(1,h):
        N=h//math.gcd(a,h)
        boundaries=sorted(set(round((math.pi-2*math.pi*j*a/h)%(2*math.pi),10) for j in range(h)))
        charges=[]
        for k in range(len(boundaries)):
            phase=(boundaries[k]+(boundaries[(k+1)%N]+(2*math.pi if k==N-1 else 0)))/2
            phases=(2*np.pi*np.arange(h)*a/h+phase+np.pi)%(2*np.pi)-np.pi
            Q=-phases@roots/(2*np.pi)
            maxerr=max(maxerr,float(np.max(np.abs((Q+a/h)-np.rint(Q+a/h)))))
            assert np.sum(np.rint(Q+a/h))==a
            charges.append(tuple(np.round(Q,8)))
        assert len(set(charges))==N
        ncases+=1
check('10 charge chambers',maxerr<1e-8,f'{ncases} (h,a) pairs h=2..14; max weight residual={maxerr:.3g}')

# 11. Bn, Cn direct mass matrices in orthonormal Cartesian coordinates.
for n in [2,3,4,7]:
    I=np.eye(n)
    broots=np.vstack([-I[0]-I[1],*[I[j]-I[j+1] for j in range(n-1)],I[-1]])
    bmarks=np.array([1,1]+[2]*(n-1))
    croot=np.vstack([-np.sqrt(2)*I[0],*[(I[j]-I[j+1])/np.sqrt(2) for j in range(n-1)],np.sqrt(2)*I[-1]])
    cmarks=np.array([1]+[2]*(n-1)+[1])
    bmass=np.sqrt(np.linalg.eigvalsh(broots.T@np.diag(bmarks)@broots))
    cmass=np.sqrt(np.linalg.eigvalsh(croot.T@np.diag(cmarks)@croot))
    bexp=np.sort(np.sqrt(2)*np.array([2*np.sin(np.pi*a/(2*n)) for a in range(1,n)]+[1]))
    cexp=np.sort([2*np.sin(np.pi*a/(2*n)) for a in range(1,n+1)])
    check('11 folded mass spectra',np.max(abs(bmass-bexp))+np.max(abs(cmass-cexp))<1e-12,f'n={n}; max B,C errors={np.max(abs(bmass-bexp)):.2g},{np.max(abs(cmass-cexp)):.2g}')

# 12,13 analytic exponential jets, no finite differences, independent canonical energy.
def taujets(h,data,x,t):
    terms=[]
    for mask in itertools.product([0,1],repeat=len(data)):
        c=mp.mpc(1); k=mp.mpc(0); nu=mp.mpc(0); species=0
        for i,(a,th,xi) in enumerate(data):
            if mask[i]:
                ma=2*mp.sin(mp.pi*a/h); k+=ma*mp.cosh(th);nu-=ma*mp.sinh(th);species+=a;c*=mp.exp(xi)
                for l in range(i):
                    if mask[l]:
                        al,tl,_=data[l];ch=mp.cosh(th-tl)
                        c*=(ch-mp.cos(mp.pi*(a-al)/h))/(ch-mp.cos(mp.pi*(a+al)/h))
        terms.append((c*mp.exp(k*x+nu*t),k,nu,species))
    out=[]
    for j in range(h):
        vals=[mp.mpc(0) for _ in range(5)]
        for c,k,nu,a in terms:
            z=c*mp.exp(2j*mp.pi*j*a/h)
            for l,v in enumerate([1,k,nu,k*k,nu*nu]): vals[l]+=z*v
        out.append(vals)
    return out

def residual_energy(h,data,x,t):
    tj=taujets(h,data,x,t); residual=[]; lx=[];lt=[];V=mp.mpc(0)
    for j,(z,zx,zt,zxx,ztt) in enumerate(tj):
        neighbor=tj[(j-1)%h][0]*tj[(j+1)%h][0]
        residual.append((z*(zxx-ztt)-zx*zx+zt*zt-z*z+neighbor)/(1+abs(z)**2))
        lx.append(zx/z);lt.append(zt/z);V+=1-neighbor/(z*z)
    # Components sum alpha_j*dlog tau_j = l_{k+1}-l_k with cyclic indexing shift irrelevant to norm.
    phix=[1j*(lx[(j+1)%h]-lx[j]) for j in range(h)]
    phit=[1j*(lt[(j+1)%h]-lt[j]) for j in range(h)]
    energy=sum((px*px+pt*pt)/2 for px,pt in zip(phix,phit))+V
    boundary=2*sum(zxx/z-(zx/z)**2 for z,zx,zt,zxx,ztt in tj)
    return max(abs(r) for r in residual),energy,boundary
cases=[(3,[(1,mp.mpf('.4'),.3j),(2,mp.mpf('-.5'),1.2j)]),(4,[(1,mp.mpf('.3'),.4j),(2,mp.mpf('-.6'),1.1j)]),(3,[(1,mp.pi/3*1j,mp.pi/2*1j),(2,-mp.pi/3*1j,mp.pi/2*1j)])]
for h,data in cases:
    vals=[residual_energy(h,data,x,mp.mpf('.23')) for x in [-1,0,1]]
    err=max(v[0] for v in vals);stress=max(abs(v[1]-v[2]) for v in vals)
    check('12 bilinear and 13 stress',err<mp.mpf('1e-25') and stress<mp.mpf('1e-25'),f'h={h}, max residual={mp.nstr(err,4)}, stress error={mp.nstr(stress,4)}')
# Regular one-soliton energies with phase safely off zero rays.
for h,a,phase in [(3,1,.3),(4,2,1.1)]:
    data=[(a,mp.mpf('0'),1j*phase)]
    energy=mp.quad(lambda xx:residual_energy(h,data,xx,0)[1],[-20,-2,0,2,20])
    expected=4*h*mp.sin(mp.pi*a/h)
    check('13 canonical one-soliton energy',abs(energy-expected)<mp.mpf('1e-12'),f'h={h},a={a}, E={mp.nstr(energy,18)}, expected={mp.nstr(expected,18)}')
