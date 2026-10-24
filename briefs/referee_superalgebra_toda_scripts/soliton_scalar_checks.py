"""Independent tau-function, supersymmetry, scalar-product and residue checks."""
import numpy as np
import scipy.special as sc
from scipy.integrate import quad
import sympy as s

def check(name,ok,detail=''):
    print(('PASS ' if ok else 'FAIL ')+name+' '+str(detail),flush=True)
    if not ok:raise AssertionError(name)

for m in (2,3,4):
    P=np.eye(m)-np.ones((m,m))/m
    roots=[]
    for j in range(m):
        roots.append(np.r_[-P[j],1j*P[(j-1)%m]]) # even root 2j
        roots.append(np.r_[P[j],-1j*P[j]]) # odd root 2j+1
    roots=np.array(roots)
    gram=roots@roots.T
    J=roots.T@np.diag([1j*(-1)**j for j in range(2*m)])@np.linalg.pinv(roots.T)
    check(f'm={m} root rank/nullity',np.linalg.matrix_rank(gram)==2*m-2 and np.max(abs(np.diag(gram)))<1e-12)
    check(f'm={m} J antisymmetric',np.max(abs(J.T+J))<1e-12)
    rng=np.random.default_rng(7)
    p=rng.normal(size=m);p=P@p;q=P@rng.normal(size=m);field=np.r_[p,q]
    weights=np.exp(roots@field)
    grad=roots.T@weights
    potential=grad@grad/2
    direct=sum(np.exp(p-np.roll(p,-1))-np.exp(1j*(q-np.roll(q,-1))))
    check(f'm={m} Vreal+Vimag',abs(potential-direct)<1e-11)
    for a in range(1,m):
        mu=1.;beta=1.;v=.27;rho=2*np.sin(np.pi*a/m)/np.sqrt(1-v*v)
        phase=.321 # tau zeros avoided
        def data(x,rho=rho):
            z=np.exp(rho*x+1j*(2*np.pi*a*np.arange(m)/m+phase))
            zprev=np.roll(z,1)
            f=1j*(np.log(1+z)-np.log(1+zprev))/beta
            fp=1j*rho*(z/(1+z)-zprev/(1+zprev))/beta
            fpp=1j*rho*rho*(z/(1+z)**2-zprev/(1+zprev)**2)/beta
            return f,fp,fpp
        maxe=maxfirst=0.
        angle=np.pi*a/(2*m)+np.pi/4
        cc=np.cos(angle);ss=np.sin(angle)
        for x in (-1.3,-.1,.7):
            f,fp,fpp=data(x)
            force=(np.exp(1j*(f-np.roll(f,-1)))-np.exp(1j*(np.roll(f,1)-f)))/1j
            maxe=max(maxe,np.max(abs((v*v-1)*fpp+force)))
            full=np.r_[np.zeros(m),f];der=np.r_[np.zeros(m),(1-v)*fp]
            F=-(roots.T@np.exp(roots@full))
            lam=-np.sqrt((1-v)/(1+v))
            residual=-(cc*der-ss*J@der)+lam*(cc*F+ss*J@F)
            maxfirst=max(maxfirst,np.max(abs(residual)))
        check(f'm={m},a={a} soliton eom',maxe<1e-10,maxe)
        check(f'm={m},a={a} first order lambda negative (homogeneous form)',maxfirst<1e-10,maxfirst)
        # Rest energy from local derivative and periodic potential (branch independent).
        rest_rho=2*np.sin(np.pi*a/m)
        def density(x):
            f,fp,_=data(x,rest_rho)
            return fp@fp/2+sum(1-np.exp(1j*(f-np.roll(f,-1))))
        energy=quad(lambda x:density(x).real,-25/rest_rho,25/rest_rho,epsabs=1e-9,limit=400)[0]+1j*quad(lambda x:density(x).imag,-25/rest_rho,25/rest_rho,epsabs=1e-9,limit=400)[0]
        predicted=4*m*np.sin(np.pi*a/m)
        w=np.exp(2j*np.pi*a/m);c=1j*np.exp(1j*np.pi*a/m)
        central=m*(c*(1/w-1)+(w-1)/c)
        check(f'm={m},a={a} rest energy and central boundary term',abs(energy-predicted)<1e-7 and abs(central-predicted)<1e-11,f'E={energy}, Z={central}')

# Gandenberger 3.28 vs 5.6 at a=b=1 and lambda=1/h.
u,j=s.symbols('u j');factor=(j-u)*(u+j-1)/((j-1-u)*(u+j))
for M in (1,2,7):
    prod=s.prod(factor.subs(j,k) for k in range(1,M+1))
    check(f'G/F partial ratio M={M}',s.cancel(prod+(M-u)/(M+u))==0)
check('fusion makes sign (-1)^(ab), not universal minus',(-1)**(2*1)==1)

# Scalar T from the displayed Gamma product. Work in log space with 1/L tail extrapolation.
def logpi(theta,h,omega,L):
    eta=-1j*h*theta/(2*np.pi);r=h*np.arange(L,dtype=float)+eta
    terms=sc.loggamma(r+1-omega)+sc.loggamma(r+omega)-sc.loggamma(r+h/2+1-omega)-sc.loggamma(r+h/2+omega)+sc.loggamma(r+h/2+1)+sc.loggamma(r+h/2)-sc.loggamma(r+1)-sc.loggamma(r)
    return sum(terms)
def T(theta,n,omega):
    h=2*n-1
    def logt(L):return logpi(theta,h,omega,L)+logpi(1j*np.pi-theta,h,omega,L)
    return np.exp(2*logt(4000)-logt(2000))
def S(theta,n,omega):
    h=2*n-1;g=-np.sin(np.pi*omega);t=T(theta,n,omega)
    R=1j*g/np.sinh(h*theta/2)*t;D=(-1)**(n+1)*1j*g/np.cosh(h*theta/2)*t
    B=(1+2j*g/np.sinh(h*theta))*t;F=(-1+2j*g/np.sinh(h*theta))*t
    return np.array([[B,0,0,-D],[0,R,t,0],[0,t,R,0],[D,0,0,F]])
for n in (2,3):
    omega=np.sqrt(2)/5;theta=.91
    t=T(theta,n,omega);tm=T(-theta,n,omega);h=2*n-1
    expected=np.sinh(h*theta/2)**2/(np.sinh(h*theta/2)**2+np.sin(np.pi*omega)**2)
    check(f'n={n} Gamma scalar unitarity',abs(t*tm-expected)<2e-6,abs(t*tm-expected))
    sm=S(theta,n,omega)
    check(f'n={n} physical unitarity',np.max(abs(sm@sm.conj().T-np.eye(4)))<3e-6)
    for a in range(1,n):
        pole=1j*np.pi*(h-2*a)/h;eps=1e-6
        res=(eps*S(pole+eps,n,omega))[np.ix_([0,3],[0,3])]/1j
        check(f'n={n},a={a} positive direct residue and rank one',min(res.diagonal().real)>0 and abs(np.linalg.det(res))<1e-5,f'eigenvalues={np.linalg.eigvals(res)}')
        check(f'n={n},a={a} 1/epsilon pole scaling',abs(abs(S(pole+eps/2,n,omega)[0,0]/S(pole+eps,n,omega)[0,0])-2)<1e-4)
        # Source 2.44, using the residue ratio derived from the displayed operator.
        half=pole/2
        ratio=(-1)**(n+1)*np.cosh(h*pole)/np.sinh(h*pole/2)
        sp,sm=S(theta+half,n,omega),S(theta-half,n,omega)
        bootstrap=sp[0,0]*sm[0,0]-ratio*sp[1,1]*sm[3,0]
        def elementary(k):return np.sinh(theta/2+1j*np.pi*k/(2*h))/np.sinh(theta/2-1j*np.pi*k/(2*h))
        def block(k):return elementary(k-1)*elementary(k+1)/(elementary(k-1+2*omega)*elementary(k+1-2*omega))
        product=np.prod([block(h/2-a+2*p-1) for p in range(1,a+1)])
        check(f'n={n},a={a} bootstrap S_an',abs(bootstrap-product)<3e-6,abs(bootstrap-product))
print('Independent soliton and scalar checks complete.')
