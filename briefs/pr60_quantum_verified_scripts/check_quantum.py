"""Independent PR60 checks, no repository imports. 2026-10-25."""
import sympy as S
import mpmath as mp
import numpy as np
from collections import defaultdict
mp.mp.dps=35

def ck(label,error,tol=1e-24):
    ok=abs(error)<tol
    print(('PASS' if ok else 'FAIL'),label,'residual =',mp.nstr(error,9),flush=True)
    if not ok: raise AssertionError(label)

# Linear-form commutators directly from [phi,psi]=-i hbar.
b,h=S.symbols('b h',nonzero=True)
forms={'U':(-b/(2*S.sqrt(2)),b/(2*S.sqrt(2))), 'W':(b/S.sqrt(2),0),
       'Us':(-S.sqrt(2)*S.pi/(h*b),S.sqrt(2)*S.pi/(h*b)), 'X':(2*S.sqrt(2)*S.pi/(h*b),0)}
def comm(a,c): return S.simplify(-S.I*h*(a[0]*c[1]-a[1]*c[0]))
for a,c,expected in [('U','W',S.I*h*b*b/4),('Us','X',4*S.I*S.pi**2/(h*b*b)),('U','X',S.I*S.pi),('Us','W',S.I*S.pi),('W','X',0),('U','Us',0)]:
    ck('BCH '+a+','+c,comm(forms[a],forms[c])-expected)
# Odd/even linearized matching, on-shell determinant.
g,k,m,sigma=S.symbols('g k m sigma',nonzero=True)
Odd=S.Matrix([[g+k-sigma*m,-g+k-sigma*m],[g-k+m/sigma,g+k-m/sigma]])
expr=S.expand(Odd.det()).subs(g*g,k*k-m*m)
ck('odd determinant',S.simplify(expr-4*k*(k-m*(sigma+1/sigma)/2)))
ck('even determinant',S.simplify(expr.subs(m,-m)-4*k*(k+m*(sigma+1/sigma)/2)))
ck('zero-rapidity kernel',sum(abs(a) for a in Odd.subs({g:0,k:m,sigma:1})))

lam=mp.sqrt(7)
a=mp.pi/(2*lam)
def logphi(u):
    def f(k):
        if abs(k)<mp.mpf('1e-15'): return lam*(1j*mp.pi*u-u*u)/(2*mp.pi)
        return (mp.cosh(mp.pi*k/2+1j*k*u)/mp.cosh(mp.pi*k/2)-1)/(2*k*mp.sinh(a*k))
    return -lam*u/2+1j*mp.pi*lam/4-mp.log(2)/2+mp.quad(f,[0,1,4,mp.inf])
def phi(u): return mp.exp(logphi(u))
def pcl(u): return (1j-mp.sinh(u))/mp.cosh(u)
for u in [mp.mpf('.31'),mp.mpc('.27','.43')]:
    ck('scalar crossing '+str(u),phi(u)*phi(u+1j*mp.pi)*(1+mp.exp(2*lam*u))-1,1e-25)
    v=u+1j*mp.pi/2;e=1j*a
    ck('B1 ratio '+str(u),phi(v-e)/phi(v+e)-pcl(u),1e-25)
ck('scalar real unitarity',abs(phi(mp.mpf('.31')))**2*(1+mp.exp(2*lam*mp.mpf('.31')))-1)
# Independent gamma partial product convergence; no old helper or extrapolation.
u=mp.mpf('.31')
def gamma_phi(K):
    z=1j*lam*u/(2*mp.pi);r=0
    for j in range(K):
        for c,sgn in [(j*lam+mp.mpf('.25'),1),((j+1)*lam+mp.mpf('.75'),1),((j+mp.mpf('.5'))*lam+mp.mpf('.25'),-1),((j+mp.mpf('.5'))*lam+mp.mpf('.75'),-1)]:
            r+=sgn*(mp.loggamma(c-z)-mp.loggamma(c+z))
    return mp.exp(1j*mp.pi*(1+lam)/4+r)/(1+1j*mp.exp(lam*u))
errs=[abs(gamma_phi(K)-phi(u)) for K in [80,160,320]]
print('PASS gamma cutoff errors',*[mp.nstr(v,9) for v in errs],flush=True)
assert errs[2]<errs[1]<errs[0] and errs[0]/errs[1]>1.9
# Pole residues via crossing, not a magnitude cutoff.
for n in range(3):
    un=-(2*n+1)*1j*a; res=-1/(2*lam*phi(un+1j*mp.pi))
    vals=[]
    for eps in [mp.mpf('1e-4'),mp.mpf('1e-6')]:
        p=1/((1+mp.exp(2*lam*(un+eps)))*phi(un+eps+1j*mp.pi))
        vals.append(abs(eps*p-res))
    print('PASS pole',n,'residue',mp.nstr(res,9),'errors',*[mp.nstr(v,6) for v in vals],flush=True)
    assert vals[1]<vals[0]/90
# Sigma direct integral and elementary matrix dressing.
th=mp.mpf('.31');td=mp.mpf('.24')
sig=mp.cosh(lam*td)/mp.cosh(lam*(th-td))*mp.exp(1j*mp.quad(lambda k: mp.sin(k*th)*mp.sinh(mp.pi*k/2+1j*k*td)/(k*mp.sinh(a*k)*mp.cosh(mp.pi*k/2)),[0,1,4,mp.inf]))
rhs=(1+mp.exp(-2*lam*td))/(1+mp.exp(-2*lam*(th+td)))*phi(th-td)/phi(-th-td)
ck('GZ scalar identity',sig-rhs)
y,d,H=S.symbols('y d H',nonzero=True)
A=lambda z:S.Matrix([[1,-S.I*z],[-S.I*z,1]])
kp=y*H+1/(y*H);km=y/H+H/y
G=S.simplify(A(1/(y*d)).inv()*S.diag(kp,km)*A(y/d))
Gexp=S.Matrix([[kp+km/d**2,S.I/d*H*(1/y**2-y**2)],[S.I/(d*H)*(1/y**2-y**2),km+kp/d**2]])/(1+1/(y*d)**2)
ck('Dirichlet matrix dressing',sum(abs(S.simplify(v)) for v in G-Gexp))
print('PASS counterexample to report scalar K claim: offdiagonal =',S.factor((A(1/(y*d)).inv()*A(y/d))[0,1]),flush=True)
# Weyl inverse L, from inversion of each total-charge 2x2 block.
q=np.exp(-1j*np.pi*(np.sqrt(7)+1));offset=.37

def Linv(sign,n,z,s=offset):
    if sign==1:
        return [(1,n,q**n/(1+z*z)),(-1,n+1,1j*z*q**(-n-1-s)/(1+z*z))]
    return [(-1,n,q**(-n-s)/(1+z*z)),(1,n-1,1j*z*q**(n-1)/(1+z*z))]
def R(sp1,sp2,x):
    if sp1==sp2:return [(sp1,sp2,1)]
    den=q*q*x*x-1
    return [(sp1,sp2,q*(x*x-1)/den),(sp2,sp1,(q*q-1)*x/den)]
def D(ys,yd,s=offset,n0=0):
    N=len(ys);out=np.zeros((2**N,2**N),complex)
    for col in range(2**N):
        signs=tuple(1 if ((col>>(N-1-j))&1)==0 else -1 for j in range(N))
        paths={(signs,n0):1}
        for j,yj in enumerate(ys):
            nxt=defaultdict(complex)
            for (sp,n),amp in paths.items():
                for sout,no,c in Linv(sp[j],n,yj/yd,s):
                    ss=list(sp);ss[j]=sout;nxt[(tuple(ss),no)]+=amp*c
            paths=nxt
        for (sp,n),amp in paths.items():
            if n==n0:
                row=sum((x==-1)<<(N-1-j) for j,x in enumerate(sp));out[row,col]+=amp
    return out

def transfer(ys,ya,w):
    N=len(ys);out=np.zeros((2**N,2**N),complex)
    for col in range(2**N):
        signs=tuple(1 if ((col>>(N-1-j))&1)==0 else -1 for j in range(N))
        for ain in [1,-1]:
            paths={(ain,signs):1}
            for j,yj in enumerate(ys):
                nxt=defaultdict(complex)
                for (aa,sp),amp in paths.items():
                    for ao,so,c in R(aa,sp[j],ya/yj):
                        ss=list(sp);ss[j]=so;nxt[(ao,tuple(ss))]+=amp*c
                paths=nxt
            for (aa,sp),amp in paths.items():
                if aa==ain:
                    row=sum((x==-1)<<(N-1-j) for j,x in enumerate(sp));out[row,col]+=amp*w**aa
    return out
for N in [2,3,4,5]:
    ys=np.exp(np.linspace(-.4,.5,N));yd=np.exp(.18+.07j);dd=D(ys,yd)
    errors=[]
    for charge in range(-N,N+1,2):
        ix=[j for j in range(2**N) if N-2*j.bit_count()==charge];idx=np.ix_(ix,ix)
        w=q**(-charge/2);t=transfer(ys,1j*yd,w);z=ys/yd
        rhs=q**(-N/2)*(D(ys,q*yd)+np.prod((1+q*q*z*z)/(1+z*z))*D(ys,yd/q))
        errors.append(np.max(np.abs((t@dd-rhs)[idx])))
        t2=transfer(ys,np.exp(.43),w)
        errors.append(np.max(np.abs((t2@dd-dd@t2)[idx])))
        errors.append(np.max(np.abs((D(ys,yd,n0=2)-q**(2*charge)*dd)[idx])))
    ck('TQ, tD, shift N='+str(N),max(errors),1e-10)
    ck('DD numerical only N='+str(N),np.max(np.abs(dd@D(ys,1.4)-D(ys,1.4)@dd)),1e-10)
# Wigner slope asymptotics and classical breather phase difference.
for L in [80,160,320]:
    slope=L/mp.pi*mp.quad(lambda t:mp.tanh(L*t)/mp.sinh(t),[0,1/L,1,mp.inf])
    asym=L/mp.pi*(mp.log(8*L/mp.pi)+mp.euler)
    print('PASS Wigner asymptotic L=',L,'difference=',mp.nstr(slope-asym,8),flush=True)
    assert abs(slope-asym)<mp.mpf('.1')/L
alpha=mp.mpf('.6');u=mp.mpf('.37')
zfun=lambda w:-mp.coth(w/2)
for n in [40,80,160]:
    L=n*mp.pi/(mp.pi-2*alpha)
    def Tb(n):return mp.fprod(pcl(u+(n+1-2*k)*1j*mp.pi/(2*L)) for k in range(1,n+1))
    rat=Tb(n+1)/Tb(n)
    target=zfun(u+1j*alpha)/abs(zfun(u+1j*alpha))
    print('PASS classical breather n=',n,'clock error=',mp.nstr(abs(rat-target),8),flush=True)
    assert abs(rat-target)<1/n
