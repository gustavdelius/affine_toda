"""Independent checks of CZ 2009b using their dimensionless real-coupling fields."""
import numpy as np
import sympy as s

p,q,l,sig,gamma=s.symbols('p q l sigma gamma',nonzero=True)
C=s.exp(q)+s.exp(-q)+gamma
f=(2*s.exp(p-l)+s.exp(-p+l)*C)/sig
g=sig*(s.exp(l)*C+2*s.exp(-l))
J=lambda f,g:s.expand((s.diff(f,q)*s.diff(g,l)-s.diff(f,l)*s.diff(g,q))/2)
Vs=s.exp(p+q)+s.exp(-p-q)-s.exp(p-q)-s.exp(-p+q)
print('SG momentum residual:',s.simplify(J(f,g)-Vs),flush=True)
Ct=s.exp(q)+s.exp(-q)
ft=(s.exp(2*p-2*l)+s.exp(-p+l)*Ct)/sig
gt=sig*(8*s.exp(-l)+s.exp(2*l)*Ct**2)/2
Vt=s.exp(2*(p+q))+2*s.exp(-(p+q))-s.exp(2*(p-q))-2*s.exp(-(p-q))
print('Tz momentum residual:',s.simplify(J(ft,gt)-Vt),flush=True)
wrong=s.exp(2*p-2*l)/sig+s.exp(-p+l)*Ct
print('Tz printed-f control:',complex((J(wrong,gt)-Vt).subs({p:.23,q:.17,l:.11,sig:.73}).evalf()),flush=True)

def data(E,z,theta,kind):
    if kind=='sg':
        def field(x):return 2*np.log((1+x)/(1-x))
        def slope(x):return 4*x/(1-x*x)
        a=np.sqrt(2)*np.cosh(theta);b=-np.sqrt(2)*np.sinh(theta)
    else:
        def field(x):return np.log((1+x)**2/(1-4*x+x*x))
        def slope(x):return 2*x/(1+x)+(4*x-2*x*x)/(1-4*x+x*x)
        a=np.sqrt(6)*np.cosh(theta);b=-np.sqrt(6)*np.sinh(theta)
    u,v=field(E),field(z*E);du,dv=slope(E),slope(z*E)
    return (u+v)/2,(u-v)/2,a*(du+dv)/2,b*(du+dv)/2,a*(du-dv)/2,b*(du-dv)/2,b

def residual(E,z,lamfunc,theta,eta,tau,kind):
    pp,qq,px,pt,qx,qt,bt=data(E,z,theta,kind)
    sig=np.exp(eta) if kind=='sg' else np.exp(eta)/np.sqrt(2)
    lam=lamfunc(E)
    h=2e-5
    # derivative in log E, with a fourth-order centered stencil
    derivative=(-lamfunc(E*np.exp(2*h))+8*lamfunc(E*np.exp(h))-8*lamfunc(E*np.exp(-h))+lamfunc(E*np.exp(-2*h)))/(12*h)
    lt=bt*derivative/lam
    c=np.exp(qq)+np.exp(-qq)
    if kind=='sg':
        c+=2*np.cosh(2*tau)
        rhs1=-(sig*lam+np.exp(-pp)*lam/sig)*(np.exp(qq)-np.exp(-qq))/2
        rhs2=-sig*(lam*c-2/lam)/2
        rhs3=(np.exp(-pp)*lam*c-2*np.exp(pp)/lam)/(2*sig)
    else:
        rhs1=-sig/2*lam**2*(np.exp(2*qq)-np.exp(-2*qq))-np.exp(-pp)*lam*(np.exp(qq)-np.exp(-qq))/(2*sig)
        rhs2=-sig/2*(lam**2*c*c-4/lam)
        rhs3=(np.exp(-pp)*lam*c-2*np.exp(2*pp)/lam**2)/(2*sig)
    return np.max(np.abs([px-pt+2*lt-rhs1,qx-qt-rhs2,qx+qt-rhs3]))

def sg_checks():
    theta=1.13;eta=.27;tau=.36;sig=np.exp(eta)
    rho=np.tanh((theta-eta+tau)/2);rhob=np.tanh((theta-eta-tau)/2)
    z1=rho*rhob
    corrected=lambda E:-(1+E)*(1+z1*E)/(np.sqrt(2)*np.cosh(tau)*(1-rho*E)*(1-rhob*E))
    printed=lambda E:(1+E)*(1+z1*E)/(np.sqrt(2)*np.cosh(tau)*(1+rho*E)*(1+rhob*E))
    worst=max(residual(E,z1,corrected,theta,eta,tau,'sg') for E in (.1j,.7j,2.1j))
    bad=max(residual(E,z1,printed,theta,eta,tau,'sg') for E in (.1j,.7j,2.1j))
    def solve_lam(E,z):
        pp,qq,px,pt,qx,qt,_=data(E,z,theta,'sg')
        return -2*(sig*sig*(qx+qt)+np.exp(pp)*(qx-qt))/(sig*(np.exp(pp)-np.exp(-pp))*(np.exp(qq)+np.exp(-qq)+2*np.cosh(2*tau)))
    plus=lambda E:solve_lam(E,1/z1)
    good2=max(residual(E,1/z1,plus,theta,eta,tau,'sg') for E in (.1j,.7j,2.1j))
    print('SG tanh sector: corrected residual=',worst,'printed control=',bad,'lambda(E=0)=',corrected(0),flush=True)
    print('SG coth sector: residual=',good2,'lambda(E~0)=',plus(1e-6j),'expected +vacuum=',1/(np.sqrt(2)*np.cosh(tau)),flush=True)
    assert max(worst,good2)<2e-8 and bad>1

def tz_checks():
    theta=1.13;eta=.27;r=np.exp(-theta+eta);d=(theta-eta)/2
    coth=lambda z:1/np.tanh(z)
    z1=coth(d-1j*np.pi/12)*coth(d+1j*np.pi/12)
    z2=coth(d-1j*np.pi/4)*np.tanh(d+1j*np.pi/12)
    z3=np.conj(z2)
    rho1=(r-np.sqrt(2))*(r+np.sqrt(2))/((r-np.exp(1j*np.pi/6))*(r-np.exp(-1j*np.pi/6)))
    rho2=(r-np.sqrt(2)*np.exp(1j*np.pi/3))*(r+np.sqrt(2)*np.exp(1j*np.pi/3))/((r+np.exp(1j*np.pi/6))*(r-np.exp(-1j*np.pi/2)))
    lam1=lambda E:(1+E)*(1+z1*E)/(1+2*rho1*E+z1*E*E)
    lam2=lambda E:np.exp(-2j*np.pi/3)*(1+E)*(1+z2*E)/(1+2*rho2*E+z2*E*E)
    lam3=lambda E:np.conj(lam2(np.conj(E)))
    errors=[]
    for z,lf in ((z1,lam1),(z2,lam2),(z3,lam3)):
        errors.append(max(residual(E,z,lf,theta,eta,0,'tz') for E in (.07+.08j,.2+.1j)))
    raw=(r-np.exp(1j*np.pi/6))*(r+np.exp(-1j*np.pi/2))/((r+np.exp(1j*np.pi/6))*(r-np.exp(-1j*np.pi/2)))
    printed=coth(d+1j*np.pi/4)*np.tanh(d-1j*np.pi/12)
    print('Tz sectors residuals=',errors,'product-1=',abs(z1*z2*z3-1),flush=True)
    print('CZ(4.9) z2 minus book z2=',abs(raw-z2),'CZ(4.11) minus z3=',abs(printed-z3),flush=True)
    assert max(errors)<2e-8 and abs(raw-z2)<1e-12 and abs(printed-z3)<1e-12

if __name__=='__main__':sg_checks();tz_checks()
