"""Independent Gamma-product, integral, pole, window and Wigner checks."""
import numpy as np
from scipy.special import loggamma, digamma
from scipy.integrate import quad


def logH(h,w,z,K=800):
    a=.5+(h-1)*w/2+np.arange(K)*h*w+0j
    return np.sum(loggamma(a+z)+loggamma(a-z)-loggamma(a+w+z)-loggamma(a+w-z)
                  -2*loggamma(a)+2*loggamma(a+w))


def extrap_H(h,w,z):
    # Each normalized summand has a complete expansion beginning at k^-2.
    vals=np.array([logH(h,w,z,K) for K in [100,200,400,800]])
    for j in range(1,4): vals=(2**j*vals[1:]-vals[:-1])/(2**j-1)
    return vals[0]


def constant(h,w):
    return np.exp(-(np.log(2)+sum(extrap_H(h,w,w*((h-1)/2-j)) for j in range(h)))/h)


def g(h,w,t,c=None):
    if c is None: c=constant(h,w)
    z=1j*h*w*(t-1j*np.pi/2)/(2*np.pi)
    return c*np.exp(-w*(t-1j*np.pi/2)/2+extrap_H(h,w,z))


def gb(h,w,t,c=None):
    return np.prod([g(h,w,t+1j*np.pi*(2*m-h+2)/h,c) for m in range(h-1)])


def phase_derivative(h,w,t):
    """From the convergent Gamma integral, not finite differences."""
    z=1j*h*w*(t-1j*np.pi/2)/(2*np.pi)
    a=.5+(h-1)*w/2
    def integrand(u):
        if u==0: return 2*z/(h*w)
        # combine the growing sinh with the decaying exponential first
        return ((np.exp((-a+z-h*w)*u/w)-np.exp((-a-z-h*w)*u/w))*(-np.expm1(-u))
                /((-np.expm1(-u/w))*(-np.expm1(-h*u))))
    val=quad(lambda u: integrand(u).real,0,120,epsabs=2e-9,limit=500)[0]
    first=digamma(a+z)-digamma(a-z)-digamma(a+w+z)+digamma(a+w-z)
    return h/(2*np.pi)*val+h*w/(2*np.pi)*first.real


def run():
    for h in [3,4,5]:
      for w in [1.37,.6]:
        q=-np.exp(-1j*np.pi*w); c=constant(h,w)
        errors=[]
        for t in [-.43,.28,.79+.09j]:
          X=np.exp(w*(t-1j*np.pi/2)); gg=g(h,w,t,c)
          errors.extend([abs(gg*gb(h,w,t-1j*np.pi,c)*(1+(-q)**(h-1)*X**h)-1),
                         abs(gb(h,w,t,c)*g(h,w,t-1j*np.pi,c)*(1-q*X**h)-1),
                         abs(g(h,w,t+2j*np.pi,c)/gg-(1+(-1/q)**(h-1)*X**h)/(1+(-1/q)**(h+1)*X**h))])
        print('scalar h,w,c,max residual',h,w,c,max(errors))
    h=3; w=1.37; z=.43+.72j
    for K in [100,1000,2000,4000]:
        AK=.5+(h-1)*w/2+h*w*K
        drift=loggamma(AK-z)-loggamma(AK)
        norm=logH(h,w,z,K)
        print('cutoff K,logH,normalized printed log,endpoint+zlogAK',K,norm,norm+drift,drift+z*np.log(AK))
    for h in [3,4]:
      for w in [10,40,160]:
        errs=[[],[]]
        for t in [-.7,.5,1.1]:
          classical=np.log(abs((np.exp(-t)-1j*np.exp(1j*np.pi/h))/(np.exp(-t)-1j*np.exp(-1j*np.pi/h))))
          vec=phase_derivative(h,w,t)
          bar=sum(phase_derivative(h,w,t+1j*np.pi*(2*m-h+2)/h) for m in range(h-1))
          errs[0].append(abs(2*np.pi/(h*(w+1))*vec-classical))
          errs[1].append(abs(2*np.pi/(h*(w+1))*bar-classical))
        print('Wigner h,w,errors g/gbar',h,w,max(errs[0]),max(errs[1]))
    for h in [3,4,5]:
      w=1.37; c=constant(h,w)
      for a,fun in [(1,g),(h-1,gb)]:
        tp=1j*np.pi/2-1j*(np.pi*(h-a)/h+np.pi/(h*w))
        residues=[eps*fun(h,w,tp+eps,c) for eps in [1e-3,1e-4,1e-5]]
        print('pole residues h,a',h,a,residues)
    for h in [3,4,5,6]:
      print('strict windows h',h,[(a,max(0,a-h/2),a) for a in range(1,h)])
    print('a2,a=2,omega=2 endpoint E/M',np.cosh(.7+1j*np.pi/6-1j*np.pi/(3*2)))
    # The sign condition is periodic unless restricted to the nearest rapidity sheet.
    w=1/6.5; h=3; a=1
    t=.7+1j*np.pi/2-1j*(np.pi*(h-a)/h+np.pi/(h*w))
    print('remote-sheet positive-Re negative-Im counterexample omega,theta,E',w,t,np.cosh(t))
    for h in [3,4]:
      for w in [1.37,10,40]:
        x=np.exp(.3*w)
        print('suppression h,w',h,w,max(min(x**j,x**(j-h+1)) for j in range(1,h-1)))
    w=1.37
    corrected=np.exp(-1j*np.pi*w/2); printed=np.exp(1j*np.pi*w/2)
    print('CZ2007 6.14 corrected residual, printed residual',abs(corrected-np.exp(w*(-1j*np.pi/2))),abs(printed-corrected))


if __name__=='__main__': run()
