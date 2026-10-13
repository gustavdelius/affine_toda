from tauvac import *
from scipy.integrate import quad
n=3; T=Tau(n); h=4
rng=np.random.default_rng(1)
# test: field equation and energy
al=roots(n)
p=np.array([np.log(0.3),0.7,np.log(0.2),0.0]) # n=3 -> k=1: 2 params + mid : nparams
print(nparams(n))
p=np.array([np.log(0.3),0.7,0.2])
cs=cs_from_p(n,p); x=np.linspace(-12,0,2001)
t,d,dd=T.tau(cs,x); print('imag',abs(t.imag).max(), 'min tau',t.real.min())
t=t.real;d=d.real;dd=dd.real
phi=-(al.T@np.log(t)); dphi=-(al.T@(d/t)); ddphi=-(al.T@(dd/t-(d/t)**2))
rhs=al.T@np.exp(al@phi)
print('eom resid',abs(ddphi-rhs).max())
dens=0.5*(dphi**2).sum(0)+(np.exp(al@phi)-1).sum(0)
Eint=np.trapz(dens,x); print('bulk energy',Eint,'W',-2*np.sum(d[:,-1]/t[:,-1]))
