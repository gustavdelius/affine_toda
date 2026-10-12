from tauvac import *
n=5; T=Tau(n); al=roots(n)
p=np.array([-10.995,5.236,1.936,10.472,0.0]); cs=cs_from_p(n,p)
x=-np.geomspace(40,1e-6,200000)[::-1]; x=np.concatenate([-np.geomspace(40,1e-6,200000),[0.0]])
t,d,dd=T.tau(cs,x); print('imag',abs(t.imag).max(),'min tau',t.real.min())
t=t.real; d=d.real; dd=dd.real
phi=-(al.T@np.log(t)); dphi=-(al.T@(d/t)); ddphi=-(al.T@(dd/t-(d/t)**2))
print('eom',abs(ddphi-al.T@np.exp(al@phi)).max())
dens=0.5*(dphi**2).sum(0)+(np.exp(al@phi)-1).sum(0)
print('bulk',np.trapz(dens,x),'W',-2*np.sum(d[:,-1]/t[:,-1]))
print('phi0',phi[:,-1])
