import numpy as np, time
from sol_boot import boot_coef
from rhonp import rho_plus, blk
om=0.37
def f(mu,T): return -blk(T/2,mu,T)**2
def g(mu,T,om): return blk(T/2-om/2,mu,T)*blk(T/2+om/2,mu,T)
for n in (2,3,4):
    T=(n+1)*om/2
    for th in (0.31+0.07j,):
        mu=-1j*T*th/np.pi
        P=boot_coef(n,om,th,1,order=-1)[0]/boot_coef(n,om,th,0,order=-1)[0]
        mus=[mu+T-k*om for k in range(1,n+1)]
        for lab,Y in (('rho',lambda m: rho_plus(m,n,om)),('rho f',lambda m: rho_plus(m,n,om)*f(m,T)),('rho f g',lambda m: rho_plus(m,n,om)*f(m,T)*g(m,T,om))):
            lhs=Y(mu); rhs=P*np.prod([Y(m) for m in mus])
            print(n,lab,'Y/(P prod Y)=',lhs/rhs)
        print(n,'f/prod f =',f(mu,T)/np.prod([f(m,T) for m in mus]),' g/prod g=',g(mu,T,om)/np.prod([g(m,T,om) for m in mus]))
