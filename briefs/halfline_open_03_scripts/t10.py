import mpmath as mp
from rho import rho_plus, blk, cproj_plus
mp.mp.dps=25
om=mp.mpf('0.37')
for n in (1,2):
    h=n+1; T=h*om/2
    def R(mu):
        ma=mu-T/2+mp.mpf(1)/2; mb=mu+T/2-mp.mpf(1)/2
        return rho_plus(ma,n,om,K=120)*rho_plus(mb,n,om,K=120)*cproj_plus(mu,om,T)
    tgt=lambda mu: blk(n*om/2-mp.mpf(1)/2,mu,T)*blk(mp.mpf(1)/2,mu,T)*blk(om/2+mp.mpf(1)/2,mu,T)/blk(mp.mpf(1)/2+h*om/2,mu,T)
    mu=mp.mpc('0.13','0.21')
    print(n,'R cross',mp.nstr(R(T/2-mu)/R(T/2+mu),10),' tgt cross',mp.nstr(tgt(T/2-mu)/tgt(T/2+mu),10))
    print(n,'R period 2T',mp.nstr(R(mu+2*T)/R(mu),10),' R period 1', mp.nstr(R(mu+1)/R(mu),10))
