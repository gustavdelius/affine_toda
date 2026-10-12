import mpmath as mp
from rho import rho_plus, blk, cproj_plus
mp.mp.dps=30
for om in (mp.mpf('0.37'),mp.mpf('1.3')):
  for n in (1,2,3):
    h=n+1; T=h*om/2
    for mu in (mp.mpc('0.13','0.21'),mp.mpc('-0.31','0.4')):
        ma=mu-T/2+mp.mpf(1)/2; mb=mu+T/2-mp.mpf(1)/2
        R=rho_plus(ma,n,om,K=150)*rho_plus(mb,n,om,K=150)*cproj_plus(mu,om,T)
        tgt=blk(n*om/2-mp.mpf(1)/2,mu,T)*blk(mp.mpf(1)/2,mu,T)*blk(om/2+mp.mpf(1)/2,mu,T)/blk(mp.mpf(1)/2+h*om/2,mu,T)
        print(float(om),n,mp.nstr(mu,4),'R/target =',mp.nstr(R/tgt,10))
