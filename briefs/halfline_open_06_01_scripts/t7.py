from linrefl import *
from t5 import blk, Kdg, pattern_of
n=2; h=3; T=Tau(n); s=np.sin(np.pi/h); d=1/(1+s)
th=0.6; Om=2*np.sin(np.pi/h)*np.cosh(th)
chiend=(np.pi-2*np.pi*2/h)%(2*np.pi)
print('chi_end',chiend)
for delta in (0.3,0.1,0.03,0.01,0.003):
  for sgn in (+1,-1):
    chi=chiend+sgn*delta
    cs=cs_from_p(n,np.array([np.log(d),chi]))
    tf=lambda x: T.tau(cs,np.array([x]))[0][:,0].real
    r,Cex=pattern_of(T,cs)
    t0=tf(0.0)
    # P = exact pattern with the sign at node 1 or 0 flipped so that product is -1, keeping C_2=-1
    for Cp in [(1,1,-1),(-1,-1,-1)]:
        K,lab=refl_matrix(n,tf,Cp,Om,rtol=1e-12,atol=1e-14)
        print('delta %+.3f exact pattern %s tau0 %s | P=%s K='%(sgn*delta,Cex,np.round(t0,4),Cp),np.round(K,4).tolist())
