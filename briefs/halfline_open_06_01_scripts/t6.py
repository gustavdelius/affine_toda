from linrefl import *
from t5 import blk, Kdg, pattern_of
from bp_table import class_b
n=3; h=4; T=Tau(n); al=roots(n)
th=0.6; Om=2*np.sin(np.pi/h)*np.cosh(th)
for d in (0.3,0.9,0.99,0.999,0.9999):
    cs=np.zeros(n,complex); cs[1]=d
    tf=lambda x: T.tau(cs,np.array([x]))[0][:,0].real
    r,C=pattern_of(T,cs)
    K,lab=refl_matrix(n,tf,C,Om,rtol=1e-12,atol=1e-14)
    print(d,C,'bc %.1e'%r,np.round(K,4).tolist())
print('DG b=2, a=1:',np.round(Kdg(1,2,h,th),4),' flat C=-1:',np.round(-(blk(1,th,h)*blk(3,th,h)),4),' flat C=+1:',np.round(-1/(blk(1,th,h)*blk(3,th,h)),4))
