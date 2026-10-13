import sys; sys.argv=['x']
exec(open('exc4b.py').read().split('from scipy.optimize import minimize_scalar')[0].split("for l1 in np.linspace")[0])
from linrefl2 import refl_matrix2
from t5 import Kdg, blk
h=5; th=0.6; Om=2*np.sin(np.pi/h)*np.cosh(th)
for eps in (1e-3,1e-4,1e-5,1e-6):
    l1=-1.355; ls=np.linspace(-1.2,0.2,300); v=[t2(l1,l)-eps for l in ls]
    r=[brentq(lambda l:t2(l1,l)-eps,ls[i],ls[i+1]) for i in range(len(ls)-1) if v[i]*v[i+1]<0][0]
    cs=cs_of(l1,r); tf=lambda x: T.tau(cs,np.array([x]))[0][:,0].real
    K,lab=refl_matrix2(n,tf,P,Om,X=14,rtol=1e-11,atol=1e-13)
    print(eps,lab,np.round(K,4).tolist())
print('DG b=2 a=1:',np.round(Kdg(1,2,h,th),4))
