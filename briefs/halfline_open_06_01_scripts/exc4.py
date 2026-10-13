from redmin2 import *
n=4; T=Tau(n); P=(1,1,-1,1,-1)
best=[]
for s1 in (0,np.pi):
  for s2 in (0,np.pi):
    for k in range(5):
      chi1=(s2-2*s1+2*np.pi*k)/5; chi2=2*chi1+s1
      for l1 in np.linspace(-3,0.5,36):
        for l2 in np.linspace(-3,0.5,36):
          cs=cs_from_p(n,np.array([l1,chi1,l2,chi2]))
          m=tau_min(T,cs)
          if m.min()<=1e-5: continue
          F,t0,eh=Fval(T,cs,P)
          best.append((F,s1,s2,k,l1,l2,tuple(np.round(t0,4))))
best.sort(key=lambda z:z[0])
for b in best[:8]: print(b)
