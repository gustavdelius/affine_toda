import numpy as np
from proj import cproj
def cf(n,s,om,th):
    T=(n+1)*om/2; mu=-1j*T*th/np.pi
    if s==1:
        return np.cos(np.pi*mu)*np.cos(np.pi*(mu+om/2))*np.cos(np.pi*(mu-T))/np.cos(np.pi*(mu-om/2))
    return np.sin(np.pi*mu)*np.sin(np.pi*(mu+om/2))*np.sin(np.pi*(mu-T))/np.sin(np.pi*(mu-om/2))
for n in (2,3,4,5):
  for s in (1,-1):
    for om in (0.37,1.83,2.6):
      T=(n+1)*om/2
      r=[cproj(n,s,om,th)[0]/cf(n,s,om,th) for th in (0.31+0.07j,-0.5+0.3j,0.9-0.2j)]
      print(n,s,om,np.round(r,8), ' sin^2(pi om/2)=',round(np.sin(np.pi*om/2)**2,6), 'cos(pi om/2)', round(np.cos(np.pi*om/2),6))
