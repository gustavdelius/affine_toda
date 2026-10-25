import numpy as np
from independent_checks import ux, phase_derivative
err=0.
for th in [.4,.8,1.3]:
 for a in [.6,1.4,-.9]:
  h=np.sinh(th); k=np.sin(a/2)/np.sqrt(h*h+np.sin(a/2)**2)
  for s in [1,-1]:
   for t in [-7.,-.3,0.,.8,5.]:
    closed=2*np.sin(a/2)+4*h*np.sqrt(1-k*k)/(s*np.cosh(h*t)+k)
    err=max(err,abs(ux(t,th,a,s)-closed))
print('PASS closed boundary derivative residual',err)
print('Elementary integral: I_s = 4*pi*s - 8*asin(k), k=sin(a/2)/sqrt(sinh(theta)^2+sin(a/2)^2)')
for L in [300+np.sqrt(2),1000+np.sqrt(2),3000+np.sqrt(2)]:
 a=1.4;t=.8
 vals=[phase_derivative((L+1)*a/2,t,L,N) for N in [4000,16000]]
 ext=(4*vals[1]-vals[0])/3
 print('large-lambda',L,'cutoff-extrapolated',ext,'target',2/np.pi*np.arctan(np.sin(a/2)/np.sinh(t)))
