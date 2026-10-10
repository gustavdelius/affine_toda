# Evaluate Hollowood's continuum integral I(a,n), eq. (3.26) of hep-th/9209024, as printed (m=1)
import numpy as np
from scipy.integrate import quad
def I(a,n):
    ms=lambda b: 2*np.sin(np.pi*b/n)
    ma=ms(a)
    def f(k):
        s=0.0
        for b in range(1,n):
            mb=ms(b); c=ma**2+mb**2-ms(a+b)**2
            s+= -4*np.sqrt(k*k+mb*mb)*c/(c*c+4*ma*ma*k*k) + 0.5*mb**2/np.sqrt(k*k+mb*mb)
        return s/(4*np.pi)
    return 2*quad(f,0,np.inf,limit=500)[0]
for n in [3,4,5,6,7,8]:
    for a in range(1,n//2+1):
        print(n,a, round(I(a,n),6), " -cot/2:",round(-0.5/np.tan(np.pi/n),6)," -cosec/2:",round(-0.5/np.sin(np.pi/n),6), " -cot/4:",round(-0.25/np.tan(np.pi/n),6))
# Hollowood bound-state sum for odd n
print("odd n bound sums / m_a, compare (cot+cosec)/2 and (cot+cosec)/4")
for n in [3,5,7,9]:
    for a in range(1,(n+1)//2):
        s=0.0
        for b in range(1,n):
            if b<a or (n/2<=b<n/2+a): s+= np.sin(np.pi*b/n)*abs(np.sin(np.pi*(a-b)/n))
        ma=2*np.sin(np.pi*a/n)
        print(n,a, round(s/ma,6), round((1/np.tan(np.pi/n)+1/np.sin(np.pi/n))/2,6), round((1/np.tan(np.pi/n)+1/np.sin(np.pi/n))/4,6))
