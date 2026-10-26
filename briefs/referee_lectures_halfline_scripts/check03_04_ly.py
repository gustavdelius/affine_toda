import cmath
from math import pi,sin

def b(x,t,h=6):return cmath.sinh(t/2+1j*pi*x/(2*h))/cmath.sinh(t/2-1j*pi*x/(2*h))
def S(t):return (cmath.sinh(t)+1j*sin(2*pi/3))/(cmath.sinh(t)-1j*sin(2*pi/3))
def K(p,t):return b(1,t)*b(3,t)/b(4,t)*b(1+p,t)*b(5-p,t)/(b(1-p,t)*b(5+p,t))
t=.71;p=.37;v=pi*(p+1)/6
exc=S(t-1j*v)*S(t+1j*v)*K(p,t)
print('3 printed shift residual',abs(exc-K(p+2,t)))
print('3 actual 4-b shift residual',abs(exc-K(4-p,t)))
assert abs(exc-K(4-p,t))<1e-12 and abs(exc-K(p+2,t))<1e-12
print('PASS 3: bootstrap shifts b to b+2, equivalently 4-b')
print('4 printed doubled S residual',abs(S(2*t)+b(2,t)*b(4,t)))
print('4 corrected doubled S residual',abs(S(2*t)+b(1,t)*b(2,t)/(b(4,t)*b(5,t))))
print('4 identity K bootstrap residual',abs(K(0,t)-S(2*t)*K(0,t+1j*pi/3)*K(0,t-1j*pi/3)))
print('4 doubling arbitrary block residual',abs(b(1.23,2*t)+b(1.23/2,t)/b(6-1.23/2,t)))
assert abs(S(2*t)+b(1,t)*b(2,t)/(b(4,t)*b(5,t)))<1e-12
print('PASS 4: corrected block doubling and positive identity K satisfy bootstrap')
