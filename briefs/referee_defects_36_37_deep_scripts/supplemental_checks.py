"""Additional independent identities for the report; no book imports."""
import numpy as np
from scipy.special import loggamma, zeta
from independent_checks import phi, log_r, bulk_scalar, six_vertex, PI


print('ITEM 3: natural noncommuting impurity candidates')
lam=1.37
t1,t2=.9,.3
s_plus=bulk_scalar(t1+t2,lam)
s_minus=bulk_scalar(t1-t2,lam)
print('Diagonal reflection with the Weyl transmission: |S0(sum)-S0(diff)| =',abs(s_plus-s_minus))
print('Anticommuting Pauli reflection/transmission: |S0(sum)+S0(diff)| =',abs(s_plus+s_minus))

print('\nITEM 5: KL/BCZ pairing and finite compact-product identity')
lam,theta=1.37,.61
z=1j*lam*theta/(2*PI)
y=np.exp(lam*theta)
for count in (50,200,800):
    k=np.arange(count)*lam
    def kl_log(zz):
        return np.sum(loggamma(k+.75-zz)+loggamma(k+lam+.25-zz)
                      -loggamma(k+lam/2+.25-zz)-loggamma(k+lam/2+.75-zz))
    kl_p=np.exp(1j*PI*(lam-1)/4+kl_log(z)-kl_log(-z))/(1-1j*y)
    k=np.arange(1,count+1)*lam
    compact=(loggamma(.25-z)-loggamma(.25+z)
             +np.sum(loggamma(.5+2*k-2*z)+loggamma(.5+2*k-lam+2*z)
                     -loggamma(.5+2*k+2*z)-loggamma(.5+2*k-lam-2*z)))
    endpoint=(-4*z*np.log(2)+loggamma(count*lam+.75-z)-loggamma(count*lam+.75+z)
              +loggamma(2*count*lam+.5+2*z)-loggamma(2*count*lam+.5-2*z))
    paired=log_r(z,lam,count)-log_r(-z,lam,count)
    print(f'N={count}: raw KL to limit={abs(kl_p-phi(theta,lam)):.6g}; '
          f'exact endpoint identity={abs(np.exp(compact+endpoint)-np.exp(paired)):.3g}; '
          f'regularized compact to paired limit={abs(np.exp(compact+2*z*np.log(count*lam))-phi(theta,lam)*(1+1j*y)*np.exp(-1j*PI*(1+lam)/4)):.6g}')

print('\nITEM 6: fusion covariance of one-kink phases')
d={a:np.sin((a+1)*PI/5)/np.sin(PI/5) for a in range(4)}
f=lambda a,b,t:(d[a]/d[b])**(1j*t/(2*PI))
a,b,c=0,1,2
u,v=.37,.29
constant=(d[a]/d[b])**(-u/(2*PI))*(d[b]/d[c])**(v/(2*PI))
for theta in (-1.2,.4,1.1):
    ratio=f(a,b,theta+1j*u)*f(b,c,theta-1j*v)/f(a,c,theta)
    print(f'theta={theta}: fusion-vertex phase is rapidity-independent: {abs(ratio-constant):.3g}')

print('\nITEM 7: independent Hirota asymptotic check')
a=np.tanh((.4-1.1)/2)**2
e1=.6
tau_ratio=lambda e2:(1+1j*e1+1j*e2-a*e1*e2)/(1-1j*e1-1j*e2-a*e1*e2)
small=(1+1j*e1)/(1-1j*e1)
large=-(1+1j*a*e1)/(1-1j*a*e1)
print('E2 -> 0 residual:',abs(tau_ratio(1e-9)-small))
print('E2 -> infinity residual:',abs(tau_ratio(1e9)-large))
print('Two even defects vs slow soliton temporal factor:',1/a)
print('A11, left-to-right spatial factor:',a)
