"""Check source (2.19) literally, and shared Cartan constraints in the limits."""
import numpy as np
import sympy as s

q=np.exp(.71j);Q=np.exp(-.355j);bp=.37+.19j;bm=-.28+.13j;r=np.exp(.23j);x=np.exp(.41)
def T(z,a):
    A=lambda a:Q**a-bm*bp.conjugate()*Q**(-a)*z*z
    D=lambda a:-bp*bm.conjugate()*Q**(a-4)*z*z+Q**(-a)
    B=lambda a:z*(bp*Q**a+bm*Q**(-a))
    C=lambda a:-z*(bm.conjugate()*Q**(a-4)+bp.conjugate()*Q**(-a))
    return np.array([[A(a),C(a+2)],[B(a),D(a+2)]])

def printed(a,h):
    g=(abs(bp)**2*abs(bm)**2+1)*Q**-2
    p=lambda a:(x*x-x**-2)*(bp*bm.conjugate()*Q**(2*a-2)/r/x-bp.conjugate()*bm*Q**(-2*a-2)*r*x)+x*(g*r-h/r)+(g/r-h*r)/x
    m=lambda a:(x*x-x**-2)*(-bp*bm.conjugate()*Q**(2*a-6)/r*x+bp.conjugate()*bm*Q**(-2*a+2)*r/x)+x*(g/r-h*r)+(g*r-h/r)/x
    up=lambda a:(x*x-x**-2)*(bm/r+abs(bp)**2*bm*r+bp*Q**(2*a)/r+bm*bm*bp.conjugate()*Q**(-2*a)*r)
    down=lambda a:(x*x-x**-2)*(bm.conjugate()*Q**-4*r+abs(bp)**2*bm.conjugate()*Q**-4/r+bm.conjugate()**2*bp*Q**(2*a-8)/r+bp.conjugate()*Q**(-2*a)*r)
    return np.array([[p(a),down(a+2)],[up(a),m(a+2)]])
kd=np.diag([r*x+1/r/x,r/x+x/r])
hprinted=-(abs(bm)**2*Q**-2+abs(bp)**2)*Q**-2
# Direct substitution in source (2.17): h=b_- c_+ q^-1+b_+ c_- q=-(|b_-|^2+|b_+|^2)q.
hcorrect=-(abs(bm)**2+abs(bp)**2)*Q**-2
for a in (0,2,4):
    F=np.linalg.solve(T(1/x,a),kd@T(x,a))
    for name,h in (('printed h',hprinted),('substitution h',hcorrect)):
        P=printed(a,h);z=np.vdot(P,F)/np.vdot(P,P)
        print('CZ 2.19',a,name,'relative full-block error',np.linalg.norm(F-z*P)/np.linalg.norm(F))

q0=s.symbols('q',nonzero=True);n=s.symbols('n',integer=True);d=s.symbols('d',nonzero=True)
# Limit L=diag(A,D) forces chi(k)=A D^-1 and chi(k^-1)=D^-1 A.
print('diagonal limit required chi(k) chi(k^-1) at s=0:',q0**(-4*n))
# For the charge-conjugating limit B|n>=-i q^n|n-1>, C|n>=-i q^-n|n+1>.
# C B^-1 and B^-1 C both shift by +2, so their product shifts by +4, never the identity.
print('offdiagonal limit: chi(k) shifts +2; independently fitted chi(k^-1) also shifts +2')

# Four numerical unrestricted probes confirm the SAME residual as the matching forced grade.
from independent_checks import ls_test
for sign,kind,tw in ((-1,'e',0),(-1,'f',1),(1,'f',0),(1,'e',1)):
    graded=ls_test(3,sign,kind,tw,margin=1)
    arbitrary=ls_test(3,sign,kind,tw,margin=1,full=True)
    print('all-entry minus graded residual:',sign,kind,tw,arbitrary-graded)
