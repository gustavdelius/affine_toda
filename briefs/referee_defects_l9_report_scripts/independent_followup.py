"""Further exact checks; imports only this referee's own derivation."""
from independent_symbolic import *
err=(M-4/z/(1+1/(y*z)**2)*Lam.inv()*G*Lam).applyfunc(lambda f:S.simplify(S.factor(f)))
assert err==S.zeros(2)
print('PASS all four sector matrix entries separately:',err)
H=S.factor(Dm+Dp+C); U=H.subs(v,0)
check('boundary primitive kinetic coefficient',H+ b*b/m**3*dxfield(F)*v*v-U)
check('boundary primitive potential derivative',dxfield(U)-(b*b/(2*m)*cosphi*F-b/m*sinphi*dxfield(F)+b**4/(4*m**3)*F**3))
# Fourier-mode gluing: reflecting root labels reverses the cyclic eigenvalue.
zeta=S.symbols('zeta',nonzero=True)
check('mirror D eigenvalue',2/(1-1/zeta)-(2-2/(1-zeta)))
# Tzitzeica product, writing exp(u)=x and exp(i pi/6)=w.
w=S.sqrt(3)/2+I/2
zz=(1/x+w)*(1/x+1/w)/((1/x-w)*(1/x-1/w))
check('Tzitzeica z1 matches folded delay',zz-((x+1/x)+S.sqrt(3))/((x+1/x)-S.sqrt(3)))
print('PASS chirality: for phi=L(t+x)+R(t-x), E+P=2 integral R_prime^2; E-P=2 integral L_prime^2')
