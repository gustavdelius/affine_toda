import numpy as np, mpmath as m
m.mp.dps=45
q=-np.exp(-1j*np.pi*np.sqrt(2)); lam=np.sqrt(2);t=.37-.12j
E0=.53-.19j;E1=-.41+.31j
R=lambda z: np.array([[1,0,0,0],[0,*np.linalg.solve([[z,1/q],[1/q,z]],[1,z/q]),0],[0,0,0,0],[0,0,0,1]],complex)
def braid(z):
 c,b=np.linalg.solve([[z,1/q],[1/q,z]],[1,z/q]);return np.array([[1,0,0,0],[0,c,b,0],[0,b,c,0],[0,0,0,1]])
def K(y):return np.array([[(q-1/q)*(E1*y+E0/y),y*y-y**-2],[y*y-y**-2,(q-1/q)*(E0*y+E1/y)]])
C=np.array([[0,1],[1,0]]);P=np.eye(4)[[0,2,1,3]]
vm=(K(np.exp(lam*(.5j*np.pi-t)))@C).T.ravel();vp=(K(np.exp(lam*(.5j*np.pi+t)))@C).T.ravel()
for label,v in [('correct S_cd',braid(np.exp(2*lam*t))@vp),('printed S_dc',braid(np.exp(2*lam*t))@P@vp)]:
 coef=np.vdot(vm,v)/np.vdot(vm,vm);er=np.linalg.norm(v-coef*vm)/np.linalg.norm(v);print(label,'matrix proportionality residual',er)
 assert (er<1e-12)==label.startswith('correct')
print('PASS crossing lower-index order is fixed by lecture ZF convention')
# Kim light/heavy a3 one-loop coefficients from original paper compared with independent derivative.
def block(x,t,h):return m.sinh(t/2+1j*m.pi*x/(2*h))/m.sinh(t/2-1j*m.pi*x/(2*h))
def KN(a,B,t,h):return m.fprod(block(c-1,t,h)*block(h-c+B/2,t,h)/(block(h-c,t,h)*block(c-1+B/2,t,h)) for c in range(1,a+1))
t=m.mpf('.71');z=m.cosh(t);v=m.sinh(t);r=m.sqrt(2)
coeff=[1j*v/16*(1/(z+1/r)-1/(z-1)),1j*v/16*(1/z-1/(z-1)-1/(z-1/r)+1/(z+1/r))]
for a in [1,2]:
 independent=m.diff(lambda B:KN(a,B,t,4),0)/(2*m.pi)
 print('Kim a3 species',a,'derivative error',m.nstr(abs(independent-coeff[a-1]),8));assert abs(independent-coeff[a-1])<m.mpf('1e-35')
print('PASS direct source Kim one-loop light/heavy coefficients match KN')
# Separate issue: proposed self-dual ansatz in lecture compared with actual Kim product.
# At h=3, Kim's [1/2][3/2] is not the printed KN*Kcl/K.
def Kim(B):return m.fprod(block(x-m.mpf('.5'),t,3)*block(x+m.mpf('.5'),t,3)/(block(x-m.mpf('.5')+B/2,t,3)*block(x+m.mpf('.5')-B/2,t,3)) for x in [m.mpf('.5'),m.mpf('1.5')])
for B in [m.mpf('.31'),m.mpf('.47')]:
 print('Kim actual a2 formula dual residual',m.nstr(abs(Kim(B)-Kim(2-B)),10))
# a1 Dirichlet limit and all-rank parity coweight calculation are symbolic below.
import sympy as s
z,L=s.symbols('z L',nonzero=True);H=s.diag(1,-1);E=s.Matrix([[0,1],[0,0]]);F=E.T
# z=e^{i betaSG phi0/2}; K_dir division specified in exercise.
KD=-s.Matrix([[0,z*L-1/(z*L)],[L/z-z/L,0]])
a=(z*(L*E-F/L)+z**-1*(L*F-E/L))/2
at=(z*(F/L-L*E)+z**-1*(E/L-L*F))/2
assert s.simplify(KD*a+at*KD)==s.zeros(2) and KD*H+H*KD==s.zeros(2)
print('PASS a1 Dirichlet K limit and arbitrary Cartan derivative')
for h in range(3,10):
 roots=[s.eye(h)[:,h-1]-s.eye(h)[:,0]]+[s.eye(h)[:,j-1]-s.eye(h)[:,j] for j in range(1,h)]
 for a in range(1,h):
  v=-sum((s.Rational(a*j,h)*roots[j] for j in range(h)),s.zeros(h,1))
  assert all((r.dot(v)).is_Integer for r in roots)
print('PASS parity vacuum shift lies in coweight lattice for h=3..9, every species; branch changes add roots')
# The actual self-dual Kim a2 product preserves the one-loop coefficient only.
d1=m.diff(Kim,0)-m.diff(lambda B:KN(1,B,t,3),0)
d2=(m.diff(Kim,0,2)-m.diff(lambda B:KN(1,B,t,3),0,2))/2
print('PASS actual Kim a2 minus KN: first derivative',m.nstr(abs(d1),8),'B^2 coefficient',m.nstr(d2,15))
assert abs(d1)<m.mpf('1e-35') and abs(d2)>m.mpf('.01')
# CDRS general product versus the dual minus branch, independent direct product.
for h in range(3,9):
 for a in range(1,h):
  B=m.sqrt(2)/4
  cdr=m.fprod(block(c,t,h)*block(h-c+B/2,t,h)/(block(c-1+B/2,t,h)*block(h-c+1,t,h)) for c in range(1,a+1))
  dualminus=-block(a,t,h)*block(h-a,t,h)*KN(a,B,t,h)
  assert abs(cdr-dualminus)<m.mpf('1e-35')
print('PASS CDRS all-rank candidate equals dual minus product, h=3..8 all a')
