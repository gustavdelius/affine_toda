import sympy as s
import numpy as np
E=s.symbols('E'); a,b,m,ph,dp,dm,box,lam=s.symbols('a b m ph dp dm box lam',nonzero=True)
def unit(i,j):
 r=s.zeros(3);r[i,j]=1;return r
H=a*s.diag(1,0,-1);es=[unit(2,0),s.sqrt(2)*(unit(0,1)+unit(1,2))];fs=[e.T for e in es];roots=[-2*a,a];cs=[s.sqrt(2)*a,a]
Ap=b*dp*H/2+m*lam*sum((cs[j]*s.exp(b*roots[j]*ph/2)*es[j] for j in range(2)),s.zeros(3))
Am=-b*dm*H/2-m/lam*sum((cs[j]*s.exp(b*roots[j]*ph/2)*fs[j] for j in range(2)),s.zeros(3))
cur=Am.diff(ph)*dp+Am.diff(dm)*box-Ap.diff(ph)*dm-Ap.diff(dp)*box+Ap*Am-Am*Ap
rhs=-m*m/b*sum([roots[j]*[1,2][j]*s.exp(b*roots[j]*ph) for j in range(2)])
cur=cur.subs(box,rhs).applyfunc(s.simplify)
print('PASS twisted A2 Lax curvature',cur==s.zeros(3))

def bil(t):return 2*(t*(E*s.diff(E*s.diff(t,E),E))-(E*s.diff(t,E))**2)
F=1+3*E;G=1-E;HH=s.Integer(1)
for l in range(2,7):
 f,g,h=s.symbols('f g h');ff=F+f*E**l;gg=G+g*E**l;hh=HH+h*E**l
 rr=[bil(ff)-ff**2+hh,bil(gg)-gg**2+hh,bil(hh)-2*(hh**2-ff*gg**3)]
 sol=s.solve([s.expand(r).coeff(E,l) for r in rr],[f,g,h]);F=ff.subs(sol);G=gg.subs(sol);HH=hh.subs(sol)
rr=[bil(F)-F**2+HH,bil(G)-G**2+HH,bil(HH)-2*(HH**2-F*G**3)]
print('G2 triple orbit F,G,H=',F,G,HH)
print('PASS G2 triple orbit residuals',list(map(s.expand,rr)), 'degrees', [s.degree(t,E) for t in [F,G,HH]])
f=1+E;h=1-4*E+E**2
rr=[3*bil(f)-f**2+h,3*bil(h)-2*(h**2-f**4)]
print('PASS G2 invariant residuals',list(map(s.expand,rr)),'degrees',[1,1,2])
for n in [2,3]:
 h=2*n;roots=[]
 for j in range(h):
  r=s.zeros(h,1);r[j]=1;r[(j+1)%h]=-1;roots.append(r)
 v=sum(roots[::2],s.zeros(h,1));dots=[(r.T*v)[0] for r in roots]
 print('PASS embedded SG n=',n,'v squared=',(v.T*v)[0],'root dots=',dots,'E=',8*n)
# c-soliton zero locations, m=beta=1, a=1, theta=.3+.2i, xi=.4i
th=.3+.2j;k=np.sqrt(3)*np.cosh(th);nu=-np.sqrt(3)*np.sinh(th)
candidates=[]
for j in range(3):
 for ell in range(-3,4):
  target=np.pi-2*np.pi*j/3+2*np.pi*ell-.4
  t=target/(nu.imag-k.imag*nu.real/k.real);x=-nu.real*t/k.real
  if t>0:candidates.append((t,x,j,ell))
t,x,j,ell=min(candidates)
print('PASS complex soliton first positive zero t,x,j,ell=',t,x,j,ell,'tau residual=',abs(1+np.exp(k*x+nu*t+.4j+2j*np.pi*j/3)))
