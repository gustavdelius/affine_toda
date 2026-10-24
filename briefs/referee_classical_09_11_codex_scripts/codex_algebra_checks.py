"""Additional independent finite-root and polynomial checks; no book imports."""
import numpy as np
import itertools
import sympy as s

def check(name,ok,result): print(('PASS ' if ok else 'FAIL ')+name+': '+str(result))
# 14: D4 reflections, Coxeter plane and full orbit partition.
I=np.eye(4)
roots=np.array([I[0]-I[1],I[1]-I[2],I[2]-I[3],I[2]+I[3]])
S=[np.eye(4)-np.outer(a,a) for a in roots]
w=S[0]@S[2]@S[3]@S[1] # outer black; central white
weights=np.linalg.inv(roots).T
orbits=[]
for lam in weights:
    gamma=(np.eye(4)-w.T)@lam
    orbit=np.array([np.linalg.matrix_power(w,k)@gamma for k in range(6)])
    orbits.append(orbit)
actual=set()
for i,j in itertools.combinations(range(4),2):
    for si,sj in itertools.product([-1,1],repeat=2):actual.add(tuple(si*I[i]+sj*I[j]))
union={tuple(np.round(v,8)) for O in orbits for v in O}
check('14 D4 orbit partition',union==actual and len(union)==24,f'{len(union)} roots, sizes {[len(set(map(tuple,O))) for O in orbits]}')
eig,vec=np.linalg.eig(w.astype(complex)); z=vec[:,np.argmin(abs(eig-np.exp(1j*np.pi/3)))]
plane=np.array([z.real/np.linalg.norm(z.real),z.imag/np.linalg.norm(z.imag)])
radii=[np.linalg.norm(O@plane.T,axis=1) for O in orbits]
ratio=np.array([r.mean() for r in radii])/radii[0].mean()
check('14 D4 radii',np.max(abs(ratio-[1,np.sqrt(3),1,1]))<1e-12,ratio)
sets=[{tuple(np.round(v,9)) for v in O@plane.T} for O in orbits]
check('14 D4 coincident outer hexagons',sets[0]==sets[2]==sets[3],len(sets[0]))
tri=next((a,b,c) for a in orbits[0] for b in orbits[0] for c in orbits[1] if np.linalg.norm(a+b+c)<1e-10)
check('14 D4 triangle',True,[list(v.astype(int)) for v in tri])
# 14: A2 explicit Coxeter orbits and cubic selection.
r=np.array([[1,-1,0],[0,1,-1]],float)
w=(np.eye(3)-np.outer(r[0],r[0]))@(np.eye(3)-np.outer(r[1],r[1]))
lams=np.array([[2,-1,-1],[1,1,-2]],float)/3
O=[np.array([np.linalg.matrix_power(w,k)@(np.eye(3)-w.T)@lam for k in range(3)]) for lam in lams]
allowed=[]
for a,b,c in itertools.combinations_with_replacement(range(2),3):
    if any(np.linalg.norm(x+y+z)<1e-10 for x in O[a] for y in O[b] for z in O[c]): allowed.append((a+1,b+1,c+1))
check('14 A2 Dorey rule',allowed==[(1,1,1),(2,2,2)],allowed)
# 16: explicit c2 polynomial from parent A3, verify each Hirota equation algebraically.
E=s.symbols('E'); taus=[1+2*E+E**2/2,1+E**2/2,1-2*E+E**2/2,1+E**2/2]
def dx(f):return s.sqrt(2)*E*s.diff(f,E)
res=[s.expand(t*dx(dx(t))-dx(t)**2-t**2+taus[(j-1)%4]*taus[(j+1)%4]) for j,t in enumerate(taus)]
check('16 c2 polynomial',res==[0]*4,res)
# 20: Dn root construction and marks.
for n in [4,5,6,7,8,10]:
    I=np.eye(n); r=np.array([*[I[j]-I[j+1] for j in range(n-1)],I[-2]+I[-1]])
    marks=np.array([1]+[2]*(n-3)+[1,1]); ar=np.vstack([-marks@r,r]); am=np.r_[1,marks]
    mass=np.sqrt(np.linalg.eigvalsh(ar.T@np.diag(am)@ar)); h=2*n-2
    target=np.sort([2*np.sqrt(2)*np.sin(np.pi*a/h) for a in range(1,n-1)]+[np.sqrt(2)]*2)
    check('20 Dn masses',np.max(abs(mass-target))<1e-12,f'n={n}, max error={np.max(abs(mass-target)):.3g}, trace={sum(mass**2):.12g}')
# E8 Bourbaki coordinates, with explicit affine marks.
I=np.eye(8)
r=np.array([(I[0]+I[7]-sum(I[j] for j in range(1,7)))/2,I[0]+I[1],*[I[j]-I[j-1] for j in range(1,7)]])
marks=np.array([2,3,4,6,5,4,3,2]); ar=np.vstack([-marks@r,r]); am=np.r_[1,marks]
check('20 E8 root lengths',np.max(abs(np.sum(ar*ar,axis=1)-2))<1e-12,np.sum(ar*ar,axis=1))
mass=np.sqrt(np.linalg.eigvalsh(ar.T@np.diag(am)@ar)); c=lambda d:np.cos(np.pi*d)
b=2*c(1/5); target=np.array([1,b,2*c(1/30),2*b*c(7/30),2*b*c(2/15),2*b*c(1/30),4*b*c(1/5)*c(7/30),4*b*c(1/5)*c(2/15)])
check('20 E8 ratios',np.max(abs(mass/mass[0]-target))<1e-12,f'max error={np.max(abs(mass/mass[0]-target)):.3g}, trace={sum(mass**2):.12g}')
# 19: energy subset sums at A3 m=beta=1.
M=np.array([8*np.sqrt(2),16,8*np.sqrt(2)])
energies=[(tuple(i+1 for i in range(3) if mask[i]),float(np.dot(mask,M))) for mask in itertools.product([0,1],repeat=3) if any(mask)]
check('19 A3 static energies',all(abs(e-16)>1e-10 or subset==(2,) for subset,e in energies),energies)
