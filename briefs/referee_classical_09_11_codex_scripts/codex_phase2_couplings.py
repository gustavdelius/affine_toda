"""Independent A4/E6 cubic tensors, charge bases, Coxeter tests."""
import numpy as np
import itertools

def cubic(roots,marks,basis,a,b,c):return np.sum(marks*(roots@basis[a])*(roots@basis[b])*(roots@basis[c]))
def area(ms):
 sem=sum(ms)/2;return np.sqrt(max(0,sem*np.prod(sem-ms)))
# A4 charge Fourier basis; m=beta=1.
h=5;eye=np.eye(h);roots=np.array([eye[-1]-eye[0]]+[eye[j]-eye[j+1] for j in range(h-1)])
marks=np.ones(h);M=roots.T@roots;w=np.exp(2j*np.pi/h)
basis=np.array([[w**(a*k)/np.sqrt(h) for k in range(h)] for a in range(1,h)])
masses=np.array([2*np.sin(np.pi*a/h) for a in range(1,h)])
err=0;bad=[];nonzero=0
for a,b,c in itertools.combinations_with_replacement(range(h-1),3):
 C=cubic(roots,marks,basis,a,b,c);allowed=(a+b+c+3)%h==0
 if (abs(C)>1e-10)!=allowed:bad.append((a,b,c))
 if allowed:nonzero+=1;err=max(err,abs(abs(C)-4/np.sqrt(h)*area(masses[[a,b,c]])))
print(f'PASS A4: {nonzero} nonzero triples; selection mismatches={bad}; max area error={err:.3g}')
# A4 real basis from real/imag Fourier modes.
rb=np.array([np.sqrt(2)*basis[0].real,np.sqrt(2)*basis[0].imag,np.sqrt(2)*basis[1].real,np.sqrt(2)*basis[1].imag]);rm=np.array([masses[0]]*2+[masses[1]]*2)
violations=[]
for a,b,c in itertools.combinations_with_replacement(range(4),3):
 C=cubic(roots,marks,rb,a,b,c);rhs=4/np.sqrt(h)*area(rm[[a,b,c]])
 if abs(C)>1e-9 and abs(abs(C)-rhs)>1e-8:violations.append((a,b,c,float(abs(C)/rhs)))
print('PASS A4 real-basis nonzero-coupling counterexamples:',violations[:4]);assert violations and err<1e-12 and not bad
# E6 Bourbaki diagram: chain 1-3-4-5-6; 2 attached to 4.
C=2*np.eye(6)
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:C[i,j]=C[j,i]=-1
r=np.linalg.cholesky(C);marks=np.array([1,2,2,3,2,1]);roots=np.vstack([-marks@r,r]);am=np.r_[1,marks];M=roots.T@np.diag(am)@roots
perm=np.array([1,6,3,5,4,2,0]) # affine diagram centre rotation
O=roots[perm[1:]].T@np.linalg.inv(r.T)
assert np.linalg.norm(O@O.T-np.eye(6))<1e-12 and np.linalg.norm(O@M-M@O)<1e-12
# Coxeter bicolouring.
black=[0,3,5];white=[1,2,4];W=np.eye(6)
for i in black+white:W=W@(np.eye(6)-np.outer(r[i],r[i]))
weights=np.linalg.inv(r).T;orbits=[np.array([np.linalg.matrix_power(W,k)@(np.eye(6)-W.T)@lam for k in range(12)]) for lam in weights]
eig,vec=np.linalg.eig(W.astype(complex));v=vec[:,np.argmin(abs(eig-np.exp(2j*np.pi/12)))];plane=np.array([v.real/np.linalg.norm(v.real),v.imag/np.linalg.norm(v.imag)])
radii=np.array([np.linalg.norm(orbit[0]@plane.T) for orbit in orbits]);target=radii*np.sqrt(np.trace(M)/sum(radii**2))
chars=np.exp(2j*np.pi*(weights@weights[0]))
vals,U=np.linalg.eigh(M);basis=[]
for ma,ch in zip(target,chars):
 sub=U[:,abs(vals-ma*ma)<1e-9];ev,evec=np.linalg.eig(sub.T@O@sub);v=sub@evec[:,np.argmin(abs(ev-ch))];v=v/np.linalg.norm(v)
 assert np.linalg.norm(O@v-ch*v)<1e-8
 basis.append(v)
basis=np.array(basis);bad=[];err=0;nonzero=[]
for a,b,c in itertools.combinations_with_replacement(range(6),3):
 value=cubic(roots,am,basis,a,b,c)
 allowed=any(np.linalg.norm(x+y+z)<1e-8 for x in orbits[a] for y in orbits[b] for z in orbits[c])
 if (abs(value)>1e-8)!=allowed:bad.append((a+1,b+1,c+1))
 if allowed:nonzero.append((a+1,b+1,c+1));err=max(err,abs(abs(value)-4/np.sqrt(12)*area(target[[a,b,c]])))
print(f'PASS E6 all 56 cubic triples: nonzero={len(nonzero)}, orbit mismatches={bad}, max area error={err:.3g}');assert not bad and err<1e-10
print('E6 node masses:',target)
# The real orthonormal eigenbasis supplied by eigh is sufficient to exhibit failure.
rb=U.T;rm=np.sqrt(vals);violations=[]
for a,b,c in itertools.combinations_with_replacement(range(6),3):
 value=cubic(roots,am,rb,a,b,c);rhs=4/np.sqrt(12)*area(rm[[a,b,c]])
 if abs(value)>1e-8 and rhs>1e-8 and abs(abs(value)-rhs)>1e-7:violations.append((a,b,c,float(abs(value)/rhs)))
print('PASS E6 real-basis nonzero-coupling counterexamples:',violations[:4]);assert violations
