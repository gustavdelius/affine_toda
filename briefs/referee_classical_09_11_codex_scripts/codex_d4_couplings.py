import numpy as np
import itertools
I=np.eye(4)
r=np.array([I[0]-I[1],I[1]-I[2],I[2]-I[3],I[2]+I[3]])
a0=-I[0]-I[1]; ar=np.vstack([a0,r]); marks=np.array([1,1,2,1,1])
S=[I-np.outer(a,a) for a in r]; w=S[0]@S[2]@S[3]@S[1]
weights=np.linalg.inv(r).T
orbits=[np.array([np.linalg.matrix_power(w,k)@(I-w.T)@lam for k in range(6)]) for lam in weights]
Z1=np.diag([-1,1,1,-1])
Z4=np.array([r[2],r[1],r[0],a0]).T@np.linalg.inv(r.T)
chars=[(1,-1),(1,1),(-1,-1),(-1,1)]
basis=[]
for c1,c4 in chars:
    projector=(I+c1*Z1)@(I+c4*Z4)/4
    _,vec=np.linalg.eigh(projector)
    basis.append(vec[:,-1])
M=ar.T@np.diag(marks)@ar
masses=np.sqrt([v@M@v for v in basis])
errors=[]; maxerr=0; nonzero=[]
for a,b,c in itertools.combinations_with_replacement(range(4),3):
    C=np.sum(marks*(ar@basis[a])*(ar@basis[b])*(ar@basis[c]))
    allowed=any(np.linalg.norm(x+y+z)<1e-10 for x in orbits[a] for y in orbits[b] for z in orbits[c])
    if (abs(C)>1e-10)!=allowed:errors.append((a+1,b+1,c+1,C,allowed))
    if allowed:
        ma,mb,mc=masses[[a,b,c]]; semi=(ma+mb+mc)/2; area=np.sqrt(semi*(semi-ma)*(semi-mb)*(semi-mc))
        maxerr=max(maxerr,abs(abs(C)-4*area/np.sqrt(6)))
        nonzero.append((a+1,b+1,c+1))
print(('PASS' if not errors and maxerr<1e-12 else 'FAIL')+f' 14 D4 all 20 cubic triples: mismatches={errors}, max area error={maxerr:.3g}')
print('D4 masses in node order:',masses)
print('D4 nonzero unordered cubic triples:',nonzero)
