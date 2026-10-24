import sympy as s
import itertools
p,q,r=s.symbols('p q r'); rt=s.sqrt(2); I=s.I
variables=[p,q,r]; masses=[rt,2,rt]; species=[1,2,3]
A={(0,1):3-2*rt,(0,2):s.Rational(1,2),(1,2):3-2*rt}
taus=[]
for j in range(4):
    tau=0
    for mask in itertools.product([0,1],repeat=3):
        term=1
        for k in range(3):
            if mask[k]:
                term*=I**(j*species[k])*variables[k]
                for l in range(k):
                    if mask[l]: term*=A[l,k]
        tau+=term
    taus.append(s.expand(tau))
def dx(f):return s.expand(sum(m*v*s.diff(f,v) for m,v in zip(masses,variables)))
res=[s.simplify(s.expand(t*dx(dx(t))-dx(t)**2-t*t+taus[(j-1)%4]*taus[(j+1)%4])) for j,t in enumerate(taus)]
print(('PASS' if res==[0,0,0,0] else 'FAIL')+' 19 exact three-species static Hirota residual: '+str(res))
# Tangent field vectors are i sum_j alpha_j i^(ja), with independent nonzero velocity multiples.
roots=s.Matrix([[-1,0,0,1],[1,-1,0,0],[0,1,-1,0],[0,0,1,-1]])
tangent=s.Matrix.hstack(*[I*roots.T*s.Matrix([I**(j*a) for j in range(4)]) for a in species])
print(('PASS' if tangent.rank()==3 else 'FAIL')+' 19 unstable-manifold tangent rank: '+str(tangent.rank()))
