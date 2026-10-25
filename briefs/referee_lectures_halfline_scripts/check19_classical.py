import sympy as s
from itertools import product
l=s.symbols('l',nonzero=True);C0,C1=s.symbols('C0 C1')
K=s.Matrix([[l*l-l**-2,-2*(C0*l-C1/l)],[-2*(C1*l-C0/l),l*l-l**-2]])
def gen(n,a,b):
    E=s.zeros(n);E[a,b]=1;F=E.T;H=E*F-F*E
    return E,F,H
for c,(a,b) in zip([C1,C0],[(0,1),(1,0)]):
    E,F,H=gen(2,a,b);res=K*(-c*H+l*E-F/l)-(c*H+l*E-F/l)*K
    assert all(s.simplify(x)==0 for x in res)
for c0,c1,c2 in product([-1,1],repeat=3):
    sig=c0*c1*c2;d=sig*l*l-1/l
    K=s.Matrix([[d,2*c1,-2*c1*c2*l],[-2*c0*c2*l,d,2*c2],[2*c0,-2*c0*c1*l,d]])
    for c,(a,b) in zip([c1,c2,c0],[(0,1),(1,2),(2,0)]):
        E,F,H=gen(3,a,b);res=K*(-c*H+l*E-F/l)-(c*H+l*E-F/l)*K
        assert all(s.simplify(x)==0 for x in res)
    assert s.simplify(K.det())!=0
print('PASS 19: exact symbolic node equations for arbitrary a1 parameters and all 8 a2 patterns')
