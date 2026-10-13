import sympy as sp
y1,y2,q=sp.symbols('y1 y2 q',nonzero=True)
N,d=sp.symbols('N d',positive=True)
k1,k2=sp.symbols('kappa1 kappa2')
def M(z):
    Z=z**N; c=(q-1/q)/(q*Z-1/q); t=(Z-1)/(q*Z-1/q)
    return sp.Matrix([[c*z**(N-d),t],[t,c*z**d]])
D=lambda k: sp.diag(1,k)
z1=y1/y2; z2=y1*y2
E=M(z1)*D(k1)*M(z2)*D(k2)-D(k2)*M(z2)*D(k1)*M(z1)
for i in range(2):
    for j in range(2):
        e=sp.factor(sp.simplify(E[i,j]))
        print(i,j,e)
