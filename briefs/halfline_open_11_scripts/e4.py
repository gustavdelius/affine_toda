import sympy as sp
y,y1,y2,xi=sp.symbols('y y1 y2 xi')
N,d=sp.symbols('N d',positive=True)
k0,k1,k2=sp.symbols('kappa0 kappa1 kappa2')
def F(k1,k2,y1,y2,N,d):
    z=y1/y2; w=y1*y2
    return k1*k2*w**d*(1-z**N)+k1*z**d*(w**N-1)+k2*z**(N-d)*(1-w**N)+w**(N-d)*(z**N-1)
for Nv,dv in [(2,1),(3,1),(3,2),(4,1),(4,2)]:
    sol=sp.solve(F(k0,k2,2,y,Nv,dv),k2)[0]
    sol=sp.factor(sp.simplify(sol))
    # reparametrize kappa0 so that solution is nice: try ansatz kappa(y)= y^(N-2d) * (xi*y^-N... )
    print(Nv,dv,sol)
    # test ansatz
    for name,ans in [('A', -y**(Nv-2*dv)*(xi - y**Nv)/(xi*y**Nv-1)), ('B', y**(-2*dv)*(xi*y**Nv - 1)/(xi- y**Nv)*y**Nv)]:
        ex=sp.simplify(F(ans.subs(y,y1),ans.subs(y,y2),y1,y2,Nv,dv))
        print('   ansatz',name,ex==0, sp.simplify(ans))
