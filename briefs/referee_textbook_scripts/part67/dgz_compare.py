# Compare Part VI closed-form spinor R-matrices with Delius-Gould-Zhang.
import sympy as sp
x,q=sp.symbols('x q')
# --- ch.24: U_q(b_n^(1)) spinor x spinor, book: <l>_+ = (x-q^l)/(1-x q^l), eigenvalue on Lambda^{n-k}
def book24(n,m):
    k=n-m
    return sp.Mul(*[(x-q**(4*k-8*i+6))/(1-x*q**(4*k-8*i+6)) for i in range(1,(k+1)//2+1)])
# DGZ 1994 (3.25),(4.21)-(4.24): <a> = (1 - x q^a)/(x - q^a)
def dgzang(a,X,Q): return (1-X*Q**a)/(X-Q**a)
def dgz94(n,m,X,Q):
    if m==n: return sp.Integer(1)
    if n%2==0:
        if m%2==0: a=m//2; return sp.Mul(*[dgzang(4*i-1,X,Q) for i in range(1,n//2-a+1)])
        a=(m-1)//2; return sp.Mul(*[dgzang(4*i-3,X,Q) for i in range(1,n//2-a+1)])
    else:
        if m%2==0: a=m//2; return sp.Mul(*[dgzang(4*i-3,X,Q) for i in range(1,(n+1)//2-a+1)])
        a=(m-1)//2; return sp.Mul(*[dgzang(4*i-1,X,Q) for i in range(1,(n-1)//2-a+1)])
print("ch24 vs DGZ1994 (4.21)-(4.24) with x_DGZ=x, q_DGZ=q^-2:")
for n in range(2,10):
    ok=all(sp.simplify(book24(n,m)-dgz94(n,m,x,q**-2))==0 for m in range(0,n+1))
    ok2=all(sp.simplify(book24(n,m)-dgz94(n,m,1/x,q**2))==0 for m in range(0,n+1))
    print(f"  n={n}: identical: {ok}   (equivalently x->1/x, q_DGZ=q^2: {ok2})")
# --- ch.23: U_q(d_{n+1}^(2)) spinor x spinor, book: <l>_s = (x - s q^l)/(1 - s x q^l), on Lambda^{n-k}: prod_{j=1}^k <2j>_{(-1)^{j+1}}
def book23(n,m):
    k=n-m
    return sp.Mul(*[(x-(-1)**(j+1)*q**(2*j))/(1-(-1)**(j+1)*x*q**(2*j)) for j in range(1,k+1)])
# DGZ 1996 (4.8),(4.31) at a=1: <a>_s = (1 + s x q^a)/(x + s q^a), rho_k (on V0(lambda_k)=Lambda^k) = prod_{i=1}^{l-k} <(a-1)/2+i>_{(-1)^i}
def dgz96(l,m,X,Q):
    return sp.Mul(*[(1+(-1)**i*X*Q**i)/(X+(-1)**i*Q**i) for i in range(1,l-m+1)])
print("ch23 vs DGZ1996 (4.31), a=1, with x_DGZ=x, q_DGZ=q^-2:")
for n in range(2,9):
    ok=all(sp.simplify(book23(n,m)-dgz96(n,m,x,q**-2))==0 for m in range(0,n+1))
    ok2=all(sp.simplify(book23(n,m)-dgz96(n,m,1/x,q**2))==0 for m in range(0,n+1))
    print(f"  n={n}: identical: {ok}   (equivalently x->1/x, q_DGZ=q^2: {ok2})")
print("n=3 ch24 eigenvalues (Lambda^3..Lambda^0):",[sp.factor(book24(3,m)) for m in (3,2,1,0)])
