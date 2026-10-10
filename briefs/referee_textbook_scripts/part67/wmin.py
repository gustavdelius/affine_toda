# Central charges and conformal weights for the table of sec-perturbed-minimal.
from fractions import Fraction as F
import itertools, numpy as np
def c_W(n,p,pp):  # W_{n+1}(p,p') minimal model, rank n
    return n*(1-F((n+1)*(n+2)*(pp-p)**2,p*pp))
def h_W(n,p,pp,lam,mu):
    # h(lam|mu) = [ (p'(lam+rho) - p(mu+rho))^2 - (p'-p)^2 rho^2 ] / (2 p p'), sl_{n+1}, weights in eps basis
    N=n+1
    def vec(dyn):  # Dynkin labels -> vector in R^N (traceless)
        v=np.zeros(N)
        for i,a in enumerate(dyn): v[:i+1]+=a
        return v-v.mean()
    rho=vec([1]*n)
    L=pp*(vec(lam)+rho)-p*(vec(mu)+rho)
    return F(int(round((L@L-(pp-p)**2*(rho@rho))*1e6)),int(2*p*pp*1e6)).limit_denominator(10**6)
adj=lambda n:[1]+[0]*(n-2)+[1] if n>1 else [2]
zero=lambda n:[0]*n
print("W3(4,5): c =",c_W(2,4,5)," h(0|adj) =",h_W(2,4,5,zero(2),adj(2)),"(3-state Potts energy: 2/5)")
print("Virasoro M(3,4): c =",c_W(1,3,4),"; h_{1,2} =",h_W(1,3,4,[0],[1]),"; h_{1,3} =",h_W(1,3,4,[0],[2]))
print("M(2,5): c =",c_W(1,2,5),"; h_{1,2} =",h_W(1,2,5,[0],[1]),"; h_{1,3} =",h_W(1,2,5,[0],[2]))
# compare h(0|adj) with the dimension of the affine-root vertex operator: (n+1) b - n, b = beta^2/4pi = p/p'
for n in (1,2,3,4):
    for p,pp in [(n+2,n+3),(n+2,n+4),(n+3,2*n+5)]:
        from math import gcd
        if gcd(p,pp)!=1: continue
        hv=(n+1)*F(p,pp)-n
        print(f"n={n} (p,p')=({p},{pp}): h(0|adj)={h_W(n,p,pp,zero(n),adj(n))}, (n+1)p/p'-n={hv}, beta^2/4pi={F(p,pp)}, beta_UV^2/4pi={F(n,n+1)}")
