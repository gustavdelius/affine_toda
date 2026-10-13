import mpmath as mp
from pr_core import Prod, K_uniform, axioms, THS, blk
def Kmin(a, n): return Prod(n + 1, {(a, 0): 1, (n + 1 - a, 0): 1}, -1)*K_uniform(a, n)
for n in range(1,5):
    print(n, axioms(lambda a, n=n: Kmin(a, n), n, 'm'))
    for a in range(1,n+1):
        th=THS[0]; print('  a',a,'K(B=0)=',mp.nstr(Kmin(a,n)(th,0),8),' Kunif(0)*(-(a)(h-a))=',mp.nstr(K_uniform(a,n)(th,0)*(-blk(a,th,n+1)*blk(n+1-a,th,n+1)),8))
