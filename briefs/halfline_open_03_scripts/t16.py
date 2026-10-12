import numpy as np
from rhonp import rho_plus, blk
from gz import A1
def order(F,x0):
    e1,e2=1e-4,1e-5
    return np.log(abs(F(x0+e2))/abs(F(x0+e1)))/np.log(e1/e2)   # >0 pole order
for om in (0.37,1.3,2.7):
  for n in (1,2,3,4):
    T=(n+1)*om/2
    sig=lambda m: blk(T/2,m,T)**2*blk(T/2-om/2,m,T)*blk(T/2+om/2,m,T)
    Ad=lambda m: A1(m,n,om,K=1500); r=lambda m: rho_plus(m,n,om,K=1500)
    out=[]
    for x0,lab in ((T/2,'th=i pi/2'),(T/2-om/2,'T/2-om/2')):
        out.append((lab, round(order(Ad,x0),2), round(order(r,x0),2), round(order(sig,x0),2)))
    print('om',om,'n',n,' [pole orders of A_d, rho, sigma]',out)
