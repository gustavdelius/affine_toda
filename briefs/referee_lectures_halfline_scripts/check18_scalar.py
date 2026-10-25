import mpmath as m
m.mp.dps=45
om=m.sqrt(2);h=3;n=2;T=h*om/2;mu=m.mpc('.29','-.08')
def ib(a,u):return m.sin(m.pi*(u+a)/(h*om))/m.sin(m.pi*(u-a)/(h*om))
def A(u,N):
    shifts=[-m.mpf('1.5')*h*om+1,-h*om/2,-m.mpf('1.5')*n*om-om/2+1,-n*om/2-m.mpf('1.5')*om]
    total=m.mpc(0)
    for k in range(1,N+1):
        for a,sign in zip(shifts,[1,1,-1,-1]):
            z=2*k*h*om+a
            total+=sign*(m.loggamma(z+2*u)-m.loggamma(z-2*u))
    return ib(-h*om/4,u)*m.exp(total)
target=ib(-n*om/2,mu)*ib(m.mpf('.5')-h*om/2,mu)*ib(-om/2-m.mpf('.5'),mu)
errs=[]
for N in [30,60,120,300]:
    d=T/2-m.mpf('.5')
    val=A(mu-d,N)*A(mu+d,N)*A(T/2-mu,N)/A(T/2+mu,N)
    err=abs(val/target-1);errs.append(err)
    print('18 scalar cutoff',N,'relative residual',m.nstr(err,12))
assert all(a>b for a,b in zip(errs,errs[1:]))
print('PASS 18: 45-digit Gamma product fusion converges to the three-block factor')
