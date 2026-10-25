"""Direct Gamma products with cutoff extrapolation, independently transcribed."""
import mpmath as mp
mp.mp.dps=38
j=mp.j; pi=mp.pi
def extrap(vals):
    # Polynomial interpolation in 1/N at zero; N doubles each entry.
    for k in range(1,len(vals)):
        vals=[(2**k*vals[i+1]-vals[i])/(2**k-1) for i in range(len(vals)-1)]
    return vals[0]
def sigma_product(lam,x,u,base=32,levels=5):
    w=lam*u/pi; X=x/pi; total=0; out=[]
    cuts=[base*2**k for k in range(levels)]
    for ell in range(1,cuts[-1]+1):
        aa=[mp.mpf('.5')+X+(2*ell-1)*lam,mp.mpf('.5')-X+(2*ell-1)*lam,
            mp.mpf('.5')-X+(2*ell-2)*lam,mp.mpf('.5')+X+2*ell*lam]
        total+=sum(sign*(mp.loggamma(a-w)-mp.loggamma(a+w)) for sign,a in zip([1,1,-1,-1],aa))
        if ell in cuts: out.append(total)
    return mp.cos(x)/mp.cos(x+lam*u)*mp.exp(extrap(out))
def phi_product(lam,theta,base=32,levels=5):
    w=j*lam*theta/(2*pi); total=0; out=[]
    cuts=[base*2**k for k in range(levels)]
    for k in range(cuts[-1]):
        aa=[k*lam+mp.mpf('.25'),(k+1)*lam+mp.mpf('.75'),(k+mp.mpf('.5'))*lam+mp.mpf('.25'),(k+mp.mpf('.5'))*lam+mp.mpf('.75')]
        total+=sum(sign*(mp.loggamma(a-w)-mp.loggamma(a+w)) for sign,a in zip([1,1,-1,-1],aa))
        if k+1 in cuts: out.append(total)
    return mp.exp(j*pi*(1+lam)/4+extrap(out))/(1+j*mp.exp(lam*theta))
def sigma_integral(lam,x,u):
    def f(t):
        if not t: return -2*lam*u/pi*(x/pi+lam/2)
        return -mp.sinh(lam*u*t/pi)*mp.sinh(x*t/pi+lam*t/2)/(t*mp.sinh(t/2)*mp.cosh(lam*t/2))
    return mp.cos(x)/mp.cos(x+lam*u)*mp.exp(mp.quad(f,[0,1,4,16,mp.inf]))
def report(name,err,tol=mp.mpf('1e-8')):
    print(('PASS ' if err<tol else 'FAIL ')+name+' residual='+mp.nstr(err,7),flush=True)
for lam in [mp.sqrt(2),mp.e/2,mp.sqrt(5)]:
    x=mp.mpc('.21','.31'); u=mp.mpc('.13','-.19')
    si=sigma_integral(lam,x,u)
    p1=sigma_product(lam,x,u); p2=sigma_product(lam,x,u,64)
    report('sigma integral lambda='+mp.nstr(lam,8),abs(si-p2))
    report('sigma cutoff stability',abs(p1-p2))
    for theta in [mp.mpf('.43'),mp.mpc('.37','.16')]:
        td=mp.mpf('.61'); yd=mp.exp(lam*td); y=mp.exp(lam*theta)
        lhs=sigma_product(lam,j*lam*td,-j*theta,64)
        rhs=(1+yd**-2)/(1+y**-2*yd**-2)*phi_product(lam,theta-td,64)/phi_product(lam,-theta-td,64)
        report('fusion Gamma identity theta='+str(theta),abs(lhs-rhs))
lam=mp.sqrt(2); theta=mp.mpf('.43')
for v in [2,4,8,12]:
    sig=sigma_integral(lam,j*v,-j*theta)
    print('Dirichlet limit vartheta='+str(v)+' abs(sigma-1)='+mp.nstr(abs(sig-1),8),flush=True)
