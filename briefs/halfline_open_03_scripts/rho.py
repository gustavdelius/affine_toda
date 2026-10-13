import mpmath as mp
mp.mp.dps=30
lg=mp.loggamma
def lpsi(c,mu): return lg(c+mu)-lg(c-mu)
def log_rho_o(mu,a,b,T,K):
    """sum_{k<2K} (-1)^k [psi_{b+(2k+1)a} psi_{1-T+(2k+1)a} / (psi_{1-b+(2k+1)a} psi_{T+(2k+1)a})]"""
    s=0
    for k in range(2*K):
        c=(2*k+1)*a
        s+=(-1)**k*(lpsi(b+c,mu)+lpsi(1-T+c,mu)-lpsi(1-b+c,mu)-lpsi(T+c,mu))
    return s
def rho_plus(mu,n,om,K=200):
    h=n+1; T=mp.mpf(h)*om/2; a=T/2; b=mp.mpf(om)/2
    f=lambda KK: log_rho_o(mu,a,b,T,KK)
    lo=2*f(2*K)-f(K)     # Richardson (tail ~ 1/K)
    return mp.sin(mp.pi*b)/mp.sin(mp.pi*(a+mu))*mp.exp(lo)
def blk(x,mu,T):  # DG block (x)_I, period 2T = h om
    return mp.sin(mp.pi*(mu+x)/(2*T))/mp.sin(mp.pi*(mu-x)/(2*T))
def cproj_plus(mu,om,T):
    return -mp.cos(mp.pi*mu)*mp.cos(mp.pi*(mu+om/2))*mp.cos(mp.pi*(mu-T))/(mp.sin(mp.pi*om/2)**2*mp.cos(mp.pi*(mu-om/2)))
if __name__=='__main__':
    om=mp.mpf('0.37'); 
    for n in (1,2,3):
        h=n+1; T=h*om/2; a=T/2; b=om/2
        mu=mp.mpc('0.13','0.21')
        r=lambda m: rho_plus(m,n,om)
        u=r(mu)*r(-mu)*mp.sin(mp.pi*(a+mu))*mp.sin(mp.pi*(a-mu))/mp.sin(mp.pi*b)**2
        c=r(a-mu)/r(a+mu)*mp.sin(mp.pi*(b-mu))/mp.sin(mp.pi*(b+mu))
        print(n,'unitarity',mp.nstr(u,12),' crossing',mp.nstr(c,12))
