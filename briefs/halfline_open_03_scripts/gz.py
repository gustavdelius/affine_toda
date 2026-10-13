import numpy as np
from scipy.special import loggamma as lg
from rhonp import rho_plus, blk
def rich(f,K): return 2*f(2*K)-f(K)
def logR0(u,lam,K=4000):
    def part(uu,KK):
        l=np.arange(1,KK+1); z=2*lam*uu/np.pi
        return np.sum(lg(4*l*lam-z)+lg(4*lam*(l-1)+1-z)-lg((4*l-3)*lam-z)-lg((4*l-1)*lam+1-z))
    return rich(lambda KK: part(u,KK)-part(-u,KK),K)
def logsigma(x,u,lam,K=4000):
    def part(uu,KK):
        l=np.arange(1,KK+1); z=lam*uu/np.pi
        return np.sum(lg(0.5+x/np.pi+(2*l-1)*lam-z)+lg(0.5-x/np.pi+(2*l-1)*lam-z)-lg(0.5-x/np.pi+(2*l-2)*lam-z)-lg(0.5+x/np.pi+2*l*lam-z))
    pre=np.log(np.cos(x)/np.cos(x+lam*u))
    return pre+rich(lambda KK: part(u,KK)-part(-u,KK),K)
def A_gz(u,lam,eta,vth=0.0):
    return np.exp(logR0(u,lam)+logsigma(eta,u,lam)+logsigma(1j*vth,u,lam))
def A1(mu,n,lam,K=4000):
    h=n+1
    G=lambda j,m,a: lg(m+j*h*lam+a)-lg(-m+j*h*lam+a)
    def part(KK):
        k=np.arange(1,KK+1)
        return np.sum(G(2*k,2*mu,-1.5*h*lam+1)+G(2*k,2*mu,-0.5*h*lam)-G(2*k,2*mu,-1.5*n*lam-0.5*lam+1)-G(2*k,2*mu,-0.5*n*lam-1.5*lam))
    T=h*lam/2
    return blk(-h*lam/4,mu,T)*np.exp(rich(part,K))
if __name__=='__main__':
    lam=0.37; n=1; T=lam
    for th in (0.3+0.2j,-0.5+0.4j,0.8-0.1j):
        u=-1j*th; mu=-1j*T*th/np.pi
        eta=np.pi*(lam+1)/2
        AG=A_gz(u,lam,eta)*np.cos(lam*u)
        Am=A1(mu,n,lam)*rho_plus(mu,n,lam)
        f=-blk(T/2,mu,T)**2; g=blk(T/2-lam/2,mu,T)*blk(T/2+lam/2,mu,T)
        print(th,'A_GZ/(A_d rho)=',AG/Am,'  1/(f g)=',1/(f*g))
        # diagonal check: DG A1 vs GZ at eta=pi/2 (matrix cos(lam u) * 1 after /cos eta -> limit)
