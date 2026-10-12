"""s = -1 scalar: rho_-(mu)rho_-(-mu) = -sin^2(pi b)/(cos pi(a+mu) cos pi(a-mu)),
rho_-(a-mu)/rho_-(a+mu) = cos pi(mu+b)/cos pi(mu-b), a = T/2, b = om/2."""
import numpy as np
from scipy.special import loggamma as lg
def lpsi(c,mu): return lg(c+mu)-lg(c-mu)
# cos pi(b+mu)/cos pi(b-mu) = gamma_{1/2+b}(mu)/gamma_{1/2-b}(mu)   [cos pi x = pi/(Gamma(1/2+x)Gamma(1/2-x))]
# rho_e = sin(pi b)/cos(pi(a+mu)) * i ? : rho_e(mu)rho_e(-mu) = sin^2/(cos cos); need extra -1 -> rho_e = i sin(pi b)/cos(pi(a+mu))
# rho_e(a+mu)/rho_e(a-mu) = cos pi(T-mu)/cos pi(T+mu) = gamma_{1/2+T}(-mu)/gamma_{1/2-T}(-mu) = gamma_{1/2-T}(mu)/gamma_{1/2+T}(mu)
# H = c(mu) * rho_e(a+mu)/rho_e(a-mu) = gamma_{1/2+b} gamma_{1/2-T} / (gamma_{1/2-b} gamma_{1/2+T})
def log_rho_o(mu,a,x0s,sgn,K):
    k=np.arange(2*K); c=(2*k+1)*a; sg=(-1.0)**k
    return np.sum(sg*sum(s*lpsi(x+c,mu) for x,s in zip(x0s,sgn)))
def rho_minus(mu,n,om,K=4000):
    T=(n+1)*om/2; a=T/2; b=om/2
    x0s=[0.5+b,0.5-T,0.5-b,0.5+T]; sgn=[1,1,-1,-1]
    lo=2*log_rho_o(mu,a,x0s,sgn,2*K)-log_rho_o(mu,a,x0s,sgn,K)
    return 1j*np.sin(np.pi*b)/np.cos(np.pi*(a+mu))*np.exp(lo)
def cproj_minus(mu,om,T):
    return np.sin(np.pi*mu)*np.sin(np.pi*(mu+om/2))*np.sin(np.pi*(mu-T))/(np.sin(np.pi*om/2)**2*np.sin(np.pi*(mu-om/2)))
if __name__=='__main__':
    om=0.37
    for n in (1,2,3):
        T=(n+1)*om/2; a=T/2; b=om/2; mu=0.13+0.21j
        r=lambda m: rho_minus(m,n,om)
        print(n,'unit',r(mu)*r(-mu)*np.cos(np.pi*(a+mu))*np.cos(np.pi*(a-mu))/(-np.sin(np.pi*b)**2),
              'cross',r(a-mu)/r(a+mu)*np.cos(np.pi*(mu-b))/np.cos(np.pi*(mu+b)))
