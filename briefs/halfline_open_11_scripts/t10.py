import numpy as np
from scipy.special import loggamma as lg
from gz import rich
def logsig(x, mu, T, K=3000):
    # GZ sigma in mu variable: z = mu, lam -> T
    def part(m, KK):
        l = np.arange(1, KK+1)
        return np.sum(lg(0.5+x/np.pi+(2*l-1)*T-m)+lg(0.5-x/np.pi+(2*l-1)*T-m)-lg(0.5-x/np.pi+(2*l-2)*T-m)-lg(0.5+x/np.pi+2*l*T-m))
    return np.log(np.cos(x)/np.cos(x+np.pi*mu)) + rich(lambda KK: part(mu, KK)-part(-mu, KK), K)
T = 0.555; x = 0.3+0.2j
for mu in (0.13+0.07j, -0.21+0.11j):
    u = np.exp(logsig(x, mu, T)+logsig(x, -mu, T))*np.cos(x+np.pi*mu)*np.cos(x-np.pi*mu)/np.cos(x)**2
    c = np.exp(logsig(x, T+mu, T)+logsig(x, T-mu, T))
    print('unit', u, ' sigma(T+mu)sigma(T-mu)=', c, ' *cos cos', c*np.cos(x+np.pi*mu)*np.cos(x-np.pi*mu), c*np.cos(x+np.pi*mu+np.pi*T)*np.cos(x-np.pi*mu+np.pi*T)/np.cos(x)**2)
