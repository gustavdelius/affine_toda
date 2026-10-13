import numpy as np
from scipy.special import loggamma as lg
from gz import rich
def logh(eps, mu, T, K=3000, pre=True):
    """'Even half' of GZ sigma: h(mu)h(-mu) ~ 1/(sin pi(eps+mu) sin pi(eps-mu)), h(T+mu)h(T-mu) ~ const."""
    def part(m, KK):
        l = np.arange(1, KK+1)
        return np.sum(-lg(1-eps+(2*l-2)*T-m)-lg(eps+2*l*T-m))
    p = np.log(np.sin(np.pi*eps)/np.sin(np.pi*(eps+mu))) if pre else 0
    return p + rich(lambda KK: part(mu, KK)-part(-mu, KK), K)
if __name__ == '__main__':
    T = 0.555; eps = 0.31+0.2j
    for mu in (0.13+0.07j, -0.21+0.11j, 0.4-0.2j):
        for pre in (True, False):
            u = np.exp(logh(eps, mu, T, pre=pre)+logh(eps, -mu, T, pre=pre))*np.sin(np.pi*(eps+mu))*np.sin(np.pi*(eps-mu))
            c = np.exp(logh(eps, T+mu, T, pre=pre)+logh(eps, T-mu, T, pre=pre))
            print(pre, 'h h(-)*S', np.round(u, 8), ' h(T+)h(T-)', np.round(c, 8))
