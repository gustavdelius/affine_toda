"""Book-style G_j form of the SP scalar factor: G_j(x, c) = Gamma(x + j h omega + c)/Gamma(-x + j h omega + c).
 Lam_2T(mu;c1)/Lam_2T(mu;c2) = Gamma(c1+mu)Gamma(c2-mu)/(Gamma(c1-mu)Gamma(c2+mu)) prod_j G_j(mu,c1)G_j(mu,1-c1)/(G_j(mu,c2)G_j(mu,1-c2))
 rho = (omega/2)_I Gamma(2mu+2T+omega)Gamma(-2mu)/(Gamma(-2mu+2T+omega)Gamma(2mu)) prod_j G_2j(2mu,2T+omega)G_2j(2mu,1-2T-omega)/(G_2j(2mu,0)G_2j(2mu,1))"""
import mpmath as mp
from sp_mp import SP, logprod, blk
mp.mp.dps = 30
G = lambda x, c, s=1: [(s, x + c), (-s, -x + c)]
err = 0
for n, l in ((1, 1), (2, 0), (2, 1), (3, 2)):
    S = SP(n, mp.mpf('0.37'), l, mp.mpc('0.41', '0.23'))
    T, om = S.T, S.om
    for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07')):
        c1, c2 = S.eps + 2*T, S.epsb + T
        ratio = (mp.gamma(c1 + mu)*mp.gamma(c2 - mu)/(mp.gamma(c1 - mu)*mp.gamma(c2 + mu))
                 * mp.exp(logprod(2*T, G(mu, c1) + G(mu, 1 - c1) + G(mu, c2, -1) + G(mu, 1 - c2, -1))))
        rho = (blk(om/2, mu, T)*mp.gamma(2*mu + 2*T + om)*mp.gamma(-2*mu)/(mp.gamma(-2*mu + 2*T + om)*mp.gamma(2*mu))
               * mp.exp(logprod(4*T, G(2*mu, 2*T + om) + G(2*mu, 1 - 2*T - om) + G(2*mu, 0, -1) + G(2*mu, 1, -1))))
        a_G = rho*mp.exp((1 - l)*om*1j*mp.pi*mu/T)*1j*ratio/mp.sin(mp.pi*(S.eps + mu))
        err = max(err, abs(a_G/S.a(mu) - 1), abs(rho/S.rho(mu) - 1))
print('G-form of a and rho vs Lam-form:', mp.nstr(err, 3))
