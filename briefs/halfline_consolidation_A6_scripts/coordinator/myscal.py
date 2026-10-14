"""Coordinator check: the proposed eq-dn-nondiagonal-scalar, transcribed from the A6 report's book text,
satisfies A(mu)A(-mu)U = 1 and A(T/2-mu)/A(T/2+mu) = F(2mu)/Psi (book A, F of dn_uniform.py; Psi of ch. 33)."""
import mpmath as mp
src = open('dn_uniform.py').read()
ns = {'mp': mp}
exec(src.split('# ---------------------------------------------------------------- 2./3. scalar factors')[1].split('\ndev2, dev3')[0],
     ns)
Ffac, Afac, blk = ns['Ffac'], ns['Afac'], ns['blk']
mp.mp.dps = 30
def logprod1(P, terms, J=60, M=40):   # j >= 1 (an_nondiagonal_scalar convention)
    head = sum(c*mp.loggamma(j*P + a) for j in range(1, J) for c, a in terms)
    tail = sum((-1)**(k + 1)*sum(c*mp.bernpoly(k + 1, a) for c, a in terms)/(k*(k + 1)*P**k)*mp.zeta(k, J)
               for k in range(2, M + 1))
    return head + tail
Gt = lambda m, a, s=1: [(s, m + a), (-s, -m + a)]
def rho(sg, mu, T, om):
    if sg == 1:
        pre = blk(-T/2, mu, T)*mp.sin(mp.pi*om/2)/mp.sin(mp.pi*(mu - T/2)); hf = 0
    else:
        pre = -1j*mp.sin(mp.pi*om/2)/mp.cos(mp.pi*(mu - T/2)); hf = mp.mpf(1)/2
    if sg == 1:
        t = Gt(mu, om/2 - 1.5*T) + Gt(mu, 1 - om/2 - T/2) + Gt(mu, om/2 - T/2, -1) + Gt(mu, 1 - om/2 - 1.5*T, -1)
    else:
        t = Gt(mu, hf + om/2 - 1.5*T) + Gt(mu, hf - om/2 - T/2) + Gt(mu, hf + om/2 - T/2, -1) + Gt(mu, hf - om/2 - 1.5*T, -1)
    return pre*mp.exp(logprod1(2*T, t))
def vterms(sg, x, s, om):          # Gamma terms of v(x)^s; x = shift (without the jP part)
    if sg == 1:
        return [(s, x + om/2), (s, 1 + x - om/2), (-s, x), (-s, 1 + x)]
    return [(s, mp.mpf(1)/2 + x + om/2), (s, mp.mpf(1)/2 + x - om/2), (-2*s, mp.mpf(1)/2 + x)]
def v(sg, x, om):
    return mp.exp(sum(c*mp.loggamma(a) for c, a in vterms(sg, x, 1, om)))
def gfac(sg, mu, T, om):
    # prod_{k>=1} [v(mu+kT)/v(-mu+kT)]^{(-1)^k}: k = 2j (+), k = 2j-1 (-), j >= 1
    t = vterms(sg, mu, 1, om) + vterms(sg, -mu, -1, om) + vterms(sg, mu - T, -1, om) + vterms(sg, -mu - T, 1, om)
    return v(sg, mu, om)*mp.exp(logprod1(2*T, t))
def U(sg, mu, T, om):
    c2 = mp.cos(2*mp.pi*mu)
    return sg*(c2 - sg*mp.cos(mp.pi*T))*(c2 - sg*mp.cos(mp.pi*om))/(2*mp.sin(mp.pi*om/2)**2*(c2 - sg))
Psi = {1: lambda mu, om: mp.sin(mp.pi*(om/2 - mu))/mp.sin(mp.pi*(om/2 + mu)),
       -1: lambda mu, om: mp.cos(mp.pi*(mu - om/2))/mp.cos(mp.pi*(mu + om/2))}
for n, om, mu in ((4, mp.mpf('0.37'), mp.mpc('0.13', '0.21')), (5, mp.mpf('0.53'), mp.mpc('-0.17', '0.11')),
                  (6, mp.mpf('0.29'), mp.mpc('0.21', '-0.07'))):
    T = (n - 1)*om
    for sg in (1, -1):
        As = lambda m: Afac(m, T, om)*rho(sg, m, T, om)*gfac(sg, m, T, om)
        e1 = abs(As(mu)*As(-mu)*U(sg, mu, T, om) - 1)
        e2 = abs(As(T/2 - mu)/As(T/2 + mu)/(Ffac(2*mu, T, om)/Psi[sg](mu, om)) - 1)
        print(n, sg, 'unitarity', mp.nstr(e1, 3), ' crossing', mp.nstr(e2, 3))
