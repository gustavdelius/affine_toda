"""Side check (outside this round's diff): is the a_3 CDD factor P_{w/2} = (T/2-w/2)_I (T/2+w/2)_I of
eq-an-nondiagonal-scalar the only one that reproduces Kim's a_3 one loop, or does P_{1/2+w/2} also?
And does P_{1/2+w/2} satisfy the a_n soliton-bootstrap condition sigma(mu) = prod_{k=1}^n sigma(mu + T - k w)?"""
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
half = mp.mpf(1)/2
def bI(x, mu, T): return mp.sin(pi*(mu + x)/(2*T))/mp.sin(pi*(mu - x)/(2*T))
n = 3; h = n + 1
def rhs(m, T, o):      # eq-an-neumann-breather RHS (with the book's CDD factor included in A^(+))
    return bI(-n*o/2, m, T)*bI(n*o/2 - half, m, T)*bI(half, m, T)
def Px(x): return lambda m, T, o: bI(T/2 - x(T, o), m, T)*bI(T/2 + x(T, o), m, T)
Pw = Px(lambda T, o: o/2); Pwh = Px(lambda T, o: half + o/2)
def alt(m, T, o):
    ma, mb = m - (T/2 - half), m + (T/2 - half)
    return rhs(m, T, o)*Pwh(ma, T, o)*Pwh(mb, T, o)/(Pw(ma, T, o)*Pw(mb, T, o))
def cont(f, th, B):
    o = -2/B; T = h*o/2; return f(-1j*T*th/pi, T, o)
def kim(th):   # Kim (33) [a_3, hep-th/9506031 appendix], d/dB with beta^2 = 2 pi B
    sh, ch = mp.sinh(th), mp.cosh(th)
    return 1j*2*pi/16*(sh/(ch + 1/mp.sqrt(2)) - sh/(ch - 1))
eps = mp.mpf('1e-14')
for name, f in (('book P_{w/2}', rhs), ('alt P_{1/2+w/2}', alt)):
    dev = 0
    for th in (mp.mpc('0.41', '0.17'), mp.mpc('-1.1', '0.33'), mp.mpc('0.7', '-0.52')):
        g = lambda B: cont(f, th, B)
        d1 = (mp.log(g(eps)) - mp.log(g(2*eps)))/(-eps)
        dev = max(dev, abs(g(eps) - 1), abs(d1 - kim(th)))
    print(f"a_3 {name}: max(|K(0)-1|, |dlogK/dB - Kim (33)|) = {mp.nstr(dev, 3)}")
# soliton bootstrap condition for sigma = P_{1/2+w/2} and P_{w/2}
for name, sig in (('P_{w/2}', Pw), ('P_{1/2+w/2}', Pwh), ('-P_{1/2+w/2}', lambda m, T, o: -Pwh(m, T, o))):
    o = mp.mpf('0.53'); T = h*o/2; mu = mp.mpc('0.17', '0.11')
    rhsb = 1
    for k in range(1, n + 1):
        rhsb *= sig(mu + T - k*o, T, o)
    print(f"soliton bootstrap sigma(mu)/prod = {mp.nstr(sig(mu, T, o)/rhsb, 8)} for sigma = {name}")
