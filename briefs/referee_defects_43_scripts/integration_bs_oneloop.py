"""Integration check, ch. 43 item 3: Bajnok-Simon 0712.4292 Sec. 4.2 one-loop factor = O(b^2) term of the exact factor
with theta_d = (1 - B/2) theta_d^cl, B/2 = b^2/(8 pi + b^2); and BS (17) => B_BS vartheta = (1 - B_BS) mu exactly."""
import sympy as sp
th, mu, t = sp.symbols('theta mu t')       # t = b^2/(8 pi)
I = sp.I
T = lambda m: -I*sp.sinh(th/2 - I*sp.pi/4 + m/2)/sp.sinh(th/2 + I*sp.pi/4 + m/2)   # BS classical / (16) with B vartheta -> mu
d = sp.simplify((sp.diff(T(mu), mu) + 1/(1 - I*sp.sinh(th + mu))).rewrite(sp.exp))
print('d/dmu T + 1/(1 - i sinh(theta + mu)) =', d)
Bhalf = t/(1 + t)                           # b^2/(8 pi + b^2)
mu_q = (1 - Bhalf)*mu                       # B_BS vartheta = -theta_d, theta_d = (1 - B/2) theta_d^cl, mu = -theta_d^cl
ser = sp.series(T(mu_q), t, 0, 2).removeO()
oneloop = T(mu) + mu*t/(1 - I*sp.sinh(th + mu))
diff = sp.simplify((ser - oneloop).rewrite(sp.exp))
print('O(b^2) term minus BS one-loop:', diff)
# BS (17): M e^{b^2 vartheta/8 pi} = M_cl e^{mu}, M -> M_cl: b^2 vartheta/(8 pi) = mu, so vartheta = mu/t,
# B_BS vartheta = t/(1+t) * mu/t = mu/(1+t) = (1 - B_BS) mu; (23): mu (1 - t + ...)
print('B_BS vartheta - (1 - B_BS) mu =', sp.simplify(Bhalf*mu/t - (1 - Bhalf)*mu), '; series', sp.series(Bhalf*mu/t, t, 0, 2))
assert d == 0 and diff == 0
