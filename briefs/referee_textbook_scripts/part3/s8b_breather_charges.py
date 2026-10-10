"""Q_2 (spin-2 current of s8, basis 1) and E on a_2^(1) solutions, with exact derivatives of the tau-functions:
single solitons with real rapidity, and a breather {(1, th0+iu), (2, th0-iu)} (Theta-invariant data)."""
import sympy as sp, mpmath as mp
exec(open('s8_spin2_charge.py').read().split('mp.mp.dps = 25')[0])
sub = {fr: (1 if i == 1 else 0) for i, fr in enumerate(free)}
Ts = sp.expand(T.subs(sol).subs(sub)); Ths = sp.expand(Th.subs(sol).subs(sub))
Tf = sp.lambdify((b, m, *p, *q, *X), Ts, 'mpmath'); Thf = sp.lambdify((b, m, *p, *q, *X), Ths, 'mpmath')
mp.mp.dps = 20
beta = mp.mpf('0.9'); h = 3; w = mp.e**(2j*mp.pi/3); m1 = 2*mp.sin(mp.pi/3)


def make(spec, ths, xis):
    A = 0
    if len(spec) == 2:
        c = mp.cosh(ths[0] - ths[1]); A = (c - mp.cos(mp.pi*(spec[0]-spec[1])/h))/(c - mp.cos(mp.pi*(spec[0]+spec[1])/h))

    def taus(x, t, op):      # op: 'p' -> d_+ = d_t + d_x, 'x', 't'; returns [(tau, d tau, d^2 tau)]_j
        rate = {'p': [m1*mp.e**(-th) for th in ths], 'x': [m1*mp.cosh(th) for th in ths], 't': [-m1*mp.sinh(th) for th in ths]}[op]
        Es = [mp.e**(m1*(x*mp.cosh(th) - t*mp.sinh(th)) + xi) for th, xi in zip(ths, xis)]
        out = []
        for j in range(3):
            terms = [(1, 0)] + [(w**(j*a)*E, r) for a, E, r in zip(spec, Es, rate)]
            if len(spec) == 2:
                terms.append((A*w**(j*(spec[0]+spec[1]))*Es[0]*Es[1], rate[0] + rate[1]))
            out.append((sum(c for c, r in terms), sum(c*r for c, r in terms), sum(c*r*r for c, r in terms)))
        return out
    return taus


def dphi(tt, k):
    d1 = 1j/beta*sum(al[j][k]*tt[j][1]/tt[j][0] for j in range(3))
    d2 = 1j/beta*sum(al[j][k]*(tt[j][2]/tt[j][0] - (tt[j][1]/tt[j][0])**2) for j in range(3))
    return d1, d2


def charges(taus, t0=mp.mpf(0)):
    def dens(x, which):
        tp = taus(x, t0, 'p')
        T_ = [v[0] for v in tp]; Xs = [T_[(j-1) % 3]*T_[(j+1) % 3]/T_[j]**2 for j in range(3)]
        if which == 0:
            P = [dphi(tp, i)[0] for i in range(2)]; Qq = [dphi(tp, i)[1] for i in range(2)]
            return Tf(1j*beta, 1, *P, *Qq, *Xs) + Thf(1j*beta, 1, *P, *Qq, *Xs)
        tx = taus(x, t0, 'x'); tt = taus(x, t0, 't')
        px = [dphi(tx, k)[0] for k in range(3)]; pt = [dphi(tt, k)[0] for k in range(3)]
        return sum(a*a for a in px)/2 + sum(a*a for a in pt)/2 - (1/beta**2)*sum(Xs[j] - 1 for j in range(3))
    return [mp.quad(lambda x: dens(x, k), mp.linspace(-25, 25, 21)) for k in range(2)]


M = 2*h*m1/beta**2
th = mp.mpf('0.3')
q2s, Es = charges(make([1], [th], [mp.mpc(0.2, 0.37)]))
print(f'single species-1 soliton, theta={float(th)}: Q2 = {mp.nstr(q2s, 10)} (static -5.7735026919 e^(2 theta) = {mp.nstr(-5.7735026919*mp.e**(2*th), 10)}); E = {mp.nstr(Es, 10)} (M cosh = {mp.nstr(M*mp.cosh(th), 10)})')
q2s2, _ = charges(make([2], [th], [mp.mpc(0.2, 0.37)]))
print(f'single species-2 soliton, theta={float(th)}: Q2 = {mp.nstr(q2s2, 10)}')
u = mp.mpf('0.5'); th0 = mp.mpf('0.2')
for xis in ([mp.mpc(0.3, 0.4), mp.mpc(-0.1, 1.9)], [mp.mpc(0.3, 0.4), mp.mpc(0.3, -0.4) + 2j*mp.pi/3]):
    tb = make([1, 2], [th0 + 1j*u, th0 - 1j*u], xis)
    mn = min(abs(v[0]) for x in mp.linspace(-25, 25, 400) for v in tb(x, 0, 'x'))
    q2b, Eb = charges(tb)
    print(f'breather (1, th0+iu),(2, th0-iu), th0={float(th0)}, u={float(u)}, xi={[mp.nstr(z_, 3) for z_ in xis]}: min|tau_j| on t=0: {mp.nstr(mn, 3)}')
    print(f'   E  = {mp.nstr(Eb, 10)};  2 M cosh(th0) cos(u) = {mp.nstr(2*M*mp.cosh(th0)*mp.cos(u), 10)}')
    pred = q2s/mp.e**(2*th)*mp.e**(2*(th0+1j*u)) + q2s2/mp.e**(2*th)*mp.e**(2*(th0-1j*u))
    print(f'   Q2 = {mp.nstr(q2b, 10)};  q(1) e^(2(th0+iu)) + q(2) e^(2(th0-iu)) = {mp.nstr(pred, 10)}')
