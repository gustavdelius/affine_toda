"""Integration A, chapter 10 items 9, 10, 11, 13, 14: own checks from the chapter's a_n^(1) formulas."""
import itertools, mpmath as mp
mp.mp.dps = 30

def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)

def mass(a, h, m=1): return 2*m*mp.sin(mp.pi*a/h)
def Acoef(a, b, h, th): return (mp.cosh(th) - mp.cos(mp.pi*(a-b)/h))/(mp.cosh(th) - mp.cos(mp.pi*(a+b)/h))
def Eexp(a, h, th, xi, x, t): return mp.e**(mass(a, h)*(x*mp.cosh(th) - t*mp.sinh(th)) + xi)
def taus(h, sol, x, t):          # sol: list of (a, theta, xi); eq-n-soliton
    om = mp.e**(2j*mp.pi/h); out = []
    for j in range(h):
        tot = 0
        for k in range(len(sol)+1):
            for S in itertools.combinations(range(len(sol)), k):
                term = 1
                for i in S: term *= om**(j*sol[i][0])*Eexp(sol[i][0], h, sol[i][1], sol[i][2], x, t)
                for i, l in itertools.combinations(S, 2): term *= Acoef(sol[i][0], sol[l][0], h, sol[i][1]-sol[l][1])
                tot += term
        out.append(tot)
    return out

# ---- item 9: c_2 from a_3^(1): coincident static pair of species 1 and 3 (parent labels)
h = 4; tc = [lambda E: 1+2*E+E**2/2, lambda E: 1+E**2/2, lambda E: 1-2*E+E**2/2, lambda E: 1+E**2/2]
ok = True
for x in (mp.mpf('-0.3'), mp.mpf('0.4')):
    T = taus(h, [(1, 0, 0), (3, 0, 0)], x, 0); E = mp.e**(mass(1, h)*x)
    ok &= max(abs(T[j] - tc[j](E)) for j in range(4)) < 1e-25
report('c_2 from a_3^(1): tau_0=1+2E+E^2/2, tau_1=tau_3=1+E^2/2, tau_2=1-2E+E^2/2 (parent labels, degree 2 = |O|)', ok)
# a_{2n}^(2) from a_{2n}^(1): the coincident pair (a, h-a) has degree 2 at every node, while n_j^vee = (2,..,2,1)
h = 5
for a in (1, 2):
    T = lambda E: [1 + (mp.e**(2j*mp.pi*j*a/h) + mp.e**(-2j*mp.pi*j*a/h))*E + Acoef(a, h-a, h, 0)*E**2 for j in range(h)]
    deg2 = all(abs(c) > 1e-20 for c in [Acoef(a, h-a, h, 0)])
report('a_4^(2) from a_4^(1): coincident pair (a, 5-a) has A(0) != 0, so degree 2 at every node, not |O| n_j^vee = 4', deg2)
kn = [2*1*2/2, 2*2*1/2, 2*2*mp.mpf(1)/2/2]   # n_j^vee = k n_j alpha_j^2/2, k=2, n=(1,2,2), alpha^2=(2,1,1/2)
report('a_4^(2): n_j^vee = k n_j alpha_j^2/2 = (2, 2, 1)', kn == [2, 2, 1])

# ---- item 10: position. |E| = 1 on x = t tanh(th) - Re(xi)/(m_a cosh th)
h, a, th, xi = 5, 2, mp.mpf('0.7'), mp.mpc('0.9', '0.4')
x0 = mp.findroot(lambda x: abs(Eexp(a, h, th, xi, x, 0)) - 1, 0)
report('position at t=0 is -Re xi/(m_a cosh theta), not -Re xi/m_a', abs(x0 + mp.re(xi)/(mass(a, h)*mp.cosh(th))) < 1e-25
       and abs(x0 + mp.re(xi)/mass(a, h)) > 0.05)

# ---- item 11: signed displacements and frame dependence of the passage-time shift
def shifts(a1, a2, h, th1, th2):
    lA = mp.log(Acoef(a1, a2, h, th1-th2))
    dx1 = -lA/(mass(a1, h)*mp.cosh(th1)); dx2 = lA/(mass(a2, h)*mp.cosh(th2))
    return dx1, dx2, -dx1/mp.tanh(th1), (-dx2/mp.tanh(th2) if th2 != 0 else None)
# numerical check of dx1 from the two-soliton tau: along soliton 1, tau_j / (single-soliton factor) at t -> -inf, +inf
h, a1, a2, th1, th2 = 6, 1, 2, mp.mpf('0.6'), mp.mpf('-0.6')
xi1, xi2 = mp.mpc('0.1', '0.3'), mp.mpc('-0.2', '0.8')
def ratio(t, shift):
    x = mp.tanh(th1)*t                                 # follow soliton 1's free trajectory through x=0 at t=0
    T = taus(h, [(a1, th1, xi1), (a2, th2, xi2)], x, t)
    single = taus(h, [(a1, th1, xi1 + shift)], x, t)
    E2 = Eexp(a2, h, th2, xi2, x, t); om = mp.e**(2j*mp.pi/h)
    return [T[j]/single[j]/(1 if shift == 0 else om**(j*a2)*E2) for j in range(h)]
lA = mp.log(Acoef(a1, a2, h, th1-th2))
r_in = ratio(mp.mpf(-40), 0); r_out = ratio(mp.mpf(40), lA)
asym = max(abs(r - 1) for r in r_in + r_out) < 1e-10
dx1, dx2, dt1, dt2 = shifts(a1, a2, h, th1, th2)
report(f'two solitons (a5, theta = +-0.6): soliton 1 emerges with xi_1 + ln A; dx1 = {mp.nstr(dx1, 6)} > 0, dx2 = {mp.nstr(dx2, 6)} < 0, both time shifts < 0',
       asym and dx1 > 0 and dx2 < 0 and dt1 < 0 and dt2 < 0)
d = shifts(1, 1, 3, mp.mpf('-0.4'), mp.mpf('-1'))
report(f'frame dependence: a2, theta = (-0.4, -1): dx1 = {mp.nstr(d[0], 8)} > 0 but dt1 = {mp.nstr(d[2], 8)} > 0', d[0] > 0 and d[2] > 0)
d = shifts(1, 2, 6, mp.mpf('0.9'), mp.mpf('0'))
report('rest frame of soliton 2: dt1 = -dx1/v1 < 0 (advance)', d[2] < 0)

# ---- item 13: fusion
ok_m = ok_p = True
for h in range(3, 9):
    for a in range(1, h):
        for b in range(1, h):
            A_, B_ = mp.pi*a/h, mp.pi*b/h
            if a + b < h:
                ok_m &= abs(mass(a, h)*mp.e**(1j*B_) + mass(b, h)*mp.e**(-1j*A_) - mass(a+b, h)) < 1e-25
                eps = mp.mpf('1e-12'); ok_p &= abs(eps*Acoef(a, b, h, 1j*mp.pi*(a+b)/h + eps)) > 1e-6    # simple pole
            elif a + b > h:
                c = a + b - h
                ok_m &= abs(mass(a, h)*mp.e**(1j*mp.pi*(h-b)/h) + mass(b, h)*mp.e**(-1j*mp.pi*(h-a)/h) - mass(c, h)) < 1e-25
                eps = mp.mpf('1e-12'); ok_p &= abs(eps*Acoef(a, b, h, 1j*mp.pi*(2*h-a-b)/h + eps)) > 1e-6
            else:
                ok_m &= abs(mass(a, h)*mp.e**(1j*B_) + mass(b, h)*mp.e**(-1j*A_)) < 1e-25            # vacuum
                eps = mp.mpf('1e-6'); ok_p &= abs(eps*Acoef(a, b, h, 1j*mp.pi + eps)) > 1e3                 # double pole
report('fusion: a+b<h -> species a+b at theta_12 = i pi(a+b)/h; a+b>h -> species a+b-h at i pi(2h-a-b)/h (simple poles); '
       'a+b=h -> zero exponent and double pole (h = 3..8)', ok_m and ok_p)
# the joint limit: theta_12 -> pole and xi_1, xi_2 -> -inf with A e^{xi_1+xi_2} fixed gives the single soliton a+b
def fused(h, a, b, th, xi, eps, x, t, upper=False):
    if not upper: t1, t2, c = th + 1j*mp.pi*b/h + eps/2, th - 1j*mp.pi*a/h - eps/2, a+b
    else: t1, t2, c = th + 1j*mp.pi*(h-b)/h + eps/2, th - 1j*mp.pi*(h-a)/h - eps/2, a+b-h
    xs = (xi - mp.log(Acoef(a, b, h, t1-t2)))/2
    T = taus(h, [(a, t1, xs), (b, t2, xs)], x, t); T0 = taus(h, [(c, th, xi)], x, t)
    return max(abs(T[j] - T0[j]) for j in range(h))
e1 = fused(3, 1, 1, mp.mpf('0.3'), mp.mpc('0.2', '0.5'), mp.mpf('1e-8'), mp.mpf('0.4'), mp.mpf('0.1'))
e2 = fused(3, 2, 2, mp.mpf('0.3'), mp.mpc('0.2', '0.5'), mp.mpf('1e-8'), mp.mpf('0.4'), mp.mpf('0.1'), upper=True)
e3 = fused(5, 1, 2, mp.mpf('-0.2'), mp.mpc('0.1', '1.1'), mp.mpf('1e-8'), mp.mpf('-0.3'), mp.mpf('0.2'))
e4 = fused(3, 1, 1, mp.mpf('0.3'), mp.mpc('0.2', '0.5'), mp.mpf('1e-12'), mp.mpf('0.4'), mp.mpf('0.1'))
report(f'joint fusion limit: a2 1+1 -> 2, a2 2+2 -> 1 (upper channel), a4 1+2 -> 3  [errors {mp.nstr(e1, 2)}, {mp.nstr(e2, 2)}, {mp.nstr(e3, 2)} at eps=1e-8; {mp.nstr(e4, 2)} at 1e-12]',
       max(e1, e2, e3) < 1e-3 and e4 < 1e-2*e1)

# ---- item 14: breathers
ok = True
for h in range(3, 9):
    for a in range(1, h):
        for u in (mp.mpf('0.3'), mp.mpf('0.9'), mp.mpf('1.4')):
            ok &= abs(Acoef(a, h-a, h, 2j*u) - (1 - mp.sin(mp.pi*a/h)**2/mp.cos(u)**2)) < 1e-25
            ok &= abs(Acoef(a, a, h, 2j*u) + mp.sin(u)**2/(mp.sin(mp.pi*a/h)**2 - mp.sin(u)**2)) < 1e-20
report('A_{a,h-a}(2iu) = 1 - sin^2(pi a/h)/cos^2 u and A_aa(2iu) = -sin^2 u/(sin^2(pi a/h) - sin^2 u)', ok)
# regular type-B family for A < 0: xi_2 = conj(xi_1) + i pi, min Re tau_j > 0
def breather_taus(h, a, u, th0, xi1, xi2, x, t):
    return taus(h, [(a, th0 + 1j*u, xi1), (h-a, th0 - 1j*u, xi2)], x, t)
h, a = 5, 1; ulo = abs(mp.pi/2 - mp.pi*a/h)
u = (ulo + mp.pi/2)/2; xi1 = mp.mpc('0.3', '0.5')
minre = min(mp.re(tj) for x in mp.linspace(-6, 6, 25) for t in mp.linspace(0, 8, 25)
            for tj in breather_taus(h, a, u, 0, xi1, mp.conj(xi1) + 1j*mp.pi, x, t))
report(f'type B, a4 species 1, u inside (|pi/2 - pi a/h|, pi/2): A < 0 and xi_2 = conj(xi_1)+i pi gives Re tau_j > 0 (min {mp.nstr(minre, 4)})',
       mp.re(Acoef(a, h-a, h, 2j*u)) < 0 and minre > 0)
# for 0 < A < 1: a zero of tau_j for every choice of xi_1, xi_2 (construction r = 1/sqrt(A))
u = ulo/2; Av = mp.re(Acoef(a, h-a, h, 2j*u)); ok0 = 0 < Av < 1
for xi1, xi2 in [(mp.mpc('0.3', '0.5'), mp.mpc('-0.4', '2.0')), (mp.mpc(0, '0.5'), mp.mpc(0, '-0.5') + 1j*mp.pi), (mp.mpc('1.0', '0'), mp.mpc('0', '1.3'))]:
    j = 1
    f = lambda x, t: breather_taus(h, a, u, 0, xi1, xi2, x, t)[j]
    best = min(((abs(f(x, t)), x, t) for x in mp.linspace(-4, 4, 41) for t in mp.linspace(0, 2*mp.pi/(mass(a, h)*mp.sin(u)), 41)), key=lambda z: z[0])
    sol = mp.findroot(lambda x, t: [mp.re(f(x, t)), mp.im(f(x, t))], (best[1], best[2]))
    ok0 &= abs(f(sol[0], sol[1])) < 1e-20
report(f'type B with 0 < A = {mp.nstr(Av, 5)} < 1: tau_1 has a zero for each of three choices of (xi_1, xi_2)', ok0)
h = 3; a = 1
report('a_2 band: regular neutral breathers have 0 < 2 M cos u < 2 M sin(pi/3) = sqrt(3) M',
       abs(2*mp.cos(abs(mp.pi/2 - mp.pi*a/h)) - mp.sqrt(3)) < 1e-25)
