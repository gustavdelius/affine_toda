"""(a) real soliton in a_3^(1) (counterexample to 10-solitons.qmd:10);
(b) spins of conserved charges for twisted algebras via twisted Coxeter elements;
(c) sine-Gordon: classical time delay vs leading semiclassical phase of the book's S_0 (eq-sg-s0);
(d) breather with xi2 = conj(xi1): tau_j real and vanishing on curves (exr-breather-real);
(e) isolated zero of tau_j: energy density has double pole with zero residue (exr-singular-points)."""
import itertools, numpy as np, mpmath as mp
def report(n, ok): print(('PASS ' if ok else 'FAIL ') + n)
mp.mp.dps = 30
# (a) species-2 soliton of a_3^(1), Im xi = pi/2
h = 4; beta = mp.mpf('0.8')
al = [[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [-1, 0, 0, 1]]
def phi_a3(x, t, th=mp.mpf('0.3'), xi=mp.mpc(0.4, mp.pi/2), a=2):
    E = mp.e**(2*mp.sin(mp.pi*a/h)*(x*mp.cosh(th) - t*mp.sinh(th)) + xi)
    T = [1 + mp.e**(2j*mp.pi*j*a/h)*E for j in range(h)]
    return [1j/beta*sum(al[j][k]*mp.log(T[j]) for j in range(h)) for k in range(h)]
im = max(abs(mp.im(c)) for x in np.linspace(-5, 5, 41) for c in phi_a3(mp.mpf(x), mp.mpf('0.2')))
far = phi_a3(mp.mpf(30), 0)
def resid(x0, t0):
    w = 0
    for k in range(h):
        ptt = mp.diff(lambda t: phi_a3(x0, t)[k], t0, 2); pxx = mp.diff(lambda x: phi_a3(x, t0)[k], x0, 2)
        P = phi_a3(x0, t0)
        F = 1j/beta*sum(al[j][k]*mp.e**(1j*beta*sum(al[j][l]*P[l] for l in range(h))) for j in range(h))
        w = max(w, abs(ptt - pxx - F))
    return w
print('   phi(+inf) =', [mp.nstr(c, 8) for c in far], ' (2pi/beta) =', mp.nstr(2*mp.pi/beta, 8))
report(f'a_3^(1) species-2 soliton, Im xi = pi/2: solves field eq (res {float(resid(mp.mpf(0.3), mp.mpf(0.1))):.0e}) and is REAL (max |Im phi| = {float(im):.0e}), interpolating 0 -> nonzero vacuum',
       im < 1e-25 and resid(mp.mpf(0.3), mp.mpf(0.1)) < 1e-15)

# (b) twisted Coxeter elements
def cartan(typ, n):
    C = 2*np.eye(n)
    if typ == 'a':
        for i in range(n-1): C[i, i+1] = C[i+1, i] = -1
    elif typ == 'd':
        for i in range(n-2): C[i, i+1] = C[i+1, i] = -1
        C[n-3, n-1] = C[n-1, n-3] = -1
    elif typ == 'e':                         # Bourbaki: 1-3-4-5-6, 2-4
        for (i, j) in [(1, 3), (3, 4), (4, 5), (5, 6), (2, 4)]: C[i-1, j-1] = C[j-1, i-1] = -1
    return C
def twisted_exponents(typ, n, sigma):
    C = cartan(typ, n); al = np.linalg.cholesky(C)
    S = np.linalg.solve(al, al[sigma]).T                   # S alpha_i = alpha_sigma(i)
    assert np.allclose(S @ al.T, al[sigma].T)
    reps, seen = [], set()
    for i in range(n):
        if i not in seen:
            orb = {i}; k = sigma[i]
            while k != i: orb.add(k); k = sigma[k]
            seen |= orb; reps.append(i)
    w = np.eye(n)
    for i in reps: w = w @ (np.eye(n) - np.outer(al[i], al[i]))
    w = w @ S
    order = next(p for p in range(1, 200) if np.allclose(np.linalg.matrix_power(w, p), np.eye(n)))
    ev = np.linalg.eigvals(w)
    ex = sorted(int(round((np.angle(e)/(2*np.pi)*order) % order)) for e in ev)
    return order, ex
cases = {'a_2^(2)': ('a', 2, [1, 0]), 'a_4^(2)': ('a', 4, [3, 2, 1, 0]), 'a_3^(2)': ('a', 3, [2, 1, 0]),
         'a_5^(2)': ('a', 5, [4, 3, 2, 1, 0]), 'd_4^(2)': ('d', 4, [0, 1, 3, 2]), 'd_5^(2)': ('d', 5, [0, 1, 2, 4, 3]),
         'e_6^(2)': ('e', 6, [5, 1, 4, 3, 2, 0]), 'd_4^(3)': ('d', 4, [2, 1, 3, 0])}
hbook = {'a_2^(2)': 3, 'a_4^(2)': 5, 'a_3^(2)': 3, 'a_5^(2)': 5, 'd_4^(2)': 4, 'd_5^(2)': 5, 'e_6^(2)': 9, 'd_4^(3)': 4}
for name, (typ, n, sig) in cases.items():
    order, ex = twisted_exponents(typ, n, sig)
    print(f'   {name}: twisted Coxeter order {order} (= k * h_book = {order//hbook[name]} * {hbook[name]}), exponents mod {order}: {ex}')
untw = {'b_n,c_n (n=3)': [1, 3, 5], 'g_2': [1, 5], 'f_4': [1, 5, 7, 11]}
print('   untwisted: c_3,b_3 exponents 1,3,5 mod 6 (all odd); g_2: 1,5 mod 6; f_4: 1,5,7,11 mod 12')
o, ex = twisted_exponents('d', 4, [2, 1, 3, 0]); s_d43 = {s for s in range(1, 200) if s % o in ex}
s_g2 = {s for s in range(1, 200) if s % 6 in (1, 5)}
o, ex = twisted_exponents('e', 6, [5, 1, 4, 3, 2, 0]); s_e62 = {s for s in range(1, 200) if s % o in ex}
s_f4 = {s for s in range(1, 200) if s % 12 in (1, 5, 7, 11)}
o, ex = twisted_exponents('d', 5, [0, 1, 2, 4, 3]); s_d52 = {s for s in range(1, 200) if s % o in ex}   # dual of c_4
s_c4 = {s for s in range(1, 200) if s % 8 in (1, 3, 5, 7)}
o, ex = twisted_exponents('a', 5, [4, 3, 2, 1, 0]); s_a52 = {s for s in range(1, 200) if s % o in ex}   # dual of b_3
s_b3 = {s for s in range(1, 200) if s % 6 in (1, 3, 5)}
report('Lie-dual pairs have identical spin sets: (g2,d4^(3)), (f4,e6^(2)), (c4,d5^(2)), (b3,a5^(2))',
       s_g2 == s_d43 and s_f4 == s_e62 and s_c4 == s_d52 and s_b3 == s_a52)
o, ex = twisted_exponents('a', 2, [1, 0])
naive = {s for s in range(1, 40) if s % 3 in {e % 3 for e in ex}}
print('   a_2^(2): "exponents + h Z" with h = 3 would give', sorted(naive)[:8], '...; true spins', sorted(s for s in range(1, 40) if s % 6 in ex)[:8])

# (c) sine-Gordon time delay vs semiclassical S_0
def lnS0(th, xi):
    f = lambda t: mp.sinh(t*(mp.pi - xi)/2)*mp.sin(th*t)/(t*mp.sinh(xi*t/2)*mp.cosh(mp.pi*t/2))
    return -mp.quad(f, [0, 1, 5, 20, 80, mp.inf])           # delta(theta) up to the constant pi
for lam in (40, 160):
    xi = mp.pi/lam
    for th in (mp.mpf('0.5'), mp.mpf('1.5')):
        dd = mp.diff(lambda s: lnS0(s, xi), th)
        print(f'   lambda={lam}: theta={float(th)}: (xi/2) d delta/d theta = {mp.nstr(xi/2*dd, 8)}, ln tanh(theta/2) = {mp.nstr(mp.log(mp.tanh(th/2)), 8)}')
# (d) a_2 breather with xi2 = conj(xi1): tau_j real, negative somewhere
h = 3; u = mp.mpf('0.6'); a = 1
th1, th2 = mp.mpc(0.2, u), mp.mpc(0.2, -u); xi1 = mp.mpc(0.3, 0.4); xi2 = mp.conj(xi1)
ma = 2*mp.sin(mp.pi*a/h)
A = (mp.cosh(th1 - th2) - mp.cos(mp.pi*(2*a - h)/h))/(mp.cosh(th1 - th2) - mp.cos(mp.pi))
def tau(j, x, t):
    E1 = mp.e**(ma*(x*mp.cosh(th1) - t*mp.sinh(th1)) + xi1); E2 = mp.e**(ma*(x*mp.cosh(th2) - t*mp.sinh(th2)) + xi2)
    w = mp.e**(2j*mp.pi/h)
    return 1 + w**(j*a)*E1 + w**(j*(h-a))*E2 + A*E1*E2
vals = [tau(j, mp.mpf(x), mp.mpf(t)) for j in range(3) for x in np.linspace(-4, 4, 33) for t in np.linspace(-4, 4, 33)]
print(f'   a_2 breather, xi2=conj(xi1): A={mp.nstr(A, 6)}, max|Im tau_j|={float(max(abs(mp.im(v)) for v in vals)):.0e}, min tau_j={mp.nstr(min(mp.re(v) for v in vals), 5)}, max tau_j={mp.nstr(max(mp.re(v) for v in vals), 5)}')
report('exr-breather-real: with xi2=conj(xi1) the tau_j are real but change sign => vanish on curves (singular field)',
       max(abs(mp.im(v)) for v in vals) < 1e-20 and min(mp.re(v) for v in vals) < 0)
# (e) residue of energy density at a zero of tau_j: a_2 two-soliton, zero of tau_0 near real axis
beta = mp.mpf(1); spec = [1, 1]; ths = [mp.mpf('0.7'), mp.mpf('-0.4')]; xis = [mp.mpc(0.1, 0.5), mp.mpc(-0.2, 2.0)]
alv = [[1, -1, 0], [0, 1, -1], [-1, 0, 1]]
def taus(x, t):
    w = mp.e**(2j*mp.pi/3); Es = [mp.e**(2*mp.sin(mp.pi/3)*(x*mp.cosh(th) - t*mp.sinh(th)) + xi) for th, xi in zip(ths, xis)]
    c = mp.cosh(ths[0] - ths[1]); A = (c - 1)/(c - mp.cos(2*mp.pi/3))
    return [1 + w**j*(Es[0] + Es[1]) + w**(2*j)*A*Es[0]*Es[1] for j in range(3)]
def T00(x, t):     # energy density 1/2 phi_t^2 + 1/2 phi_x^2 + V with phi = (i/beta) sum alpha_j ln tau_j, as meromorphic fn of complex x
    T = taus(x, t)
    dx = [mp.diff(lambda y: taus(y, t)[j], x) for j in range(3)]; dt = [mp.diff(lambda s: taus(x, s)[j], t) for j in range(3)]
    px = [1j/beta*sum(alv[j][k]*dx[j]/T[j] for j in range(3)) for k in range(3)]
    pt = [1j/beta*sum(alv[j][k]*dt[j]/T[j] for j in range(3)) for k in range(3)]
    V = -(1/beta**2)*sum(T[(j-1) % 3]*T[(j+1) % 3]/T[j]**2 - 1 for j in range(3))
    return sum(p*p for p in px)/2 + sum(p*p for p in pt)/2 + V
t0 = mp.mpf('0.3')
z = mp.findroot(lambda x: taus(x, t0)[0], mp.mpc(0.5, 0.3))
r = mp.mpf('0.05')
res = mp.quad(lambda s: T00(z + r*mp.e**(1j*s), t0)*1j*r*mp.e**(1j*s), [0, mp.pi, 2*mp.pi])/(2j*mp.pi)
dbl = mp.quad(lambda s: T00(z + r*mp.e**(1j*s), t0)*(r*mp.e**(1j*s))*1j*r*mp.e**(1j*s), [0, mp.pi, 2*mp.pi])/(2j*mp.pi)
print(f'   zero of tau_0 at x*={mp.nstr(z, 8)} (t={float(t0)}): residue of T00 = {mp.nstr(res, 5)}, double-pole coefficient = {mp.nstr(dbl, 8)}')
report('exr-singular-points: T00 has a double pole with zero residue at a zero of tau_j', abs(res) < 1e-15 and abs(dbl) > 1e-3)
