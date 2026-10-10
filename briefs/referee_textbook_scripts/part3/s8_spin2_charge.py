"""a_2^(1): construct the spin-2 conserved currents d_-T + d_+Theta = 0 (real coupling, d_pm = d_t pm d_x), continue
beta -> i beta, evaluate Q_2 = int (T + Theta) dx on static single solitons of species 1 and 2 (imaginary coupling).
Tests 11-reality.qmd:43 ('real rapidities ... => all conserved charges real')."""
import sympy as sp, mpmath as mp
b, m = sp.symbols('beta m', positive=True)
p = sp.symbols('p1 p2'); q = sp.symbols('q1 q2'); X = sp.symbols('X0 X1 X2')
al = [(-1, 0, 1), (1, -1, 0), (0, 1, -1)]          # alpha_j = e_j - e_{j+1} (indices mod 3), j = 0,1,2
p3 = [p[0], p[1], -p[0]-p[1]]
adot = lambda j, v: sum(al[j][k]*v[k] for k in range(3))
F = [-(m**2/b)*sum(al[j][i]*X[j] for j in range(3)) for i in range(2)]                  # d+d- u_i (real coupling)
dF = [-(m**2/b)*sum(al[j][i]*X[j]*b*adot(j, p3) for j in range(3)) for i in range(2)]   # d+ F_i
mons3 = [p[0]**3, p[0]**2*p[1], p[0]*p[1]**2, p[1]**3]
c = sp.symbols('c0:4'); d = sp.symbols('d0:4'); f = sp.symbols('f0:6')
T = sum(ci*mi for ci, mi in zip(c, mons3)) + b*sum(d[2*i+k]*p[i]*q[k] for i in range(2) for k in range(2))
Th = m**2*sum(f[2*j+i]*X[j]*p[i] for j in range(3) for i in range(2))
dT = sum(sp.diff(T, p[i])*F[i] + sp.diff(T, q[i])*dF[i] for i in range(2))
dTh = sum(sp.diff(Th, X[j])*b*adot(j, p3)*X[j] for j in range(3)) + sum(sp.diff(Th, p[i])*q[i] for i in range(2))
eqs = sp.Poly(sp.expand(dT + dTh), *p, *q, *X).coeffs()
unk = list(c) + list(d) + list(f)
sol = sp.solve(eqs, unk, dict=True)[0]
free = [s for s in unk if s not in sol]
mp.mp.dps = 25
beta = mp.mpf('0.9'); h = 3; w = mp.e**(2j*mp.pi/3)
def single(a, xi=mp.mpc(0.2, 0.37)):
    ma = 2*mp.sin(mp.pi*a/h)
    def f_(x):
        E = mp.e**(ma*x + xi)
        return ([1 + w**(j*a)*E for j in range(3)], [ma*w**(j*a)*E for j in range(3)], [ma**2*w**(j*a)*E for j in range(3)])
    return f_
def Q(taus_static, Tf, Thf):
    def dens(x):
        T_, d1, d2 = taus_static(x)
        L1 = [d1[j]/T_[j] for j in range(3)]; L2 = [d2[j]/T_[j] - (d1[j]/T_[j])**2 for j in range(3)]
        P = [1j/beta*sum(al[j][i]*L1[j] for j in range(3)) for i in range(2)]          # phi = (i/beta) sum alpha_j ln tau_j
        Qq = [1j/beta*sum(al[j][i]*L2[j] for j in range(3)) for i in range(2)]
        Xs = [T_[(j-1) % 3]*T_[(j+1) % 3]/T_[j]**2 for j in range(3)]                 # e^{i beta alpha_j . phi}
        return Tf(1j*beta, 1, *P, *Qq, *Xs) + Thf(1j*beta, 1, *P, *Qq, *Xs)
    return mp.quad(dens, mp.linspace(-30, 30, 13))
for k in range(len(free)):
    sub = {fr: (1 if i == k else 0) for i, fr in enumerate(free)}
    Ts = sp.expand(T.subs(sol).subs(sub)); Ths = sp.expand(Th.subs(sol).subs(sub))
    assert sp.expand((dT + dTh).subs(sol).subs(sub)) == 0
    Tf = sp.lambdify((b, m, *p, *q, *X), Ts, 'mpmath'); Thf = sp.lambdify((b, m, *p, *q, *X), Ths, 'mpmath')
    q1 = Q(single(1), Tf, Thf); q2 = Q(single(2), Tf, Thf)
    print(f'basis {k}: T = {Ts};  Theta = {Ths}')
    print(f'     Q_2(species-1 soliton) = {mp.nstr(q1, 12)},  Q_2(species-2 soliton) = {mp.nstr(q2, 12)}')
