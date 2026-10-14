"""A1 side claim (ch. 30, after prp-no-dirichlet; ch. 35 open problem 7): for a_2^(1) no combination Q_2 + c Q_{-2}
of the spin-+-2 charges survives the full Dirichlet condition phi(0,t) = phi_0.
Real coupling, m = beta = 1, V = sum_j (e^{alpha_j.phi} - 1), d_+ d_- phi = -grad V, d_+- = d_t +- d_x.
1. Solve D_- T = D_+ Theta for T a general spin-3 polynomial in u_k = d_+^{k+1} phi, Theta = sum_j e^{alpha_j.phi} f_j.u_0.
   Parity gives Tbar(v), Thetabar(v), v_k = d_-^{k+1} phi.  d/dt Q = (T + Theta)|_0, d/dt Qbar = -(Tbar + Thetabar)|_0.
2. At x = 0 with phi_t = 0: u_0 = p, u_1 = 2 p' + gradV, u_2 = 4 p'' + HessV p; v_0 = -p, v_1 = -2 p' + gradV,
   v_2 = -4 p'' - HessV p, with p = phi_x(0, t) free (' = d/dt).  The flux F_c = T + Theta - c (Tbar + Thetabar) must be
   a total time derivative: its Euler-Lagrange expression in p must vanish identically.  Find all c.
"""
import itertools
import sympy as sp

r2, r3 = sp.sqrt(2), sp.sqrt(3)
al = [sp.Matrix([r2, 0]), sp.Matrix([-1/r2, r3/r2])]
al = [-(al[0] + al[1])] + al                                 # alpha_0, alpha_1, alpha_2, |alpha|^2 = 2
E = sp.symbols('E0:3')                                        # e^{alpha_j . phi}, independent symbols
u = [sp.Matrix(sp.symbols(f'u{k}_1:3')) for k in range(3)]

# ansatz for T (spin 3): cubic in u0, u0*u1, u2
monos = []
for i, j, k in itertools.combinations_with_replacement(range(2), 3):
    monos.append(u[0][i]*u[0][j]*u[0][k])
for i in range(2):
    for j in range(2):
        monos.append(u[0][i]*u[1][j])
for i in range(2):
    monos.append(u[2][i])
cT = sp.symbols(f'cT0:{len(monos)}')
T = sum(c*m for c, m in zip(cT, monos))
fs = [sp.Matrix(sp.symbols(f'f{j}_1:3')) for j in range(3)]
Theta = sum(E[j]*(fs[j].T*u[0])[0] for j in range(3))

# derivations on the jet variables
gradV = sum((al[j]*E[j] for j in range(3)), sp.zeros(2, 1))
Dm_u = [-gradV,
        -sum((al[j]*E[j]*(al[j].T*u[0])[0] for j in range(3)), sp.zeros(2, 1)),
        -sum((al[j]*E[j]*((al[j].T*u[0])[0]**2 + (al[j].T*u[1])[0]) for j in range(3)), sp.zeros(2, 1))]


def Dminus(expr):
    return sp.expand(sum(sp.diff(expr, u[k][i])*Dm_u[k][i] for k in range(3) for i in range(2)))


def Dplus(expr):
    """d_+ on functions of E_j (phi) and u0, u1 (enough for Theta)."""
    out = sum(sp.diff(expr, E[j])*E[j]*(al[j].T*u[0])[0] for j in range(3))
    out += sum(sp.diff(expr, u[k][i])*u[k + 1][i] for k in range(2) for i in range(2))
    return sp.expand(out)


eq = sp.expand(Dminus(T) - Dplus(Theta))
unknowns = list(cT) + [x for f in fs for x in f]
coeffs = sp.Poly(eq, *(list(E) + [x for k in range(3) for x in u[k]])).coeffs()
sol = sp.solve(coeffs, unknowns, dict=True)[0]
Tsol, Thsol = sp.expand(T.subs(sol)), sp.expand(Theta.subs(sol))
free = sorted(set().union(Tsol.free_symbols, Thsol.free_symbols) & set(unknowns), key=str)
print('free parameters of the spin-3 current:', free)

# trivial currents T = D_+ X, Theta = D_- X (X spin 2) do not change the flux test; keep the general solution.
# boundary substitution
t = sp.symbols('t')
p = sp.Matrix([sp.Function('p1')(t), sp.Function('p2')(t)])
e0 = [sp.Rational(3, 2), sp.Rational(2, 5)]
Ev = {E[1]: e0[0], E[2]: e0[1], E[0]: 1/(e0[0]*e0[1])}       # e^{alpha_j.phi_0}, product 1
gV = gradV.subs(Ev)
HV = sum((al[j]*al[j].T*Ev[E[j]] for j in range(3)), sp.zeros(2, 2))
U = [p, 2*p.diff(t) + gV, 4*p.diff(t, 2) + HV*p]
Vb = [-p, -2*p.diff(t) + gV, -4*p.diff(t, 2) - HV*p]
sub_u = {u[k][i]: U[k][i] for k in range(3) for i in range(2)}
sub_v = {u[k][i]: Vb[k][i] for k in range(3) for i in range(2)}
flux_plus = sp.expand((Tsol + Thsol).subs(Ev).subs(sub_u, simultaneous=True))
flux_minus = sp.expand((Tsol + Thsol).subs(Ev).subs(sub_v, simultaneous=True))     # parity image
c = sp.symbols('c')
F = flux_plus - c*flux_minus
EL = [sp.expand(sp.euler_equations(F, [p[0], p[1]], t)[i].lhs) for i in range(2)]
P = sp.symbols('P1:3'); Pd = sp.symbols('Pd1:3'); Pdd = sp.symbols('Pdd1:3'); P3 = sp.symbols('P3_1:3'); P4 = sp.symbols('P4_1:3')
rep = {}
for i in range(2):
    rep[p[i].diff(t, 4)] = P4[i]; rep[p[i].diff(t, 3)] = P3[i]; rep[p[i].diff(t, 2)] = Pdd[i]; rep[p[i].diff(t)] = Pd[i]; rep[p[i]] = P[i]
ELp = [sp.expand(e.subs(rep)) for e in EL]
conds = []
for e in ELp:
    conds += sp.Poly(e, *(list(P) + list(Pd) + list(Pdd) + list(P3) + list(P4))).coeffs()
conds = [sp.simplify(x) for x in conds]
# conditions are polynomials in c and the free parameters of the current
sols = sp.solve(conds, [c] + free, dict=True)
print('solutions (c, current parameters) making the flux a total time derivative:', sols)
# report the flux for the nontrivial current at c = 1 restricted to p only (energy-like): the cubic part
Fp = sp.expand(F.subs(c, 1).subs(rep))
cubic = sp.expand(Fp.subs({Pd[0]: 0, Pd[1]: 0, Pdd[0]: 0, Pdd[1]: 0}))
print('flux at c = 1 with dp/dt = 0:', sp.factor(cubic))
print('cubic coefficients of T in each solution:', [[sp.simplify(Tsol.subs(s).coeff(u[0][0], 3)), sp.simplify(Tsol.subs(s).coeff(u[0][1], 3))] for s in sols])
