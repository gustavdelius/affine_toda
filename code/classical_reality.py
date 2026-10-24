"""Chapter 11 (sec-classical): reality of classical charges, singular solutions, breather band.

Conventions of chapters 9-10, with m = beta = 1: d_pm = d_t pm d_x, Hirota ansatz
phi = (i/beta) sum_j (2 alpha_j/alpha_j^2) ln tau_j, a_n^(1) N-soliton tau_j of eq-n-soliton with
E_k = exp(m_a (x cosh th - t sinh th) + xi), m_a = 2 sin(pi a/h), M_a = 2 h m_a/beta^2.
psi = i beta phi = -sum_j alpha_j ln tau_j (alpha_j = e_j - e_{j+1} in R^h) solves the real-coefficient
equation d_+ d_- psi = -sum_j alpha_j e^{alpha_j.psi}, and Theta: phi -> -conj(phi) is psi -> conj(psi).
Energy and momentum are the canonical integrals int (1/2)(phi_t^2 + phi_x^2) + V_imag and -int phi_t.phi_x,
evaluated from the tau-functions by the trapezoidal rule on [-L, L] (no boundary-term formula is used).

1. prp-hirota-energy: the canonical energy density equals (2/beta^2) sum_j (2/alpha_j^2) d_x^2 ln tau_j;
   E, P additive for complex rapidities with Re cosh > 0; sine-Gordon at theta = i pi,
   xi = i pi/2 is a regular real antikink with E = +8, not 8 cosh(i pi) = -8; a_2 soliton with Re cosh < 0
   has E = -M cosh theta, alone or in a two-soliton, being the conjugate species at theta - i pi up to
   a constant shift of phi.
2. prp-theta-reality: breather E, P real; the sine-Gordon kink has Theta-invariant data and real
   topological charge of parity -1; a field fixed by Theta has E_imag = -E_real(-i phi) < 0.
3. prp-breather-charges: the a_2^(1) currents of spin 2 and 4, constructed here from the field equation,
   have q_2 = -q_1 on single solitons, q_a e^{s theta} at complex rapidity, and on breathers
   Q_s = q_1 (e^{s(th0+iu)} - e^{s(th0-iu)}) = 2 i q_1 e^{s th0} sin(s u), for several phase choices.
4. sec-singular-solutions: near a zero of tau on a slice the energy integral diverges
   (8(1 - coth 2 eps) for static sine-Gordon); slices that avoid the zeros give the additive value;
   a zero crossing the real line changes the winding of tau_j by 2 pi, so the charge jumps by a root.
5. sec-breathers / exr-breather-real: A_{a abar}(2iu) = 1 - sin^2(pi a/h)/cos^2 u for a_5^(1); A < 0 gives a
   regular breather, 0 < A < 1 a zero of tau; xi_2 = conj(xi_1) gives real tau_j that change sign.
"""
import itertools
import warnings
import numpy as np
import sympy as sp
from scipy.optimize import fsolve
from sympy.polys.matrices import DomainMatrix

TOL = 1e-9
warnings.simplefilter('ignore', RuntimeWarning)   # fsolve trial points far out overflow harmlessly


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ------------------------------------------------------------------ Hirota solutions of a_n^(1)
def A_int(h, a, b, th):
    c = np.cosh(th)
    return (c - np.cos(np.pi*(a-b)/h))/(c - np.cos(np.pi*(a+b)/h))


def tau_terms(h, sols):
    """sols = [(a, theta, xi), ...]; tau_j = sum over subsets S of c_{S,j} exp(k_S x + w_S t)."""
    om = np.exp(2j*np.pi/h); out = []
    for r in range(len(sols)+1):
        for S in itertools.combinations(range(len(sols)), r):
            k = sum(2*np.sin(np.pi*sols[i][0]/h)*np.cosh(sols[i][1]) for i in S)
            w = -sum(2*np.sin(np.pi*sols[i][0]/h)*np.sinh(sols[i][1]) for i in S)
            c = np.prod([np.exp(sols[i][2]) for i in S]) if S else 1.0
            for i1, i2 in itertools.combinations(S, 2):
                c = c*A_int(h, sols[i1][0], sols[i2][0], sols[i1][1]-sols[i2][1])
            out.append((np.array([c*om**(j*sum(sols[i][0] for i in S)) for j in range(h)]), k, w))
    return out


def taus(h, sols, x, t):
    return sum(np.outer(c, np.exp(k*x + w*t)) for c, k, w in tau_terms(h, sols))


def fields(h, sols, x, t, nmax=5):
    """D^n psi (D = d_-), d_+ psi, d_+ d_- psi and E_j = e^{alpha_j.psi}, from exact tau derivatives."""
    Dn = np.zeros((nmax+1, h, len(x)), complex); P = np.zeros((h, len(x)), complex); PM = np.zeros_like(P)
    for c, k, w in tau_terms(h, sols):
        ex = np.exp(k*x + w*t)
        for n in range(nmax+1):
            Dn[n] += np.outer(c, (w-k)**n*ex)
        P += np.outer(c, (w+k)*ex); PM += np.outer(c, (w+k)*(w-k)*ex)
    tau = Dn[0]; r = [None] + [Dn[n]/tau for n in range(1, nmax+1)]
    L = {1: r[1], 2: r[2] - r[1]**2, 3: r[3] - 3*r[1]*r[2] + 2*r[1]**3,
         4: r[4] - 4*r[1]*r[3] - 3*r[2]**2 + 12*r[1]**2*r[2] - 6*r[1]**4,
         5: r[5] - 5*r[1]*r[4] - 10*r[2]*r[3] + 20*r[1]**2*r[3] + 30*r[1]*r[2]**2 - 60*r[1]**3*r[2] + 24*r[1]**5}
    comp = lambda arr: np.array([arr[(c-1) % h] - arr[c] for c in range(h)])   # psi_c = ln tau_{c-1} - ln tau_c
    Ej = np.array([tau[(j-1) % h]*tau[(j+1) % h]/tau[j]**2 for j in range(h)])
    return dict(tau=tau, Ej=Ej, q={n: comp(L[n]) for n in L}, qp=comp(P/tau),
                qpm=comp(PM/tau - P*Dn[1]/tau**2))


def energy_momentum(h, sols, t=0.0, L=40.0, N=40001):
    x = np.linspace(-L, L, N); F = fields(h, sols, x, t)
    sq = lambda v: np.sum(v*v, axis=0)
    eps = -0.25*(sq(F['qp']) + sq(F['q'][1])) + np.sum(1 - F['Ej'], axis=0)   # (1/2)(phi_t^2+phi_x^2)+V
    mom = 0.25*(sq(F['qp']) - sq(F['q'][1]))                                     # -phi_t.phi_x
    eom = np.max(np.abs(F['qpm'] + np.array([F['Ej'][c] - F['Ej'][(c-1) % h] for c in range(h)])))
    return np.trapz(eps, x), np.trapz(mom, x), eom, np.min(np.abs(F['tau']))


def mass(h, a):
    return 2*h*2*np.sin(np.pi*a/h)


# ------------------------------------------------------------------ 1. prp-hirota-energy
h = 3; M = mass(3, 1)
sols = [(1, 0.23+1j*np.pi/3, 0.1+0.2j), (2, 0.23-1j*np.pi/3, 0.1+0.8j)]
x = np.linspace(-6, 6, 13); F = fields(h, sols, x, 0.4)
eps = -0.25*np.sum(F['qp']**2 + F['q'][1]**2, axis=0) + np.sum(1 - F['Ej'], axis=0)
T0 = sum(np.outer(c, np.exp(k*x + w*0.4)) for c, k, w in tau_terms(h, sols))
T1 = sum(np.outer(c, k*np.exp(k*x + w*0.4)) for c, k, w in tau_terms(h, sols))
T2 = sum(np.outer(c, k*k*np.exp(k*x + w*0.4)) for c, k, w in tau_terms(h, sols))
err = np.max(np.abs(eps - 2*np.sum(T2/T0 - (T1/T0)**2, axis=0)))
report(f'a_2 breather: canonical energy density = (2/beta^2) sum_j (2/alpha_j^2) d_x^2 ln tau_j pointwise [err {err:.1e}]',
       err < 1e-12)
cases = [[(1, 0.3+0.4j, 0.2+0.5j)], [(2, -0.6+0.25j, -0.4+1.3j)],
         [(1, 0.4+0.3j, 0.1+0.7j), (2, -0.5+0.2j, -0.3+0.4j)]]
worst = 0; eom_w = 0
for sols in cases:
    E, P, eom, mt = energy_momentum(h, sols, t=0.05)
    worst = max(worst, abs(E - sum(M*np.cosh(s[1]) for s in sols)), abs(P - sum(M*np.sinh(s[1]) for s in sols)))
    eom_w = max(eom_w, eom)
report(f'a_2: canonical E, P = sum M cosh, sum M sinh for complex rapidities with Re cosh > 0 '
       f'(one soliton of each species, a non-conjugate pair) [err {worst:.1e}, field eq. {eom_w:.1e}]',
       worst < TOL and eom_w < 1e-10)

E, P, eom, mt = energy_momentum(2, [(1, 1j*np.pi, 1j*np.pi/2)])
x = np.linspace(-5, 5, 11); T = taus(2, [(1, 1j*np.pi, 1j*np.pi/2)], x, 0.3)
phi = np.sqrt(2)*1j*(np.log(T[1]) - np.log(T[0]))           # phi.alpha_1/sqrt2, beta = 1
report(f'sine-Gordon theta = i pi, xi = i pi/2: regular (min|tau| = {mt:.2f}) real antikink '
       f'2 sqrt2 arctan e^(-2x) with E = {E.real:.10f}, while 8 cosh(i pi) = -8',
       abs(E - 8) < TOL and mt > 0.99 and np.allclose(phi, 2*np.sqrt(2)*np.arctan(np.exp(-2*x)), atol=1e-12))

th = 0.3 + 1j*(np.pi - 0.4)
E, P, eom, mt = energy_momentum(3, [(1, th, 0.2+0.5j)], t=0.1)
x = np.linspace(-3, 3, 13); ratio = 0
for t in (0.0, 0.7):
    TA = taus(3, [(1, th, 0.2+0.5j)], x, t); TB = taus(3, [(2, th - 1j*np.pi, -0.2-0.5j)], x, t)
    Ep = np.exp(2*np.sin(np.pi/3)*(x*np.cosh(th - 1j*np.pi) - t*np.sinh(th - 1j*np.pi)) - 0.2-0.5j)
    om = np.exp(2j*np.pi/3)
    ratio = max(ratio, max(np.max(np.abs(TA[j]*Ep/(om**j*TB[j]) - 1)) for j in range(3)))
th2 = -0.2 + 0.3j
E2, P2, _, mt2 = energy_momentum(3, [(1, th, 0.2+0.5j), (2, th2, -0.3+0.9j)], t=0.1)
err2 = max(abs(E2 + M*np.cosh(th) - M*np.cosh(th2)), abs(P2 + M*np.sinh(th) - M*np.sinh(th2)))
report(f'a_2, Re cosh theta < 0: E = -M cosh theta, P = -M sinh theta (also as one constituent of a two-soliton), and (a, theta, xi) is (abar, theta - i pi, -xi) '
       f'(tau_j agree up to the factors omega^j/E, i.e. a constant shift of phi) [err {max(abs(E + M*np.cosh(th)), abs(P + M*np.sinh(th))):.1e}]',
       abs(E + M*np.cosh(th)) < TOL and abs(P + M*np.sinh(th)) < TOL and ratio < 1e-12 and err2 < TOL and mt2 > 0.05)

idn = max(abs(A_int(hh, a, b, z + 1j*np.pi)*A_int(hh, hh - a, b, z) - 1)
          for hh in (3, 4, 6) for a in range(1, hh) for b in range(1, hh) for z in (0.3+0.2j, -1.1+0.7j, 2.0))
report(f'A_ab(theta + i pi) A_(h-a)b(theta) = 1 (relabelling inside multi-solitons) [err {idn:.1e}]', idn < 1e-12)

# ------------------------------------------------------------------ 2. prp-theta-reality
th0, u = 0.23, np.pi/3; worst = 0
for xis, t in [((1j*np.pi/2, 1j*np.pi/2), 0.17), ((0.2+1j*np.pi/2, -0.2+1j*np.pi/2), 0.43), ((0.1+0.2j, 0.1+0.8j), 1.1)]:
    E, P, eom, mt = energy_momentum(3, [(1, th0+1j*u, xis[0]), (2, th0-1j*u, xis[1])], t=t)
    worst = max(worst, abs(E - 2*M*np.cosh(th0)*np.cos(u)), abs(P - 2*M*np.sinh(th0)*np.cos(u)), mt < 0.05)
report(f'a_2 breather (1, th0+iu), (2, th0-iu): E = 2M cosh th0 cos u, P = 2M sinh th0 cos u, real, '
       f'for three phase choices [err {worst:.1e}]', worst < TOL)

x = np.linspace(-30, 30, 6001); T = taus(2, [(1, 0.0, 1j*np.pi/2)], x, 0.0)
lnT = np.log(np.abs(T)) + 1j*np.unwrap(np.angle(T), axis=1)
dphi = (1j*np.sqrt(2)*(lnT[1] - lnT[0]))                    # alpha_1 = -alpha_0 = sqrt2 in R^1
Qtop = (dphi[-1] - dphi[0])*np.sqrt(2)/(2*np.pi)            # (beta/2pi) Delta phi in units of alpha_1/2
report('sine-Gordon kink: data (1, 0) Theta-invariant, topological charge 1 (real), but Theta reverses it '
       '(parity -1), so "Q = eps_Q conj Q" fails for topological charges', abs(Qtop - 1) < 1e-10)

rng = np.random.default_rng(1)
x = np.linspace(-12, 12, 4801)
chi = np.array([c*np.exp(-(x - x0)**2) for c, x0 in rng.normal(size=(3, 2))]); chi -= chi.mean(axis=0)
chit = np.array([c*np.exp(-(x - x0)**2/2) for c, x0 in rng.normal(size=(3, 2))]); chit -= chit.mean(axis=0)
chix = np.gradient(chi, x, axis=1)
al = lambda v, j: v[j] - v[(j+1) % 3]
E_real = np.trapz(0.5*np.sum(chit**2 + chix**2, axis=0) + sum(np.exp(al(chi, j)) - 1 for j in range(3)), x)
ph, pht, phx = -1j*chi, -1j*chit, -1j*chix                  # phi = -i chi is fixed by Theta
E_imag = np.trapz(0.5*np.sum(pht**2 + phx**2, axis=0) + sum(1 - np.exp(1j*al(ph, j)) for j in range(3)), x)
report(f'a_2 field fixed by Theta, phi = -i chi: E_imag[phi] = -E_real[chi] = {E_imag.real:.6f} < 0',
       abs(E_imag + E_real) < 1e-12 and E_real > 0)

# ------------------------------------------------------------------ 3. a_2^(1) higher charges
K = 6
q = [[sp.Symbol(f'q{k}_{c}') for c in range(2)] for k in range(K + 2)]
Es = sp.symbols('E0 E1 E2'); GENS = [v for row in q for v in row] + list(Es)
cmp = lambda k, c: q[k][c] if c < 2 else -q[k][0] - q[k][1]


def Dm(e):              # total d_- on the jet; d_- E_j = E_j alpha_j.d_- psi
    out = sum(sp.diff(e, q[k][c])*q[k+1][c] for k in range(1, K+1) for c in range(2))
    return sp.expand(out + sum(sp.diff(e, Es[j])*Es[j]*(cmp(1, j % 3) - cmp(1, (j+1) % 3)) for j in range(3)))


DF = [[sp.expand(-(Es[c] - Es[(c-1) % 3]))] for c in range(2)]   # d_+ d_- psi_c
for c in range(2):
    for k in range(1, K):
        DF[c].append(Dm(DF[c][-1]))
Dp = lambda T: sp.expand(sum(sp.diff(T, q[k][c])*DF[c][k-1] for k in range(1, K+1) for c in range(2)))


def monomials(w):
    vs = [(q[k][c], k) for k in range(1, w+1) for c in range(2)]; res = set()
    def rec(i0, left, cur):
        if left == 0:
            res.add(sp.Mul(*cur)); return
        for i in range(i0, len(vs)):
            if vs[i][1] <= left:
                rec(i, left - vs[i][1], cur + [vs[i][0]])
    rec(0, w, [])
    return sorted(res, key=sp.default_sort_key)


def qq(v):
    return sp.QQ.from_sympy(sp.Rational(v)) if isinstance(v, sp.Basic) else v


def rank(vecs):
    return DomainMatrix([[qq(a) for a in v] for v in vecs], (len(vecs), len(vecs[0])), sp.QQ).rank() if vecs else 0


def current(s):
    """Solve d_+ T = d_- X, T of weight s+1 in the q, X = sum_j E_j (weight s-1). Q = int (T - X) dx."""
    mT = monomials(s+1); mX = [Es[j]*m for j in range(3) for m in (monomials(s-1) if s > 1 else [sp.Integer(1)])]
    cols = [dict(sp.Poly(Dp(m), *GENS).terms()) for m in mT] + [dict(sp.Poly(-Dm(m), *GENS).terms()) for m in mX]
    keys = sorted({k for c in cols for k in c}); idx = {k: i for i, k in enumerate(keys)}
    rows = [[sp.QQ(0)]*len(cols) for _ in keys]
    for j, c in enumerate(cols):
        for k, v in c.items():
            rows[idx[k]][j] = qq(v)
    ns = [list(r) for r in DomainMatrix(rows, (len(keys), len(cols)), sp.QQ).nullspace().to_Matrix().tolist()]
    mono = lambda m: sp.Poly(m, *GENS).monoms()[0]
    triv = []                                                    # T = D G, X = d_+ G
    for g in monomials(s):
        dT = dict(sp.Poly(Dm(g), *GENS).terms()); dX = dict(sp.Poly(Dp(g), *GENS).terms())
        triv.append([dT.get(mono(m), 0) for m in mT] + [dX.get(mono(m), 0) for m in mX])
    rt = rank(triv); new = rank(triv + ns) - rt
    pick = next((v for v in ns if rank(triv + [v]) > rt), None)
    if pick is None:
        return new, None
    T = sum(sp.QQ.to_sympy(qq(c))*m for c, m in zip(pick[:len(mT)], mT))
    X = sum(sp.QQ.to_sympy(qq(c))*m for c, m in zip(pick[len(mT):], mX))
    assert sp.expand(Dp(T) - Dm(X)) == 0
    return new, sp.lambdify([q[k][c] for k in range(1, 6) for c in range(2)] + list(Es), T - X, 'numpy')


cur = {s: current(s) for s in (2, 3, 4)}
report('a_2^(1): one primitive local current of spin 2 and of spin 4, none of spin 3 (from the field equation; '
       'real coefficients in psi, so Theta-even before normalization)',
       cur[2][0] == 1 and cur[4][0] == 1 and cur[3][0] == 0)


def charge(s, sols, t=0.0, L=40.0, N=40001):
    x = np.linspace(-L, L, N); F = fields(3, sols, x, t)
    args = [F['q'][k][c] for k in range(1, 6) for c in range(2)] + [F['Ej'][j] for j in range(3)]
    return np.trapz(cur[s][1](*args)*np.ones_like(x), x)


for s in (2, 4):
    q1 = charge(s, [(1, 0.0, 0.3j)]); q2 = charge(s, [(2, 0.0, 0.3j)])
    qc = charge(s, [(1, 0.3+0.4j, 0.2+0.5j)], t=0.05)
    report(f'spin {s}: q_1 imaginary in the real-coefficient normalization ({q1:.6f}), q_2/q_1 = -1 '
           f'= (-1)^(s-1), complex-rapidity soliton gives q_1 e^(s theta)',
           abs(q1.real) < TOL and abs(q2/q1 + 1) < TOL and abs(qc/(q1*np.exp(s*(0.3+0.4j))) - 1) < TOL)
    worst = 0
    for xis, t in [((1j*np.pi/2, 1j*np.pi/2), 0.17), ((0.2+1j*np.pi/2, -0.2+1j*np.pi/2), 0.43),
                   ((0.1+0.2j, 0.1+0.8j), 1.1)]:
        Qb = charge(s, [(1, th0+1j*u, xis[0]), (2, th0-1j*u, xis[1])], t=t)
        worst = max(worst, abs(Qb/q1 - 2j*np.exp(s*th0)*np.sin(s*u)))
    report(f'spin {s}: breather charge Q/q_1 = 2i e^(s th0) sin(s u), imaginary, for three phase choices '
           f'[err {worst:.1e}]', worst < TOL)

# ------------------------------------------------------------------ 4. singular solutions
x0 = np.array([0.05, 0.2, 0.7]); T = taus(2, [(1, 0.0, 0.0)], x0, 0.0)          # tau_0,1 = 1 +- e^{2x}
F = fields(2, [(1, 0.0, 0.0)], x0, 0.0)
dens = -0.25*(np.sum(F['qp']**2 + F['q'][1]**2, axis=0)) + np.sum(1 - F['Ej'], axis=0)
ok_d = np.allclose(dens, -8/np.sinh(2*x0)**2, rtol=1e-10)
ok_i = True
for e in (1e-2, 1e-3):
    xs = np.concatenate([-np.geomspace(20, e, 200001), np.geomspace(e, 20, 200001)])
    Fx = fields(2, [(1, 0.0, 0.0)], xs, 0.0)
    de = -0.25*(np.sum(Fx['qp']**2 + Fx['q'][1]**2, axis=0)) + np.sum(1 - Fx['Ej'], axis=0)
    val = np.trapz(de[:200001].real, xs[:200001]) + np.trapz(de[200001:].real, xs[200001:])
    ok_i &= abs(val - 8*(1 - 1/np.tanh(2*e))) < 1e-6*abs(val)
report('static sine-Gordon tau = 1 +- e^(2x): energy density -8 csch^2 2x (double pole, no residue); '
       'with (-eps, eps) removed the integral is 8(1 - coth 2 eps) ~ -4/eps, divergent', ok_d and ok_i)

us = 0.2; M3 = mass(3, 1); worst = 0
for t in (0.0, 1.9, 5.1):
    E, P, eom, mt = energy_momentum(3, [(1, 0.1+1j*us, 0.3+0.2j), (2, 0.1-1j*us, -0.1+1.1j)], t=t, N=80001)
    worst = max(worst, abs(E - 2*M3*np.cosh(0.1)*np.cos(us)), 1.0 if mt < 0.04 else 0.0)
report(f'a_2 breather with u = 0.2 (A = {A_int(3, 1, 2, 2j*us).real:.3f} > 0, tau has zeros): slices avoiding '
       f'the zeros give E = 2M cosh th0 cos u [err {worst:.1e}]', worst < 1e-8)

sol = [(1, 0.3+0.4j, 0.2+0.5j)]
f = lambda v: [taus(3, sol, np.array([v[0]]), v[1])[1][0].real, taus(3, sol, np.array([v[0]]), v[1])[1][0].imag]
xz, tz = fsolve(f, [0.3, 0.5], xtol=1e-13)
jac = np.array([[taus(3, sol, np.array([xz+d*(i == 0)]), tz+d*(i == 1))[1][0] for i in (0, 1)] for d in (1e-6,)])[0]
jac = (jac - taus(3, sol, np.array([xz]), tz)[1][0])/1e-6
xl = np.linspace(-25, 25, 50001)
wind = [np.unwrap(np.angle(taus(3, sol, xl, tz + d)[1]))[[0, -1]] @ [-1, 1] for d in (-1e-5, 1e-5)]
report(f'complex-rapidity a_2 soliton: tau_1 has an isolated zero at (x, t) = ({xz:.4f}, {tz:.4f}) '
       f'(Jacobian det {abs(np.imag(np.conj(jac[0])*jac[1])):.2f} != 0); as it crosses the real line the winding of '
       f'tau_1 changes by 2 pi, so the topological charge jumps by a root',
       max(abs(np.array(f([xz, tz])))) < 1e-12 and abs(np.imag(np.conj(jac[0])*jac[1])) > 0.1
       and abs(abs(wind[1] - wind[0]) - 2*np.pi) < 1e-3)

# ------------------------------------------------------------------ 5. breather band, a_5^(1)
h = 6; ok_A = ok_reg = ok_zero = ok_real = True
for a in (1, 2, 3):
    for uu in np.linspace(0.05, np.pi/2 - 0.05, 25):
        Aab = A_int(h, a, h - a, 2j*uu).real
        ok_A &= abs(Aab/(1 - np.sin(np.pi*a/h)**2/np.cos(uu)**2) - 1) < 1e-12
        ok_A &= (Aab < 0) == (uu > abs(np.pi/2 - np.pi*a/h))
    uu = (abs(np.pi/2 - np.pi*a/h) + np.pi/2)/2                # inside the regular band
    xg, tg = np.meshgrid(np.linspace(-12, 12, 241), np.linspace(0, 2*np.pi/(2*np.sin(np.pi*a/h)*np.sin(uu)), 121))
    T = taus(h, [(a, 1j*uu, 0.4+0.3j), (h - a, -1j*uu, 0.4-0.3j+1j*np.pi)], xg.ravel(), tg.ravel())
    ok_reg &= np.min(T.real) > 0.999
    T2 = taus(h, [(a, 1j*uu, 0.4+0.3j), (h - a, -1j*uu, 0.4-0.3j)], xg.ravel(), tg.ravel())   # xi_2 = conj xi_1
    ok_real &= np.max(np.abs(T2.imag)) < 1e-8*np.max(np.abs(T2)) and np.min(T2.real) < 0
for a in (1, 2):                                               # 0 < A < 1 below the band
    uu = 0.5*abs(np.pi/2 - np.pi*a/h)
    for xi2 in (0.3+0.9j, -0.7+2.0j, 1.1-0.4j):
        sols = [(a, 1j*uu, 0.2+0.1j), (h - a, -1j*uu, xi2)]
        found = False
        for x0 in np.linspace(-6, 6, 13):
            for t0 in np.linspace(0, 12, 13):
                g = lambda v: [taus(h, sols, np.array([v[0]]), v[1])[0][0].real, taus(h, sols, np.array([v[0]]), v[1])[0][0].imag]
                v, info, ier, msg = fsolve(g, [x0, t0], full_output=True, xtol=1e-13)
                if ier == 1 and max(abs(np.array(g(v)))) < 1e-10:
                    found = True; break
            if found:
                break
        ok_zero &= found
report('a_5^(1): A_{a,h-a}(2iu) = 1 - sin^2(pi a/h)/cos^2 u, negative exactly for |pi/2 - pi a/h| < u < pi/2',
       ok_A)
report('a_5^(1): in that band xi_2 = conj(xi_1) + i pi gives Re tau_j > 0 (regular breather, mass in (0, 2 M_a sin(pi a/h)))',
       ok_reg)
report('a_5^(1): for 0 < A < 1 tau_0 vanishes somewhere, for each of three choices of xi_2', ok_zero)
report('exr-breather-real: xi_2 = conj(xi_1) gives real tau_j that take negative values (Theta-fixed, singular field)',
       ok_real)
