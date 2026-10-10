"""Method of images for a_n^(1) solitons on the half-line x <= 0 at imaginary coupling (Part VII, ch. 31).

Book conventions: phi = (i/beta) sum_j alpha_j ln tau_j, tau_j = 1 + omega^{ja} E, E = exp(m_a(x cosh th - t sinh th) + xi),
m_a = 2 m sin(pi a/h); boundary condition d_x phi = -grad B_imag = (i m/beta) sum_j C_j alpha_j e^{i beta alpha_j.phi/2}.
Mirror pair: soliton (a, th, xi_1) + antisoliton (h-a, -th, xi_2), s = sin(pi a/h)/cosh th,
   Neumann (C = 0):     xi_1 + xi_2 = -ln A,          A = 1 - s^2 (interaction coefficient of the pair)
   C_j = C = +-1:        xi_1 + xi_2 = -2 ln(1 + C s)   (Delius hep-th/9807189 (3.19) with eps = -C, xi_D = -xi)
Checks:
1. a_2, a_3 (a=1, 2): Neumann and C = +-1 pairs satisfy the boundary condition at x=0 for all t; the same-species
   mirror cannot be tuned to satisfy it; the C=+1 shift fails for C=-1.
2. Time delay and virtual collision point fitted from the numerical solution:
   Delta t = (2/(m_a sinh th)) ln(1 + C s) (C=+-1),  (1/(m_a sinh th)) ln(1 - s^2) (Neumann);
   d = (C/(2 m_a cosh th)) ln((1+s)/(1-s)) (behind the boundary for C=+1);  hep-th/9807189 (1.4), (1.7), (1.8).
3. Coweight shifts: C_j -> C e^{i pi alpha_j.lambda}; prod_j C_j = C^h is invariant (all patterns for h odd, half for h even).
4. Soliton-preserving: the Theta P mirror (a, -th, xi_2 = -conj xi_1 - ln A) satisfies phi(x) = 2v - conj phi(-x),
   v = Re phi(0) in (2 pi/beta) x coweights, and d_x Im phi(0) = 0.
5. Boundary breathers (th = i eta, u = tan eta): Neumann regular on x<=0 iff A < 0 (u^2 > cot^2(pi a/h)) and
   |rho_-| < arcsinh sqrt|A|; C = +1 with rho_- = 0 regular on x<0 with zeros of tau exactly at x = 0; C = -1 singular.
6. Two pairs (4 solitons): xi_r + xi_rbar = [pair value] - sum_{p != r} (ln A_{rp} + ln A_{r pbar}) (hep-th/9807189 (3.26)).
"""
import itertools
import numpy as np
from scipy.optimize import least_squares, minimize


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


beta, m = 0.9, 1.0


def roots_an(n):
    h = n + 1
    return [np.eye(h)[j] - np.eye(h)[(j + 1) % h] for j in range(h)]


def interaction(a, b, th1, th2, h):
    """book eq-interaction-coefficient: A_ab = (cosh th - cos(A-B))/(cosh th - cos(A+B)), A = pi a/h."""
    A_, B_ = np.pi * a / h, np.pi * b / h
    c = np.cosh(th1 - th2)
    return (c - np.cos(A_ - B_)) / (c - np.cos(A_ + B_))


class Sol:
    """N-soliton tau functions of a_n^(1) (book eq-n-soliton), vectorized in x."""
    def __init__(self, n, species, thetas, xis):
        self.n, self.h = n, n + 1
        self.sp, self.th, self.xi = list(species), np.array(thetas, complex), np.array(xis, complex)
        self.ma = 2 * m * np.sin(np.pi * np.array(self.sp) / self.h)
        self.om = np.exp(2j * np.pi / self.h)
        N = len(self.sp)
        self.subsets = [S for k in range(N + 1) for S in itertools.combinations(range(N), k)]
        self.coef = []
        for S in self.subsets:
            c = 1.0 + 0j
            for i, k in itertools.combinations(S, 2):
                c *= interaction(self.sp[i], self.sp[k], self.th[i], self.th[k], self.h)
            self.coef.append(c)

    def taus(self, x, t):
        x = np.atleast_1d(np.asarray(x, float))
        k = self.ma[:, None] * np.cosh(self.th)[:, None]
        w = self.ma[:, None] * np.sinh(self.th)[:, None]
        expo = k * x[None, :] - w * t + self.xi[:, None]          # (N, nx)
        T = np.zeros((self.h, x.size), complex); dT = np.zeros_like(T)
        for S, c in zip(self.subsets, self.coef):
            if not S:
                T += 1; continue
            e = c * np.exp(expo[list(S)].sum(0)); dk = k[list(S)].sum(0)
            for j in range(self.h):
                ph = self.om ** (j * sum(self.sp[i] for i in S))
                T[j] += ph * e; dT[j] += ph * dk * e
        return T, dT

    def phi(self, x, t):
        """phi on an increasing grid x (starting far left where tau ~ 1), with continuous logarithms."""
        T, dT = self.taus(x, t)
        L = np.log(np.abs(T)) + 1j * np.unwrap(np.angle(T), axis=1)
        al = roots_an(self.n)
        ph = (1j / beta) * sum(al[j][:, None] * L[j][None, :] for j in range(self.h))
        dph = (1j / beta) * sum(al[j][:, None] * (dT[j] / T[j])[None, :] for j in range(self.h))
        return ph, dph, T


def vac_round(v, n):
    """nearest vacuum (2 pi/beta) x coweight lattice to v (coordinates beta alpha_k.v/(2pi) rounded)."""
    al = roots_an(n)
    c = np.array([np.real(beta * (al[k] @ v) / (2 * np.pi)) for k in range(1, n + 1)])
    A = np.array([al[k] for k in range(1, n + 1)] + [np.ones(n + 1)])
    return np.linalg.solve(A, np.concatenate([2 * np.pi / beta * np.round(c), [0]])), np.max(np.abs(c - np.round(c)))


L_box = 30.0
xgrid = np.linspace(-L_box, 0, 1501)


def bc_residual(sol, C, ts, t0):
    """max_t |d_x phi(0) - (i m/beta) sum_j C_j alpha_j e^{i beta alpha_j.(phi(0) - v0)/2}|, v0 = vacuum at x=0 at t0."""
    al = roots_an(sol.n)
    ph0, _, _ = sol.phi(xgrid, t0)
    v0, err = vac_round(ph0[:, -1], sol.n)
    assert err < 1e-6
    worst = 0
    for t in ts:
        ph, dph, T = sol.phi(xgrid, t)
        p0, dp0 = ph[:, -1] - v0, dph[:, -1]
        rhs = (1j * m / beta) * sum(C[j] * al[j] * np.exp(1j * beta * (al[j] @ p0) / 2) for j in range(sol.h))
        worst = max(worst, np.max(np.abs(dp0 - rhs)))
    return worst


def pair(n, a, th, xi1, C):
    h = n + 1; s = np.sin(np.pi * a / h) / np.cosh(th)
    S = -np.log(1 - s**2 + 0j) if C == 0 else -2 * np.log(1 + C * s + 0j)
    return Sol(n, [a, h - a], [th, -th], [xi1, S - xi1])


# ---------------------------------------------------------------- 1. boundary conditions of the mirror pairs
th = 0.6
ts = np.linspace(-6, 6, 13)
for n, a, xi1 in [(2, 1, 0.4 + 0.25j), (2, 2, -0.3 + 0.1j), (3, 1, 0.2 + 1j * np.pi / 4), (3, 2, 0.1 + 0.6j)]:
    xi1 = complex(xi1)
    res = {C: bc_residual(pair(n, a, th, xi1, C), [C] * (n + 1), ts, -40.0) for C in (0, 1, -1)}
    wrong = bc_residual(pair(n, a, th, xi1, 1), [-1] * (n + 1), ts, -40.0)
    report(f'a_{n}^(1), a={a}: soliton + mirror antisoliton satisfies the BC for Neumann, C=+1, C=-1 at 13 times '
           f'[res {max(res.values()):.1e}]; C=+1 shift with C=-1 BC fails [res {wrong:.2f}]',
           max(res.values()) < 1e-8 and wrong > 1e-2)

# same-species mirror: no xi_2 makes it work
for n, a in [(2, 1), (3, 1)]:
    h = n + 1; xi1 = 0.4 + 0.25j if n == 2 else 0.2 + 1j * np.pi / 4
    for C in (0, 1, -1):
        def F(p):
            sol = Sol(n, [a, a], [th, -th], [xi1, p[0] + 1j * p[1]])
            try:
                return [bc_residual(sol, [C] * h, np.linspace(-3, 3, 5), -40.0)]
            except AssertionError:
                return [10.0]
        best = min((minimize(lambda p: F(p)[0], x0, method='Nelder-Mead', options={'xatol': 1e-4, 'fatol': 1e-6, 'maxiter': 200})
                    for x0 in [(0.0, 0.0), (-1.0, 1.0), (1.0, -2.0)]), key=lambda r: r.fun)
        report(f'a_{n}^(1), a={a}, C={C:+d}: no same-species mirror (a, -th, xi_2) satisfies the BC  [min over xi_2 of residual {best.fun:.2f}]',
               best.fun > 1e-2)

# ---------------------------------------------------------------- 2. time delay and virtual collision point (fitted)
def fit_xi(phdata, xs, t, n, a, rap, offset, guess):
    """fit the complex xi of a single soliton (a, rap) to phdata on the window xs (field = offset + soliton)."""
    def r(p):
        sol = Sol(n, [a], [rap], [p[0] + 1j * p[1]])
        xx = np.concatenate([np.linspace(-L_box - 40, xs[0], 400)[:-1], xs])
        ph, _, _ = sol.phi(xx, t)
        d = (ph[:, -len(xs):] + offset[:, None]) - phdata
        return np.concatenate([d.real.ravel(), d.imag.ravel()])
    out = least_squares(r, [guess.real, guess.imag], xtol=1e-14, ftol=1e-14)
    return out.x[0] + 1j * out.x[1], np.max(np.abs(out.fun))


ok_dt = ok_d = True
for n, a in [(2, 1), (3, 1), (3, 2)]:
    h = n + 1; ma = 2 * m * np.sin(np.pi * a / h); s = np.sin(np.pi * a / h) / np.cosh(th)
    xi1 = 0.4 + 0.25j if n == 2 else 0.2 + 1j * np.pi / 4
    for C in (0, 1, -1):
        sol = pair(n, a, th, xi1, C)
        t_in, t_out = -25.0, 25.0
        # incoming soliton on x<0 at t_in
        xc = np.tanh(th) * t_in
        xs = np.linspace(xc - 8, xc + 8, 321)
        xx = np.concatenate([np.linspace(-L_box - 60, xs[0], 4000)[:-1], xs])
        ph, _, _ = sol.phi(xx, t_in)
        xi_in, r1 = fit_xi(ph[:, -321:], xs, t_in, n, a, th, np.zeros(h), xi1 + 0.7 - 0.3j)
        # outgoing antisoliton on x<0 at t_out
        xc = -np.tanh(th) * t_out
        xs = np.linspace(xc - 8, xc + 8, 321)
        xx = np.concatenate([np.linspace(-L_box - 60, xs[0], 4000)[:-1], xs])
        ph, _, _ = sol.phi(xx, t_out)
        xi_out, r2 = fit_xi(ph[:, -321:], xs, t_out, n, h - a, -th, np.zeros(h), sol.xi[1] + 0.5 + 0.2j)
        # incoming mirror on x>0 at t_in (offset = vacuum between the two)
        xc = -np.tanh(th) * t_in
        xs = np.linspace(xc - 8, xc + 8, 321)
        xx = np.concatenate([np.linspace(-L_box - 60, xs[0], 8000)[:-1], xs])
        ph, _, _ = sol.phi(xx, t_in)
        v0, _ = vac_round(ph[:, np.argmin(np.abs(xx))], n)
        xi_m, r3 = fit_xi(ph[:, -321:], xs, t_in, n, h - a, -th, v0, sol.xi[1] - 0.6)
        dt_num = -(xi_in.real + xi_out.real) / (ma * np.sinh(th))
        d_num = -(xi_in.real + xi_m.real) / (2 * ma * np.cosh(th))
        dt_th = np.log(1 - s**2) / (ma * np.sinh(th)) if C == 0 else 2 * np.log(1 + C * s) / (ma * np.sinh(th))
        d_th = C / (2 * ma * np.cosh(th)) * np.log((1 + s) / (1 - s))
        # Delius (1.4), (1.8) with m=1, eps=-C, v = tanh th, his m_a = 2 sin(pi a/h)
        v = np.tanh(th); mD = 2 * np.sin(np.pi * a / h)
        dt_D = (np.sqrt(1 - v**2) / (mD * v) * np.log(1 - mD**2 / 4 * (1 - v**2)) if C == 0 else
                2 * np.sqrt(1 - v**2) / (mD * v) * np.log(1 + C * mD / 2 * np.sqrt(1 - v**2)))
        d_D = C * np.sqrt(1 - v**2) / (2 * mD) * np.log((1 + mD / 2 * np.sqrt(1 - v**2)) / (1 - mD / 2 * np.sqrt(1 - v**2)))
        print(f'   a_{n} a={a} C={C:+d}: Delta t fit {dt_num:+.10f}  formula {dt_th:+.10f}  Delius(1.4/1.8) {dt_D:+.10f};'
              f'  d fit {d_num:+.10f}  formula {d_th:+.10f}  Delius(1.7) {d_D:+.10f}  [fit res {max(r1, r2, r3):.0e}]')
        ok_dt &= abs(dt_num - dt_th) < 1e-8 and abs(dt_th - dt_D) < 1e-12 and max(r1, r2, r3) < 1e-8
        ok_d &= abs(d_num - d_th) < 1e-8 and abs(d_th - d_D) < 1e-12
report('time delay fitted from the numerical solution = (2/(m_a sinh th)) ln(1 + C s); Neumann (1/(m_a sinh th)) ln(1-s^2): '
       'C=+1 delay, C=-1 and Neumann advance (= Delius (1.4), (1.8) with eps = -C)', ok_dt)
report('virtual collision point d = (C/(2 m_a cosh th)) ln((1+s)/(1-s)): behind the boundary for C=+1, in front for C=-1, at it for Neumann (= Delius (1.7))', ok_d)

# ---------------------------------------------------------------- 3. coweight shifts of the boundary condition
def patterns(n, C):
    al = roots_an(n); h = n + 1
    A = np.array([al[k] for k in range(1, n + 1)] + [np.ones(h)])
    cow = [np.linalg.solve(A, np.concatenate([np.eye(n)[k], [0]])) for k in range(n)]       # fundamental coweights
    pats = set()
    for ks in itertools.product([0, 1], repeat=n):
        lam = sum(k * w for k, w in zip(ks, cow))
        pats.add(tuple(int(round((C * np.exp(1j * np.pi * (al[j] @ lam))).real)) for j in range(h)))
    return pats
p2 = patterns(2, 1) | patterns(2, -1); p3 = patterns(3, 1) | patterns(3, -1)
ok = len(p2) == 8 and len(p3) == 8 and all(np.prod(p) == 1 for p in p3)
# numerical: a_2 pair for C=+1 shifted by -(2 pi/beta) w_1 satisfies the BC with C_j = e^{i pi alpha_j.w_1}
al = roots_an(2); A_ = np.array([al[1], al[2], np.ones(3)]); w1 = np.linalg.solve(A_, [1, 0, 0])
Cs = [np.exp(1j * np.pi * (al[j] @ w1)).real for j in range(3)]
sol = pair(2, 1, th, 0.4 + 0.25j, 1)
ph0, _, _ = sol.phi(xgrid, -40.0); v0, _ = vac_round(ph0[:, -1], 2)
worst = 0
for t in ts:
    ph, dph, _ = sol.phi(xgrid, t)
    p0 = ph[:, -1] - v0 - 2 * np.pi / beta * w1
    rhs = (1j * m / beta) * sum(Cs[j] * al[j] * np.exp(1j * beta * (al[j] @ p0) / 2) for j in range(3))
    worst = max(worst, np.max(np.abs(dph[:, -1] - rhs)))
report(f'coweight shifts: C_j -> C e^(i pi alpha_j.lambda) give all 8 sign patterns for a_2, only the 8 with prod C_j = +1 for a_3; '
       f'a_2 C=(+1,+1,+1) solution shifted by -(2pi/beta)w_1 obeys C = {tuple(int(c) for c in Cs)}  [res {worst:.1e}]', ok and worst < 1e-8)

# ---------------------------------------------------------------- 4. soliton-preserving (Theta P) mirrors
xs_full = np.linspace(-L_box, L_box, 12001)
ok = True; lat = 0; resP = 0; resN = 0
for n, a, xi1 in [(2, 1, 0.4 + 0.25j), (3, 1, 0.2 + 1j * np.pi / 4), (3, 2, 0.1 + 0.6j)]:
    h = n + 1
    Aaa = interaction(a, a, th, -th, h).real
    sol = Sol(n, [a, a], [th, -th], [xi1, -np.conj(xi1) - np.log(Aaa)])
    al = roots_an(n)
    vals = []
    for t in np.linspace(-6, 6, 7):
        ph, dph, _ = sol.phi(xs_full, t)
        i0 = len(xs_full) // 2
        c = ph + np.conj(ph[:, ::-1])                         # phi(x) + conj phi(-x) should be constant = 2v
        resP = max(resP, np.max(np.abs(c - c[:, [i0]])))
        vals.append(ph[:, i0].real)
        resN = max(resN, np.max(np.abs(dph[:, i0].imag)))     # Neumann on Im phi
        v, err = vac_round(ph[:, i0].real, n)
        lat = max(lat, err)
        resP = max(resP, np.max(np.abs(c[:, i0] - 2 * v)))
    ok &= np.ptp(np.array(vals), axis=0).max() < 1e-8
report(f'Theta P mirror (a, -th, xi_2 = -conj xi_1 - ln A_aa): phi(x) = 2v - conj phi(-x), Re phi(0) = v in (2pi/beta) coweights, '
       f'd_x Im phi(0) = 0 (a_2 a=1; a_3 a=1,2)  [res {max(resP, resN, lat):.1e}]', ok and max(resP, resN, lat) < 1e-8)
# and the soliton-conjugating pair does not satisfy it
sol = pair(2, 1, th, 0.4 + 0.25j, 0)
ph, dph, _ = sol.phi(xs_full, 1.0); i0 = len(xs_full) // 2
report('the Neumann pair (soliton + antisoliton) violates the soliton-preserving condition', np.max(np.abs(dph[:, i0].imag)) > 1e-2 or vac_round(ph[:, i0].real, 2)[1] > 1e-2)

# ---------------------------------------------------------------- 5. boundary breathers
def breather(n, a, u, rho_m, C, zeta_m=0.3):
    h = n + 1; eta = np.arctan(u); s = np.sin(np.pi * a / h) / np.cos(eta)
    A = np.cos(np.pi * a / h)**2 - u**2 * np.sin(np.pi * a / h)**2
    S = -np.log(A + 0j) if C == 0 else -2 * np.log(1 + C * s + 0j)
    xm = rho_m + 1j * zeta_m
    return Sol(n, [a, h - a], [1j * eta, -1j * eta], [S / 2 + xm, S / 2 - xm]), A, s


def min_tau(sol, xmax, period, nx=400, nt=400):
    """min over j, x in [-12, xmax], t in one period of |tau_j|, refined by local minimization."""
    xs = np.linspace(-12, xmax, nx); best = (np.inf, None)
    for t in np.linspace(0, period, nt, endpoint=False):
        T, _ = sol.taus(xs, t)
        a_ = np.abs(T); k = np.unravel_index(np.argmin(a_), a_.shape)
        if a_[k] < best[0]:
            best = (a_[k], (k[0], xs[k[1]], t))
    j, x0, t0 = best[1]
    f = lambda p: abs(sol.taus([min(p[0], xmax)], p[1])[0][j][0])
    r = minimize(f, [x0, t0], method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 2000})
    return min(best[0], r.fun), (j, min(r.x[0], xmax), r.x[1])


n, a = 2, 1; h = 3
per = lambda u: 2 * np.pi / (2 * m * np.sin(np.pi * a / h) * np.sin(np.arctan(u)))
cot2, tan2 = 1 / np.tan(np.pi * a / h)**2, np.tan(np.pi * a / h)**2
# Neumann, A > 0 (u^2 < cot^2): singular on x <= 0
sol, A, s = breather(n, a, 0.4, 0.0, 0)
mt, _ = min_tau(sol, 0.0, per(0.4))
report(f'Neumann boundary breather with u^2 = 0.16 < cot^2(pi/3) = {cot2:.3f} (A = {A:.3f} > 0): tau_j has a zero on x <= 0  [min|tau| {mt:.1e}]', mt < 1e-6)
# Neumann, A < 0 with cot^2 < u^2 < tan^2: Delius (3.36) "u^2 > tan^2" would exclude it
u = 1.0
sol, A, s = breather(n, a, u, 0.0, 0)
thr = np.arcsinh(np.sqrt(abs(A)))
res = []
for rm in [0.0, 0.5 * thr, 0.95 * thr, 1.05 * thr, 1.5 * thr]:
    sol, A, s = breather(n, a, u, rm, 0)
    res.append(min_tau(sol, 0.0, per(u))[0])
print(f'   Neumann, u^2 = 1 (cot^2 = {cot2:.3f} < u^2 < tan^2 = {tan2:.3f}), A = {A:.3f}: arcsinh sqrt|A| = {thr:.4f}, '
      f'min|tau| on x<=0 for rho_- = (0, .5, .95, 1.05, 1.5) x threshold: ' + ', '.join(f'{r:.1e}' for r in res))
report('Neumann boundary breathers: regular on x<=0 iff A<0 (u^2 > cot^2(pi a/h), not tan^2 as in hep-th/9807189 (3.36)) '
       'and |rho_-| < arcsinh sqrt|A| (|sinh rho_-|, not |cosh rho_-|, in (3.37))', all(r > 1e-3 for r in res[:3]) and all(r < 1e-6 for r in res[3:]))
# C = +1 (Delius eps = -1), rho_- = 0: regular on x < 0, zeros of tau at x = 0
ok = True; info = []
for u in [0.3, 1.0, 2.0]:
    sol, A, s = breather(n, a, u, 0.0, 1)
    m_in, _ = min_tau(sol, -0.05, per(u))
    m0, (j, x0, t0) = min_tau(sol, 0.0, per(u))
    info.append(f'u={u}: min|tau| x<=-0.05 {m_in:.2e}, x<=0 {m0:.1e} at x={x0:.1e}')
    ok &= m_in > 1e-3 and m0 < 1e-6 and abs(x0) < 1e-6
    sol2, _, _ = breather(n, a, u, 0.3, 1)
    ok &= min_tau(sol2, 0.0, per(u))[0] < 1e-6        # rho_- != 0: zero moves into x < 0
print('   C=+1 breathers, rho_-=0: ' + '; '.join(info))
report('C=+1 boundary breathers (rho_- = 0): no zero of tau on x<0, zeros exactly at x=0 (hep-th/9807189 sec. 3.6.1); rho_- = 0.3 singular', ok)
ok = True
for u in [0.3, 1.0, 2.0]:
    for rm in [0.0, 0.4]:
        sol, A, s = breather(n, a, u, rm, -1)
        ok &= min_tau(sol, 0.0, per(u))[0] < 1e-6
report('C=-1 boundary breathers: singular on x<=0 for all tested u, rho_-', ok)

# ---------------------------------------------------------------- 6. two pairs (four solitons)
def two_pairs(n, specs, C):
    """specs = [(a_r, th_r, xi_r)], mirrors (h-a_r, -th_r, xi_rbar) with the shift of hep-th/9807189 (3.26)."""
    h = n + 1
    sp_ = [a for a, _, _ in specs] + [h - a for a, _, _ in specs]
    ths = [t for _, t, _ in specs] + [-t for _, t, _ in specs]
    xis = [x for _, _, x in specs]
    for r, (a, t, x) in enumerate(specs):
        s = np.sin(np.pi * a / h) / np.cosh(t)
        S = -np.log(1 - s**2 + 0j) if C == 0 else -2 * np.log(1 + C * s + 0j)
        for p, (b, tp, _) in enumerate(specs):
            if p != r:
                S -= np.log(interaction(a, b, t, tp, h) + 0j) + np.log(interaction(a, h - b, t, -tp, h) + 0j)
        xis.append(S - x)
    return Sol(n, sp_, ths, xis)


ok = True; worst = 0
for C in (0, 1, -1):
    sol = two_pairs(2, [(1, 0.6, 0.4 + 0.25j), (2, 1.1, -0.5 + 0.1j)], C)
    r = bc_residual(sol, [C] * 3, np.linspace(-5, 5, 11), -60.0)
    worst = max(worst, r)
    bad = Sol(2, sol.sp, sol.th, list(sol.xi[:3]) + [sol.xi[3] + 0.3])
    ok &= bc_residual(bad, [C] * 3, np.linspace(-5, 5, 11), -60.0) > 1e-3
report(f'two soliton-antisoliton pairs (a_2: species 1 and 2, th = 0.6, 1.1) satisfy the BC for C = 0, +1, -1 with the '
       f'shifts of (3.26)  [res {worst:.1e}]; detuned mirror fails', ok and worst < 1e-8)
