"""Classical (linearized) reflection matrices about the real-coupling vacua of a_n^(1) on the half-line x <= 0, against
the B -> 0 limit of Delius-Gandenberger's solitonic amplitudes (Part VII, ch. 31: sec-classical-reflection-factors,
prp-solitonic-vacua, prp-self-conjugate-vacuum, prp-dg-weights, sec-boundary-missing-charges, tbl-boundary-types,
tbl-missing-vacua; ch. 35: prp-solitonic-amplitude, sec-solitonic-boundaries).
DG = Delius-Gandenberger hep-th/9904002, Bow = Bowcock hep-th/9609233, BP = Bowcock-Perkins hep-th/9909174.

Conventions (book, sec-boundary-potential; m = 1): alpha_j = e_j - e_{j+1 mod h}, j = 0..n, h = n+1,
phi = -(1/b) sum_j alpha_j ln tau_j, e^{b alpha_j.phi} = tau_{j-1} tau_{j+1}/tau_j^2 =: e^{psi_j}.  The fluctuation
eps (rescaled by b) about a static vacuum obeys
   eps'' = M^2(x) eps - Omega^2 eps,   M^2(x) = sum_j alpha_j alpha_j^T e^{psi_j},
   eps'(0) + H eps(0) = 0,             H = (1/2) sum_j C_j alpha_j alpha_j^T e^{psi_j(0)/2},
and K[c', c] is defined by eps -> rho_c e^{ipx} + sum_c' K[c', c] rho_c' e^{-ip'x} + (decaying), with the Fourier
polarizations rho_c = (w^{jc})_j, w = e^{2 pi i/h} (particle c; m_c = 2 sin(pi c/h)).  Bow (4.1)-(4.4).
Method (Bow sec. 4, BP sec. 5; exact, no ODE): a linear wave of species c is one more Hirota soliton,
E = e^{sigma x - i Omega t}, sigma = m_c cosh(theta_k), added to the static background solitons to first order,
   dtau_j = w^{jc} e^{sigma x} sum_mu prod A_ik^{mu_i mu_k} prod_i (A(c, a_i) c_i w^{j a_i} e^{m_{a_i} x})^{mu_i},
   A(a, b) = (cosh theta - cos(pi(a-b)/h))/(cosh theta - cos(pi(a+b)/h)),   eps = -sum_j alpha_j dtau_j/tau_j;
open channels have sigma = -+ i p (in/out), closed ones sigma = kappa > 0.  Everything is a finite sum of exponentials,
differentiated exactly; the boundary condition gives an n x n linear system (mpmath).
DG amplitudes: K_a^(b)(theta) at B = 0 is K_uniform(a)*M_sol(a, b) of particle_reflection.py (DG (3.16), (5.4)-(5.5)).

Checks:
 1. The dressed waves solve the linearized field equation (residual at 30 and 60 digits) for pair, middle-soliton
    and multi-soliton backgrounds.
 2. Flat vacuum phi = 0, uniform C = +-1: K = -((a)(h-a))^{-C} (eq-linear-reflection), diagonal, a_2..a_5.
 3. Continued stationary pair (eq-vacuum-pair), every regular chi-interval of every species a < h/2, a_2..a_5, two
    chi per interval: the tau_j are those of boundary_vacua.pair_taus; exactly one sign pattern C is satisfied, the
    pattern of the soliton charge; its class by DG (4.8)-(4.9) (prp-dg-weights) is b = a.
 4. Main check (DG (5.7) <-> linearization): for every background of check 3 and for the self-conjugate family
    tau_j = 1 + (-1)^j d e^{2x} (C = -1, b = h/2), the linearized K is diagonal and K_cc = DG K_c^(b)(B -> 0) for
    every open channel, at Omega = m_1 cosh 0.6 (channels 1, h-1), m_2 cosh 0.6 and 2.3 (all open); flux-normalized
    K is unitary.
 5. Bow (4.14) with (4.21)-(4.22), -1/((b)(h-b)) (h/2+a-b)(h/2-a+b)/((h/2+a+b)(h/2-a-b)), is DG (5.7); the diagonal
    entries of BP Table 2 for the regular (pair and flat) rows, r = 2..5, are DG's amplitudes of the class.
 6. Boundary-singular endpoint of the a_5 species-2 pair, chi -> 5 pi/3 (tau_2(0) = tau_5(0) = 0), boundary Hessian
    of the target pattern (+,+,-,-,-,-) (class b = 2; BP Table 2 lists it as diagonal): along the regular solutions
    chi = 5pi/3 +- delta the limit exists, is the same from both sides, and is NOT diagonal: K_11 -> 0.155633 + 0.677155i,
    |K_15| -> 0.719193 at theta = 0.6.  The limit depends on the path: along the non-solution configurations
    d = (1-e)/(1+s), chi = 5pi/3 (boundary residual -> 0) it is |K_15| -> 0.937437; for (+,+,-,-,+,-) it depends on the
    side (|K_15| -> 0.0913 / 0.2721); the regular pattern's own Hessian gives DG's value from its own side only.
 7. Middle-soliton endpoint d -> 1 (a_3 (+,-,+,-), class b = 2; a_5 (+,-,+,-,-,-), b = 3): the off-diagonal entries
    vanish like (1-d), K_11 -> DG K_1^(h/2)(B -> 0) = C = -1 value; the same limit along a non-solution path with a
    species-1 pair admixture of amplitude (1-d)^2.
 8. a_4 exceptional background of BP (sec. 4.4.1; species 1, 4, 2, 3, phases chi_1 = 9pi/5, chi_2 = 18pi/5): the
    regular-on-x<0 static solution with tau_2(0) = tau_4(0) = 0 and tau_i'/tau_i(0) = m_1 + m_2 - sqrt 5 (BP's C) for
    pattern (+,+,-,+,-) (class b = 2); the linearized K approached along two different paths of regular configurations
    has the same non-diagonal limit, K_11 -> 0.32535 + 0.93408i, |K_14| -> 0.14713, not DG's 0.34282 + 0.93940i.
Runtime: about two minutes.
"""
import contextlib
import io
import itertools
import os
import sys

import mpmath as mp
import numpy as np

from boundary_common import report

with contextlib.redirect_stdout(io.StringIO()):
    import boundary_vacua as BV                      # runs its own checks (a few seconds); output suppressed

_here = next(d for d in [os.path.dirname(os.path.abspath(__file__))] + sys.path
             if os.path.exists(os.path.join(d or '.', 'particle_reflection.py')))
_PR = {}
exec(open(os.path.join(_here or '.', 'particle_reflection.py')).read()
     .split('\n# ---------------------------------------------------------------- 1. blocks')[0], _PR)
blk, K_uniform, M_sol = _PR['blk'], _PR['K_uniform'], _PR['M_sol']
I = mp.mpc(0, 1)


def K_dg(c, b, n, th):
    """DG (5.5) at B = 0: K_c(theta) M_c^(b)(theta)."""
    return (K_uniform(c, n) * M_sol(c, b, n))(th, 0)


# ---------------------------------------------------------------- Hirota backgrounds and dressed linear waves
class ExpSum:
    def __init__(self, terms):
        self.t = terms

    def __call__(self, x, k=0):
        return mp.fsum(c * r**k * mp.exp(r * x) for c, r in self.t)


class Background:
    """static real-coupling Hirota background of a_n^(1): sols = [(species, coefficient)], rapidity 0."""
    def __init__(self, n, sols):
        self.n, self.h, self.sols = n, n + 1, sols
        self.w = mp.exp(2 * I * mp.pi / self.h)
        N = len(sols)
        self.Ast = {(i, k): self.A(sols[i][0], sols[k][0], mp.mpf(1)) for i in range(N) for k in range(i + 1, N)}
        self.tau = [self._sum(j, mp.mpc(1), mp.mpf(0), [mp.mpf(1)] * N) for j in range(self.h)]

    def mass(self, c):
        return 2 * mp.sin(mp.pi * c / self.h)

    def A(self, a, b, ch):
        return (ch - mp.cos(mp.pi * (a - b) / self.h)) / (ch - mp.cos(mp.pi * (a + b) / self.h))

    def _sum(self, j, c0, r0, Aw):
        terms = []
        for mu in itertools.product([0, 1], repeat=len(self.sols)):
            coef, rate = c0, r0
            for i, (a, ci) in enumerate(self.sols):
                if mu[i]:
                    coef *= Aw[i] * ci * self.w**(j * a); rate += self.mass(a)
            for (i, k), Aik in self.Ast.items():
                if mu[i] and mu[k]:
                    coef *= Aik
            terms.append((coef, rate))
        return ExpSum(terms)

    def wave(self, c, sigma):
        Aw = [self.A(c, a, sigma / self.mass(c)) for a, _ in self.sols]
        return [self._sum(j, self.w**(j * c), sigma, Aw) for j in range(self.h)]

    def epsi(self, x):
        """e^{psi_j(x)}"""
        h = self.h; t = [tt(x) for tt in self.tau]
        return [t[(j - 1) % h] * t[(j + 1) % h] / t[j]**2 for j in range(h)]


def eps_of(bg, dtau, x):
    """eps, eps', eps'' at x of the wave dtau: eps_k = q_{k-1} - q_k, q_j = dtau_j/tau_j."""
    h = bg.h; q = []
    for j in range(h):
        t0, t1, t2 = (bg.tau[j](x, k) for k in range(3))
        d0, d1, d2 = (dtau[j](x, k) for k in range(3))
        q.append((d0 / t0, d1 / t0 - d0 * t1 / t0**2,
                  d2 / t0 - 2 * d1 * t1 / t0**2 - d0 * t2 / t0**2 + 2 * d0 * t1**2 / t0**3))
    return [[q[(k - 1) % h][r] - q[k][r] for k in range(h)] for r in range(3)]


def alpha_form(h, coeffs, v):
    """sum_j coeffs_j alpha_j (alpha_j . v)"""
    out = [mp.mpc(0)] * h
    for j in range(h):
        s = v[j] - v[(j + 1) % h]
        out[j] += coeffs[j] * s; out[(j + 1) % h] -= coeffs[j] * s
    return out


def pde_residual(bg, dtau, x, Om):
    e0, _, e2 = eps_of(bg, dtau, x)
    Mv = alpha_form(bg.h, bg.epsi(x), e0)
    return max(abs(e2[k] - Mv[k] + Om**2 * e0[k]) for k in range(bg.h)) / max(abs(v) for v in e0 + e2)


def bc_vec(bg, dtau, C):
    e0, e1, _ = eps_of(bg, dtau, mp.mpf(0))
    Hv = alpha_form(bg.h, [C[j] * mp.sqrt(v) / 2 for j, v in enumerate(bg.epsi(mp.mpf(0)))], e0)
    return [e1[k] + Hv[k] for k in range(bg.h)]


def kmatrix(bg, C, Om):
    """open channels, K[(c', c)] (polarization basis rho_c), flux-unitarity defect, max BC residual (all h rows)."""
    n, h = bg.n, bg.h
    opn = [c for c in range(1, h) if bg.mass(c) < Om]
    p = {c: mp.sqrt(Om**2 - bg.mass(c)**2) for c in opn}
    unk = [(c, -I * p[c]) for c in opn] + [(c, mp.sqrt(bg.mass(c)**2 - Om**2)) for c in range(1, h) if c not in opn]
    cols = [bc_vec(bg, bg.wave(c, s), C) for c, s in unk]
    M = mp.matrix([[cols[i][k] for i in range(n)] for k in range(n)])
    nf = lambda c: -(1 - bg.w**(-c))                  # eps -> nf(c) rho_c e^{sigma x} for unit Hirota amplitude
    K, res = {}, 0
    for c in opn:
        r = bc_vec(bg, bg.wave(c, I * p[c]), C)
        sol = mp.lu_solve(M, mp.matrix([-v for v in r[:n]]))
        res = max(res, max(abs(r[k] + mp.fsum(cols[i][k] * sol[i] for i in range(n))) for k in range(h)))
        for i, c2 in enumerate(opn):
            K[(c2, c)] = sol[i] * nf(c2) / nf(c)
    U = mp.matrix([[K[(a, b)] * mp.sqrt(p[a] / p[b]) for b in opn] for a in opn])
    return opn, K, mp.mnorm(U.H * U - mp.eye(len(opn)), 1), res


def offdiag(opn, K):
    return max([abs(K[(a, b)]) for a in opn for b in opn if a != b] + [mp.mpf(0)])


def bc_patterns(bg, tol):
    """sign patterns C whose (zeroth-order) boundary condition the background satisfies at x = 0."""
    h = bg.h; t0 = [bg.tau[j](0) for j in range(h)]; t1 = [bg.tau[j](0, 1) for j in range(h)]
    eh = [mp.sqrt(v) for v in bg.epsi(mp.mpf(0))]
    dphi = [t1[(k - 1) % h] / t0[(k - 1) % h] - t1[k] / t0[k] for k in range(h)]
    return [C for C in itertools.product([1, -1], repeat=h)
            if max(abs(dphi[k] + sum(C[j] * eh[j] * ((k == j) - (k == (j + 1) % h)) for j in range(h))) for k in range(h)) < tol]


def dg_class(C):
    """DG (4.8)-(4.9) (prp-dg-weights) in DG labels (book alpha_j = DG alpha_{j+1}, book alpha_n = DG alpha_0)."""
    h = len(C); Cdg = [C[h - 1]] + list(C[:h - 1])
    ks = [k for k in range(1, h) if Cdg[k] == -1]
    b = sum((-1)**(len(ks) - i) * k for i, k in enumerate(ks, 1))
    return min(b % h, h - b % h)


def pair(n, a, chi, d=None):
    h = n + 1; s = mp.sin(mp.pi * a / h); d = 1 / (1 + s) if d is None else d
    return Background(n, [(a, d * mp.exp(I * chi)), (h - a, d * mp.exp(-I * chi))])


def regular(bg, xs):
    return min(bg.tau[j](x).real for j in range(bg.h) for x in xs) > 0


TH = mp.mpf('0.6')
XS = [mp.mpf(0)] + [-mp.mpf(k) / 40 for k in range(1, 480)]

# ---------------------------------------------------------------- checks 1-4: regular vacua, 30 digits
mp.mp.dps = 30
w = {'pde': 0, 'flat': 0, 'bv': 0, 'off': 0, 'diag': 0, 'unit': 0, 'res': 0}
ok_bg = True
ncases = 0
for n in (2, 3, 4, 5):
    h = n + 1
    for C in (1, -1):                                                       # check 2
        bg = Background(n, [])
        for Om in (bg.mass(1) * mp.cosh(mp.mpf('0.4')), mp.mpf('2.3')):
            opn, K, u, res = kmatrix(bg, (C,) * h, Om)
            for c in opn:
                thc = mp.acosh(Om / bg.mass(c))
                w['flat'] = max(w['flat'], abs(K[(c, c)] + (blk(c, thc, h) * blk(h - c, thc, h))**(-C)), offdiag(opn, K))
    for a in range(1, h // 2 + 1):
        if 2 * a == h:
            cases = [('d', mp.mpf(d)) for d in ('0.4', '-0.7')]
        else:
            ends = sorted({float((mp.pi - 2 * mp.pi * j * a / h) % (2 * mp.pi)) for j in range(h)})
            gaps = [(ends[k], ends[k + 1] if k + 1 < len(ends) else ends[0] + 2 * np.pi) for k in range(len(ends))]
            cases = [('chi', mp.mpf(lo) + mp.mpf(fr) * (mp.mpf(hi) - mp.mpf(lo))) for lo, hi in gaps for fr in ('0.23', '0.71')]
        for kind, par in cases:
            ncases += 1
            if kind == 'chi':
                bg = pair(n, a, par); chi = float(par)
            else:
                bg = Background(n, [(a, par)]); chi = float(mp.acos(par))
            ref = BV.pair_taus(n, a, chi, np.array([-0.7, -0.2, 0.0]))
            w['bv'] = max(w['bv'], max(abs(complex(bg.tau[j](mp.mpf(x))) - ref[j][0][i]) for j in range(h) for i, x in enumerate((-0.7, -0.2, 0.0))))
            pats = bc_patterns(bg, mp.mpf(10)**-20)
            ok_bg &= len(pats) == 1 and regular(bg, XS)
            C = pats[0]
            if kind == 'chi':
                lam = BV.as_weight(BV.soliton_charge(n, a, chi), h)
                ok_bg &= BV.pattern(np.array(lam, float), BV.book_roots(n)) == C and sum(lam) == a
            else:
                ok_bg &= C == (-1,) * h
            b = dg_class(C)
            ok_bg &= b == min(a, h - a)
            for Om in (bg.mass(1) * mp.cosh(TH), bg.mass(2) * mp.cosh(TH), mp.mpf('2.3')):
                opn, K, u, res = kmatrix(bg, C, Om)
                w['off'] = max(w['off'], offdiag(opn, K))
                w['diag'] = max(w['diag'], max(abs(K[(c, c)] - K_dg(c, b, n, mp.acosh(Om / bg.mass(c)))) for c in opn))
                w['unit'] = max(w['unit'], u); w['res'] = max(w['res'], res)
                c0 = opn[0]; p0 = mp.sqrt(Om**2 - bg.mass(c0)**2)
                w['pde'] = max(w['pde'], max(pde_residual(bg, bg.wave(c0, sg), mp.mpf(x), Om)
                                            for sg in (I * p0, -I * p0) for x in ('-0.37', '-1.9')))
print(f'   {ncases} backgrounds (all regular chi-intervals of the pairs and the self-conjugate family, a_2..a_5), 3 frequencies each; '
      + ', '.join(f'max {k} = {mp.nstr(v, 2)}' for k, v in w.items()))
report('1. the dressed Hirota waves solve the linearized field equation eps\'\' = M^2(x) eps - Omega^2 eps (residual < 1e-25)', w['pde'] < 1e-25)
report('2. flat vacuum, C = +-1: K_a = -((a)(h-a))^{-C} (eq-linear-reflection), diagonal, a_2..a_5', w['flat'] < 1e-25)
report('3. continued pair = boundary_vacua.pair_taus; regular on x <= 0; exactly one pattern, that of the soliton charge; '
       'DG (4.8)-(4.9) class b = a (b = h/2, C = -1 for the self-conjugate family)', ok_bg and w['bv'] < 1e-12)
report('4. DG (5.7): about every regular continued vacuum of a_2..a_5 the linearized K is diagonal and equals DG K_c^(b) at B -> 0 '
       'in every open channel (Omega = m_1 cosh 0.6, m_2 cosh 0.6, 2.3); flux-normalized K unitary',
       w['off'] < 1e-24 and w['diag'] < 1e-24 and w['unit'] < 1e-24 and w['res'] < 1e-24)

# ---------------------------------------------------------------- check 5: Bowcock (4.14)-(4.22) and BP Table 2
ok = True
for n in range(2, 7):
    h = n + 1; hh = mp.mpf(h) / 2
    for a in range(1, h):
        for b in range(1, h):
            for th in (mp.mpc('0.37', '0.11'), mp.mpf('1.4')):
                bow = -1 / (blk(b, th, h) * blk(h - b, th, h)) * (blk(hh + a - b, th, h) * blk(hh - a + b, th, h)) / (
                    blk(hh + a + b, th, h) * blk(hh - a - b, th, h))
                ok &= abs(bow - K_dg(b, a, n, th)) < 1e-25
B_ = lambda th, h, *xs: mp.fprod(blk(mp.mpf(x), th, h)**e for x, e in xs)
bp2 = [  # (BP pattern A_i = -C_i, particle, value)
    ('-++', 1, lambda t: -B_(t, 3, (0.5, 1), (1.5, 2), (2.5, 1), (1, -1), (2, -1))),
    ('++++', 1, lambda t: -B_(t, 4, (1, 1), (3, 1))), ('++++', 2, lambda t: -B_(t, 4, (2, 2))),
    ('--++', 1, lambda t: B_(t, 4, (2, 2), (1, -1), (3, -1))), ('--++', 2, lambda t: -B_(t, 4, (1, 2), (3, 2), (2, -2))),
    ('-++++', 1, lambda t: -B_(t, 5, (0.5, 1), (1.5, 1), (3.5, 1), (4.5, 1), (1, -1), (4, -1))),
    ('-++++', 2, lambda t: -B_(t, 5, (1.5, 1), (2.5, 2), (3.5, 1), (2, -1), (3, -1))),
    ('---++', 1, lambda t: -B_(t, 5, (2.5, 2), (0.5, -1), (1, -1), (4, -1), (4.5, -1))),
    ('---++', 2, lambda t: -B_(t, 5, (0.5, 1), (1.5, 1), (3.5, 1), (4.5, 1), (2, -1), (3, -1))),
    ('-++-++', 1, lambda t: B_(t, 6, (2, 1), (4, 1), (1, -1), (5, -1))), ('-++-++', 2, lambda t: -B_(t, 6, (1, 1), (3, 2), (5, 1), (2, -1), (4, -1))),
    ('-++-++', 3, lambda t: -B_(t, 6, (2, 2), (4, 2), (3, -2))),
    ('----++', 1, lambda t: -B_(t, 6, (3, 2), (1, -2), (5, -2))), ('----++', 2, lambda t: mp.mpf(1)),
    ('----++', 3, lambda t: -B_(t, 6, (1, 2), (5, 2), (3, -2))),
    ('++++++', 1, lambda t: -B_(t, 6, (1, 1), (5, 1))), ('++++++', 2, lambda t: -B_(t, 6, (2, 1), (4, 1))),
    ('++++++', 3, lambda t: -B_(t, 6, (3, 2)))]
for pat, c, f in bp2:
    C = tuple(-1 if ch == '+' else 1 for ch in pat); n = len(C) - 1
    for th in (mp.mpf('0.6'), mp.mpc('1.1', '0.3')):
        ok &= abs(f(th) - K_dg(c, dg_class(C), n, th)) < 1e-25
report('5. Bow (4.14) with (4.21)-(4.22) equals DG (5.7) (n = 2..6, all a, b); BP Table 2 diagonal entries of the regular rows '
       '(-++; ++++, --++; -++++, ---++; ++++++, -++-++, ----++) are DG K_c^(b) at B -> 0', ok)

# ---------------------------------------------------------------- check 6: a_5 pair endpoint, 60 digits
mp.mp.dps = 60
n, h = 5, 6
s = mp.sin(mp.pi / 3); d0 = 1 / (1 + s); ce = 5 * mp.pi / 3
Om = 2 * mp.sin(mp.pi / 6) * mp.cosh(TH)
P = (1, 1, -1, -1, -1, -1); Q = (1, 1, -1, -1, 1, -1)
lim = {}
ok_pde = pde_residual(pair(n, 2, ce + mp.mpf('1e-3')), pair(n, 2, ce + mp.mpf('1e-3')).wave(1, I * mp.sqrt(Om**2 - 1)), mp.mpf('-0.01'), Om) < 1e-50
for side in (1, -1):
    for dl in ('1e-4', '1e-6'):
        bg = pair(n, 2, ce + side * mp.mpf(dl))
        reg = bc_patterns(bg, mp.mpf(10)**-14)
        ok_pde &= len(reg) == 1 and regular(bg, XS) and reg[0] == ((1, -1, -1, 1, -1, -1) if side == 1 else (-1, 1, -1, -1, 1, -1))
        for C in (P, Q, reg[0], ((-1, 1, -1, -1, 1, -1) if side == 1 else (1, -1, -1, 1, -1, -1))):
            opn, K, u, res = kmatrix(bg, C, Om)
            ok_pde &= u < 1e-20 and res < 1e-20
            lim[(side, dl, C if C in (P, Q) else ('own' if C == reg[0] else 'other'))] = (K[(1, 1)], abs(K[(5, 1)]))
Kp = lim[(1, '1e-6', P)]
conv = max(abs(lim[(sd, dl, P)][0] - Kp[0]) + abs(lim[(sd, dl, P)][1] - Kp[1]) for sd in (1, -1) for dl in ('1e-4', '1e-6'))
dgv = K_dg(1, dg_class(P), 5, TH)
print(f'   a_5 (+,+,-,-,-,-), chi -> 5pi/3: K_11 = {mp.nstr(Kp[0], 9)}, |K_15| = {mp.nstr(Kp[1], 9)} (spread over both sides, '
      f'delta = 1e-4, 1e-6: {mp.nstr(conv, 2)}); DG b = {dg_class(P)} = BP Table 2: {mp.nstr(dgv, 9)}')
ok6a = ok_pde and conv < 1e-6 and abs(Kp[0] - mp.mpc('0.15563283', '0.67715549')) < 1e-7 and abs(Kp[1] - mp.mpf('0.71919316')) < 2e-7 \
    and abs(Kp[0] - dgv) > 0.5 and dg_class(P) == 2
report('6a. a_5 (+,+,-,-,-,-) about the species-2 pair endpoint, as the limit along the regular solutions chi = 5pi/3 +- delta: the same '
       'from both sides and not diagonal, |K_15| -> 0.71919, K_11 -> 0.15563 + 0.67716i, not BP/DG\'s (2)(4)/((1)(5))', ok6a)
dpath = []
for e in ('1e-7', '1e-9'):
    bg = pair(n, 2, ce, d=d0 * (1 - mp.mpf(e)))
    opn, K, u, res = kmatrix(bg, P, Om)
    dpath.append((K[(1, 1)], abs(K[(5, 1)]), u, regular(bg, XS)))
q1, q2 = lim[(1, '1e-6', Q)], lim[(-1, '1e-6', Q)]
own, other = lim[(1, '1e-6', 'own')], lim[(1, '1e-6', 'other')]
print(f'   path dependence: d-path (non-solutions) |K_15| -> {mp.nstr(dpath[-1][1], 7)}, K_11 -> {mp.nstr(dpath[-1][0], 7)}; '
      f'(+,+,-,-,+,-): |K_15| -> {mp.nstr(q1[1], 5)} (chi > 5pi/3), {mp.nstr(q2[1], 5)} (chi < 5pi/3); regular pattern\'s own Hessian: '
      f'|K_15| = {mp.nstr(own[1], 2)} on its side, {mp.nstr(other[1], 5)} from the other side')
report('6b. the limit depends on the approach: non-solution path d -> 1/(1+s) gives |K_15| -> 0.93744 for (+,+,-,-,-,-); for (+,+,-,-,+,-) the two '
       'sides give 0.0913 and 0.2721; the regular pattern\'s own Hessian gives DG only from its own side',
       all(x[2] < 1e-30 and x[3] for x in dpath) and abs(dpath[-1][1] - mp.mpf('0.937437')) < 2e-6 and abs(dpath[0][1] - dpath[1][1]) < 1e-6
       and abs(q1[1] - mp.mpf('0.0913009')) < 2e-6 and abs(q2[1] - mp.mpf('0.272069')) < 2e-6
       and own[1] < 1e-20 and abs(own[0] - dgv) < 1e-20 and other[1] > 0.5)

# ---------------------------------------------------------------- check 7: middle-soliton endpoints
ok = True; lines = []
for n, C in ((3, (1, -1, 1, -1)), (5, (1, -1, 1, -1, -1, -1))):
    h = n + 1; b = dg_class(C); Om = 2 * mp.sin(mp.pi / h) * mp.cosh(TH)
    dgv = K_dg(1, b, n, TH)
    ok &= b == h // 2 and abs(dgv + blk(1, TH, h) * blk(h - 1, TH, h)) < 1e-40       # = C = -1 flat value
    offs = []
    for e in ('1e-4', '1e-6'):
        ee = mp.mpf(e)
        for sols in ([(h // 2, 1 - ee)], [(h // 2, 1 - ee), (1, ee**2 * mp.exp(I * mp.mpf('0.3'))), (h - 1, ee**2 * mp.exp(-I * mp.mpf('0.3')))]):
            bg = Background(n, sols)
            opn, K, u, res = kmatrix(bg, C, Om)
            ok &= regular(bg, XS) and u < 1e-30
            offs.append((offdiag(opn, K), abs(K[(1, 1)] - dgv)))
    ok &= all(o[0] < 3e-6 and o[1] < 3e-6 for o in offs[2:]) and offs[0][0] / offs[2][0] > 50
    lines.append(f'a_{n} {"".join("+" if x > 0 else "-" for x in C)}: off-diagonal {mp.nstr(offs[0][0], 2)}, {mp.nstr(offs[2][0], 2)} at 1-d = 1e-4, 1e-6')
print('   ' + '; '.join(lines))
report('7. middle-soliton endpoints d -> 1 (a_3 (+,-,+,-), a_5 (+,-,+,-,-,-)): K -> diagonal DG K_1^(h/2) = -(1)(h-1), '
       'off-diagonal O(1-d), also along a path with a species-1 admixture', ok)

# ---------------------------------------------------------------- check 8: a_4 exceptional background (BP sec. 4.4.1)
n, h = 4, 5
Pe = (1, 1, -1, 1, -1)
def exc(d1, d2):
    c1, c2 = d1 * mp.exp(I * 9 * mp.pi / 5), d2 * mp.exp(I * 18 * mp.pi / 5)
    return Background(n, [(1, c1), (4, mp.conj(c1)), (2, c2), (3, mp.conj(c2))])
d2_of = lambda d1: mp.findroot(lambda d2: exc(d1, d2).tau[2](0).real, mp.mpf('0.5637'))
def bcfin(d1, eta=mp.mpf('1e-30')):
    bg = exc(d1, d2_of(d1)); x = -eta
    t = [(bg.tau[j](x), bg.tau[j](x, 1)) for j in range(h)]
    eh = [mp.sqrt(v) for v in bg.epsi(x)]
    return (t[h - 1][1] / t[h - 1][0] - t[0][1] / t[0][0] + Pe[0] * eh[0] - Pe[h - 1] * eh[h - 1]).real
# the book BC is tau_j'/tau_j - C_j e^{psi_j/2} = kappa for all j; on the nodes 0, 1, 3 next to the zeros e^{psi_j(0)/2} = 0
def lograt(d1, i):
    bg = exc(d1, d2_of(d1)); return (bg.tau[i](0, 1) / bg.tau[i](0)).real
d1s = mp.findroot(lambda d1: lograt(d1, 0) - lograt(d1, 3), mp.mpf('0.2579')); d2s = d2_of(d1s)
bgs = exc(d1s, d2s)
Cbp = 2 * mp.sin(mp.pi / 5) + 2 * mp.sin(2 * mp.pi / 5) - mp.sqrt(5)
ok = (abs(bcfin(d1s)) < 1e-12 and abs(bgs.tau[2](0)) < 1e-40 and abs(bgs.tau[4](0)) < 1e-40 and regular(bgs, XS[1:])
      and all(abs(bgs.tau[i](0, 1) / bgs.tau[i](0) - Cbp) < 1e-25 for i in (0, 1, 3)) and dg_class(Pe) == 2)
Om = 2 * mp.sin(mp.pi / 5) * mp.cosh(TH)
res8 = []
for e in ('1e-5', '1e-7'):
    ee = mp.mpf(e)
    for bg in (exc(d1s, d2s * (1 - ee)), exc(d1s * (1 - ee), d2s)):
        opn, K, u, res = kmatrix(bg, Pe, Om)
        ok &= regular(bg, XS) and u < 1e-30
        res8.append((K[(1, 1)], abs(K[(4, 1)])))
dgv = K_dg(1, 2, 4, TH)
print(f'   a_4 exceptional: d_1 = {mp.nstr(d1s, 12)}, d_2 = {mp.nstr(d2s, 12)}; K_11 -> {mp.nstr(res8[-1][0], 8)}, |K_14| -> {mp.nstr(res8[-1][1], 7)} '
      f'(two paths at 1e-7: spread {mp.nstr(abs(res8[-1][0] - res8[-2][0]), 2)}); DG b = 2: {mp.nstr(dgv, 8)}')
print(f'   a_4 exceptional: finite part of the boundary condition at x = 0 (evaluated at x = -1e-30): {mp.nstr(abs(bcfin(d1s)), 2)}; '
      f'tau_i\'/tau_i(0) - (m_1 + m_2 - sqrt 5) = {mp.nstr(max(abs(bgs.tau[i](0, 1) / bgs.tau[i](0) - Cbp) for i in (0, 1, 3)), 2)}')
ok &= abs(res8[-1][0] - res8[-2][0]) < 1e-5 and abs(res8[-1][1] - mp.mpf('0.147135')) < 1e-5 and abs(res8[-1][0] - dgv) > 0.01
report('8. a_4 exceptional background (BP): tau_2(0) = tau_4(0) = 0, tau_i\'/tau_i(0) = m_1 + m_2 - sqrt 5 (BP\'s C), boundary condition of (+,+,-,+,-); the linearized K for (+,+,-,+,-), '
       'approached along two families of regular configurations, has one non-diagonal limit |K_14| = 0.1471, not DG\'s diagonal b = 2 value', ok)
