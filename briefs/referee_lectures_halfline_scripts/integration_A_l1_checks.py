#!/usr/bin/env python3
"""Integration agent A: independent checks of the findings on the half-line Lecture 1 (and index).

Written from the lecture's and the book's definitions; imports only the book helper
code/boundary_common.py (which prints nothing).  Each check prints PASS/FAIL plus the numbers.
"""
import sys
from pathlib import Path
import numpy as np
import sympy as sp
import mpmath as mp

REPO = Path('/home/gustav/Git/affine_toda/.claude/worktrees/referee-halfline-lectures')
sys.path.insert(0, str(REPO / 'code'))

mp.mp.dps = 40
I = mp.mpc(0, 1)
PI = mp.pi


def rep(name, ok, extra=''):
    print(('PASS ' if ok else 'FAIL ') + name + (f'  [{extra}]' if extra else ''))


# ---------------------------------------------------------------- L1-1 continuity equations
t, x = sp.symbols('t x', real=True)
T = sp.Function('T')(t, x)
X = sp.Function('X')(t, x)
dm = lambda f: sp.diff(f, t) - sp.diff(f, x)   # partial_-  (lecture: d_t - d_x)
dp = lambda f: sp.diff(f, t) + sp.diff(f, x)   # partial_+
ident = sp.simplify((dm(T) - dp(X)) - (sp.diff(T - X, t) - sp.diff(T + X, x)))
rep('L1-1: d_-T - d_+X == d_t(T-X) - d_x(T+X) identically', ident == 0)
ident2 = sp.simplify((dm(T) + dp(X)) - (sp.diff(T + X, t) - sp.diff(T - X, x)))
rep('L1-1: printed d_t(T+X) = d_x(T-X) is d_-T = -d_+X (opposite sign of X)', ident2 == 0)
# counterexample to the printed form (the referee's T=t+x, X=t-x satisfies both forms)
for (TT, XX, lab) in [(t - x, t + x, 'T=t-x,X=t+x'), (t + x, t - x, 'referee T=t+x,X=t-x')]:
    holds = sp.simplify(dm(TT) - dp(XX)) == 0
    printed = sp.simplify(sp.diff(TT + XX, t) - sp.diff(TT - XX, x))
    corrected = sp.simplify(sp.diff(TT - XX, t) - sp.diff(TT + XX, x))
    print(f'     {lab}: d_-T=d_+X holds: {holds}; printed residual {printed}; corrected residual {corrected}')

# explicit s = 1 current of a scalar field with potential V, x^pm = (t pm x)/2, d_pm = d_t pm d_x
phi = sp.Function('phi')(t, x)
V = sp.Function('V')
eom = {sp.diff(phi, t, 2): sp.diff(phi, x, 2) - V(phi).diff(phi)}
T2 = dp(phi)**2
X0 = -2*V(phi)
cur = sp.expand(dm(T2) - dp(X0))
cur = sp.simplify(cur.subs(eom))
rep('L1-1: (T_2, X_0) = ((d_+phi)^2, -2V) obeys d_-T_2 = d_+X_0 on shell', cur == 0)
dens = sp.expand(T2 - X0)
Edens = sp.Rational(1, 2)*sp.diff(phi, t)**2 + sp.Rational(1, 2)*sp.diff(phi, x)**2 + V(phi)
Pdens = -sp.diff(phi, t)*sp.diff(phi, x)          # momentum density: P > 0 for a right-mover f(x - v t)
rep('L1-1: its density T_2 - X_0 = 2(E-density - P-density): spin -1 (book: (T_{s+1},X_{s-1}) has spin -s)',
    sp.simplify(dens - 2*(Edens - Pdens)) == 0)
T2b = dm(phi)**2
densb = sp.expand(T2b - X0)
rep('L1-1: parity partner (d_-phi)^2 + 2V has density 2(E + P)', sp.simplify(densb - 2*(Edens + Pdens)) == 0)
# book normalization Q_{+-1} = E +- P: Q_{-1} = int (T-X)/2, Q_1 = int (Tb-Xb)/2; flux of their sum at x=0
flux = sp.expand(((T2 + X0) - (T2b + X0))/2)   # d/dt(Q_1+Q_-1) = [(T+X) - (Tb+Xb)]/2 at x=0
rep('L1-1: d/dt(Q_1+Q_-1) = 2 phi_t phi_x|_0 = -2 dB/dt with phi_x = -B\'(phi): Sigma_1 = 2B',
    sp.simplify(flux - 2*sp.diff(phi, t)*sp.diff(phi, x)) == 0)

# ---------------------------------------------------------------- L1-4 braid <-> S/K conversion
rng = np.random.default_rng(7)
P = np.zeros((4, 4)); P[0, 0] = P[3, 3] = 1; P[1, 2] = P[2, 1] = 1
I2 = np.eye(2)
S_d = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))     # S_12(theta_1 - theta_2), arbitrary
S_s = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))     # S_12(theta_1 + theta_2)
Ka = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))     # K(theta_1)
Kb = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))     # K(theta_2)
K1a, K2a, K2b = np.kron(Ka, I2), np.kron(I2, Ka), np.kron(I2, Kb)
S21 = lambda S: P @ S @ P
lhs_sk = S_d @ K1a @ S21(S_s) @ K2b
rhs_sk = K2b @ S_s @ K1a @ S21(S_d)
for name, Sc in [('S12 P', lambda S: S @ P), ('P S21', lambda S: P @ S21(S))]:
    lhs_br = Sc(S_d) @ K2a @ Sc(S_s) @ K2b
    rhs_br = K2b @ Sc(S_s) @ K2a @ Sc(S_d)
    e = max(np.abs(lhs_br - lhs_sk).max(), np.abs(rhs_br - rhs_sk).max())
    rep(f'L1-4: Scheck = {name}, K in slot 2: braid sides equal S/K sides term by term (arbitrary matrices)',
        e < 1e-12, f'{e:.1e}')
Sc = lambda S: P @ S                                   # the lecture's P S_12
e = max(np.abs(Sc(S_d) @ K2a @ Sc(S_s) @ K2b - lhs_sk).max(), np.abs(K2b @ Sc(S_s) @ K2a @ Sc(S_d) - rhs_sk).max())
rep('L1-4: with the printed Scheck = P S_12 the terms do not match (12 <-> 21 exchanged)', e > 1e-3, f'{e:.2f}')

# ---------------------------------------------------------------- L1-5 domains of the twisted equation
def braid(a, b):        # Rcheck_{a,b}: (a, b) -> (b, a)
    return ('R', a, b)


def Kslot2(mu):         # 1 (x) K_mu: (*, mu) -> (*, mubar)
    return ('K', mu)


def apply(ops, state):
    """apply ops right to left; state = (left, right) labels; returns final state or None."""
    for op in reversed(ops):
        if op[0] == 'R':
            if state != (op[1], op[2]):
                return None
            state = (op[2], op[1])
        else:
            if state[1] != op[1]:
                return None
            state = (state[0], bar(op[1]))
    return state


def bar(l):
    return l[:-1] if l.endswith('~') else l + '~'


mu, nu = 'm', 'n'
lhs = [braid(bar(nu), bar(mu)), Kslot2(mu), braid(mu, bar(nu)), Kslot2(nu)]
rhs_printed = [Kslot2(nu), braid(bar(mu), nu), Kslot2(mu), braid(mu, nu)]
rhs_fixed = [Kslot2(nu), braid(nu, bar(mu)), Kslot2(mu), braid(mu, nu)]
target = (bar(mu), bar(nu))
rep('L1-5: LHS maps V_mu (x) V_nu to V_mubar (x) V_nubar', apply(lhs, (mu, nu)) == target)
rep('L1-5: printed RHS with Rcheck_{mubar,nu} cannot be applied', apply(rhs_printed, (mu, nu)) is None)
rep('L1-5: RHS with Rcheck_{nu,mubar} maps V_mu (x) V_nu to V_mubar (x) V_nubar', apply(rhs_fixed, (mu, nu)) == target)

# ---------------------------------------------------------------- L1-6 rotation and rapidity
Mrot = np.array([[0, 1], [-1, 0]])        # (x, t_E) -> (t_E, -x)
Mpr = np.array([[0, -1], [-1, 0]])        # printed (x, t_E) -> (-t_E, -x)
Mind = np.array([[0, -1], [1, 0]])        # index (x, t_E) -> (-t_E, x)
rep('L1-6: det of (t_E,-x) map = +1 (rotation), printed (-t_E,-x) = -1 (reflection), index (-t_E,x) = +1',
    round(np.linalg.det(Mrot)) == 1 and round(np.linalg.det(Mpr)) == -1 and round(np.linalg.det(Mind)) == 1)
pt = np.array([-0.7, 0.3])                 # a point with x < 0
rep('L1-6: (t_E,-x) sends x<0 to new time > 0; index map (-t_E,x) sends it to new time < 0',
    (Mrot @ pt)[1] > 0 and (Mind @ pt)[1] < 0)
th = mp.mpc('0.43', '0.0')
rep('L1-6: tanh(i pi/2 - theta) = -coth(theta) (reflected world line x = -t tanh theta -> velocity -coth theta)',
    abs(mp.tanh(I*PI/2 - th) + 1/mp.tanh(th)) < 1e-30)
tt = mp.mpf('0.61')                       # theta-tilde
thv = I*PI/2 - tt
E, p = mp.cosh(thv), mp.sinh(thv)
rep('L1-6: printed eq-lec1-rapidity-rotation inconsistent: with theta = i pi/2 - tt, -iE = -sinh(tt), not +sinh(tt)',
    abs(-I*p - mp.cosh(tt)) < 1e-30 and abs(-I*E + mp.sinh(tt)) < 1e-30)

# ---------------------------------------------------------------- L1-8 crossing index order (SG, book helpers)
from boundary_common import q_phys, sl2_spin, rcheck, k_intertwiner   # noqa: E402

LAM = np.sqrt(2.0)
q = q_phys(LAM)
Cc = np.array([[0, 1], [1, 0]])


def Rsg(z):
    return rcheck(sl2_spin(1, q, z), sl2_spin(1, q, 1.0))[0]


def Kco(y, e0, e1):
    null, _ = k_intertwiner(sl2_spin(1, q, y), sl2_spin(1, q, 1/y), [1/q, 1/q], [e0, e1])
    assert len(null) == 1
    return null[0]


def kv(theta, e0, e1):
    """k^{ab}(theta) = K_{abar}^b(i pi/2 - theta) = K[b, abar]; vector index a*2+b."""
    K = Kco(np.exp(LAM*(1j*np.pi/2 - theta)), e0, e1)
    return np.array([K[b, 1 - a] for a in range(2) for b in range(2)])


def propres(u, v):
    c = np.vdot(v, u)/np.vdot(v, v)
    return np.linalg.norm(u - c*v)/np.linalg.norm(u)


res_fix, res_pr = [], []
for (e0, e1) in [(0.37 + 0.21j, -0.52 + 0.11j), (1.3 - 0.4j, 0.2 + 0.9j)]:
    for th_ in (0.29, 0.83):
        kp, km = kv(th_, e0, e1), kv(-th_, e0, e1)
        Rm = Rsg(np.exp(2*LAM*th_))             # Rcheck(2 theta): Rm[(a,b),(c,d)] = S_cd^{ba}
        res_fix.append(propres(kp, Rm @ km))     # sum_cd S_cd^{ba} k^{cd}(-theta)
        res_pr.append(propres(kp, Rm @ P @ km))  # sum_cd S_dc^{ba} k^{cd}(-theta)  (printed)
rep('L1-8: k^{ab}(t) ~ S_{cd}^{ba}(2t) k^{cd}(-t) (corrected order) for the SG coideal K',
    max(res_fix) < 1e-10, f'{max(res_fix):.1e}')
rep('L1-8: printed order S_{dc}^{ba} fails', min(res_pr) > 1e-3, f'{min(res_pr):.3f}')

# ---------------------------------------------------------------- blocks
def blk(xx, th_, h):
    return mp.sinh(th_/2 + I*PI*xx/(2*h))/mp.sinh(th_/2 - I*PI*xx/(2*h))


# ---------------------------------------------------------------- L1-10 diagonal bootstrap channels, a_2
B = mp.sqrt(2)/4
h3 = 3
curly = lambda p, th_: (blk(p - 1, th_, h3)*blk(p + 1, th_, h3)
                        / (blk(p - 1 + B, th_, h3)*blk(p + 1 - B, th_, h3)))
S11 = lambda th_: curly(1, th_)
S12 = lambda th_: curly(2, th_)
S22, S21 = S11, S12
Sab = {(1, 1): S11, (1, 2): S12, (2, 1): S21, (2, 2): S22}
K1 = lambda th_: 1/(blk(2, th_, h3)*blk(1 - B/2, th_, h3)*blk(3 + B/2, th_, h3))   # DG (3.16), a = 1
K2 = K1                                                                               # K_a = K_{h-a}
th = mp.mpc('0.71', '0.13')
bb_ok = abs(K2(th) - K1(th + I*PI/3)*K1(th - I*PI/3)*S11(2*th))
bb_bad = abs(K2(th) - K1(th + I*PI/3)*K1(th - I*PI/3)*S12(2*th))
rep('L1-10: a_2 bulk bootstrap K_2 = K_1(+i pi/3) S_11(2 theta) K_1(-i pi/3) holds (S_ab channel)', bb_ok < 1e-30,
    f'{mp.nstr(bb_ok, 3)}')
rep('L1-10: the printed channel S_{b abar} = S_12 fails', bb_bad > 1e-3, f'{mp.nstr(bb_bad, 6)}')
v = PI/3 * mp.mpf('0.77')       # generic binding angle; test the dressing factors
uok, ubad, cok, cbad = 0, 0, 0, 0
for b in (1, 2):
    a = 1
    D = lambda th_, b=b: Sab[(b, a)](th_ - I*v)*Sab[(b, a)](th_ + I*v)
    Dbar = lambda th_, b=b: Sab[(3 - b, a)](th_ - I*v)*Sab[(3 - b, a)](th_ + I*v)
    Dp = lambda th_, b=b: Sab[(b, a)](th_ - I*v)*Sab[(b, 3 - a)](th_ + I*v)
    Dpbar = lambda th_, b=b: Sab[(3 - b, a)](th_ - I*v)*Sab[(3 - b, 3 - a)](th_ + I*v)
    uok = max(uok, abs(D(th)*D(-th) - 1)); ubad = max(ubad, abs(Dp(th)*Dp(-th) - 1))
    cok = max(cok, abs(D(I*PI/2 - th) - Dbar(I*PI/2 + th)))
    cbad = max(cbad, abs(Dp(I*PI/2 - th) - Dpbar(I*PI/2 + th)))
rep('L1-10: excited-boundary dressing S_ba(th-iv) S_ba(th+iv) is unitary and crossing-compatible (a_2)',
    uok < 1e-30 and cok < 1e-30, f'{mp.nstr(uok, 3)}, {mp.nstr(cok, 3)}')
rep('L1-10: printed S_ba(th-iv) S_{b abar}(th+iv) violates unitarity and crossing', ubad > 1e-3 and cbad > 1e-3,
    f'{mp.nstr(ubad, 4)}, {mp.nstr(cbad, 4)}')

# ---------------------------------------------------------------- Lee-Yang (h = 6 blocks)
b6 = lambda xx, th_: blk(xx, th_, 6)
SLY = lambda th_: -b6(2, th_)*b6(4, th_)
fx = lambda xx, th_: (mp.sinh(th_) + I*mp.sin(PI*xx))/(mp.sinh(th_) - I*mp.sin(PI*xx))
K1LY = lambda th_, eps=1: eps*b6(1, th_)*b6(3, th_)/b6(4, th_)
Kb = lambda th_, b: K1LY(th_)*b6(1 + b, th_)*b6(5 - b, th_)/(b6(1 - b, th_)*b6(5 + b, th_))
th = mp.mpc('0.71', '0.17')
xg = mp.mpf('1.37') + mp.sqrt(2)/10
rep('L1-13: S_LY = f_{2/3} = -(2)(4)', abs(SLY(th) - fx(mp.mpf(2)/3, th)) < 1e-30)
rep('L1-13: doubling identity (x)(2 theta) = -(x/2)/(6 - x/2) at generic x',
    abs(b6(xg, 2*th) + b6(xg/2, th)/b6(6 - xg/2, th)) < 1e-30)
rep('L1-13: (x)(2 theta) = -(x/2)(x/2 + 6) (same identity, period 12)',
    abs(b6(xg, 2*th) + b6(xg/2, th)*b6(xg/2 + 6, th)) < 1e-30)
rep('L1-13: S_LY(2 theta) = -(1)(2)/((4)(5))', abs(SLY(2*th) + b6(1, th)*b6(2, th)/(b6(4, th)*b6(5, th))) < 1e-30)
err_pr = abs(SLY(2*th) + b6(2, th)*b6(4, th))
rep('L1-13: printed eq-lec1-s2theta S_LY(2 theta) = -(2)(4) is false', err_pr > 1e-2, mp.nstr(err_pr, 6))
rep('L1-13: shift rule (x)(th + i pi a/6)(x)(th - i pi a/6) = (x+a)(x-a)',
    abs(b6(xg, th + I*PI/3)*b6(xg, th - I*PI/3) - b6(xg + 2, th)*b6(xg - 2, th)) < 1e-30)
rep('L1-13: (6) = -1 and (-1)(1) = 1', abs(b6(6, th) + 1) < 1e-30 and abs(b6(-1, th)*b6(1, th) - 1) < 1e-30)
for eps in (1, -1):
    shifted = K1LY(th + I*PI/3, eps)*K1LY(th - I*PI/3, eps)
    rep(f'L1-13: eps={eps}: K(th+i pi/3)K(th-i pi/3) = -eps^2 (3)(5)/(2)',
        abs(shifted + eps**2*b6(3, th)*b6(5, th)/b6(2, th)) < 1e-30)
    rhs = SLY(2*th)*shifted
    rep(f'L1-13: eps={eps}: S(2th) K K = eps^2 (1)(3)/(4); bootstrap {"holds" if eps == 1 else "fails"}',
        abs(rhs - b6(1, th)*b6(3, th)/b6(4, th)) < 1e-30 and ((abs(rhs - K1LY(th, eps)) < 1e-30) == (eps == 1)))

# residues (1/eps scaling)
def res(f, z0):
    vals = [e*f(z0 + e) for e in (mp.mpf('1e-15'), mp.mpf('1e-25'))]
    return vals[1], abs(vals[0] - vals[1])


Rbulk, d1 = res(SLY, 2*I*PI/3)
Gam2 = Rbulk/I
rep('L1-14: Res_{2 pi i/3} S_LY = -2 i sqrt3, so Gamma^2 = -2 sqrt3 (convention Res S = i Gamma^2)',
    abs(Gam2 + 2*mp.sqrt(3)) < 1e-12, f'Gamma^2 = {mp.nstr(Gam2, 12)}')
rep('L1-14: printed Gamma^2 = -2 sqrt3 (2 + sqrt3) is not this', abs(Gam2 + 2*mp.sqrt(3)*(2 + mp.sqrt(3))) > 1)
R2, _ = res(K1LY, I*PI/2)
g2 = R2/(I/2)
rep('L1-14: Res_{i pi/2} K_(1) = -2i(2 sqrt3 - 3) = (i/2) g^2, g^2 = -4(2 sqrt3 - 3)',
    abs(R2 + 2*I*(2*mp.sqrt(3) - 3)) < 1e-12, f'g^2 = {mp.nstr(g2, 12)}')
R6, _ = res(K1LY, I*PI/6)
rep('L1-12: Res_{i pi/6} K_(1) = i(3 - sqrt3)', abs(R6 - I*(3 - mp.sqrt(3))) < 1e-12, mp.nstr(R6, 15))
rep('L1-12: (Res_{i pi/6})^2 = -Gamma^2 g^2/4, consistent with +-(i/2) Gamma g',
    abs(R6**2 + Gam2*g2/4) < 1e-12)
Gam = I*mp.sqrt(2*mp.sqrt(3)); g = 2*I*mp.sqrt(2*mp.sqrt(3) - 3)
rep('L1-12: with g = 2i sqrt(2 sqrt3 - 3) and Gamma = i sqrt(2 sqrt3): Res = -(i/2) Gamma g',
    abs(R6 + I/2*Gam*g) < 1e-12)
# Jacobian: pair amplitude k(t) = K(i pi/2 - t): Res_t k at t0 = - Res_theta K at i pi/2 - t0
kLY = lambda tt_: K1LY(I*PI/2 - tt_)
Rk6, _ = res(kLY, I*PI/3)                 # t0 = i u/2 = i pi/3
Rk0, _ = res(kLY, mp.mpf(0))
rep('L1-12: pair-variable residues: k ~ -(R6)/(t - i pi/3), k ~ -(i/2) g^2/t (GZ (3.45) form)',
    abs(Rk6 + R6) < 1e-12 and abs(Rk0 + I/2*g2) < 1e-12)

# ---------------------------------------------------------------- L1-15 CDD remark and b-family
xx = mp.mpf('1.23')
lhs_ = b6(xx, I*PI/2 - th); rhs_ = b6(6 - xx, I*PI/2 + th)
rep('L1-15: (x) at i pi/2 - theta equals -(6-x) at i pi/2 + theta', abs(lhs_ + rhs_) < 1e-30)
rep('L1-15: ... and not (6-x)^{-1} at i pi/2 + theta (printed parenthesis)', abs(lhs_ - 1/rhs_) > 1e-3)
bv = mp.mpf('0.37') + mp.sqrt(2)/100
Phi = lambda th_: b6(1 + bv, th_)*b6(5 - bv, th_)/(b6(1 - bv, th_)*b6(5 + bv, th_))
rep('L1-15: Phi_b obeys unitarity, CDD crossing and Phi(th) = Phi(th+i pi/3) Phi(th-i pi/3)',
    abs(Phi(th)*Phi(-th) - 1) < 1e-30 and abs(Phi(I*PI/2 - th) - Phi(I*PI/2 + th)) < 1e-30
    and abs(Phi(th) - Phi(th + I*PI/3)*Phi(th - I*PI/3)) < 1e-30)
rep('L1-15: K_b satisfies unitarity, centred crossing and the bootstrap (b generic)',
    abs(Kb(th, bv)*Kb(-th, bv) - 1) < 1e-30
    and abs(Kb(I*PI/2 - th, bv) - SLY(2*th)*Kb(I*PI/2 + th, bv)) < 1e-30
    and abs(Kb(th, bv) - SLY(2*th)*Kb(th + I*PI/3, bv)*Kb(th - I*PI/3, bv)) < 1e-30)
vb = PI*(bv + 1)/6
rep('L1-15: excited factor S(th - i v) S(th + i v) K_b = K_{b+2} = K_{4-b}',
    abs(SLY(th - I*vb)*SLY(th + I*vb)*Kb(th, bv) - Kb(th, bv + 2)) < 1e-30
    and abs(Kb(th, bv + 2) - Kb(th, 4 - bv)) < 1e-30)
rep('L1-17: K_b at b = -2 is GZ second minimal (3)^-1(4)^-1(5)^-1: no pole at i pi/6',
    abs(Kb(th, -2) - 1/(b6(3, th)*b6(4, th)*b6(5, th))) < 1e-30)
rep('L1-17: b = 0 gives K_(1) (equal amplitudes)', abs(Kb(th, 0) - K1LY(th)) < 1e-30)

# ---------------------------------------------------------------- L1-16 f_x in h = 6 blocks
rep('L1-16: f_x = -(6x)(6-6x)', abs(fx(xx/6, th) + b6(xx, th)*b6(6 - xx, th)) < 1e-30)
rep('L1-16: printed f_x = (2x)(6-2x) is false', abs(fx(xx/6, th) - b6(2*xx/6, th)*b6(6 - 2*xx/6, th)) > 1e-3)
rep('L1-16: K_b/K_(1) = f_{(1+b)/6}/f_{(1-b)/6}', abs(Kb(th, bv)/K1LY(th) - fx((1 + bv)/6, th)/fx((1 - bv)/6, th)) < 1e-30)
rep('L1-16: f_x(2 theta) unitary but not CDD-crossing invariant',
    abs(fx(xx/6, 2*th)*fx(xx/6, -2*th) - 1) < 1e-30
    and abs(fx(xx/6, 2*(I*PI/2 - th)) - fx(xx/6, 2*(I*PI/2 + th))) > 1e-3)

# ---------------------------------------------------------------- L1-3 elasticity exercise, illustration
th1 = mp.mpf('0.8')
for r in (1, 2):
    s = 301
    lhs_l = mp.cosh(s*th1)/(mp.exp(s*th1)/2)
    rhs_l = r*mp.cosh(s*th1)/(mp.exp(s*th1)/2) + mp.cosh(s*mp.mpf('0.3'))/(mp.exp(s*th1)/2)
    print(f'     L1-3: r = {r} maximal outgoing rapidities = theta_1: LHS -> {mp.nstr(lhs_l, 6)}, RHS -> {mp.nstr(rhs_l, 6)}')
