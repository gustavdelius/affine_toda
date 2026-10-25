"""Chirality of the non-local charges of sec-nonlocal-charges (ch. 21), in the book's conventions.

Conventions: x^pm = t pm x, d_pm = d_t pm d_x, E + P has Lorentz spin +1, T = (beta/2pi)(phi(+oo) - phi(-oo)),
[phi(x), pi(y)] = i delta(x - y), and e^{i beta phi} has weights h = hbar = beta^2/8pi.

1. Light cone and spin (classical, on a grid). A boost by rapidity vartheta, applied to a wave packet,
   multiplies E + P = (1/2) int (d_- phi)^2 by e^vartheta and E - P = (1/2) int (d_+ phi)^2 by e^-vartheta.
   A conserved density of weight h in x^- (here (d_- phi)^h) gives a charge of spin h - 1, one in x^+ gives
   spin 1 - h. So right-movers phi_R(x^-) carry positive spin and left-movers phi_L(x^+) negative spin.
2. phi_R(x^-) = (phi - int_{-oo}^x pi)/2 and phi_L(x^+) = (phi + int_{-oo}^x pi)/2 up to constants.
3. Topological charge and vacuum phase (lattice canonical commutators): e^{-i a phi_R} and e^{+i a phi_L},
   a = 8pi/beta, raise T by 2; e^{+i a phi_R} lowers it. The phase of J_+ = e^{-i a phi_R} on a vacuum of
   label n (phi = 2 pi n/beta) is q^n, q = e^{-i pi h_J}; the braiding of J_+(x) with a right-moving
   soliton vertex operator at y < x is q^T. For e^{+i a phi_L} both are inverted.
4. Coproduct bookkeeping: with the density at the second soliton seeing the vacuum label left of it,
   Q_+ = int e^{-i a phi_R} has Delta(Q) = Q x 1 + q^T x Q (eq-q-physical); int e^{+i a phi_L}, of the
   same topological charge, would give q^{-T}; int e^{-i a phi_L} q^{-T} gives Delta(f) = f x k^-1 + 1 x f.
5. Conservation for a_n^{(1)}: with V_imag ~ sum_j e^{i beta alpha_j . phi}, the current of node j must be
   e^{-i (4pi/beta) alpha_j^vee . phi_chi} for either chirality (double pole with node j, regular with the
   others); the opposite sign has simple poles at neighbouring nodes for n >= 2. The right-moving one raises
   T by alpha_j^vee (e_j, spin +), the left-moving one lowers it (f_j, spin -).
6. Conservation law: d_+ A + d_- B = 0 gives the conserved charge int (A + B) dx, so d_+ J = d_- H gives
   Q = int (J - H).
7. Boundary sine-Gordon node pairing (sec-bsg-parameters): J_+ = e^{-i a phi_R} and Jbar_- = e^{-i a phi_L},
   which build b_1, have their double poles with e^{i beta phi} and with the Neumann boundary operator
   e^{i beta phi/2}, the C_1 terms; cnj-bsg-uv-ir with eq-bsg-eta-theta makes eps-hat_1 a function of
   mu_+ ~ C_0 alone; charge conjugation exchanges eps-hat_0, eps-hat_1 in eq-bsg-k-matrix and C_0, C_1.
"""
import numpy as np
import sympy as sp


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


rng = np.random.default_rng(7)
beta = np.e / 2                                 # generic coupling (book's beta_SG); not a root of unity
a = 8 * np.pi / beta
hJ = a**2 / (8 * np.pi)
q = np.exp(-1j * np.pi * hJ)

# ------------------------------------------------------------------ 1. boost of chiral wave packets
x = np.linspace(-60, 60, 24001)
dx = x[1] - x[0]


def packet(s):                                  # smooth profile with nonzero h-th moments
    return np.exp(-s**2 / 4) * (1 + 0.3 * np.sin(1.7 * s)) + 0.5 * np.exp(-(s - 2)**2)


def fprime(s, d=1e-5):
    return (packet(s + d) - packet(s - d)) / (2 * d)


def charges(vt, chir):
    """phi = packet(x^-) (chir='R') or packet(x^+) (chir='L'), boosted by rapidity vt, at t' = 0.
    The boosted field is phi'(t', x') = phi(t' cosh vt - x' sinh vt, x' cosh vt - t' sinh vt)."""
    t0 = -x * np.sinh(vt)                       # pre-image of (t' = 0, x')
    x0 = x * np.cosh(vt)
    # derivatives of phi' at t' = 0 by the chain rule
    if chir == 'R':
        g = fprime(t0 - x0)                     # d phi / d x^- at the pre-image
        dt_, dx_ = g * (np.cosh(vt) + np.sinh(vt)), -g * (np.cosh(vt) + np.sinh(vt))
    else:
        g = fprime(t0 + x0)
        dt_, dx_ = g * (np.cosh(vt) - np.sinh(vt)), g * (np.cosh(vt) - np.sinh(vt))
    dm, dp = dt_ - dx_, dt_ + dx_
    E = 0.5 * np.sum(dt_**2 + dx_**2) * dx
    P = -np.sum(dt_ * dx_) * dx
    Qh = {h: np.sum((dm if chir == 'R' else dp)**h) * dx for h in (2, 3, 4)}
    return E, P, Qh


vt = 0.37
ok_EP = ok_spin = True
for chir in ('R', 'L'):
    E0, P0, Q0 = charges(0.0, chir)
    E1, P1, Q1 = charges(vt, chir)
    if chir == 'R':
        ok_EP &= abs((E1 + P1) / (E0 + P0) - np.exp(vt)) < 1e-6 and abs(E0 - P0) < 1e-8 * E0
    else:
        ok_EP &= abs((E1 - P1) / (E0 - P0) - np.exp(-vt)) < 1e-6 and abs(E0 + P0) < 1e-8 * E0
    sgn = 1 if chir == 'R' else -1
    for h in (2, 3, 4):
        ok_spin &= abs(Q1[h] / Q0[h] - np.exp(sgn * (h - 1) * vt)) < 1e-6
report('a boost that raises rapidities multiplies E+P of a right-mover by e^vt and E-P of a left-mover by e^-vt; '
       'a right-mover has E = P, a left-mover E = -P', ok_EP)
report('a density of weight h in x^- gives a charge of spin h-1, one in x^+ spin 1-h (h = 2, 3, 4, numerical boost)',
       ok_spin)

# ------------------------------------------------------------------ 2. chiral halves as functionals
t_, x_ = sp.symbols('t x', real=True)
f, g = sp.Function('f'), sp.Function('g')
phi = f(t_ - x_) + g(t_ + x_)
pi_ = sp.diff(phi, t_)
dR = sp.simplify(sp.diff((phi - sp.Integral(pi_, x_).doit()) / 2, x_) - sp.diff(f(t_ - x_), x_))
dL = sp.simplify(sp.diff((phi + sp.Integral(pi_, x_).doit()) / 2, x_) - sp.diff(g(t_ + x_), x_))
dpm = [sp.simplify(sp.diff(f(t_ - x_), t_) + sp.diff(f(t_ - x_), x_)),
       sp.simplify(sp.diff(g(t_ + x_), t_) - sp.diff(g(t_ + x_), x_))]
report('d_x (phi -+ int^x pi)/2 = d_x phi_{R,L}; d_+ phi_R(x^-) = d_- phi_L(x^+) = 0',
       dR == 0 and dL == 0 and dpm == [0, 0])

# ------------------------------------------------------------------ 3. lattice commutators
N = 400
h_ = 0.05


def lin(phi_c=None, pi_c=None):
    """A linear functional sum_j (c_j phi_j + d_j pi_j) as a pair of coefficient vectors."""
    return (np.zeros(N) if phi_c is None else phi_c, np.zeros(N) if pi_c is None else pi_c)


def comm(A, B):
    """[A, B] for [phi_j, pi_k] = i delta_jk / h_ (lattice delta function)."""
    return 1j * (A[0] @ B[1] - A[1] @ B[0]) / h_


def phi_at(m):
    c = np.zeros(N); c[m] = 1
    return lin(phi_c=c)


def phitilde_at(m):                             # int_{-oo}^{x_m} pi dy
    d = np.zeros(N); d[:m + 1] = h_
    return lin(pi_c=d)


def chiral(m, chir):
    p, pt = phi_at(m), phitilde_at(m)
    s = -1 if chir == 'R' else 1
    return (0.5 * (p[0] + s * pt[0]), 0.5 * (p[1] + s * pt[1]))


Dphi = lin(phi_c=np.eye(N)[N - 1] - np.eye(N)[0])          # phi(+oo) - phi(-oo)
mid = N // 2


def dT(c, chir):
    """Change of T when e^{i c phi_chi(x)} acts: <Dphi>' = <Dphi> - i c [phi_chi, Dphi]."""
    return (beta / (2 * np.pi)) * (-1j * c * comm(chiral(mid, chir), Dphi)).real


ok_T = (abs(dT(-a, 'R') - 2) < 1e-9 and abs(dT(+a, 'L') - 2) < 1e-9
        and abs(dT(+a, 'R') + 2) < 1e-9 and abs(dT(-a, 'L') + 2) < 1e-9)
report('topological charge: e^{-i a phi_R} and e^{+i a phi_L} raise T by 2, e^{+i a phi_R} and e^{-i a phi_L} lower it',
       ok_T)

# local half: c phi_chi = (c/2) phi + ..., so on a vacuum of label n the phase is e^{i (c/2) 2 pi n / beta}
coef_phi_R = chiral(mid, 'R')[0][mid]
coef_phi_L = chiral(mid, 'L')[0][mid]
n_ = 3
phase_R = np.exp(1j * (-a) * coef_phi_R * 2 * np.pi * n_ / beta)
phase_L = np.exp(1j * (+a) * coef_phi_L * 2 * np.pi * n_ / beta)
report('on a vacuum of label n the phase of J_+ = e^{-i a phi_R} is q^n, that of e^{+i a phi_L} is q^-n, q = e^{-i pi h_J}',
       abs(phase_R - q**n_) < 1e-9 and abs(phase_L - q**(-n_)) < 1e-9)


def braid(cJ, chirJ, cP, chirP, mx, my):
    """e^{A} e^{B} = e^{B} e^{A} e^{[A, B]} for A = i cJ phi_chirJ(x), B = i cP phi_chirP(y)."""
    return np.exp(-cJ * cP * comm(chiral(mx, chirJ), chiral(my, chirP)))


ok_br = True
for T in (1, -1):
    ok_br &= abs(braid(-a, 'R', -T * 4 * np.pi / beta, 'R', mid + 20, mid - 20) - q**T) < 1e-9
    ok_br &= abs(braid(+a, 'L', +T * 4 * np.pi / beta, 'L', mid + 20, mid - 20) - q**(-T)) < 1e-9
# the soliton vertex operators carry T = +-1
ok_br &= abs(dT(-4 * np.pi / beta, 'R') - 1) < 1e-9 and abs(dT(4 * np.pi / beta, 'L') - 1) < 1e-9
report('braiding: J_+(x) = e^{-i a phi_R(x)} passes a soliton vertex operator of charge T at y < x with q^T; '
       'e^{+i a phi_L} with q^-T', ok_br)
self_R = braid(-a, 'R', -a, 'R', mid + 20, mid - 20)
report('self-braiding J_+(x) J_+(y) = q^2 J_+(y) J_+(x) for x > y (right-movers)', abs(self_R - q**2) < 1e-9)

# ------------------------------------------------------------------ 4. coproduct bookkeeping
# Two solitons with charges T1, T2 in {+1, -1}; vacuum label n0 at -oo. The density at soliton k sees the
# label to its left through its local half, phase p^{label}; the local action A raises T_k by 2 (A v_- = v_+).
# Left normalization (strip p^{n0}) gives Delta(Q) = Q x 1 + K x Q; read off K.
basis = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
idx = {b: i for i, b in enumerate(basis)}
A1 = np.array([[0, 1], [0, 0]])                 # v_- -> v_+ in basis (v_+, v_-)


def coproduct_K(p, dTsign=+1):
    n0 = 5
    M = np.zeros((4, 4), complex)
    for (T1, T2) in basis:
        # density at soliton 1: label n0 to its left
        if (T1 == -1 and dTsign == 1) or (T1 == 1 and dTsign == -1):
            M[idx[(T1 + 2 * dTsign, T2)], idx[(T1, T2)]] += p**n0
        # density at soliton 2: label n0 + T1 to its left
        if (T2 == -1 and dTsign == 1) or (T2 == 1 and dTsign == -1):
            M[idx[(T1, T2 + 2 * dTsign)], idx[(T1, T2)]] += p**(n0 + T1)
    return M / p**n0


lhs_R = coproduct_K(q)                          # local half e^{-i(a/2)phi}: phase q^n
rhs_R = np.kron(A1, np.eye(2)) + np.kron(np.diag(q**np.array([1., -1.])), A1)
lhs_L = coproduct_K(1 / q)                      # local half e^{+i(a/2)phi}: phase q^-n
rhs_L = np.kron(A1, np.eye(2)) + np.kron(np.diag(q**-np.array([1., -1.])), A1)
report('Q_+ = int e^{-i a phi_R}: Delta(Q_+) = Q_+ x 1 + q^T x Q_+ with q = e^{-i pi h_J} (eq-q-physical)',
       np.allclose(lhs_R, rhs_R))
report('int e^{+i a phi_L} (same T) would give q^{-T} x Q: the old labels contradict eq-q-physical',
       np.allclose(lhs_L, rhs_L) and not np.allclose(lhs_L, rhs_R))
# Qbar_- = int e^{-i a phi_L} lowers T; local half e^{-i(a/2)phi}: phase q^n; f_1 = Qbar_- q^{-T}
Qb = coproduct_K(q, dTsign=-1)
kinv = np.diag(q**-np.array([1., -1.]))
F = Qb @ np.kron(kinv, kinv)
Fsingle = A1.T @ kinv
report('f_1 = Qbar_- q^{-T}, Qbar_- = int e^{-i a phi_L}: Delta(f_1) = f_1 x k_1^-1 + 1 x f_1 with k_1 = q^T',
       np.allclose(F, np.kron(Fsingle, kinv) + np.kron(np.eye(2), Fsingle)))

# ------------------------------------------------------------------ 5. conservation for a_n^(1)
def an_roots(n):
    e = np.eye(n + 1)
    simple = [e[i] - e[i + 1] for i in range(n)]
    return [-(e[0] - e[n])] + simple             # alpha_0 = -theta, alpha_1..alpha_n


sR, sL = dT(1.0, 'R'), dT(1.0, 'L')               # T shift per unit coefficient, from the lattice
ok_cons = ok_neg = ok_Tsign = True
for n in (1, 2, 3, 4):
    R = an_roots(n)
    for j, aj in enumerate(R):
        av = 2 * aj / (aj @ aj)
        for s in (-1, +1):
            # e^{i s (4pi/beta) av . phi_chi} times e^{i beta alpha_k . phi}: exponent s av . alpha_k
            ex = np.array([s * av @ ak for ak in R])
            if s == -1:
                ok_cons &= abs(ex[j] + 2) < 1e-12 and all(e_ > -1e-12 for k, e_ in enumerate(ex) if k != j)
            elif n >= 2:
                ok_neg &= any(abs(e_ + 1) < 1e-12 for e_ in ex)
        # each component of the vector field is an independent copy of the lattice field
        cvec = -(4 * np.pi / beta) * av
        ok_Tsign &= np.allclose(sR * cvec, av) and np.allclose(sL * cvec, -av)
report('a_n^(1), n = 1..4: e^{-i (4pi/beta) alpha_j^vee . phi_chi} has a double pole with node j and none with the others',
       ok_cons)
report('a_n^(1), n = 2..4: the opposite sign has simple poles with the neighbouring nodes (not conserved)', ok_neg)
report('the right-moving node-j current raises T by alpha_j^vee (e_j, spin +), the left-moving one lowers it (f_j, spin -)',
       ok_Tsign)

# ------------------------------------------------------------------ 6. conservation law
A_, B_ = sp.Function('A')(t_, x_), sp.Function('B')(t_, x_)
dp_ = lambda F: sp.diff(F, t_) + sp.diff(F, x_)
dm_ = lambda F: sp.diff(F, t_) - sp.diff(F, x_)
ident = sp.simplify(sp.diff(A_ + B_, t_) + sp.diff(A_ - B_, x_) - (dp_(A_) + dm_(B_)))
report('d_t (A + B) + d_x (A - B) = d_+ A + d_- B, so d_+ J = d_- H conserves int (J - H) dx', ident == 0)

# ------------------------------------------------------------------ 7. boundary sine-Gordon pairing
# OPE exponent of e^{i c1 phi_chi} with e^{i c2 phi_chi}: c1 c2 / 4pi (weights c^2/8pi). The bulk e^{i b beta phi}
# contains e^{i b beta phi_chi}. At a Neumann boundary d_x phi = phi_L' - phi_R' = 0, so phi_L = phi_R + const at
# x = 0 and the boundary operator e^{i b beta phi(0)/2} is e^{i b beta phi_chi(0)} for either chirality.
tt = sp.symbols('tau', real=True)
neu = sp.simplify(sp.diff(f(t_ - x_) + g(t_ + x_), x_).subs(x_, 0)
                  - (sp.diff(g(tt), tt) - sp.diff(f(tt), tt)).subs(tt, t_))
cB_bulk, cB_bdy = beta, 2 * (beta / 2)          # chiral exponents of e^{i beta phi} and e^{i beta phi(0)/2}
ok_b = neu == 0
for cJ in (-a,):                                # J_+ = e^{-i a phi_R}, Jbar_- = e^{-i a phi_L}: same exponent
    for b, want in ((+1, -2), (-1, +2)):
        ok_b &= abs(cJ * b * cB_bulk / (4 * np.pi) - want) < 1e-12           # bulk e^{+- i beta phi}
        ok_b &= abs(cJ * b * cB_bdy / (4 * np.pi) - want) < 1e-12            # boundary e^{+- i beta phi/2}
report('b_1 is built from J_+ = e^{-i a phi_R} and Jbar_- = e^{-i a phi_L}: double poles with e^{i beta phi} and the '
       'boundary term e^{i beta phi/2} (coefficient C_1), so the charges pair node 1 with C_1', ok_b)

lam = 8 * np.pi / beta**2 - 1


def eps_hat(eta, vt):
    return (-np.cos(eta - 1j * vt) / np.sin(np.pi * lam), -np.cos(eta + 1j * vt) / np.sin(np.pi * lam))


def eta_vt(mup, mum):
    """cos((eta -+ i vt)/(lam + 1)) = mu_+-  (cnj-bsg-uv-ir)."""
    zp, zm = np.arccos(mup + 0j), np.arccos(mum + 0j)
    return (lam + 1) * (zp + zm) / 2, (lam + 1) * (zm - zp) / (2j)


mup = 0.4 + 0.3j
e1a, _ = eps_hat(*eta_vt(mup, 0.4 - 0.3j))
e1b, _ = eps_hat(*eta_vt(mup, 0.7 + 0.1j))
report('cnj-bsg-uv-ir with eq-bsg-eta-theta: eps-hat_1 depends on mu_+ ~ C_0 alone (the crossed pairing)',
       abs(e1a - e1b) < 1e-9)


def Kc(y, e0, e1, qq):
    return np.array([[(qq - 1 / qq) * (e1 * y + e0 / y), y**2 - y**-2],
                     [y**2 - y**-2, (qq - 1 / qq) * (e0 * y + e1 / y)]])


sx = np.array([[0, 1], [1, 0]])
y0, e0_, e1_ = np.exp(0.3 + 0.2j), 0.7 - 0.2j, -0.4 + 0.9j
qph = -np.exp(-1j * np.pi * lam)
C0 = lambda phi0: np.exp(1j * beta * phi0 / 2)
C1 = lambda phi0: np.exp(-1j * beta * phi0 / 2)
report('charge conjugation s <-> sbar exchanges eps-hat_0 and eps-hat_1 in eq-bsg-k-matrix, and phi_0 -> -phi_0 '
       'exchanges C_0 and C_1 in eq-bsg-c-dictionary',
       np.allclose(sx @ Kc(y0, e0_, e1_, qph) @ sx, Kc(y0, e1_, e0_, qph)) and abs(C0(-0.3) - C1(0.3)) < 1e-15)
