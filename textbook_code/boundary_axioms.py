"""Axioms of boundary scattering (Part VII, ch. 28: sec-reflection-equation, sec-boundary-unitarity-crossing,
sec-boundary-bound-states, sec-boundary-lee-yang).

Checks:
 1. Reflection equation, sine-Gordon K of boundary_common (book conventions, physical q): the braid (R-check)
    form (eq-reflection-braid) and the S/K form (eq-reflection-sk) with S_12 = R-check P hold, and their
    two sides agree term by term; a random constant K fails.
 2. Diagonal boundary crossing-unitarity (eq-boundary-crossing-diagonal): the book form
    K_a(i pi/2 - theta) = S_{a abar}(2 theta) K_abar(i pi/2 + theta), the form
    K_a(theta + i pi/2) K_abar(theta - i pi/2) S_{a abar}(2 theta) = 1 (Corrigan-Delius (2.7)), and
    Delius-Gandenberger (3.17) K_a(theta) K_abar(theta + i pi) = S_aa(2 theta) are equivalent: all three hold
    for the sinh-Gordon Dirichlet factor and for a CDD-modified factor, and all three fail together for a
    factor violating crossing-unitarity.
 3. Boundary scaling Lee-Yang (eq-ly-reflection, Dorey-Pocklington-Tateo-Watts hep-th/9712197 (2.8)-(2.9)):
    R_(1) and R_b satisfy unitarity, crossing-unitarity with S = f_{2/3}, and the boundary bootstrap
    R(theta) = S(2 theta) R(theta + i pi/3) R(theta - i pi/3); with the overall sign of Ghoshal-Zamolodchikov
    (3.51) R_(1) violates the bootstrap; R_{b=0} = R_(1), R_{b=-2} = R_(2) (GZ (3.52)); the residues of R_(1)
    at i pi/2 and i pi/6 agree with GZ (3.41), (3.45) and g = 2i sqrt(2 sqrt3 - 3) (GZ (3.54)).
 4. Boundary bound states of R_b: pole at theta = i(b+1)pi/6, boundary bootstrap
    R'(theta) = R_b(theta) S(theta - iv) S(theta + iv) gives R' = R_{b+2}; energy gap M cos((b+1)pi/6).
 5. Boundary Bethe-Yang (eq-boundary-bethe-yang): free Dirichlet boson (S = 1, R = -1) gives p L = pi k.
Runtime: a few seconds.
"""
import numpy as np
import mpmath as mp
from boundary_common import report, q_phys, sl2_spin, rcheck, k_intertwiner, prop_residual, one_k

mp.mp.dps = 30

# ---------------------------------------------------------------- 1. reflection equation: braid form vs S/K form
LAM = 1.37
q = q_phys(LAM)
e0, e1 = 0.43 - 0.2j, -0.71 + 0.35j


def Rch(z):
    return rcheck(sl2_spin(1, q, z), sl2_spin(1, q, 1.0))[0]


def Kf(y):
    null, _ = k_intertwiner(sl2_spin(1, q, y), sl2_spin(1, q, 1/y), [1/q, 1/q], [e0, e1])
    return null[0]


P = np.zeros((4, 4))
for a in range(2):
    for b in range(2):
        P[2*b + a, 2*a + b] = 1
y1, y2 = np.exp(LAM*(0.31 + 0.2j)), np.exp(LAM*(-0.45 + 0.1j))
K1, K2 = Kf(y1), Kf(y2)
Kslot1 = lambda K: np.kron(K, np.eye(2))
Lb = Rch(y1/y2) @ one_k(K1, 2) @ Rch(y1*y2) @ one_k(K2, 2)
Rb = one_k(K2, 2) @ Rch(y1*y2) @ one_k(K1, 2) @ Rch(y1/y2)
S12 = lambda z: Rch(z) @ P            # S-matrix on fixed tensor slots
S21 = lambda z: P @ Rch(z)            # = P S12 P
Ls = S12(y1/y2) @ Kslot1(K1) @ S21(y1*y2) @ one_k(K2, 2)
Rs = one_k(K2, 2) @ S12(y1*y2) @ Kslot1(K1) @ S21(y1/y2)
rb, _ = prop_residual(Lb, Rb)
rs, _ = prop_residual(Ls, Rs)
ident = np.abs(Lb - Ls).max() + np.abs(Rb - Rs).max()
report(f'reflection equation: braid form [{rb:.0e}] and S/K form with S_12 = R-check P [{rs:.0e}] hold; '
       f'the two sides agree term by term [{ident:.0e}]', rb < 1e-12 and rs < 1e-12 and ident < 1e-12)
rng = np.random.default_rng(1)
Kr = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
Lr = Rch(y1/y2) @ one_k(Kr, 2) @ Rch(y1*y2) @ one_k(Kr, 2)
Rr = one_k(Kr, 2) @ Rch(y1*y2) @ one_k(Kr, 2) @ Rch(y1/y2)
report(f'control: a random constant K fails the reflection equation [{prop_residual(Lr, Rr)[0]:.2f}]',
       prop_residual(Lr, Rr)[0] > 1e-3)


# ---------------------------------------------------------------- 2. diagonal crossing-unitarity: three printed forms
def blk(x, h):
    """eq-blocks: (x) = sinh(theta/2 + i pi x/2h)/sinh(theta/2 - i pi x/2h)."""
    return lambda th: mp.sinh(th/2 + 1j*mp.pi*x/(2*h))/mp.sinh(th/2 - 1j*mp.pi*x/(2*h))


Bc = mp.mpf('0.63')                   # sinh-Gordon B
E0 = mp.mpf('1.41')
S_shg = lambda th: -1/(blk(Bc, 2)(th)*blk(2 - Bc, 2)(th))                       # eq-sinh-gordon-smatrix
KD = lambda th: blk(2 - Bc/2, 2)(th)*blk(1 + Bc/2, 2)(th)/blk(1, 2)(th)         # Dirichlet (CD (2.6))
K_E = lambda th: KD(th)/(blk(1 - E0, 2)(th)*blk(1 + E0, 2)(th))                 # with a CDD-type factor
K_bad = lambda th: KD(th)*blk(mp.mpf('0.7'), 2)(th)                             # unitary, but not crossing-unitary


def forms(K):
    pts = [mp.mpc('0.37', '0.11'), mp.mpc('-0.8', '0.3'), mp.mpc('1.3', '-0.2')]
    ip2 = 1j*mp.pi/2
    f_book = max(abs(K(ip2 - t) - S_shg(2*t)*K(ip2 + t)) for t in pts)
    f_cd = max(abs(K(t + ip2)*K(t - ip2)*S_shg(2*t) - 1) for t in pts)
    f_dg = max(abs(K(t)*K(t + 1j*mp.pi) - S_shg(2*t)) for t in pts)
    f_unit = max(abs(K(t)*K(-t) - 1) for t in pts)
    return f_book, f_cd, f_dg, f_unit


ok = True
for K in (KD, K_E):
    fb, fc, fd, fu = forms(K)
    ok &= max(fb, fc, fd, fu) < 1e-25
fb, fc, fd, fu = forms(K_bad)
report('diagonal crossing-unitarity: the book form, CD (2.7) and DG (3.17) all hold for the sinh-Gordon K_D and '
       'K_D/((1-E)(1+E)), and all fail for a unitary K_D (0.7) (S self-conjugate, S(theta+2 pi i) = S(theta))',
       ok and min(fb, fc, fd) > 1e-3 and fu < 1e-25)

# ---------------------------------------------------------------- 3. boundary scaling Lee-Yang
br = lambda x: blk(x, 6)              # (x) with h = 6: sinh(theta/2 + i pi x/12)/sinh(theta/2 - i pi x/12)
fx = lambda x: (lambda th: (mp.sinh(th) + 1j*mp.sin(mp.pi*x))/(mp.sinh(th) - 1j*mp.sin(mp.pi*x)))
S_ly = fx(mp.mpf(2)/3)                # eq-lee-yang-smatrix
pts = [mp.mpc('0.37', '0.11'), mp.mpc('-0.8', '0.3'), mp.mpc('1.3', '-0.2')]
report('S_LY = f_{2/3} = -(2)(4) in the blocks (x) = sinh(theta/2 + i pi x/12)/sinh(theta/2 - i pi x/12)',
       max(abs(S_ly(t) + br(2)(t)*br(4)(t)) for t in pts) < 1e-25)

R1 = lambda th: br(1)(th)*br(3)(th)/br(4)(th)                                      # DPTW (2.8)
R2 = lambda th: 1/(br(3)(th)*br(4)(th)*br(5)(th))                                  # DPTW (2.8) = GZ (3.52)
Rb = lambda b: (lambda th: R1(th)*br(1 + b)(th)*br(5 - b)(th)/(br(1 - b)(th)*br(5 + b)(th)))   # DPTW (2.9)
R1_gz = lambda th: -R1(th)                                                         # GZ (3.51) as printed


def axioms(R):
    ip2 = 1j*mp.pi/2; ip3 = 1j*mp.pi/3
    u = max(abs(R(t)*R(-t) - 1) for t in pts)
    c = max(abs(R(ip2 - t) - S_ly(2*t)*R(ip2 + t)) for t in pts)
    bs = max(abs(R(t) - S_ly(2*t)*R(t + ip3)*R(t - ip3)) for t in pts)
    return u, c, bs


ok = True
for R in (R1, R2, Rb(mp.mpf('0.4')), Rb(mp.mpf('-1.3')), Rb(mp.mpf('1.7'))):
    ok &= max(axioms(R)) < 1e-25
report('Lee-Yang: R_(1), R_(2), R_b (b = 0.4, -1.3, 1.7) satisfy unitarity, crossing-unitarity with f_{2/3} and the '
       'boundary bootstrap R(theta) = S(2theta) R(theta + i pi/3) R(theta - i pi/3)', ok)
u, c, bs = axioms(R1_gz)
report(f'GZ (3.51) with its overall minus sign satisfies unitarity and crossing-unitarity but violates the '
       f'bootstrap [{float(bs):.2f}]: the sign must be +', u < 1e-25 and c < 1e-25 and bs > 0.1)
report('R_b at b = 0 is R_(1); at b = -2 it is R_(2) = GZ (3.52)',
       max(abs(Rb(0)(t) - R1(t)) + abs(Rb(-2)(t) - R2(t)) for t in pts) < 1e-25)


def residue(f, z0, eps=mp.mpf('1e-12')):
    return (f(z0 + eps) - f(z0 - eps))*eps/2


g2 = -4*(2*mp.sqrt(3) - 3)                             # g = 2i sqrt(2 sqrt3 - 3), GZ (3.54)
rho = residue(R1, 1j*mp.pi/2)
f2 = residue(S_ly, 2j*mp.pi/3)/1j                      # eq-residue: S ~ i Gamma^2/(theta - i u)
r6 = residue(R1, 1j*mp.pi/6)
# K(theta) = R(i pi/2 - theta); GZ (3.45): K ~ -(i/2) g^2/theta  <=>  Res_{i pi/2} R = (i/2) g^2
# GZ (3.41): K ~ (i/2) f g/(theta - i u/2)          <=>  Res_{i pi/6} R = -(i/2) f g, so (Res)^2 = -f^2 g^2/4
report(f'R_(1): Res at i pi/2 = (i/2) g^2 with g^2 = -4(2 sqrt3 - 3) [{float(abs(rho - 1j*g2/2)):.0e}]; '
       f'Res at i pi/6 squared = -Gamma^2 g^2/4 [{float(abs(r6**2 + f2*g2/4)):.0e}] (GZ (3.41), (3.45), (3.54))',
       abs(rho - 1j*g2/2) < 1e-10 and abs(r6**2 + f2*g2/4) < 1e-10)
report(f'R_(1) with GZ sign: Res at i pi/2 = -(i/2) g^2, so GZ (3.54) belongs to the + sign '
       f'[{float(abs(residue(R1_gz, 1j*mp.pi/2) + 1j*g2/2)):.0e}]', abs(residue(R1_gz, 1j*mp.pi/2) + 1j*g2/2) < 1e-10)

# ---------------------------------------------------------------- 4. boundary bound states of R_b
ok_pole, ok_boot = True, True
for b in (mp.mpf('-0.5'), mp.mpf('0.6'), mp.mpf('1.4')):
    v = (b + 1)*mp.pi/6
    ok_pole &= abs(1/Rb(b)(1j*v + mp.mpf('1e-20'))) < 1e-15
    Rex = lambda th, b=b, v=v: Rb(b)(th)*S_ly(th - 1j*v)*S_ly(th + 1j*v)
    ok_boot &= max(abs(Rex(t) - Rb(b + 2)(t)) for t in pts) < 1e-25
report('R_b has a pole at theta = i(b+1)pi/6 (b = -0.5, 0.6, 1.4); the excited boundary has R_{b+2} = '
       'R_b(theta) S(theta - iv) S(theta + iv)', ok_pole and ok_boot)
b = mp.mpf('0.6'); v = (b + 1)*mp.pi/6
gap = mp.cos(v)
report(f'boundary bound-state energy e_1 - e_0 = M cos v = M cos((b+1)pi/6) = {float(gap):.6f} M at b = 0.6 '
       '(DPTW sec. 6)', abs(gap - mp.cos((b + 1)*mp.pi/6)) < 1e-30 and 0 < gap < 1)
# second pole: i(b-1)pi/6 enters the physical strip at b = 1
report('R_b has a second physical-strip pole at i(b-1)pi/6 for 1 < b < 4 (b = 1.4)',
       abs(1/Rb(mp.mpf('1.4'))(1j*mp.mpf('0.4')*mp.pi/6 + mp.mpf('1e-20'))) < 1e-15)

# ---------------------------------------------------------------- 5. boundary Bethe-Yang, free Dirichlet boson
Lint, k = 3.7, 5
p = np.pi*k/Lint
report('boundary Bethe-Yang e^{2ipL} R_alpha R_beta = 1 with R = -1 (Dirichlet free boson) gives p L = pi k',
       abs(np.exp(2j*p*Lint)*(-1)*(-1) - 1) < 1e-12)
