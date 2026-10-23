"""The massive N = 2 family sl(m|m)^(1) of Evans-Hollowood (1992) and its
classical solitons (sec-super-n2).

Conventions of Evans-Hollowood 1992: x^+- = (t +- x)/2, d_+- = d_t +- d_x,
L = (i/2) D_+Phi.D_-Phi - beta^{-2} sum_j mu_j exp(beta alpha_j.Phi), and the
root parametrization (5.1)

    alpha_{2j-1}.Phi = Phi^(1)_j - i Phi^(2)_j,   alpha_{2j-2}.Phi = i Phi^(2)_{j-1} - Phi^(1)_j,

j = 1..m, indices mod m, with traceless Phi^(1), Phi^(2) and the bilinear
(not hermitian) Euclidean inner product.

Checks (m = 3 and m = 4):
 1. The 2m roots are null and form a cycle with alternating products +-1; the
    Cartan matrix has rank 2m - 2, with the two null vectors (1,0,1,0,...),
    (0,1,0,1,...) (Evans-Hollowood (2.8)).
 2. At psi = 0 and mu_j = mu the bosonic potential is V_real(phi^(1)) + V_imag(phi^(2))
    of eq-imaginary-potential for a_{m-1}^(1), with the book's m equal to mu and the
    same beta (the constants cancel).
 3. The one-soliton (5.3) solves the complex Toda equation (5.2).
 4. Its energy, integrated numerically with the Lagrangian of check 2, is
    M_a = 4 m mu sin(pi a/m)/beta^2 (5.4) = 2 h m_a/beta^2 of sec-soliton-masses.
 5. The central-charge formula (5.9), evaluated from the asymptotic values of the
    fields, gives the same M_a.
 6. The soliton satisfies a first-order equation of the form (5.7),
    -(1 - sigma J) d_+phi + lambda (1 + sigma J) F = 0, F = -(mu/beta) sum_j alpha_j e^{beta alpha_j.phi},
    taken componentwise with u_j = alpha_j.phi and J alpha_j = i(-1)^j alpha_j (Evans-Hollowood
    (2.9)), with sigma = tan(pi a/2m + pi/4) as in their (5.8) but lambda = -sqrt((1-v)/(1+v)).
    With lambda = +sqrt((1-v)/(1+v)) it fails for either sign of J; the sign of lambda depends
    on the sign conventions of the supersymmetry variations (5.5), which are not checked here.
    For a = m/2, sigma is infinite and the check is skipped.
"""
import itertools
import os
import sys
import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from superspace import components, potential   # noqa: E402

FAILS = []


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        FAILS.append(name)


def roots_slmm(m):
    """alpha_0 .. alpha_{2m-1} as coefficient vectors on (phi^(1)_1..m, phi^(2)_1..m)."""
    al = []
    for k in range(2 * m):
        v = [0] * (2 * m)
        if k % 2 == 1:                      # k = 2j-1
            j = (k + 1) // 2
            v[(j - 1) % m] += 1
            v[m + (j - 1) % m] += -sp.I
        else:                               # k = 2j-2
            j = k // 2 + 1
            v[m + (j - 2) % m] += sp.I
            v[(j - 1) % m] += -1
        al.append(v)
    return al


for m in (3, 4):
    al = roots_slmm(m)
    C = sp.Matrix(2 * m, 2 * m, lambda i, j: sum(a * b for a, b in zip(al[i], al[j])))
    cyc = all(C[i, i] == 0 for i in range(2 * m)) and all(
        abs(C[i, (i + 1) % (2 * m)]) == 1 for i in range(2 * m))
    alt = all(C[i, (i + 1) % (2 * m)] == -C[(i + 1) % (2 * m), (i + 2) % (2 * m)] for i in range(2 * m))
    nulls = [sp.Matrix([1 - k % 2 for k in range(2 * m)]), sp.Matrix([k % 2 for k in range(2 * m)])]
    report(f'm = {m}: null roots on a cycle with alternating products +-1, rank {C.rank()} = 2m-2, '
           f'null vectors (1,0,1,0,..), (0,1,0,1,..)',
           cyc and alt and C.rank() == 2 * m - 2 and all(C * n == sp.zeros(2 * m, 1) for n in nulls))

# ======================================================= 2. bosonic potential
for m in (3,):
    beta, mu = sp.symbols('beta mu', positive=True)
    al = roots_slmm(m)
    Gram = sp.eye(2 * m)
    rts = [tuple(beta * x for x in a) for a in al]
    # Evans-Hollowood sign: -beta^-2 mu exp(...)  <=>  mu_j = -mu/beta^2 in the convention of components.py
    J, Lc, _ = components(Gram, rts, [True] * (2 * m), [-mu / beta**2] * (2 * m), name='n')
    U = potential(J, Lc, 2 * m)
    p1 = [J.phi(j) for j in range(m)]
    p2 = [J.phi(m + j) for j in range(m)]
    Vr = (mu**2 / beta**2) * sum(sp.exp(beta * (p1[j] - p1[(j + 1) % m])) - 1 for j in range(m))
    Vi = -(mu**2 / beta**2) * sum(sp.exp(sp.I * beta * (p2[(j - 1) % m] - p2[j])) - 1 for j in range(m))
    report(f'm = {m}: at psi = 0 the potential is V_real(phi^(1)) + V_imag(phi^(2)) of a_{m - 1}^(1) '
           '(eq-imaginary-potential) with the book m equal to mu',
           sp.simplify(sp.expand(U - Vr - Vi)) == 0)

# ============================================================== 3.-6. solitons
mp.mp.dps = 30


def soliton(m, a, beta, mu, v, xi, x, t):
    rho = 2 * mu * mp.sin(mp.pi * a / m) / mp.sqrt(1 - v ** 2)
    w = mp.exp(2j * mp.pi * a / m)
    E = mp.exp(rho * (x - v * t) + xi)
    # separate logarithms keep phi continuous along the soliton
    return [-1 / (1j * beta) * (mp.log(1 + E * w ** j) - mp.log(1 + E * w ** (j - 1))) for j in range(1, m + 1)]


def phi_d(m, a, beta, mu, v, xi, x, t, kt, kx):
    f = lambda tt, xx, j: soliton(m, a, beta, mu, v, xi, xx, tt)[j]
    return [mp.diff(lambda tt, xx: f(tt, xx, j), (t, x), (kt, kx)) for j in range(m)]


ok3 = ok4 = ok5 = True
ok6, ok6bad = True, True
for m, a in ((3, 1), (3, 2), (4, 1), (4, 2)):
    # a complex xi keeps 1 + E w^j away from zero (for real xi the solitons with
    # w^j = -1 are singular)
    beta, mu, v, xi = mp.mpf('0.8'), mp.mpf('1.3'), mp.mpf('0.35'), mp.mpc('0.2', '0.3')
    for (x, t) in ((mp.mpf('0.31'), mp.mpf('0.17')), (mp.mpf('-0.6'), mp.mpf('0.4'))):
        ph = soliton(m, a, beta, mu, v, xi, x, t)
        tt = phi_d(m, a, beta, mu, v, xi, x, t, 2, 0)
        xx = phi_d(m, a, beta, mu, v, xi, x, t, 0, 2)
        for j in range(m):
            lhs = (tt[j] - xx[j]) + mu ** 2 / (1j * beta) * (
                mp.exp(1j * beta * (ph[j] - ph[(j + 1) % m])) - mp.exp(1j * beta * (ph[(j - 1) % m] - ph[j])))
            ok3 &= abs(lhs) < 1e-15
    # 4. energy at rest
    v0 = mp.mpf(0)

    def edens(x):
        ph = soliton(m, a, beta, mu, v0, xi, x, 0)
        dx = [mp.diff(lambda xx: soliton(m, a, beta, mu, v0, xi, xx, 0)[j], x) for j in range(m)]
        V = -(mu ** 2 / beta ** 2) * sum(mp.exp(1j * beta * (ph[(j - 1) % m] - ph[j])) - 1 for j in range(m))
        return sum(d ** 2 for d in dx) / 2 + V
    mp.mp.dps = 20
    E = mp.quad(edens, [-25, -8, -2, 0, 2, 8, 25])
    mp.mp.dps = 30
    Ma = 4 * m * mu * mp.sin(mp.pi * a / m) / beta ** 2
    ok4 &= abs(E - Ma) < 1e-8
    # 5. central-charge formula (5.9) from the asymptotic values
    sig = mp.tan(mp.pi * a / (2 * m) + mp.pi / 4)
    al = roots_slmm(m)

    def u_vals(x):
        ph = soliton(m, a, beta, mu, v0, xi, x, 0)
        fld = [0] * m + ph                  # phi^(1) = 0
        return [sum(complex(c) * f for c, f in zip(al[k], fld)) for k in range(2 * m)]

    def boundary(x):
        u = u_vals(x)
        ev = sum(mp.exp(beta * u[k]) for k in range(0, 2 * m, 2))
        od = sum(mp.exp(beta * u[k]) for k in range(1, 2 * m, 2))
        return mu / beta ** 2 * ((1 + 1j * sig) / (1 - 1j * sig) * ev + (1 - 1j * sig) / (1 + 1j * sig) * od)
    M59 = boundary(mp.mpf(40)) - boundary(mp.mpf(-40))
    ok5 &= abs(M59 - Ma) < 1e-10
    # 6. first-order equation (5.7), componentwise
    lam = mp.sqrt((1 - v) / (1 + v))
    Cm = [[sum(complex(p) * complex(q) for p, q in zip(al[i], al[j])) for j in range(2 * m)] for i in range(2 * m)]
    if 2 * a == m:
        continue
    for sgn, store, lamv in ((+1, 'good', -lam), (+1, 'bad', lam), (-1, 'bad', lam)):
        lj = [sgn * 1j * (-1) ** k for k in range(2 * m)]
        res = 0
        for (x, t) in ((mp.mpf('0.31'), mp.mpf('0.17')), (mp.mpf('-0.6'), mp.mpf('0.4'))):
            ph = soliton(m, a, beta, mu, v, xi, x, t)
            dt = phi_d(m, a, beta, mu, v, xi, x, t, 1, 0)
            dx = phi_d(m, a, beta, mu, v, xi, x, t, 0, 1)
            fld = [0] * m + ph
            dpl = [0] * m + [dt[j] + dx[j] for j in range(m)]        # d_+ = d_t + d_x
            u = [sum(complex(c) * f for c, f in zip(al[k], fld)) for k in range(2 * m)]
            du = [sum(complex(c) * f for c, f in zip(al[k], dpl)) for k in range(2 * m)]
            for k in range(2 * m):
                aF = -(mu / beta) * sum(Cm[k][l] * mp.exp(beta * u[l]) for l in range(2 * m))
                res = max(res, abs(-(1 + sig * lj[k]) * du[k] + lamv * (1 - sig * lj[k]) * aF))
        if store == 'good':
            ok6 &= res < 1e-15
        else:
            ok6bad &= res > 1e-3

report('the one-soliton (5.3) solves the complex Toda equation (5.2) (m = 3, 4; a = 1, 2)', ok3)
report('its energy at rest is M_a = 4 m mu sin(pi a/m)/beta^2 (5.4) = 2 h m_a/beta^2 with m_a = 2 mu sin(pi a/h)',
       ok4)
report('the central-charge formula (5.9) from the asymptotic values gives the same M_a', ok5)
report('for a != m/2 the soliton satisfies the first-order equation (5.7) with J alpha_j = i(-1)^j alpha_j, '
       'sigma = tan(pi a/2m + pi/4) and lambda = -sqrt((1-v)/(1+v)); with lambda = +sqrt((1-v)/(1+v)) '
       'it fails for either sign of J', ok6 and ok6bad)

sys.exit(1 if FAILS else 0)
