"""Sine-Gordon on the half-line: soliton reflection matrix in book conventions (Part VII, ch. 29/32).

Part VI conventions (boundary_common.py): coproduct Delta(e)=e(x)1+k(x)e, physical q = -exp(-i pi lambda),
principal gradation y = exp(lambda theta); spin-1/2 basis (soliton, antisoliton).  Coideal generators
b_j = e_j + q^-1 f_j k_j + eps_j k_j.  GZ = Ghoshal-Zamolodchikov hep-th/9306002, BPT = Bajnok-Palla-
Takacs-Toth hep-th/0106070 (their (3.2) parametrization of GZ).

Checks:
 1. K: V(y) -> V(1/y) from the coideal: unique; solves the reflection equation (braid form) with the
    book's R-check; the reflection equation is invariant under q -> -q (S_T -> -S_T) for the same K,
    and K(eps; -q) = K(-eps; q): the sign of q is invisible in the reflection equation.
 2. Dictionary to GZ (5.12): eps_1 = eps e^{i xi}, eps_0 = eps e^{-i xi}, k = 1/(eps sin(pi lambda));
    to BPT (3.2): eps_1 = -cos(eta - i vartheta)/sin(pi lambda), eps_0 = -cos(eta + i vartheta)/sin(pi lambda);
    consistency with GZ (5.25).  Neumann (GZ xi=0, k=1/sin(pi lambda/2); BPT eta=pi(lambda+1)/2) is
    eps_0 = eps_1 = 1/(2cos(pi lambda/2)) = (2-q-q^-1)^(-1/2), the a_n special value; eps = 0 is eta = pi/2.
 3. Boundary crossing-unitarity in book form  k(theta) = S(2 theta) k(-theta),
    k = sum_ab K_{b abar}(i pi/2 - theta) v_a (x) v_b, S = S_0 R-check (eq-sg-s0, eq-sg-smatrix):
    the matrix part reduces it to a scalar equation at the physical q with C = sigma_x; at -q only with
    C' = [[0,1],[-1,0]] (the bulk kappa = -1 conjugation of sec-crossing-sign).
 4. Scalar factor: GZ's a(u) = sin(lambda(pi-u)) rho(u) equals the book S_0; GZ/BPT R_0 and sigma
    satisfy their functional equations; the full BPT K satisfies boundary unitarity and boundary
    crossing-unitarity with the book S_0 (numerically, two parameter points).
 5. Breather reflection factors: the B_1, B_2 amplitudes obtained from the soliton K by the boundary
    bootstrap equal Ghoshal (3.5)-(3.10); B_1 factor unitary and crossing-unitary with S_11 (eq-sg-s11).
 6. Continuation lambda -> -2/B: Ghoshal B_1 factor = Corrigan-Delius (2.5) with E = B eta/pi, F = i B vt/pi;
    Dirichlet and Neumann values.
 7. Excited boundary states (Bajnok-Palla-Takacs-Toth): bootstrap identity for Q on |0> (eta -> etabar),
    poles nu_0, nu_1 of sigma(eta,u), K on |0> (s <-> sbar, eta -> etabar) unitary and crossing-unitary.
 8. Limits of Al. B. Zamolodchikov's UV-IR relation (BPT hep-th/0108157 (3.2)): Neumann, Dirichlet, the eps-hat = 0
    point, and the classical value of M_crit (the C = 1 boundary).
Runtime: about 40 seconds.
"""
import numpy as np
import mpmath as mp
from boundary_common import report, q_phys, sl2_spin, rcheck, k_intertwiner, re_preserving

mp.mp.dps = 20
LAM = 1.37                         # lambda = 8 pi/beta_SG^2 - 1 (attractive regime, two breathers)
q = q_phys(LAM)
C = np.array([[0, 1], [1, 0]])     # charge conjugation s <-> sbar


def Rsg(z, qq=q):
    return rcheck(sl2_spin(1, qq, z), sl2_spin(1, qq, 1.0))[0]


def K_coideal(y, e0, e1, qq=q):
    null, _ = k_intertwiner(sl2_spin(1, qq, y), sl2_spin(1, qq, 1/y), [1/qq, 1/qq], [e0, e1])
    assert len(null) == 1
    return null[0]


def K_closed(y, e0, e1, qq=q):
    d = y**2 - y**-2
    return np.array([[(qq - 1/qq)*(e1*y + e0/y), d], [d, (qq - 1/qq)*(e0*y + e1/y)]])/d


def prop(A, B):
    c = np.vdot(B.ravel(), A.ravel())/np.vdot(B.ravel(), B.ravel())
    return np.linalg.norm(A - c*B)/np.linalg.norm(A)


# ---------------------------------------------------------------- 1. K from the coideal, reflection equation
e0, e1 = 0.43 - 0.2j, -0.71 + 0.35j
report('K(theta) from the coideal equals the closed form (unique)',
       max(prop(K_coideal(y, e0, e1), K_closed(y, e0, e1)) for y in (0.8+0.3j, 1.9)) < 1e-12)
r, c = re_preserving(Rsg, lambda y: K_coideal(y, e0, e1), 1.4 + 0.3j, 0.6 - 0.2j)
report(f'reflection equation (braid form) at the physical q [{r:.0e}]', r < 1e-10)
r2, _ = re_preserving(lambda z: Rsg(z, -q), lambda y: K_coideal(y, e0, e1), 1.4 + 0.3j, 0.6 - 0.2j)
report(f'the same K also solves it with the R-check at -q [{r2:.0e}]: the RE does not see the sign of q', r2 < 1e-10)
report('K(eps_0, eps_1; -q) = K(-eps_0, -eps_1; q)',
       prop(K_coideal(0.8+0.3j, e0, e1, -q), K_coideal(0.8+0.3j, -e0, -e1, q)) < 1e-12)
R_q, R_mq = Rsg(1.7+0.2j), Rsg(1.7+0.2j, -q)
Pi = np.diag([1, -1])
report('R-check(-q) = (Pi (x) 1) R-check(q) (Pi (x) 1), Pi = diag(1,-1): S_T -> -S_T',
       np.allclose(R_mq, np.kron(Pi, np.eye(2)) @ R_q @ np.kron(Pi, np.eye(2))))

# ---------------------------------------------------------------- 2. dictionary to GZ and BPT
ok = True
for xi, kk in ((0.4, 0.7), (1.1 + 0.2j, 2.3)):
    eps = 1/(kk*np.sin(np.pi*LAM))
    for th in (0.37, -0.81 + 0.1j):
        u = -1j*th; y = np.exp(LAM*th)
        GZ = np.array([[np.cos(xi + LAM*u), kk/2*np.sin(2*LAM*u)], [kk/2*np.sin(2*LAM*u), np.cos(xi - LAM*u)]])
        ok &= prop(K_closed(y, eps*np.exp(-1j*xi), eps*np.exp(1j*xi)), GZ) < 1e-12
report('GZ (5.12): K = R(u)[[cos(xi+lambda u), k/2 sin 2lambda u],[.., cos(xi-lambda u)]] with '
       'eps_1 = eps e^{i xi}, eps_0 = eps e^{-i xi}, k = 1/(eps sin pi lambda)', ok)
ok = True
for eta, vt in ((0.6, 0.3), (2.0, 1.1)):
    E1 = -np.cos(eta - 1j*vt)/np.sin(np.pi*LAM); E0 = -np.cos(eta + 1j*vt)/np.sin(np.pi*LAM)
    for th in (0.37, -0.81 + 0.1j):
        u = -1j*th; y = np.exp(LAM*th)
        Pp = np.cos(LAM*u)*np.cos(eta)*np.cosh(vt) - np.sin(LAM*u)*np.sin(eta)*np.sinh(vt)
        Pm = np.cos(LAM*u)*np.cos(eta)*np.cosh(vt) + np.sin(LAM*u)*np.sin(eta)*np.sinh(vt)
        Q0 = -np.sin(LAM*u)*np.cos(LAM*u)
        ok &= prop(K_closed(y, E0, E1), np.array([[Pp, Q0], [Q0, Pm]])) < 1e-12
    # GZ (5.25): cos eta cosh vt = -cos(xi)/k, cos^2 eta + cosh^2 vt = 1 + 1/k^2
    eps = np.sqrt(E0*E1); xi = -1j*np.log(E1/eps); kk = 1/(eps*np.sin(np.pi*LAM))
    ok &= abs(np.cos(eta)*np.cosh(vt) + np.cos(xi)/kk) < 1e-12 and abs(np.cos(eta)**2 + np.cosh(vt)**2 - 1 - 1/kk**2) < 1e-12
report('BPT (3.2): eps_1 = -cos(eta - i vt)/sin(pi lambda), eps_0 = -cos(eta + i vt)/sin(pi lambda); '
       'consistent with GZ (5.25)', ok)
eps_N = 1/(2*np.cos(np.pi*LAM/2))
report('Neumann (GZ xi=0, k=1/sin(pi lambda/2); BPT eta=pi(lambda+1)/2, vt=0): eps_0 = eps_1 = 1/(2cos(pi lambda/2)) '
       '= (2-q-1/q)^(-1/2)',
       abs(1/((1/np.sin(np.pi*LAM/2))*np.sin(np.pi*LAM)) - eps_N) < 1e-14 and
       abs(-np.cos(np.pi*(LAM+1)/2)/np.sin(np.pi*LAM) - eps_N) < 1e-14 and abs(eps_N**2*(2 - q - 1/q) - 1) < 1e-12)
eps_bpt = lambda eta, vt: (-np.cos(eta + 1j*vt)/np.sin(np.pi*LAM), -np.cos(eta - 1j*vt)/np.sin(np.pi*LAM))
report('dictionary: eps_0 = eps_1 = 0 <-> BPT (eta, vt) = (pi/2, 0); BPT (0, 0) <-> eps_0 = eps_1 = -1/sin(pi lambda)',
       np.allclose(eps_bpt(np.pi/2, 0), 0) and np.allclose(eps_bpt(0, 0), -1/np.sin(np.pi*LAM)))


# ---------------------------------------------------------------- 3./4. scalar factor and crossing-unitarity
def lg(z):
    return mp.loggamma(z)


L = mp.mpf(LAM)


def log_rho(u):
    """GZ (5.7)."""
    u = mp.mpmathify(u)
    def logF(l, uu):
        v = L*uu/mp.pi
        return lg(2*l*L - v) + lg(1 + 2*l*L - v) - lg((2*l+1)*L - v) - lg(1 + (2*l-1)*L - v)
    s = mp.nsum(lambda l: logF(l, u) + logF(l, mp.pi - u) - logF(l, 0) - logF(l, mp.pi), [1, mp.inf])
    return mp.log(-1/mp.pi) + lg(L) + lg(1 - L*u/mp.pi) + lg(1 - L + L*u/mp.pi) + s


def a_gz(theta):
    u = -1j*mp.mpmathify(theta)
    return mp.sin(L*(mp.pi - u))*mp.exp(log_rho(u))


def S0_book(theta):
    """eq-sg-s0, xi_SG = pi/lambda, real theta."""
    xs = mp.pi/L
    f = lambda t: mp.sinh((mp.pi - xs)*t/2)*mp.sin(theta*t)/(t*mp.sinh(xs*t/2)*mp.cosh(mp.pi*t/2))
    return -mp.exp(-1j*mp.quad(f, [0, 1, 5, 20, mp.inf]))


def log_R0(u):
    """BPT (3.2) R_0 (= GZ (5.21))."""
    u = mp.mpmathify(u)
    def t(l, uu):
        v = 2*L*uu/mp.pi
        return lg(4*l*L - v) + lg(4*L*(l-1) + 1 - v) - lg((4*l-3)*L - v) - lg((4*l-1)*L + 1 - v)
    return mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


def log_sigma(x, u):
    """BPT form of GZ's sigma: sigma(x,u) sigma(x,-u) = cos^2 x/(cos(x+lambda u)cos(x-lambda u))."""
    x, u = mp.mpmathify(x), mp.mpmathify(u)
    def t(l, uu):
        v = L*uu/mp.pi; a = x/mp.pi
        return (lg(0.5 + a + (2*l-1)*L - v) + lg(0.5 - a + (2*l-1)*L - v)
                - lg(0.5 - a + (2*l-2)*L - v) - lg(0.5 + a + 2*l*L - v))
    return mp.log(mp.cos(x)/mp.cos(x + L*u)) + mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


# GZ (5.23) as printed: the normalizer Pi^2(x, pi/2) Pi^2(-x, pi/2) diverges, log Pi(x, pi/2) ~ (lambda/4) log(l_max)
def logPi_partial(x, u, lmax):
    return sum(lg(0.5 + (2*l+0.5)*L + x/mp.pi - L*u/mp.pi) + lg(0.5 + (2*l+1.5)*L + x/mp.pi)
               - lg(0.5 + (2*l+1.5)*L + x/mp.pi - L*u/mp.pi) - lg(0.5 + (2*l+0.5)*L + x/mp.pi) for l in range(lmax))
d1 = mp.re(logPi_partial(0.7, mp.pi/2, 400) - logPi_partial(0.7, mp.pi/2, 100))
report(f'GZ (5.23) normalizer diverges: log Pi(x,pi/2) grows by (lambda/4) log 4 per factor 4 in l_max '
       f'[{float(d1):.4f} vs {float(L/4*mp.log(4)):.4f}]; BPT form used instead', abs(d1 - L/4*mp.log(4)) < 1e-3)
th_list = [0.23, 0.61]
dev = max(abs(a_gz(th) - S0_book(th)) for th in th_list)
report(f"GZ a(theta) = sin(lambda(pi-u)) rho(u) equals the book S_0 (eq-sg-s0) [{float(dev):.0e}]", dev < 1e-10)
x0 = mp.mpf('0.7')
dev1 = max(abs(mp.exp(log_sigma(x0, u) + log_sigma(x0, -u))*mp.cos(x0 + L*u)*mp.cos(x0 - L*u)/mp.cos(x0)**2 - 1)
           for u in (0.3j, 0.2 + 0.1j))
dev2 = max(abs(mp.exp(log_sigma(x0, mp.pi/2 - u) - log_sigma(x0, mp.pi/2 + u)) - 1) for u in (0.3j, 0.2 + 0.1j))
report(f'sigma(x,u)sigma(x,-u) = cos^2x/(cos(x+lambda u)cos(x-lambda u)) and sigma(x,pi/2-u) = sigma(x,pi/2+u) '
       f'(GZ (5.24), BPT normalization) [{float(max(dev1, dev2)):.0e}]', max(dev1, dev2) < 1e-10)


def K_full(theta, eta, vt):
    """BPT (3.2): K = R_0 sigma(eta,u) sigma(i vt,u)/(cos eta cosh vt) [[P+, Q0],[Q0, P-]], u = -i theta."""
    u = -1j*mp.mpmathify(theta)
    P0 = lambda s: mp.cos(L*u)*mp.cos(eta)*mp.cosh(vt) - s*mp.sin(L*u)*mp.sin(eta)*mp.sinh(vt)
    Q0 = -mp.sin(L*u)*mp.cos(L*u)
    pref = mp.exp(log_R0(u) + log_sigma(eta, u) + log_sigma(1j*vt, u))/(mp.cos(eta)*mp.cosh(vt))
    M = np.array([[complex(P0(1)), complex(Q0)], [complex(Q0), complex(P0(-1))]])
    return complex(pref)*M


def kvec(Kmat):
    """k[(a,b)] = K_{b abar}: K^{ab} of GZ (3.32)-(3.33) with the plain conjugation abar."""
    Kc = Kmat @ C            # (K C)_{b a} = K_{b abar}
    return Kc.T.reshape(4)   # index a*2 + b


# matrix part of crossing-unitarity: at physical q it reduces to a scalar with the plain conjugation
# C = sigma_x (the C for which the bulk S is crossing symmetric, kappa = +1); at -q it fails with sigma_x
# and holds only with C' = [[0,1],[-1,0]], the conjugation of the kappa = -1 bulk amplitude (sec-crossing-sign)
def cu_res(qq, Cm):
    res = []
    for eta, vt in ((0.6, 0.3), (2.0, 1.1)):
        E1 = -np.cos(eta - 1j*vt)/np.sin(np.pi*LAM); E0 = -np.cos(eta + 1j*vt)/np.sin(np.pi*LAM)
        for th in (0.29, 0.77):
            km = (K_closed(np.exp(LAM*(1j*np.pi/2 - th)), E0, E1) @ Cm).T.reshape(4)
            kp = (K_closed(np.exp(LAM*(1j*np.pi/2 + th)), E0, E1) @ Cm).T.reshape(4)
            v = Rsg(np.exp(2*LAM*th), qq) @ kp
            res.append(np.linalg.norm(np.outer(km, v) - np.outer(v, km))/(np.linalg.norm(km)*np.linalg.norm(v)))
    return max(res)


Cp = np.array([[0, 1], [-1, 0]])
r_q, r_mq, r_mq2 = cu_res(q, C), cu_res(-q, C), cu_res(-q, Cp)
report(f'crossing-unitarity k(theta) = f(theta) R(2theta) k(-theta): matrix part consistent at q with C = sigma_x '
       f'[{r_q:.0e}], inconsistent at -q [{r_mq:.2f}], restored at -q only by C = [[0,1],[-1,0]] [{r_mq2:.0e}]',
       r_q < 1e-12 and r_mq > 1e-3 and r_mq2 < 1e-12)

# full scalar factor: boundary unitarity and crossing-unitarity with the book S_0
devU, devC = 0.0, 0.0
for eta, vt in ((mp.mpf('0.6'), mp.mpf('0.3')), (mp.pi*(L+1)/2, mp.mpf(0))):     # generic; Neumann
    for th in (0.29,):
        devU = max(devU, np.abs(K_full(th, eta, vt) @ K_full(-th, eta, vt) - np.eye(2)).max())
        km = kvec(K_full(1j*np.pi/2 - th, eta, vt)); kp = kvec(K_full(1j*np.pi/2 + th, eta, vt))
        S2 = complex(S0_book(2*th))*Rsg(np.exp(2*LAM*th))
        devC = max(devC, np.abs(km - S2 @ kp).max()/np.abs(km).max())
report(f'boundary unitarity K(theta)K(-theta) = 1 with R_0 sigma(eta) sigma(i vt)/(cos eta cosh vt) [{devU:.0e}]',
       devU < 1e-9)
report(f'boundary crossing-unitarity k(theta) = S_0(2theta) R(2theta) k(-theta) with the book S_0 [{devC:.0e}] '
       '(generic point and Neumann)', devC < 1e-9)

# ---------------------------------------------------------------- 5. breather reflection factors by the bootstrap (ch. 29)
# The breather B_n sits in V(theta - i u_n/2) (x) V(theta + i u_n/2) (u_n = pi - n pi/lambda) as the U_q-invariant
# vector w (in Part VI conventions the singlet exists in this ordering only).  Reflecting the two constituents,
# (1 (x) K(theta_a)) S(theta_a + theta_b) (1 (x) K(theta_b)) w = R_B(theta) w, with the full scalar factors.
from boundary_common import cop


def inv_vec(ya, yb):
    A, B = sl2_spin(1, q, ya), sl2_spin(1, q, yb)
    M = np.vstack([cop(A, B, kd, j) - (np.eye(4) if kd == 'k' else 0) for j in (0, 1) for kd in 'efk'])
    _, s, Vh = np.linalg.svd(M)
    w = Vh[-1].conj()
    return s[-1]/s[0], w/w[1]                       # normalize the s(x)sbar component to 1


def ghoshal_breather(n, u, eta, vt, lam=None):
    """Ghoshal hep-th/9310188 (3.5)-(3.10) for n = 1, 2."""
    lam = L if lam is None else lam
    pre = ((-1)**(n+1)*mp.cos(u/2 + n*mp.pi/(4*lam))*mp.cos(u/2 - mp.pi/4 - n*mp.pi/(4*lam))*mp.sin(u/2 + mp.pi/4)
           / (mp.cos(u/2 - n*mp.pi/(4*lam))*mp.cos(u/2 + mp.pi/4 + n*mp.pi/(4*lam))*mp.sin(u/2 - mp.pi/4)))
    for l in range(1, n):
        pre *= (mp.sin(u + l*mp.pi/(2*lam))*mp.cos(u/2 - mp.pi/4 - l*mp.pi/(4*lam))**2
                / (mp.sin(u - l*mp.pi/(2*lam))*mp.cos(u/2 + mp.pi/4 + l*mp.pi/(4*lam))**2))

    def Sx(x):
        if n == 1:
            return (mp.cos(x/lam) - mp.sin(u))/(mp.cos(x/lam) + mp.sin(u))
        c1, c2 = mp.cos(x/lam - mp.pi/(2*lam)), mp.cos(x/lam + mp.pi/(2*lam))
        return (mp.sin(u) - c1)*(mp.sin(u) - c2)/((mp.sin(u) + c1)*(mp.sin(u) + c2))
    return pre*Sx(eta)*Sx(1j*vt)


eta, vt = mp.mpf('0.6'), mp.mpf('0.3')
devs, ok_order = [], True
for n in (1, 2):
    un = np.pi - n*np.pi/LAM
    for th in (0.3, -0.55):
        ta, tb = th - 1j*un/2, th + 1j*un/2
        r_in, w = inv_vec(np.exp(LAM*ta), np.exp(LAM*tb))
        ok_order &= r_in < 1e-12 and inv_vec(np.exp(LAM*tb), np.exp(LAM*ta))[0] > 1e-3
        S2 = complex(a_gz(ta + tb))*Rsg(np.exp(LAM*(ta + tb)))
        v = np.kron(np.eye(2), K_full(ta, eta, vt)) @ S2 @ np.kron(np.eye(2), K_full(tb, eta, vt)) @ w
        RB = v[1]
        devs.append(np.linalg.norm(v - RB*w)/np.linalg.norm(v))
        devs.append(abs(RB/complex(ghoshal_breather(n, -1j*mp.mpf(th), eta, vt)) - 1))
report(f'breather bootstrap from the BPT soliton K: B_1, B_2 reflect diagonally with exactly Ghoshal (3.5)-(3.10) '
       f'[{max(devs):.0e}]; singlet in V(theta - i u_n/2) (x) V(theta + i u_n/2) only', max(devs) < 1e-10 and ok_order)
S11 = lambda th: (mp.sinh(th) + 1j*mp.sin(mp.pi/L))/(mp.sinh(th) - 1j*mp.sin(mp.pi/L))     # eq-sg-s11, xi = pi/lambda
R1B = lambda th: ghoshal_breather(1, -1j*th, eta, vt)
pts = [mp.mpc('0.31', '0.07'), mp.mpc('-0.6', '0.2')]
dev = max(max(abs(R1B(t)*R1B(-t) - 1), abs(R1B(1j*mp.pi/2 - t) - S11(2*t)*R1B(1j*mp.pi/2 + t))) for t in pts)
report(f'B_1 reflection factor: unitarity and crossing-unitarity with S_11 = f_(xi/pi) [{float(dev):.0e}]', dev < 1e-20)

# ---------------------------------------------------------------- 6. continuation to sinh-Gordon (Corrigan-Delius (2.5))
Bs = mp.mpf('0.71')
lam_shg = -2/Bs                                         # lambda -> -2/B under beta_SG^2 -> -beta_SG^2
cdb = lambda x, th: mp.sinh(th/2 + 1j*mp.pi*x/4)/mp.sinh(th/2 - 1j*mp.pi*x/4)        # CD (2.3) = eq-blocks, h = 2
dev = 0
for eta_, vt_ in ((mp.mpf('0.9'), mp.mpf('0.4')), (mp.mpf('2.3'), mp.mpf('0'))):
    E, F = Bs*eta_/mp.pi, 1j*Bs*vt_/mp.pi
    for th in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
        cd = (cdb(1, th)*cdb(2 - Bs/2, th)*cdb(1 + Bs/2, th)
              / (cdb(1 - E, th)*cdb(1 + E, th)*cdb(1 - F, th)*cdb(1 + F, th)))
        dev = max(dev, abs(ghoshal_breather(1, -1j*th, eta_, vt_, lam_shg)/cd - 1))
report(f'Ghoshal B_1 factor at lambda = -2/B equals Corrigan-Delius (2.5) with E = B eta/pi, F = i B vt/pi [{float(dev):.0e}]',
       dev < 1e-20)
th = mp.mpc('0.4', '0.1')
KD_shg = cdb(2 - Bs/2, th)*cdb(1 + Bs/2, th)/cdb(1, th)
lim = ghoshal_breather(1, -1j*th, mp.mpf(0), mp.mpf(60), lam_shg)
EN = Bs*(mp.pi*(lam_shg + 1)/2)/mp.pi
report(f'Dirichlet phi_0 = 0 (eta = 0, vt -> infinity) gives K_D = (2-B/2)(1+B/2)/(1) [{float(abs(lim/KD_shg - 1)):.0e}]; '
       f'Neumann eta = pi(lambda+1)/2 gives E = B/2 - 1, i.e. |E| = 1 - B/2', abs(lim/KD_shg - 1) < 1e-8 and
       abs(EN - (Bs/2 - 1)) < 1e-25)

# ---------------------------------------------------------------- 7. excited boundary states (BPTT hep-th/0106070)
L = mp.mpf('3.3')                                       # lambda large enough for two boundary states nu_0, nu_1


eta = mp.mpf('5.5'); etab = mp.pi*(L + 1) - eta
nu = lambda n: eta/L - (2*n + 1)*mp.pi/(2*L)
a_u = lambda u: a_gz(1j*u)
b_u = lambda u: mp.sin(L*u)/mp.sin(L*(mp.pi - u))*a_gz(1j*u)
dev = 0
for u in (mp.mpc('0.3', '0.2'), mp.mpc('-0.5', '0.1')):
    lhs = a_u(u - nu(0))*b_u(u + nu(0))
    rhs = mp.exp(log_sigma(etab, u) - log_sigma(eta, u))*mp.cos(eta)/mp.cos(etab)
    dev = max(dev, abs(lhs/rhs - 1))
ok_ph = 0 < nu(1) < nu(0) < mp.pi/2 and nu(2) < 0 and eta <= mp.pi*(L + 1)/2
report(f'BPTT bootstrap on the first excited state: a(u - nu_0) Q(eta, vt, u) b(u + nu_0) = Q(etabar, vt, u), '
       f'etabar = pi(lambda+1) - eta [{float(dev):.0e}]; nu_n = eta/lambda - (2n+1)pi/(2lambda)', dev < 1e-15 and ok_ph)
# pole of the ground-state P^+ at u = nu_n (from sigma(eta, u)); simple pole: (u - nu) sigma finite and nonzero
for n in (0, 1):
    d = mp.mpf('1e-8')
    r1 = abs(mp.exp(log_sigma(eta, nu(n) + d))*d); r2 = abs(mp.exp(log_sigma(eta, nu(n) + 2*d))*2*d)
    ok_ph &= abs(r1/r2 - 1) < 1e-5 and r1 > 1e-12
report('sigma(eta, u) has simple poles at u = nu_0, nu_1 in the physical strip (boundary bound states |0>, |1>)', ok_ph)
L = mp.mpf(LAM)
Kx = np.array([[0, 1], [1, 0]])
eta0, vt0 = mp.mpf('2.9'), mp.mpf('0.3')
etab0 = mp.pi*(L + 1) - eta0
K0 = lambda th: Kx @ K_full(th, etab0, vt0) @ Kx       # K on |0>: s <-> sbar and eta -> etabar (BPTT (3.6))
th = 0.29
devU = np.abs(K0(th) @ K0(-th) - np.eye(2)).max()
km = kvec(K0(1j*np.pi/2 - th)); kp = kvec(K0(1j*np.pi/2 + th))
devC = np.abs(km - complex(S0_book(2*th))*Rsg(np.exp(2*LAM*th)) @ kp).max()/np.abs(km).max()
report(f'K on the first excited boundary |0> satisfies unitarity [{devU:.0e}] and crossing-unitarity [{devC:.0e}]',
       devU < 1e-9 and devC < 1e-9)

# ---------------------------------------------------------------- 8. limits of Al. B. Zamolodchikov's UV-IR relation (BPT hep-th/0108157 (3.2))
# cos(eta/(lambda+1)) cosh(vt/(lambda+1)) = r cos(b phi0/2),  sin(eta/(lambda+1)) sinh(vt/(lambda+1)) = r sin(b phi0/2),
# r = M_0/M_crit, beta_SG^2/8pi = 1/(lambda+1), M_crit = sqrt(2 mu/sin(beta_SG^2/8)).
def uvir_lhs(eta, vt, lam=LAM):
    return mp.cos(eta/(lam+1))*mp.cosh(vt/(lam+1)), mp.sin(eta/(lam+1))*mp.sinh(vt/(lam+1))


ok = True
c, s_ = uvir_lhs(mp.pi*(LAM + 1)/2, mp.mpf(0))
ok &= abs(c) < 1e-20 and abs(s_) < 1e-20                                               # Neumann <-> M_0 = 0
bphi = mp.mpf('0.4'); V = mp.mpf(40)
c, s_ = uvir_lhs((LAM + 1)*bphi, (LAM + 1)*V)
ok &= abs(s_/c - mp.tan(bphi)) < 1e-20 and c > 1e15                                     # Dirichlet: r -> oo, eta = 4pi phi_0/beta
c, s_ = uvir_lhs(mp.pi/2, mp.mpf(0))
ok &= abs(c - mp.cos(mp.pi/(2*(LAM + 1)))) < 1e-20 and abs(s_) < 1e-20                  # eps-hat = 0 point, phi_0 = 0
bs2 = mp.mpf('1e-4'); m0 = mp.mpf(1); mu = m0**2/bs2                                     # classical limit: mu = m_0^2/beta^2
Mcrit = mp.sqrt(2*mu/mp.sin(bs2/8))
ok &= abs(Mcrit/(4*m0/bs2) - 1) < 1e-8                                                  # = M_0 of the C = 1 boundary
report('Zamolodchikov UV-IR relation: M_0 -> 0 gives Neumann (eta, vt) = (pi(lambda+1)/2, 0); M_0 -> infinity gives '
       'the Dirichlet eta = 4 pi phi_0/beta_SG; (pi/2, 0) <-> M_0/M_crit = cos(pi/(2(lambda+1))); classically '
       'M_crit -> 4 m_0/beta_SG^2, the C_0 = C_1 = 1 boundary', ok)
