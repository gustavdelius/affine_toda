"""The a_r^(1) defect matrix as a gauge transformation in both light-cone directions (Part VIII, ch. 38: sec-defect-backlund-an).

Written from the book's definitions, without defect_common.py.  Conventions: roots alpha_j = e_j - e_{j+1} in R^h (indices mod h),
fields in the hyperplane sum = 0, H.v = diag(v_0, ..., v_r), E_alpha_j = e_{j,j+1}, E_-alpha_j = e_{j+1,j}, S = sum_j E_alpha_j,
pi e_k = e_{k-1}, D = 2(1 - pi)^-1, E = 1 - D (type I_-), (D, E) -> (-D, -E) and sigma <-> 1/sigma in B_d (type I_+), b = beta or
i beta.  Lax pair eq-lax-pair (c_j = 1), a_t, a_x of eq-boundary-lax-at, A_+- = a_t +- a_x.  Defect conditions eq-defect-conditions,
potential eq-defect-potential, defect matrix eq-defect-lax-matrix: L_d = exp(-b H.(E phi + D psi)/2) Ltilde exp(b H.(D^T phi - E psi)/2),
Ltilde_I- = 1 + (sigma/lambda) S^T, Ltilde_I+ = sum_k (lambda/sigma)^k S^k.

Checks (book ids they support):
 1. Components (sympy, h = 3): L_d^{I-} has (L_d)_kk = exp(b(phi_k - psi_k)/2) and (L_d)_{j+1,j} = (sigma/lambda) exp(b(psi_j - phi_{j+1})/2),
    all other entries 0; det L_d = 1 - (-sigma/lambda)^h for h = 2, ..., 5.
 2. Hand proof (sympy, h = 3 and 4, b and sigma symbolic): d_+ L_d - L_d A_+(psi) + A_+(phi) L_d vanishes at order lambda identically,
    its diagonal at order lambda^0 is b exp(b(phi_k - psi_k)/2) times [d_+(phi_k - psi_k) + (m sigma/b)(f_k - f_{k-1})], and its
    order lambda^-1 vanishes identically; d_- L_d - L_d A_-(psi) + A_-(phi) L_d vanishes at orders lambda^0 and lambda^-2 identically,
    and its entry (j+1, j) at order lambda^-1 is sigma b exp(b(psi_j - phi_{j+1})/2) times
    [d_-(psi_j - phi_{j+1}) - (m/(b sigma))(g_j - g_{j+1})].  So each light-cone equation is one half of the component Baecklund
    transformation of sec-defect-backlund-an.
 3. Both types, h = 2, ..., 5, real and imaginary coupling, complex lambda (numerical): with phi_x, psi_x from the defect conditions
    imposed at a point and arbitrary phi, psi, phi_t, psi_t there, d_x L_d = L_d a_x(psi) - a_x(phi) L_d and
    d_t L_d = L_d a_t(psi) - a_t(phi) L_d (chain rule, relative residual < 1e-12).  Controls: sigma -> 1.2 sigma in L_d, and
    E -> 0.6 E in the conditions, fail.
Runtime: about 10 seconds.
"""
import numpy as np
import sympy as sp


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ---------------- numerical set-up for general h
def setup(h):
    I = np.eye(h)
    Pm = np.zeros((h, h))
    for k in range(h):
        Pm[(k - 1) % h, k] = 1.0          # pi e_k = e_{k-1}
    D = 2*np.linalg.pinv(I - Pm)
    P = np.linalg.pinv(I - Pm) @ (I - Pm)  # projector onto the hyperplane
    E = P - D
    alphas = [I[j] - I[(j + 1) % h] for j in range(h)]
    Ep = [np.outer(I[j], I[(j + 1) % h]) for j in range(h)]   # E_alpha_j = e_{j,j+1}
    Em = [M.T for M in Ep]
    S = sum(Ep)
    return dict(h=h, D=D, E=E, P=P, alphas=alphas, Ep=Ep, Em=Em, S=S)


def Hdot(v):
    return np.diag(v)


def grads(d, phi, psi, b, m, sig, typ):
    """Gradients of B_d (eq-defect-potential) with the type's sigma; D in f_j, g_j is the type I_- matrix."""
    D = d['D']
    sf, sg = (sig, 1/sig) if typ == '-' else (1/sig, sig)
    gphi, gpsi = np.zeros(d['h'], complex), np.zeros(d['h'], complex)
    for a in d['alphas']:
        fj = np.exp(b*a @ (D.T @ phi + D @ psi)/2)
        gj = np.exp(b*a @ (D @ (phi - psi))/2)
        gphi += m/b**2*(sf*(b/2)*(D @ a)*fj + sg*(b/2)*(D.T @ a)*gj)
        gpsi += m/b**2*(sf*(b/2)*(D.T @ a)*fj - sg*(b/2)*(D.T @ a)*gj)
    return gphi, gpsi


def kinetic(d, typ, Escale=1.0):
    D, E = d['D'], Escale*d['E']
    return (D, E) if typ == '-' else (-D, -E)


def conditions(d, phi, psi, phit, psit, b, m, sig, typ, Escale=1.0):
    """phi_x, psi_x from eq-defect-conditions (F = E)."""
    Dk, Ek = kinetic(d, typ, Escale)
    gphi, gpsi = grads(d, phi, psi, b, m, sig, typ)
    phix = Ek @ phit + Dk @ psit - gphi
    psix = Dk.T @ phit - Ek @ psit + gpsi
    return phix, psix


def Ltilde(d, lam, sig, typ):
    if typ == '-':
        return np.eye(d['h']) + sig/lam*d['S'].T
    return sum(np.linalg.matrix_power((lam/sig)*d['S'], k) for k in range(d['h']))


def Ld(d, phi, psi, lam, b, sig, typ):
    Dk, Ek = kinetic(d, typ)
    u, w = Ek @ phi + Dk @ psi, Dk.T @ phi - Ek @ psi
    return np.diag(np.exp(-b*u/2)) @ Ltilde(d, lam, sig, typ) @ np.diag(np.exp(b*w/2))


def dLd(d, phi, psi, dphi, dpsi, lam, b, sig, typ):
    Dk, Ek = kinetic(d, typ)
    L = Ld(d, phi, psi, lam, b, sig, typ)
    du, dw = Ek @ dphi + Dk @ dpsi, Dk.T @ dphi - Ek @ dpsi
    return -b/2*Hdot(du) @ L + L @ (b/2*Hdot(dw))


def lax(d, f, ft, fx, lam, b, m):
    """a_t, a_x of eq-boundary-lax-at (c_j = 1)."""
    pos = sum(np.exp(b*a @ f/2)*Ep for a, Ep in zip(d['alphas'], d['Ep']))
    neg = sum(np.exp(b*a @ f/2)*Em for a, Em in zip(d['alphas'], d['Em']))
    at = b/2*Hdot(fx) + m/2*(lam*pos - neg/lam)
    ax = b/2*Hdot(ft) + m/2*(lam*pos + neg/lam)
    return at, ax


rng = np.random.default_rng(38)


def hyper(h):
    v = rng.normal(size=h) + 1j*rng.normal(size=h)
    return v - v.mean()


# ---------------- 1. components and determinant (sympy for h = 3)
def sym_setup(h):
    phis = sp.symbols('phi0:%d' % h)
    psis = sp.symbols('psi0:%d' % h)
    # fields in the hyperplane: last component minus the sum of the others
    phi = sp.Matrix(list(phis[:-1]) + [-sum(phis[:-1])])
    psi = sp.Matrix(list(psis[:-1]) + [-sum(psis[:-1])])
    Pm = sp.zeros(h, h)
    for k in range(h):
        Pm[(k - 1) % h, k] = 1
    I = sp.eye(h)
    D = 2*(I - Pm).pinv()
    P = (I - Pm).pinv()*(I - Pm)
    E = P - D
    S = sp.zeros(h, h)
    for j in range(h):
        S[j, (j + 1) % h] = 1
    return phi, psi, D, E, S


b, sig, lam, m = sp.symbols('b sigma lambda m')
ok1 = True
for h in [3]:
    phi, psi, D, E, S = sym_setup(h)
    u, w = E*phi + D*psi, D.T*phi - E*psi
    L = sp.diag(*[sp.exp(-b*x/2) for x in u])*(sp.eye(h) + sig/lam*S.T)*sp.diag(*[sp.exp(b*x/2) for x in w])
    target = sp.zeros(h, h)
    for k in range(h):
        target[k, k] = sp.exp(b*(phi[k] - psi[k])/2)
        target[(k + 1) % h, k] = sig/lam*sp.exp(b*(psi[k] - phi[(k + 1) % h])/2)
    ok1 &= all(sp.simplify(x) == 0 for x in (L - target))
for h in [2, 3, 4, 5]:
    d = setup(h)
    for _ in range(3):
        lv, sv, bv = rng.normal() + 1j*rng.normal(), 0.6 + 0.3j, 0.9
        Lv = Ld(d, hyper(h), hyper(h), lv, bv, sv, '-')
        ok1 &= abs(np.linalg.det(Lv) - (1 - (-sv/lv)**h)) < 1e-12
report('1. L_d^{I-} in components: diagonal exp(b(phi_k - psi_k)/2), subdiagonal (sigma/lambda) exp(b(psi_j - phi_{j+1})/2) '
       '(sympy, h = 3); det L_d = 1 - (-sigma/lambda)^h (h = 2, ..., 5)', ok1)

# ---------------- 2. the hand proof, order by order (sympy)
ok2 = True
for h in [3, 4]:
    phi, psi, D, E, S = sym_setup(h)
    dpp = sp.Matrix(sp.symbols('dpphi0:%d' % h))   # d_+ phi (arbitrary components)
    dps = sp.Matrix(sp.symbols('dppsi0:%d' % h))
    dmp = sp.Matrix(sp.symbols('dmphi0:%d' % h))
    dms = sp.Matrix(sp.symbols('dmpsi0:%d' % h))
    target = sp.zeros(h, h)
    for k in range(h):
        target[k, k] = sp.exp(b*(phi[k] - psi[k])/2)
        target[(k + 1) % h, k] = sig/lam*sp.exp(b*(psi[k] - phi[(k + 1) % h])/2)
    L = target

    def dL(dph, dps_):
        out = sp.zeros(h, h)
        for k in range(h):
            out[k, k] = b/2*(dph[k] - dps_[k])*L[k, k]
            kp = (k + 1) % h
            out[kp, k] = b/2*(dps_[k] - dph[kp])*L[kp, k]
        return out

    def Apm(f, df, sign):
        A = sp.zeros(h, h)
        for k in range(h):
            A[k, k] = sign*b/2*df[k]
        for j in range(h):
            jp = (j + 1) % h
            ex = sp.exp(b*(f[j] - f[jp])/2)
            if sign > 0:
                A[j, jp] += m*lam*ex
            else:
                A[jp, j] += -m/lam*ex
        return A
    fj = [sp.exp(b*(psi[j] - phi[(j + 1) % h])) for j in range(h)]
    gj = [sp.exp(b*(phi[j] - psi[j])) for j in range(h)]
    Rp = sp.expand(dL(dpp, dps) - L*Apm(psi, dps, 1) + Apm(phi, dpp, 1)*L)
    Rm = sp.expand(dL(dmp, dms) - L*Apm(psi, dms, -1) + Apm(phi, dmp, -1)*L)

    def order(M, n):
        """Coefficient of lambda^n of a matrix that is a Laurent polynomial in lambda."""
        return M.applyfunc(lambda x: sp.Poly(sp.expand(x*lam**3), lam).coeff_monomial(lam**(n + 3)))

    zero = lambda M: all(sp.simplify(x) == 0 for x in M)
    ok2 &= zero(order(Rp, 1)) and zero(order(Rp, -1))
    R0 = order(Rp, 0)
    for k in range(h):
        pred = b*L[k, k]*((dpp[k] - dps[k]) + m*sig/b*(fj[k] - fj[(k - 1) % h]))
        ok2 &= sp.simplify(R0[k, k] - pred) == 0
    ok2 &= all(sp.simplify(R0[i, j]) == 0 for i in range(h) for j in range(h) if i != j)
    ok2 &= zero(order(Rm, 0)) and zero(order(Rm, -2))
    Rm1 = order(Rm, -1)
    for j in range(h):
        jp = (j + 1) % h
        pred = sig*b*sp.exp(b*(psi[j] - phi[jp])/2)*((dms[j] - dmp[jp]) - m/(b*sig)*(gj[j] - gj[jp]))
        ok2 &= sp.simplify(Rm1[jp, j] - pred) == 0
    ok2 &= all(sp.simplify(Rm1[i, j]) == 0 for i in range(h) for j in range(h) if i != (j + 1) % h)
report('2. order by order in lambda, the d_+ equation for L_d is d_+(phi_k - psi_k) = -(m sigma/b)(f_k - f_{k-1}) and the d_- '
       'equation is d_-(psi_j - phi_{j+1}) = (m/(b sigma))(g_j - g_{j+1}); all other orders vanish identically (h = 3, 4)', ok2)

# ---------------- 3. both directions, both types, numerically
def gauge_resid(h, typ, b_, sig_scale=1.0, Escale=1.0, npts=4):
    d = setup(h)
    worst = 0.0
    for _ in range(npts):
        phi, psi, phit, psit = hyper(h), hyper(h), hyper(h), hyper(h)
        lv, sv, mv = rng.normal() + 1j*rng.normal(), 0.7 + 0.4j, 1.3
        phix, psix = conditions(d, phi, psi, phit, psit, b_, mv, sv, typ, Escale)
        L = Ld(d, phi, psi, lv, b_, sig_scale*sv, typ)
        at_f, ax_f = lax(d, phi, phit, phix, lv, b_, mv)
        at_s, ax_s = lax(d, psi, psit, psix, lv, b_, mv)
        dLt = dLd(d, phi, psi, phit, psit, lv, b_, sig_scale*sv, typ)
        dLx = dLd(d, phi, psi, phix, psix, lv, b_, sig_scale*sv, typ)
        rt = np.linalg.norm(dLt - (L @ at_s - at_f @ L))/np.linalg.norm(dLt)
        rx = np.linalg.norm(dLx - (L @ ax_s - ax_f @ L))/np.linalg.norm(dLx)
        worst = max(worst, rt, rx)
    return worst


ok3, ctrl3 = True, True
for h in [2, 3, 4, 5]:
    for typ in '-+':
        for b_ in [0.8, 0.8j]:
            ok3 &= gauge_resid(h, typ, b_) < 1e-12
            ctrl3 &= gauge_resid(h, typ, b_, sig_scale=1.2) > 1e-3
            if h > 2:      # for h = 2, E = 0
                ctrl3 &= gauge_resid(h, typ, b_, Escale=0.6) > 1e-3
report('3. d_x L_d = L_d a_x(psi) - a_x(phi) L_d and d_t L_d = L_d a_t(psi) - a_t(phi) L_d on the defect conditions, '
       'h = 2..5, both types, both couplings (controls: 1.2 sigma in L_d and 0.6 E in the conditions fail)', ok3 and ctrl3)
