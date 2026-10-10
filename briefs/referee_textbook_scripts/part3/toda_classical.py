"""Classical affine Toda checks for Part III (chapters 9-11).

1. Lax pair (sec-lax-pair): zero curvature for a_2^(1) <=> field equation, with
   coefficients sqrt(n_j alpha_j^2/2).
2. Spin-3 conserved current of sinh-Gordon (sec-conserved-charges).
3. Masses: eigenvalues of m^2 sum n_j alpha_j alpha_j^T are proportional to the
   Perron-Frobenius vector (simply laced), sum rule sum m_a^2 = 2h m^2.
4. Cubic couplings |C_abc| proportional to the area of the mass triangle, and
   Dorey's fusing rule with gamma_a=(1-w^{-1})lambda_a, w=s_B s_W, for e7, e8, a_n.
5. Hirota (sec-hirota): bilinear equations, one- and two-soliton tau functions of
   a_n^(1) solve the time-dependent field equation; mass M_a = 2 h m_a / beta^2.
"""
import itertools
import numpy as np
import sympy as sp
import mpmath as mp


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ---------------------------------------------------------------- 1. Lax pair
xp, xm, lam, m, b = sp.symbols('xp xm lambda m beta', positive=True)
u1, u2 = sp.Function('u1')(xp, xm), sp.Function('u2')(xp, xm)
u = [u1, u2, -u1-u2]                       # phi as traceless diag(u1,u2,u3)
Phi = sp.diag(*u)
eij = lambda i, j: sp.Matrix(3, 3, lambda r, c: 1 if (r, c) == (i, j) else 0)
Ep = [eij(0, 1), eij(1, 2), eij(2, 0)]      # E_{alpha_1}, E_{alpha_2}, E_{alpha_0}
al = [lambda v: v[0]-v[1], lambda v: v[1]-v[2], lambda v: v[2]-v[0]]
coef = 1                                   # sqrt(n_j alpha_j^2/2) = 1 for a_n
Ap = b/2*sp.diff(Phi, xp) + m*lam*sum((coef*sp.exp(b*al[j](u)/2)*Ep[j] for j in range(3)), sp.zeros(3))
Am = -b/2*sp.diff(Phi, xm) - m/lam*sum((coef*sp.exp(b*al[j](u)/2)*Ep[j].T for j in range(3)), sp.zeros(3))
Fc = sp.diff(Am, xp) - sp.diff(Ap, xm) + Ap*Am - Am*Ap
# field equation d+d- phi = -(m^2/beta) sum_j n_j alpha_j e^{beta alpha_j.phi}, in diag components
alvec = [sp.Matrix([1, -1, 0]), sp.Matrix([0, 1, -1]), sp.Matrix([-1, 0, 1])]
rhs = -(m**2/b)*sum((alvec[j]*sp.exp(b*al[j](u)) for j in range(3)), sp.zeros(3, 1))
sub = {sp.Derivative(u1, xp, xm): rhs[0], sp.Derivative(u2, xp, xm): rhs[1]}
Fs = sp.simplify(Fc.subs(sub))
report('Lax pair zero curvature <=> a_2^(1) field equation', Fs == sp.zeros(3))

# ---------------------------------------------------------------- 2. sinh-Gordon spin 3
f = sp.Function('f')(xp, xm)
T4 = sp.diff(f, xp, 2)**2 + b**2/4*sp.diff(f, xp)**4
Th2 = m**2*sp.cosh(b*f)*sp.diff(f, xp)**2
expr = sp.diff(T4, xm) + sp.diff(Th2, xp)
eom = -(m**2/b)*sp.sinh(b*f)           # d+d- f for V=(m^2/b^2)(cosh b f - 1), d+-=d_t+-d_x
expr = expr.subs(sp.Derivative(f, xp, xp, xm), sp.diff(eom, xp)).subs(sp.Derivative(f, xp, xm), eom)
report('sinh-Gordon: d_-[(d_+^2 f)^2+(b^2/4)(d_+ f)^4] + d_+[m^2 cosh(bf)(d_+ f)^2] = 0', sp.simplify(expr) == 0)

# ---------------------------------------------------------------- root data
def cartan(typ, n):
    C = 2*np.eye(n)
    if typ == 'a':
        for i in range(n-1): C[i, i+1] = C[i+1, i] = -1
    elif typ == 'd':
        for i in range(n-2): C[i, i+1] = C[i+1, i] = -1
        C[n-3, n-1] = C[n-1, n-3] = -1
    elif typ == 'e':
        for i in range(n-2): C[i, i+1] = C[i+1, i] = -1   # chain 1..n-1, node n attached to node 3
        C[2, n-1] = C[n-1, 2] = -1
    return C


def simple_roots(C):
    L = np.linalg.cholesky(C)                 # rows: simple roots with alpha_i.alpha_j = C_ij
    return L


def all_roots(al):
    r = len(al); roots = {tuple(np.round(a, 12)) for a in al}; frontier = list(al)
    while frontier:
        new = []
        for v in frontier:
            for a in al:
                for s in (1, -1):
                    w_ = v - (v @ a)*a if s == 1 else v
                    t = tuple(np.round(w_, 12))
                    if t not in roots:
                        roots.add(t); new.append(w_)
        frontier = new
    roots |= {tuple(-np.array(t)) for t in roots}
    return [np.array(t) for t in roots]


def data(typ, n):
    C = cartan(typ, n); al = simple_roots(C)
    roots = all_roots(list(al))
    theta = max(roots, key=lambda v: np.linalg.solve(al.T, v).sum())
    marks = np.round(np.linalg.solve(al.T, theta)).astype(int)
    h = 1 + marks.sum()
    a0 = -theta
    A = np.vstack([a0, al]); nj = np.concatenate([[1], marks])
    M2 = sum(nj[j]*np.outer(A[j], A[j]) for j in range(n+1))
    return C, al, roots, h, A, nj, M2


ok_pf = True
for typ, n in [('a', 2), ('a', 5), ('d', 4), ('d', 6), ('e', 6), ('e', 7), ('e', 8)]:
    C, al, roots, h, A, nj, M2 = data(typ, n)
    ev = np.sort(np.linalg.eigvalsh(M2))
    pf = np.linalg.eigh(C)[1][:, 0]; pf = np.abs(pf)
    # masses sqrt(ev) should be proportional to PF components (as multisets)
    ratio = np.sort(np.sqrt(ev))/np.sqrt(ev).min()
    pfr = np.sort(pf)/pf.min()
    ok_pf &= np.allclose(ratio, pfr, atol=1e-9) and abs(ev.sum() - 2*h) < 1e-9 and len(roots) == n*h
report('masses = Perron-Frobenius vector; sum m_a^2 = 2h m^2; |Phi| = r h (a2,a5,d4,d6,e6,e7,e8)', ok_pf)

# ---------------------------------------------------------------- 4. couplings and Dorey's rule
def dorey(typ, n):
    C, al, roots, h, A, nj, M2 = data(typ, n)
    lamw = np.linalg.inv(C) @ al                   # fundamental weights (rows)
    col = {0: 1}; stack = [0]
    while stack:
        i = stack.pop()
        for j in range(n):
            if j != i and C[i, j] < 0 and j not in col:
                col[j] = -col[i]; stack.append(j)
    refl = lambda v: np.eye(n) - np.outer(v, v)
    sW = np.eye(n); sB = np.eye(n)
    for i in range(n):
        if col[i] == 1: sW = refl(al[i]) @ sW
        else: sB = refl(al[i]) @ sB
    w = sB @ sW; wi = np.linalg.inv(w)
    gam = [(np.eye(n) - wi) @ lamw[a] for a in range(n)]
    orbits = [[np.linalg.matrix_power(w, p) @ gam[a] for p in range(h)] for a in range(n)]
    allorb = {tuple(np.round(v, 6)) for o in orbits for v in o}
    partition = len(allorb) == n*h and all(any(np.allclose(g, rt) for rt in roots) for g in gam)
    # masses from Coxeter-plane projections: PF eigenvector of C gives the plane
    ev, V = np.linalg.eigh(M2)
    # match species a with mass eigenvectors via the PF vector ordering
    return C, al, h, A, nj, M2, orbits, partition


def cubic(A, nj, vecs):
    return lambda a, b_, c: sum(nj[j]*(A[j] @ vecs[a])*(A[j] @ vecs[b_])*(A[j] @ vecs[c]) for j in range(len(nj)))


ok_rule = True; ok_area = True; ok_part = True
for typ, n in [('e', 7), ('e', 8)]:
    C, al, h, A, nj, M2, orbits, partition = dorey(typ, n)
    ok_part &= partition
    ev, V = np.linalg.eigh(M2); mass = np.sqrt(ev)
    # assign species: mass of orbit a = length of projection of gamma_a onto the PF plane (up to scale)
    pf = np.abs(np.linalg.eigh(C)[1][:, 0])
    # species index a <-> PF component; mass eigenvectors sorted by mass, PF sorted likewise
    order_pf = np.argsort(pf); order_m = np.argsort(mass)
    spec_vec = {order_pf[k]: V[:, order_m[k]] for k in range(n)}
    spec_mass = {order_pf[k]: mass[order_m[k]] for k in range(n)}
    Cf = cubic(A, nj, spec_vec)
    consts = []
    for a, b_, c in itertools.combinations_with_replacement(range(n), 3):
        val = Cf(a, b_, c)
        ma, mb, mc = spec_mass[a], spec_mass[b_], spec_mass[c]
        s = (ma+mb+mc)/2; area2 = s*(s-ma)*(s-mb)*(s-mc)
        tri = area2 > 1e-9
        rule = any(np.allclose(x+y+z, 0) for x in orbits[a] for y in orbits[b_] for z in orbits[c])
        if (abs(val) > 1e-8) != rule: ok_rule = False
        if abs(val) > 1e-8:
            if not tri: ok_area = False
            else: consts.append(abs(val)/np.sqrt(area2))
    ok_area &= np.ptp(consts) < 1e-9*max(consts) and abs(consts[0]-4/np.sqrt(h)) < 1e-9
    print(f'   {typ}{n}: {len(consts)} nonzero couplings, |C|/area = {consts[0]:.6f} (expected 4/sqrt(h) = {4/np.sqrt(h):.6f})')
report('Coxeter orbits of gamma_a=(1-w^-1)lambda_a partition the roots (e7, e8)', ok_part)
report("Dorey's rule: C_abc != 0 iff roots in Omega_a, Omega_b, Omega_c sum to zero (e7, e8)", ok_rule)
report('|C_abc| = (4 beta/sqrt h) * area of mass triangle for every triple (e7, e8; m=beta=1)', ok_area)

# ---------------------------------------------------------------- 5. Hirota, a_n^(1)
mp.mp.dps = 30


def an_alpha(n):
    h = n+1; vecs = []
    for j in range(h):
        v = [0]*h; v[j] = 1; v[(j+1) % h] = -1; vecs.append(v)
    return vecs                                   # alpha_j in R^{h}, sum zero; j=0..n


def tau_funcs(n, species, thetas, xis, mm=1):
    h = n+1; om = mp.e**(2j*mp.pi/h)
    def taus(x, t):
        Es = [mp.e**(2*mm*mp.sin(mp.pi*a/h)*(x*mp.cosh(th) - t*mp.sinh(th)) + xi) for a, th, xi in zip(species, thetas, xis)]
        res = []
        for j in range(h):
            tj = 0
            for k in range(len(species)+1):
                for sub in itertools.combinations(range(len(species)), k):
                    term = 1
                    for i in sub: term *= om**(j*species[i])*Es[i]
                    for i1, i2 in itertools.combinations(sub, 2):
                        A_, B_ = mp.pi*species[i1]/h, mp.pi*species[i2]/h
                        c = mp.cosh(thetas[i1]-thetas[i2])
                        term *= (c - mp.cos(A_-B_))/(c - mp.cos(A_+B_))
                    tj += term
            res.append(tj)
        return res
    return taus


def field_residual(n, species, thetas, xis, beta=mp.mpf('0.7'), mm=1, pts=4):
    h = n+1; alv = an_alpha(n); taus = tau_funcs(n, species, thetas, xis, mm)
    def phi(x, t):
        T = taus(x, t)
        return [1j/beta*sum(alv[j][k]*mp.log(T[j]) for j in range(h)) for k in range(h)]
    worst = 0
    rng = np.random.default_rng(3)
    for _ in range(pts):
        x0, t0 = mp.mpf(rng.uniform(-1, 1)), mp.mpf(rng.uniform(-1, 1))
        for k in range(h):
            ptt = mp.diff(lambda t: phi(x0, t)[k], t0, 2)
            pxx = mp.diff(lambda x: phi(x, t0)[k], x0, 2)
            P = phi(x0, t0)
            force = 1j*mm**2/beta*sum(alv[j][k]*mp.e**(1j*beta*sum(alv[j][l]*P[l] for l in range(h))) for j in range(h))
            worst = max(worst, abs(ptt - pxx - force))
    return worst


r1 = max(field_residual(n, [a], [mp.mpf('0.4')], [mp.mpc('0.3', '0.7')]) for n, a in [(2, 1), (3, 2), (4, 1)])
report(f'one-soliton tau_j=1+omega^(ja)E solves (d_t^2-d_x^2)phi = (i m^2/beta) sum alpha_j e^(i beta alpha_j.phi)  [res {float(r1):.1e}]', r1 < 1e-15)
r2 = max(field_residual(n, sp_, [mp.mpf('0.5'), mp.mpf('-0.3')], [mp.mpc('0.2', '0.9'), mp.mpc('-0.1', '0.4')])
         for n, sp_ in [(2, [1, 2]), (3, [1, 2]), (3, [1, 1])])
report(f'two-soliton tau with A_ab=(cosh th-cos(A-B))/(cosh th-cos(A+B)) solves field equation  [res {float(r2):.1e}]', r2 < 1e-15)


def soliton_mass(n, a, beta=mp.mpf('0.9'), mm=1):
    """Energy int (1/2) phi'.phi' + V over [-L, L] with phi' = (i/beta) sum alpha_j tau_j'/tau_j."""
    h = n+1; alv = an_alpha(n); om = mp.e**(2j*mp.pi/h); ma = 2*mm*mp.sin(mp.pi*a/h); xi = mp.mpc(0, '0.37')
    def dens(x):
        E = mp.e**(ma*x + xi)
        T = [1 + om**(j*a)*E for j in range(h)]; dT = [ma*om**(j*a)*E for j in range(h)]
        dP = [1j/beta*sum(alv[j][k]*dT[j]/T[j] for j in range(h)) for k in range(h)]
        # e^{i beta alpha_j.phi} = prod_l tau_l^{-a_lj} = tau_{j-1} tau_{j+1} / tau_j^2 for a_n
        V = -(mm**2/beta**2)*sum(T[(j-1) % h]*T[(j+1) % h]/T[j]**2 - 1 for j in range(h))
        return sum(d*d for d in dP)/2 + V
    L = 60/ma
    E = mp.quad(dens, mp.linspace(-L, L, 13))
    return E, 2*h*ma/beta**2


ok_mass = True
for n, a in [(2, 1), (3, 1), (3, 2)]:
    E, pred = soliton_mass(n, a)
    ok_mass &= abs(E - pred) < 1e-12
    print(f'   a_{n}^(1) species {a}: E = {mp.nstr(E, 15)}, 2 h m_a/beta^2 = {mp.nstr(pred, 15)}')
report('static soliton energy = 2 h m_a / beta^2 (real) for complex xi', ok_mass)

# Theta on Hirota data: conj(tau_j(a, theta, xi)) = tau_j(h-a, conj theta, conj xi)
n, a = 3, 1; h = 4
t1 = tau_funcs(n, [a], [mp.mpc('0.2', '0.1')], [mp.mpc('0.3', '0.5')])(mp.mpf('0.7'), mp.mpf('0.2'))
t2 = tau_funcs(n, [h-a], [mp.mpc('0.2', '-0.1')], [mp.mpc('0.3', '-0.5')])(mp.mpf('0.7'), mp.mpf('0.2'))
report('Theta: conj tau_j(a,theta,xi) = tau_j(h-a, conj theta, conj xi)', max(abs(mp.conj(x)-y) for x, y in zip(t1, t2)) < 1e-25)

# Sine-Gordon breather data (exr-breather-real): xi2 = conj(xi1)+i pi gives a real field,
# xi2 = conj(xi1) a purely imaginary one (mod the branch of the logarithm).
def sg_two(x, t, u, xi1, xi2):
    th1, th2 = 1j*u, -1j*u
    E1 = mp.e**(2*(x*mp.cosh(th1) - t*mp.sinh(th1)) + xi1); E2 = mp.e**(2*(x*mp.cosh(th2) - t*mp.sinh(th2)) + xi2)
    A = (mp.cosh(th1 - th2) - 1)/(mp.cosh(th1 - th2) + 1)
    tau = [1 + (-1)**j*(E1 + E2) + A*E1*E2 for j in (0, 1)]
    return 1j*mp.log(tau[1]/tau[0])


pts = [(0.3, 0.1), (-0.5, 0.9), (1.2, -0.4)]
xi1 = mp.mpc(0.3, 0.4)
real_ok = all(abs(mp.im(sg_two(x, t, 0.7, xi1, mp.conj(xi1) + 1j*mp.pi))) < 1e-14 for x, t in pts)
imag_ok = all(min(abs(mp.re(sg_two(x, t, 0.7, xi1, mp.conj(xi1))) - k*mp.pi) for k in (-1, 0, 1)) < 1e-14 for x, t in pts)
report('sine-Gordon breather: xi2=conj(xi1)+i pi real field; xi2=conj(xi1) imaginary field (mod pi)', real_ok and imag_ok)
