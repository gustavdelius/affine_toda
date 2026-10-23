"""The Sklyanin relation of the a_1^(1) defect matrix by hand (Part VIII, ch. 38: prp-defect-involution).

Written from the book's definitions, without defect_common.py or defect_lax_rmatrix.py.
Conventions: a_1^(1) in the two-dimensional representation, H = diag(1, -1)/sqrt2, sigma_x = e + f, Pi = e (x) f + f (x) e,
b = beta (real coupling) or i beta (imaginary coupling), fields phi (x < 0) and psi (x > 0) in the normalization
L = 1/2 (d phi)^2 - V of sec-lagrangians.  Type I_-: L_d = exp(-b H psi/2) (1 + (sigma/lambda) sigma_x) exp(b H phi/2)
= J + (sigma/lambda) O, with J = exp(b H (phi - psi)/2) and O = exp(-b H (phi + psi)/2) sigma_x (eq-defect-lax-matrix with
D = 1, E = 0).  Defect bracket {phi(0), psi(0)} = eps, with eps = -1 for type I_- and +1 for type I_+ (eq-dc-defect-bracket).
Book r-matrix (prp-sklyanin) for h = 2: r_book(z) = b^2/(4(z^2 - 1)) [(z^2 + 1) H (x) H + 2 z Pi].

Checks (book ids they support):
 1. The general formula eq-defect-lax-matrix in h = 2 components (D = 2(1 - pi)^-1 = 1, E = 0, Ltilde = 1 + (sigma/lambda) S^T)
    equals J + (sigma/lambda) O; r_book built from the general formula of prp-sklyanin equals the h = 2 form above.
 2. Brackets (chain rule, sympy): {J_1, O_2} = k H(x)H J_1 O_2, {O_1, J_2} = -k H(x)H O_1 J_2, {J_1, J_2} = {O_1, O_2} = 0,
    with k = kappa = -b^2 eps/2; hence {L_1, L_2} = k sigma H(x)H (X/l2 - X'/l1), X = J (x) O, X' = O (x) J.
 3. Commutators (sympy): H(x)H and Pi commute with J (x) J and O (x) O; [H(x)H, X] = 2 H(x)H X, [Pi, X] = -2 H(x)H X',
    [Pi, X'] = -2 H(x)H X.  Control: [Pi, X] + 2 H(x)H X' is not zero with the sign flipped.
 4. Matching (sympy): for r = A H(x)H + G Pi the relation {L_1, L_2} = [r(l1/l2), L_1 L_2] holds iff 2A - 2G/z = k and
    2A - 2Gz = -k; the solution is r = k/(2(z^2 - 1)) [(z^2 + 1) H(x)H + 2z Pi] = -eps r_book, so for eps = -1 the bracket of the
    defect term gives exactly the bulk r-matrix.  Control: eps = +1 gives -r_book.
 5. Uniqueness without an ansatz (numerical, both couplings, complex l1, l2): the 16 entries of a matrix r with
    {L_1, L_2} = [r, L_1 L_2] at random defect values form r_book + C 1; the commutant of L_d(l1) (x) L_d(l2) over all defect
    values is one-dimensional (scalars) for l1^2 != l2^2.  Controls: two-dimensional for l2 = l1 and for l2 = -l1.
 6. Type I_+ (eq-defect-lax-matrix with D -> -D, Ltilde_I+ = 1 + (lambda/sigma) S): L_d^{I+}(phi, psi) sigma_x
    = (lambda/sigma) L_d^{I-}(phi, -psi) (sympy); with eps = +1 it satisfies the relation with r_book (numerical, both couplings);
    with eps = -1 it does not.
 7. Bulk (sympy): {a_x(l1; x) (x) a_x(l2; y)} = delta(x - y) [r_book(l1/l2), a_x(l1) (x) 1 + 1 (x) a_x(l2)] with
    {phi(x), phi_t(y)} = delta(x - y) and a_x of eq-boundary-lax-at; so the one-site transport T = 1 - eps a_x with
    {phi, phi_t} = 1/eps has {T_1, T_2} = -[r_book, T_1 T_2] + O(eps^2).  Control: -r_book fails.
Runtime: a few seconds.
"""
import numpy as np
import sympy as sp


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def kron(A, B):
    A, B = sp.Matrix(A), sp.Matrix(B)
    return sp.Matrix(A.rows*B.rows, A.cols*B.cols,
                     lambda i, j: A[i // B.rows, j // B.cols]*B[i % B.rows, j % B.cols])


def zero(M):
    return all(sp.simplify(x) == 0 for x in M)


phi, psi, b, sig, l1, l2, c, z, m = sp.symbols('phi psi b sigma lambda1 lambda2 epsilon z m')   # c = eps
s2 = sp.sqrt(2)
h = sp.diag(1, -1)
H = h/s2
e = sp.Matrix([[0, 1], [0, 0]])
f = sp.Matrix([[0, 0], [1, 0]])
sx = e + f
I2 = sp.eye(2)
HH = kron(H, H)
Pi = kron(e, f) + kron(f, e)


def expH(t):
    """exp(t H) for the diagonal H."""
    return sp.diag(sp.exp(t/s2), sp.exp(-t/s2))


J = expH(b*(phi - psi)/2)
O = expH(-b*(phi + psi)/2)*sx


def Lm(lam, ph=phi, ps=psi):
    return expH(-b*ps/2)*(I2 + sig/lam*sx)*expH(b*ph/2)


def r_book(zz):
    return b**2/(4*(zz**2 - 1))*((zz**2 + 1)*HH + 2*zz*Pi)


def pb(A, B, cc):
    """{A (x) B} for matrices of functions of (phi, psi) with {phi, psi} = cc."""
    return cc*(kron(A.diff(phi), B.diff(psi)) - kron(A.diff(psi), B.diff(phi)))


# 1. the general formula in h = 2 components, and the general r-matrix for h = 2
# components: phi_vec = phi mu with mu = (1, -1)/sqrt2, H.v = diag(v_0, v_1), E_alpha_j = e_{j, j+1} (indices mod 2)
mu = sp.Matrix([1, -1])/s2
Pm = sp.Matrix([[0, 1], [1, 0]])            # pi e_k = e_{k-1}: on R^2 it swaps e_0 and e_1
Dmat = 2*(sp.eye(2) - Pm).pinv()             # 2 (1 - pi)^-1 on the hyperplane
Proj = (sp.eye(2) - Pm).pinv()*(sp.eye(2) - Pm)   # projector onto the hyperplane
Emat = Proj - Dmat                           # E = 1 - D on the hyperplane
S = sp.Matrix([[0, 1], [1, 0]])              # sum_j E_alpha_j = e_01 + e_10
u = Emat*(phi*mu) + Dmat*(psi*mu)
w = Dmat.T*(phi*mu) - Emat*(psi*mu)
Lgen = sp.diag(*[sp.exp(-b*x/2) for x in u]) * (I2 + sig/l1*S.T) * sp.diag(*[sp.exp(b*x/2) for x in w])
ok1a = zero(Lgen - (J + sig/l1*O)) and zero(Lm(l1) - (J + sig/l1*O)) and zero(Dmat*mu - mu)
# r from the general formula: orthonormal Cartan basis of sl_2 (one element, H), positive root e_01 of height 1, alpha^2 = 2
rgen = b**2/(4*(z**2 - 1))*((z**2 + 1)*HH + 2*(z*kron(e, f) + z**(2 - 1)*kron(f, e)))
ok1b = zero(rgen - r_book(z))
report('1. L_d of eq-defect-lax-matrix for h = 2 is J + (sigma/lambda) O; r of prp-sklyanin for h = 2 is '
       'b^2/(4(z^2-1))[(z^2+1) H(x)H + 2z Pi]', ok1a and ok1b)

# 2. brackets of D and O
k = -b**2*c/2
X = kron(J, O)
Xp = kron(O, J)
ok2 = (zero(pb(J, O, c) - k*HH*X) and zero(pb(O, J, c) + k*HH*Xp)
       and zero(pb(J, J, c)) and zero(pb(O, O, c)))
L1, L2 = Lm(l1), Lm(l2)
ok2b = zero(pb(L1, L2, c) - k*sig*HH*(X/l2 - Xp/l1))
report('2. {J_1,O_2} = k H(x)H J_1O_2, {O_1,J_2} = -k H(x)H O_1J_2, k = -b^2 eps/2; {L_1,L_2} = k sigma H(x)H (X/l2 - X\'/l1)',
       ok2 and ok2b)

# 3. commutators
def comm(A, B):
    return A*B - B*A


DD, OO = kron(J, J), kron(O, O)
ok3 = (zero(comm(HH, DD)) and zero(comm(Pi, DD)) and zero(comm(HH, OO)) and zero(comm(Pi, OO))
       and zero(comm(HH, X) - 2*HH*X) and zero(comm(HH, Xp) - 2*HH*Xp)
       and zero(comm(Pi, X) + 2*HH*Xp) and zero(comm(Pi, Xp) + 2*HH*X))
ctrl3 = not zero(comm(Pi, X) - 2*HH*Xp)
report('3. H(x)H, Pi commute with J(x)J, O(x)O; [H(x)H, X] = 2 H(x)H X; [Pi, X] = -2 H(x)H X\', [Pi, X\'] = -2 H(x)H X '
       '(control: the opposite sign fails)', ok3 and ctrl3)

# 4. matching
A, G = sp.symbols('A G')
rr = A*HH + G*Pi
resid = pb(L1, L2, c) - comm(rr, kron(L1, L2))
# by checks 2-3 the residual is sigma H(x)H [(k/l2 - 2A/l2 + 2G/l1) X + (-k/l1 - 2A/l1 + 2G/l2) X']
pred = sig*HH*((k/l2 - 2*A/l2 + 2*G/l1)*X + (-k/l1 - 2*A/l1 + 2*G/l2)*Xp)
ok4a = zero(resid - pred)
sol = sp.solve([2*A - 2*G/z - k, 2*A - 2*G*z + k], [A, G], dict=True)
ok4b = len(sol) == 1
rsol = (sol[0][A]*HH + sol[0][G]*Pi)
ok4c = zero(rsol - (-c)*r_book(z))
resid_m1 = (pb(L1, L2, -1) - comm(r_book(l1/l2), kron(L1, L2)))
resid_p1 = (pb(L1, L2, 1) - comm(r_book(l1/l2), kron(L1, L2)))
ok4d = zero(resid_m1) and not zero(resid_p1)
report('4. matching gives 2A - 2G/z = k, 2A - 2Gz = -k with the unique solution -eps r_book; for eps = -1 {L_1,L_2} = '
       '[r_book, L_1 L_2] exactly (control: eps = +1 fails)', ok4a and ok4b and ok4c and ok4d)

# 5. uniqueness without an ansatz, numerically
rng = np.random.default_rng(5)
fL = sp.lambdify((phi, psi, b, sig, l1), Lm(l1), 'numpy')
fPB = sp.lambdify((phi, psi, b, sig, l1, l2, c), pb(L1, L2, c), 'numpy')
fr = sp.lambdify((b, z), r_book(z), 'numpy')


def comm_map(M):
    """Matrix of r -> [r, M] acting on the row-major vec(r)."""
    I4 = np.eye(4)
    return np.kron(I4, M.T) - np.kron(M, I4)


def stack(bb, s, la1, la2, cc, npts=6, Lfun=None):
    rows, rhs = [], []
    for _ in range(npts):
        ph, ps = rng.normal(size=2) + 1j*rng.normal(size=2)
        LL = np.kron(Lfun(ph, ps, bb, s, la1), Lfun(ph, ps, bb, s, la2))
        rows.append(comm_map(LL))
        rhs.append(np.array(fPB(ph, ps, bb, s, la1, la2, cc), dtype=complex).ravel())
    return np.vstack(rows), np.concatenate(rhs)


ok5 = True
for bb in [0.83, 0.83j]:
    s = 0.7 + 0.2j
    la1, la2 = 1.3 + 0.4j, -0.6 + 0.9j
    Mt, rhs = stack(bb, s, la1, la2, -1, Lfun=fL)
    sol5, *_ = np.linalg.lstsq(Mt, rhs, rcond=None)
    sv = np.linalg.svd(Mt, compute_uv=False)
    nullity = int(np.sum(sv < 1e-10*sv[0]))
    diff = sol5.reshape(4, 4) - fr(bb, la1/la2)
    # the difference must be a multiple of the identity
    ok5 &= (np.linalg.norm(Mt @ sol5 - rhs) < 1e-10*np.linalg.norm(rhs) and nullity == 1
            and np.linalg.norm(diff - np.trace(diff)/4*np.eye(4)) < 1e-10)
ctrl5 = True
for la2 in [1.3 + 0.4j, -(1.3 + 0.4j)]:
    Mt, _ = stack(0.83, 0.7 + 0.2j, 1.3 + 0.4j, la2, -1, Lfun=fL)
    sv = np.linalg.svd(Mt, compute_uv=False)
    ctrl5 &= int(np.sum(sv < 1e-10*sv[0])) == 2
report('5. without ansatz, {L_1,L_2} = [r, L_1 L_2] fixes r = r_book + C 1; the commutant of L_d(l1)(x)L_d(l2) is the scalars '
       '(controls: dimension 2 at l2 = +-l1), both couplings', ok5 and ctrl5)

# 6. type I_+
def Lp(lam, ph=phi, ps=psi):
    # eq-defect-lax-matrix with (D, E) -> (-1, 0): exp(+b H psi/2) Ltilde_I+ exp(-b H phi/2), Ltilde_I+ = 1 + (lambda/sigma) S
    return expH(b*ps/2)*(I2 + lam/sig*S)*expH(-b*ph/2)


ok6a = zero(Lp(l1)*sx - l1/sig*Lm(l1, phi, -psi))
fLp = sp.lambdify((phi, psi, b, sig, l1), Lp(l1), 'numpy')
P1, P2 = Lp(l1), Lp(l2)
fPBp = sp.lambdify((phi, psi, b, sig, l1, l2, c), pb(P1, P2, c), 'numpy')


def resid6(bb, cc):
    s, la1, la2 = 0.7 + 0.2j, 1.3 + 0.4j, -0.6 + 0.9j
    res = 0.0
    for _ in range(4):
        ph, ps = rng.normal(size=2) + 1j*rng.normal(size=2)
        LL = np.kron(fLp(ph, ps, bb, s, la1), fLp(ph, ps, bb, s, la2))
        R = fr(bb, la1/la2)
        lhs = np.array(fPBp(ph, ps, bb, s, la1, la2, cc), dtype=complex)
        res = max(res, np.linalg.norm(lhs - (R @ LL - LL @ R))/np.linalg.norm(lhs))
    return res


ok6b = all(resid6(bb, 1) < 1e-12 for bb in [0.83, 0.83j]) and resid6(0.83, -1) > 0.1
report('6. type I_+: L_d^{I+}(phi,psi) sigma_x = (lambda/sigma) L_d^{I-}(phi,-psi); with {phi,psi} = +1 it satisfies the '
       'relation with r_book (control: {phi,psi} = -1 fails)', ok6a and ok6b)

# 7. the bulk bracket of a_x for a_1 in the book's normalization
phit = sp.symbols('phi_t')


def a_x(lam):
    # eq-boundary-lax-at with c_j = 1: alpha_1 = sqrt2 (E_alpha1 = e), alpha_0 = -sqrt2 (E_alpha0 = f), loop variable 1
    return (b/2*phit*H + m/2*(sp.exp(b*s2*phi/2)*(lam*e + f/lam) + sp.exp(-b*s2*phi/2)*(lam*f + e/lam)))


A1, A2 = a_x(l1), a_x(l2)
bulk_pb = kron(A1.diff(phi), A2.diff(phit)) - kron(A1.diff(phit), A2.diff(phi))   # {phi(x), phi_t(y)} = delta
rhs7 = comm(r_book(l1/l2), kron(A1, I2) + kron(I2, A2))
ok7 = zero(bulk_pb - rhs7) and not zero(bulk_pb + rhs7)
report('7. {a_x(l1) (x) a_x(l2)} = delta [r_book(l1/l2), a_x(l1)(x)1 + 1(x)a_x(l2)] for a_1 (control: -r_book fails); '
       'hence {T_1,T_2} = -[r_book, T_1T_2] for T = 1 - eps a_x', ok7)
