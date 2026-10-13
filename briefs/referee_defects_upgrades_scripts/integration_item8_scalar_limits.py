"""Item 8 (ch. 42, prp-diagonal-borel (5)): the scalar factor Phi of eq-sg-transmission-matrix does not cancel
from the limit statement.  T = Phi * M with M the matrix part (ch. 37); L = P T.  Built from the chapter's
definitions; imports only the book helper sg_transmission_scalar (eq-sgd-scalar).
Checks: (a) |Phi(D)|^2 (1 + Y^2) = 1, Y = e^{lam D}; (b) as theta_d -> +inf (D -> -inf) Phi -> 1, so the full
L tends to diag(q^-n, q^(n+s)) as stated; (c) as theta_d -> -inf (D -> +inf) the entries of L / Y are
Phi*q^-n/Y (diagonal) and -i Phi q^-n (converting): both tend to 0, so 'L/(y/y_d) -> conjugation' holds only for
the matrix part P M; (d) the matrix part's converting entries / Y are theta-independent (-i q^-n, -i q^(n+s));
(e) the intertwining equations L X_V = X_V' L are homogeneous of degree 1 in L, so Phi cancels (sympy check
on the e_1 equation of the untwisted diagonal with symbolic scalar)."""
import sys
import numpy as np
import sympy as sp
sys.path.insert(0, '/home/gustav/Git/affine_toda/.claude/worktrees/integrable-defects-framework-a54a84/textbook_code')
from defect_common import sg_transmission_scalar

lam = 1.37
q = -np.exp(-1j*np.pi*lam)
ok = True
print('(a) |Phi|^2 (1+Y^2):')
for D in (-3.0, -1.0, 0.0, 0.7, 2.0, 4.0):
    P = sg_transmission_scalar(D, lam)
    v = abs(P)**2*(1 + np.exp(2*lam*D))
    print(f'   D = {D:5.1f}: {v:.12f}')
    ok &= abs(v - 1) < 1e-8
print('(b) Phi -> 1 as D -> -inf:')
for D in (-2.0, -4.0, -6.0):
    P = sg_transmission_scalar(D, lam)
    print(f'   D = {D:5.1f}: |Phi - 1| = {abs(P - 1):.2e}')
ok &= abs(sg_transmission_scalar(-6.0, lam) - 1) < 1.2*np.exp(-6.0)   # corrections of order e^{min(1,lam) D}
print('(c) entries of L/Y as D -> +inf (n = 0, s = 0.3):')
s, n = 0.3, 0
prev = None
for D in (1.0, 2.0, 4.0, 6.0):
    P = sg_transmission_scalar(D, lam)
    Y = np.exp(lam*D)
    diag = abs(P*q**(-n)/Y)
    conv = abs(-1j*P*q**(-n))
    print(f'   D = {D:4.1f}: |diag/Y| = {diag:.3e}, |converting/Y| = {conv:.3e}  (1/sqrt(1+Y^2) = {1/np.sqrt(1+Y*Y):.3e})')
    prev = conv
ok &= prev < 1e-3
print('   while the matrix part P M / Y has converting entries -i q^-n, -i q^(n+s) of modulus', abs(q**(-n)), abs(q**(n+s)))
# (e) homogeneity: e_1 equation of (1 (x) t)Delta(U) on a 2x2 block model with a scalar phi
A, B, C, Dd, Yop, K = sp.symbols('A B C D Y kappa', commutative=False)
phi, y, t, qs = sp.symbols('phi y t q', nonzero=True)
Ep = sp.Matrix([[0, 1], [0, 0]])
L = phi*(sp.Matrix([[A, 0], [0, Dd]]) + y*sp.Matrix([[0, B], [C, 0]]))
kV = sp.Matrix([[qs, 0], [0, 1/qs]])
XV = y*Ep*sp.eye(2) + kV*Yop        # pi(e_1) (x) 1 + pi(k_1) (x) Y, operators on W as commuting-with-V symbols
XVp = Yop*sp.eye(2) + t*y*Ep*K
lhs = sp.expand(L*XV - XVp*L)
lhs0 = sp.expand((L/phi)*XV - XVp*(L/phi))
ok_e = all(sp.expand(lhs[i, j] - phi*lhs0[i, j]) == 0 for i in range(2) for j in range(2))
print('(e) the e_1 equation with T = phi*M equals phi times the equation for M:', ok_e)
ok &= ok_e
print('PASS' if ok else 'FAIL', ': item 8 - the limit in (5) holds for the matrix part P M; with Phi, L/(y/y_d) -> 0')
