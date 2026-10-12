"""Open problem 3 (Neumann from the soliton side): consolidated checks.  `python3 check_all.py [full]`
(full adds the n = 5 soliton bootstrap, ~3 min).  Book conventions throughout (boundary_common.py).

Results checked:
 1. Phi_n:  k(th) = Phi(th) R_{Vbar Vbar}(2th) k(-th) with Phi_+ = sinh(i pi om/2 - T th)/sinh(i pi om/2 + T th),
            Phi_- = cosh(T th - i pi om/2)/cosh(T th + i pi om/2)                                   (n = 1..5)
 2. U_n:    K(-th)^T K(th) = U(th) 1,  U = (s cosh 2T th - cos pi T)/(2 sin^2(pi om/2));  K_W = J K J = K^T is the
            antisoliton intertwiner Vbar(th) -> V(-th)                                                 (n = 1..5)
 3. rho_+-: unitarity and crossing (soliton-level functional equations)                              (n = 1..4)
 4. singlet projection: (1(x)K)R(1(x)K^T) maps the B_1 singlet to c_s(mu) times itself,
            c_+ = -cos(pi mu)cos(pi(mu+om/2))cos(pi(mu-T))/(sin^2(pi om/2) cos(pi(mu-om/2))),
            c_- = sin(pi mu)sin(pi(mu+om/2))sin(pi(mu-T))/(sin^2(pi om/2) sin(pi(mu-om/2)))           (n = 2..5)
 5. n = 1: A_d rho_+ sigma = GZ/BPTT scalar at the Neumann point; -A_d rho_- = GZ at (eta, vartheta) = (pi/2, i pi om/2)
 6. n = 2: A_d rho_+ sigma = Gandenberger (4.14) x sigma^+ (5.3)
 7. soliton bootstrap V^(x n) -> Vbar: Y(mu) = P(mu) prod_k Y(mu + T - k om) for Y = rho_+ sigma, -rho_- (n = 2..4[5])
 8. breather: A(mu_a) F11(2mu) A(mu_b) c_s(mu) = (-n om/2)(n om/2-1/2)(1/2) [s=+1], DG(3.11) sigma_1 [s=-1] (n=1..5)
 9. semiclassics of A^(+): (beta^2/4h) d(-i ln A)/dth -> (1/2) ln(1 - sin^2(pi/h)/cosh^2 th)
"""
import sys
import numpy as np
from boundary_common import q_phys, an_vector, k_intertwiner, report
from nk_common import k_book, fast_rcheck
from proj import cproj
from rhonp import rho_plus, blk, cproj_plus
from rhom import rho_minus, cproj_minus
from gz import A1, logR0, logsigma, rich
from scipy.special import loggamma as lg

FULL = len(sys.argv) > 1 and sys.argv[1] == 'full'


def prop(A, B):
    c = np.vdot(B.ravel(), A.ravel())/np.vdot(B.ravel(), B.ravel())
    return np.linalg.norm(A - c*B)/np.linalg.norm(A)


def Phif(n, s, om, th):
    T = (n + 1)*om/2; a = 0.5j*np.pi*om
    return -np.sinh(T*th - a)/np.sinh(T*th + a) if s == 1 else np.cosh(T*th - a)/np.cosh(T*th + a)


def Uf(n, s, om, th):
    T = (n + 1)*om/2
    return (s*np.cosh(2*T*th) - np.cos(np.pi*T))/(2*np.sin(np.pi*om/2)**2)


def sigma(mu, n, om):
    T = (n + 1)*om/2
    return blk(T/2, mu, T)**2*blk(T/2 - om/2, mu, T)*blk(T/2 + om/2, mu, T)


def A_plus(mu, n, om):
    return A1(mu, n, om)*rho_plus(mu, n, om)*sigma(mu, n, om)


def A_minus(mu, n, om):
    return -A1(mu, n, om)*rho_minus(mu, n, om)


def F11(mu, h, lam, K=4000):
    def part(KK):
        j = np.arange(1, KK + 1); a = h*j*lam
        return np.sum(lg(mu + a - (h-1)*lam + 1) + lg(mu + a - lam) + lg(-mu + a - h*lam + 1) + lg(-mu + a)
                      - lg(-mu + a - (h-1)*lam + 1) - lg(-mu + a - lam) - lg(mu + a - h*lam + 1) - lg(mu + a))
    return -np.exp(rich(part, K))


# 1, 2 -------------------------------------------------------------------------------------------
wP = wU = wI = 0
for n in (1, 2, 3, 4, 5):
    N = n + 1; J = np.eye(N)[::-1]
    for s in (1, -1):
        for om in (0.37, 1.83, 2.6):
            q = q_phys(om); w = np.exp(-0.5j*np.pi*om); e = np.ones(N)/(w + 1/w); e[0] *= s
            for th in (0.31 + 0.2j, -0.7 + 0.1j):
                km = k_book(n, om, 0.5j*np.pi - th, s).T.reshape(N*N)
                kp = k_book(n, om, 0.5j*np.pi + th, s).T.reshape(N*N)
                R2 = fast_rcheck(an_vector(n, q, np.exp(2*om*th), True), an_vector(n, q, 1.0, True))
                v = R2 @ kp; c = np.vdot(v, km)/np.vdot(v, v)
                wP = max(wP, abs(c/Phif(n, s, om, th) - 1), np.linalg.norm(km - c*v)/np.linalg.norm(km))
                P = k_book(n, om, -th, s).T @ k_book(n, om, th, s)
                wU = max(wU, np.abs(P - Uf(n, s, om, th)*np.eye(N)).max()/abs(Uf(n, s, om, th)))
                y = np.exp(om*th)
                null, _ = k_intertwiner(an_vector(n, q, y, True), an_vector(n, q, 1/y), [1/q]*N, e)
                wI = max(wI, prop(null[0], J @ k_book(n, om, th, s) @ J), prop(J @ k_book(n, om, th, s) @ J,
                                                                               k_book(n, om, th, s).T))
report(f'1. Phi_n closed form, n = 1..5, s = +-1, three omega [{wP:.1e}]', wP < 1e-10)
report(f'2. U_n closed form for K(-th)^T K(th) [{wU:.1e}]; K_W = J K J = K^T is the Vbar intertwiner [{wI:.1e}]',
       wU < 1e-12 and wI < 1e-12)

# 3 ----------------------------------------------------------------------------------------------
w3 = 0
for n in (1, 2, 3, 4):
    om = 0.37; T = (n + 1)*om/2; a, b = T/2, om/2; mu = 0.13 + 0.21j
    r = lambda m: rho_plus(m, n, om)
    w3 = max(w3, abs(r(mu)*r(-mu)*np.sin(np.pi*(a + mu))*np.sin(np.pi*(a - mu))/np.sin(np.pi*b)**2 - 1),
             abs(r(a - mu)/r(a + mu)*np.sin(np.pi*(b - mu))/np.sin(np.pi*(b + mu)) - 1))
    r = lambda m: rho_minus(m, n, om)
    w3 = max(w3, abs(-r(mu)*r(-mu)*np.cos(np.pi*(a + mu))*np.cos(np.pi*(a - mu))/np.sin(np.pi*b)**2 - 1),
             abs(r(a - mu)/r(a + mu)*np.cos(np.pi*(mu - b))/np.cos(np.pi*(mu + b)) - 1))
report(f'3. rho_+-: rho(mu)rho(-mu) = 1/U, rho(T/2-mu)/rho(T/2+mu) = 1/Phi, n = 1..4 [{w3:.1e}]', w3 < 1e-7)

# 4 ----------------------------------------------------------------------------------------------
w4 = 0
for n in (2, 3, 4, 5):
    for s in (1, -1):
        for om in (0.37, 1.83):
            T = (n + 1)*om/2
            for th in (0.31 + 0.07j, -0.5 + 0.3j):
                mu = -1j*T*th/np.pi
                for ca in (False, True):
                    c, r, rs = cproj(n, s, om, th, ca=ca)
                    cf = cproj_plus(mu, om, T) if s == 1 else cproj_minus(mu, om, T)
                    w4 = max(w4, r, abs(c/cf - 1))
report(f'4. singlet -> c_s(mu) x singlet, both orderings, closed forms c_+-, n = 2..5 [{w4:.1e}]', w4 < 1e-10)

# 5 ----------------------------------------------------------------------------------------------
lam = 0.37; w5 = 0
for th in (0.3 + 0.2j, -0.5 + 0.4j):
    u = -1j*th; mu = -1j*lam*th/np.pi
    AG = np.exp(logR0(u, lam) + logsigma(np.pi*(lam + 1)/2, u, lam) + logsigma(0, u, lam))*np.cos(lam*u)
    w5 = max(w5, abs(AG/A_plus(mu, 1, lam) - 1))
    eta = np.pi/2 - 1e-8; vt = 0.5j*np.pi*lam
    AGm = (np.exp(logR0(u, lam) + logsigma(eta, u, lam) + logsigma(1j*vt, u, lam))/(np.cos(eta)*np.cosh(vt))
           * 1j*np.sin(lam*u)*np.sin(np.pi*lam/2))
    w5 = max(w5, abs(AGm/A_minus(mu, 1, lam) - 1))
report(f'5. n = 1: A^(+) = GZ/BPTT at the Neumann point, A^(-) = GZ at (pi/2, i pi om/2) [{w5:.1e}]', w5 < 1e-6)

# 6 ----------------------------------------------------------------------------------------------
lam = 0.37


def G(j, m, a):
    return lg(m + 3*j*lam + a) - lg(-m + 3*j*lam + a)


def A_G(mu, K=4000):
    def part(KK):
        j = np.arange(1, KK + 1)
        return (np.sum(G(j, mu, -7*lam/4) + G(j, mu, -5*lam/4 + 1) - G(j, mu, -lam/4) - G(j, mu, -11*lam/4 + 1))
                + np.sum(G(2*j, 2*mu, -1.5*lam) + G(2*j, 2*mu, -4.5*lam + 1) - G(2*j, 2*mu, -2.5*lam)
                         - G(2*j, 2*mu, -3.5*lam + 1)))
    return np.sin(np.pi*lam/2)/np.sin(np.pi*(mu - 0.75*lam))*np.exp(rich(part, K))


T2 = 1.5*lam
w6 = max(abs(A_G(mu)*blk(lam/4, mu, T2)*blk(5*lam/4, mu, T2)/A_plus(mu, 2, lam) - 1)
         for mu in (0.13 + 0.21j, -0.31 + 0.1j))
report(f"6. n = 2: A^(+) = Gandenberger's A^+ (4.14) times sigma^+ = (lam/4)(5lam/4) (5.3) [{w6:.1e}]", w6 < 1e-7)

# 7 ----------------------------------------------------------------------------------------------
from sol_boot import boot_coef
w7 = 0
for n in ((2, 3, 4, 5) if FULL else (2, 3, 4)):
    om = 0.37; T = (n + 1)*om/2; th = 0.31 + 0.07j; mu = -1j*T*th/np.pi
    mus = [mu + T - k*om for k in range(1, n + 1)]
    cd, rd = boot_coef(n, om, th, 0, order=-1)
    for s, Y in ((1, lambda m: rho_plus(m, n, om)*sigma(m, n, om)), (-1, lambda m: -rho_minus(m, n, om))):
        cs, rs = boot_coef(n, om, th, s, order=-1)
        w7 = max(w7, rd, rs, abs(Y(mu)/(cs/cd*np.prod([Y(m) for m in mus])) - 1))
report(f'7. soliton bootstrap V^(x n) -> Vbar (image = iota_out K^T), A^(+) = A_d rho_+ sigma and '
       f'A^(-) = -A_d rho_-, n = 2..{5 if FULL else 4} [{w7:.1e}]', w7 < 1e-7)

# 8 ----------------------------------------------------------------------------------------------
w8 = 0
for om in (0.37, 0.83):
    for n in (1, 2, 3, 4, 5):
        h = n + 1; T = h*om/2
        for mu in (0.13 + 0.21j, -0.31 + 0.1j):
            ma, mb = mu - T/2 + 0.5, mu + T/2 - 0.5
            KBp = A_plus(ma, n, om)*F11(2*mu, h, om)*A_plus(mb, n, om)*cproj_plus(mu, om, T)
            KBm = A_minus(ma, n, om)*F11(2*mu, h, om)*A_minus(mb, n, om)*cproj_minus(mu, om, T)
            dg311 = blk(-n*om/2, mu, T)*blk(0.5 - T, mu, T)*blk(-om/2 - 0.5, mu, T)
            w8 = max(w8, abs(KBp/(blk(-n*om/2, mu, T)*blk(n*om/2 - 0.5, mu, T)*blk(0.5, mu, T)) - 1),
                     abs(KBm/(dg311*(-blk(om/2, mu, T)*blk(T - om/2, mu, T))) - 1))
report(f'8. breather: K_B1^(+) = (-n om/2)(n om/2 - 1/2)(1/2), K_B1^(-) = DG(3.11) x sigma_1, n = 1..5, '
       f'two omega [{w8:.1e}]', w8 < 1e-6)

# 9 ----------------------------------------------------------------------------------------------
errs = []
for n in (2, 3):
    h = n + 1; th = 0.6; s_ = np.sin(np.pi/h)/np.cosh(th); tgt = 0.5*np.log(1 - s_**2); e_n = []
    for L in (40., 80., 160.):
        T = h*L/2
        lnA = lambda t: np.log(A1(-1j*T*t/np.pi, n, L, K=3000)*rho_plus(-1j*T*t/np.pi, n, L, K=3000)
                               * sigma(-1j*T*t/np.pi + 0j, n, L))
        d = (lnA(th + 1e-5) - lnA(th - 1e-5))/2e-5
        e_n.append(abs(np.real(-1j*d)*4*np.pi/(L + 1)/(4*h) - tgt))
    errs.append(e_n)
report('9. semiclassics: (beta^2/4h) d(-i ln A^(+))/dth -> (1/2) ln(1 - s^2) (Neumann time advance), errors '
       + str([[f'{x:.4f}' for x in e] for e in errs]) + ' at om = 40, 80, 160',
       all(e[2] < e[1] < e[0] and e[2] < 0.03 for e in errs))
