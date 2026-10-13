"""Referee check (2026-10-12), item 1 of halfline_open_06_01_opus_prompt.md.
a_5 (h=6), species-2 stationary pair: tau_j = 1 + 2d cos(chi + 2 pi j/3) e^{sqrt3 x} + (d^2/4) e^{2 sqrt3 x}, d = 4-2 sqrt3.
At chi_end = 5pi/3, tau_2(0) = tau_5(0) = 0 (simple zeros in x), the singular 2 2bar background of BP.
Approach chi_end from both sides with regular members of the family, and compute the linearized reflection matrix
with the boundary Hessian H = 1/2 sum_j C_j e^{psi_j(0)/2} alpha_j alpha_j^T of a TARGET pattern C:
  (i) my own solver (written from scratch: own background, own basis, own channel/BC assembly; DOP853 and Radau),
  (ii) the worker's linrefl2.refl_matrix2 (imported read-only from the repository).
Book conventions: phi = -sum alpha_j ln tau_j, psi_j = alpha_j.phi, M^2 = sum e^{psi_j} alpha_j alpha_j^T,
BC phi'(0) = -sum C_j alpha_j e^{psi_j(0)/2}; linearized y'(0) = -H y(0).  K[b,a] = out-amplitude in channel b
for unit in-wave in channel a (channel a = Fourier vector omega^{ka}).  m=beta=1."""
import sys, itertools, numpy as np
from scipy.integrate import solve_ivp
h = 6; n = 5
s3 = np.sqrt(3.0); d = 4 - 2*s3; m2 = s3
ALPHA = np.array([np.eye(h)[j] - np.eye(h)[(j+1) % h] for j in range(h)])      # alpha_j, j=0..5
ma = np.array([2*np.sin(np.pi*a/h) for a in range(1, h)])

def tau_all(x, chi):
    j = np.arange(h)
    e = np.exp(m2*x)
    t = 1 + 2*d*np.cos(chi + 2*np.pi*j/3)*e + d*d/4*e*e
    tp = 2*d*np.cos(chi + 2*np.pi*j/3)*m2*e + d*d/4*2*m2*e*e
    tpp = 2*d*np.cos(chi + 2*np.pi*j/3)*m2**2*e + d*d/4*4*m2**2*e*e
    return t, tp, tpp

def psi(x, chi):
    t = tau_all(x, chi)[0]
    return np.log(np.roll(t, 1)*np.roll(t, -1)) - 2*np.log(t)

def bc_residual(chi, C):
    t, tp, _ = tau_all(0.0, chi)
    phid = -(ALPHA.T @ (tp/t))
    eh = np.sqrt(np.roll(t, 1)*np.roll(t, -1))/t
    return np.abs(phid + ALPHA.T @ (np.asarray(C)*eh)).max()

def exact_pattern(chi):
    best = min((bc_residual(chi, C), C) for C in itertools.product([1, -1], repeat=h))
    return best

def toda_residual(chi, x=-0.37, eps=1e-4):
    # static Toda: psi_j'' = 2 e^{psi_j} - e^{psi_{j-1}} - e^{psi_{j+1}}
    p = psi(x, chi); pp = (psi(x+eps, chi) - 2*p + psi(x-eps, chi))/eps**2
    e = np.exp(p)
    return np.abs(pp - (2*e - np.roll(e, 1) - np.roll(e, -1))).max()

# orthonormal basis of the sum-zero subspace: real Fourier basis (my own construction)
def real_basis():
    k = np.arange(h); cols = []
    for a in (1, 2):
        cols.append(np.cos(2*np.pi*a*k/h)); cols.append(np.sin(2*np.pi*a*k/h))
    cols.append((-1.0)**k)
    B = np.array(cols).T
    return B/np.linalg.norm(B, axis=0)
QB = real_basis()                                   # 6 x 5
AQ = ALPHA @ QB                                     # rows: alpha_j in the basis

def K_mine(chi, C, Om, X=16.0, method='DOP853', rtol=1e-12, atol=1e-14):
    w = np.exp(2j*np.pi/h); k = np.arange(h)
    V = np.array([QB.T @ (w**(k*a)/np.sqrt(h)) for a in range(1, h)]).T      # 5 x 5, column a-1 = channel a
    e0 = np.exp(psi(0.0, chi)/2)
    H = 0.5*(AQ.T*(np.asarray(C)*e0)) @ AQ
    def rhs(x, Y):
        y = Y[:n]; yp = Y[n:]
        M2 = (AQ.T*np.exp(psi(x, chi))) @ AQ
        return np.concatenate([yp, (M2 - Om**2*np.eye(n)) @ y])
    open_ = [a for a in range(n) if Om > ma[a]]; closed = [a for a in range(n) if Om <= ma[a]]
    res = {}
    for kind, lst in (('in', open_), ('out', open_), ('cl', closed)):
        for a in lst:
            if kind == 'cl':
                kap = np.sqrt(ma[a]**2 - Om**2); f, fp = np.exp(-kap*X), kap*np.exp(-kap*X)
            else:
                p = np.sqrt(Om**2 - ma[a]**2); sgn = 1 if kind == 'in' else -1
                f = np.exp(-1j*sgn*p*X); fp = 1j*sgn*p*f
            Y0 = np.concatenate([V[:, a]*f, V[:, a]*fp]).astype(complex)
            s = solve_ivp(rhs, [-X, 0.0], Y0, method=method, rtol=rtol, atol=atol)
            Y = s.y[:, -1]
            res[(kind, a)] = Y[n:] + H @ Y[:n]
    Rin = np.array([res[('in', a)] for a in open_]).T
    Rrest = np.array([res[('out', a)] for a in open_] + [res[('cl', a)] for a in closed]).T
    sol = np.linalg.solve(Rrest, -Rin)
    return sol[:len(open_), :], [a+1 for a in open_]

def blk(x, th): return np.sinh(th/2 + 1j*np.pi*x/(2*h))/np.sinh(th/2 - 1j*np.pi*x/(2*h))
def K_bp_diag(a, th):
    """BP Table 2, r=5, row {-++-++, --++++}: K_1=(2)(4)/((1)(5)), K_2=-(1)(3)^2(5)/((2)(4)), K_3=-(2)^2(4)^2/(3)^2"""
    if a in (1, 5): return blk(2, th)*blk(4, th)/(blk(1, th)*blk(5, th))
    if a in (2, 4): return -blk(1, th)*blk(3, th)**2*blk(5, th)/(blk(2, th)*blk(4, th))
    return -blk(2, th)**2*blk(4, th)**2/blk(3, th)**2

def s(C): return ''.join('+' if c > 0 else '-' for c in C)


"""(lib) Referee (2026-10-12): outward-integration solver for the linearized reflection matrix (independent of linrefl2):
start at the wall with the n solutions y(0)=u_k, y'(0)=-H u_k (BC imposed exactly), integrate to x=-X, decompose
z_c = v_c^dagger y into in/out (open) or growing/decaying (closed) modes, impose zero growing amplitude.
Background, basis, channels and BP values are taken from ref0601_a5_K.py (same directory)."""
import sys, time, numpy as np
from scipy.integrate import solve_ivp
def K_out(chi, C, Om, X=10.0, rtol=1e-12, atol=1e-14):
    w = np.exp(2j*np.pi/h); k = np.arange(h)
    V = np.array([QB.T @ (w**(k*a)/np.sqrt(h)) for a in range(1, h)]).T
    e0 = np.exp(psi(0.0, chi)/2); H = 0.5*(AQ.T*(np.asarray(C)*e0)) @ AQ
    def rhs(x, Y):
        y = Y[:n]; yp = Y[n:]
        return np.concatenate([yp, ((AQ.T*np.exp(psi(x, chi))) @ AQ - Om**2*np.eye(n)) @ y])
    sols = []
    for kk in range(n):
        u = np.eye(n)[kk]
        s = solve_ivp(rhs, [0.0, -X], np.concatenate([u, -H @ u]), method='DOP853', rtol=rtol, atol=atol)
        sols.append(s.y[:, -1])
    S = np.array(sols).T                      # 2n x n
    Z = V.conj().T @ S[:n]; Zp = V.conj().T @ S[n:]   # channel amplitudes at x=-X
    xX = -X; rows_in, rows_out, rows_grow = [], [], []
    open_ = [a for a in range(n) if Om > ma[a]]
    for c in range(n):
        if Om > ma[c]:
            p = np.sqrt(Om**2 - ma[c]**2)
            rows_in.append((Zp[c] + 1j*p*Z[c])/(2j*p)*np.exp(-1j*p*xX))      # coefficient of e^{ipx}
            rows_out.append((-Zp[c] + 1j*p*Z[c])/(2j*p)*np.exp(1j*p*xX))     # coefficient of e^{-ipx}
        else:
            kap = np.sqrt(ma[c]**2 - Om**2)
            rows_grow.append((Z[c] - Zp[c]/kap)/2*np.exp(kap*xX))               # coefficient of e^{-kappa x}
    Ain = np.array(rows_in); Aout = np.array(rows_out)
    Ag = np.array(rows_grow) if rows_grow else np.zeros((0, n))
    K = np.zeros((len(open_), len(open_)), complex)
    for i, a in enumerate(open_):
        rhs_v = np.zeros(len(open_) + Ag.shape[0], complex); rhs_v[i] = 1.0
        beta = np.linalg.lstsq(np.vstack([Ain, Ag]), rhs_v, rcond=None)[0]
        K[:, i] = Aout @ beta
    return K, [a+1 for a in open_]
