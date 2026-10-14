"""Referee (2026-10-12): outward-integration solver for the linearized reflection matrix (independent of linrefl2):
start at the wall with the n solutions y(0)=u_k, y'(0)=-H u_k (BC imposed exactly), integrate to x=-X, decompose
z_c = v_c^dagger y into in/out (open) or growing/decaying (closed) modes, impose zero growing amplitude.
Background, basis, channels and BP values are taken from ref0601_a5_K.py (same directory)."""
import sys, time, numpy as np
from scipy.integrate import solve_ivp
exec(open(sys.argv[1]).read().split("if __name__ == '__main__':")[0])
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
if __name__ == '__main__':
    chi_end = 5*np.pi/3
    targets = [(1, 1, -1, -1, -1, -1), (1, 1, -1, -1, 1, -1), (1, 1, -1, 1, -1, -1), (1, 1, -1, 1, 1, -1)]
    deltas = [float(x) for x in sys.argv[2].split(',')]
    for Om in [float(x) for x in sys.argv[3].split(',')]:
        ths = {a: np.arccosh(Om/ma[a-1]) for a in range(1, 6) if Om > ma[a-1]}
        print('\n== Omega=%.4f open channels %s; BP diagonal values %s' % (Om, list(ths), {a: np.round(K_bp_diag(a, ths[a]), 6) for a in ths}), flush=True)
        for side in (+1, -1):
            regC = exact_pattern(chi_end + side*0.2)[1]
            for C in [regC] + targets:
                for delta in deltas:
                    t = time.time(); K, lab = K_out(chi_end + side*delta, C, Om)
                    off = max(abs(K[i, j]) for i in range(len(lab)) for j in range(len(lab)) if i != j)
                    dg = max(abs(K[i, i] - K_bp_diag(lab[i], ths[lab[i]])) for i in range(len(lab)))
                    uni = np.abs(K.conj().T @ K - np.eye(len(lab))).max()
                    print('side %+d %s%s delta=%.0e tau2(0)=%.1e max|offdiag|=%.2e max|Kdiag-BP|=%.2e unitarity %.0e  (%.0fs)  K11=%s' % (
                        side, s(C), '(reg)' if C == regC else '     ', delta, tau_all(0.0, chi_end+side*delta)[0][2], off, dg, uni, time.time()-t, np.round(K[0, 0], 5)), flush=True)
