"""Coordinator's spot check of A5 task 1 by direct ODE integration (scipy), a_3 species-1 pair vacua."""
import sys, contextlib, io, itertools
sys.path.insert(0, '/home/gustav/Git/affine_toda/code')
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
with contextlib.redirect_stdout(io.StringIO()):
    import boundary_vacua as BV
src = open('/home/gustav/Git/affine_toda/code/particle_reflection.py').read()
PR = {}; exec(src.split('\n# ---------------------------------------------------------------- 1. blocks')[0], PR)

n, h, a = 3, 4, 1
al = np.array([[1 if k == j else -1 if k == (j + 1) % h else 0 for k in range(h)] for j in range(h)], float)
mass = lambda c: 2*np.sin(np.pi*c/h)
print('BV.m =', BV.m)

def psi(chi, x):
    t = BV.pair_taus(n, a, chi, x)
    t0 = [tt[0] for tt in t]; t1 = [tt[1] for tt in t]
    e = np.array([t0[(j - 1) % h]*t0[(j + 1) % h]/t0[j]**2 for j in range(h)])
    return e, t0, t1

def patterns(chi):
    e, t0, t1 = psi(chi, 0.0)
    # phi = -sum_j alpha_j ln tau_j (units b = 1): phi_k = -(ln tau_k - ln tau_{k-1}) ... compute phi' directly
    dphi = -sum(al[j]*t1[j]/t0[j] for j in range(h))
    eh = np.sqrt(e)
    good = []
    for C in itertools.product([1, -1], repeat=h):
        r = dphi + sum(C[j]*eh[j]*al[j] for j in range(h))
        good.append((np.abs(r).max(), C))
    return sorted(good)[:2]

def kmat(chi, C, Om, X=10.0):
    e0 = psi(chi, 0.0)[0]
    H = 0.5*sum(C[j]*np.sqrt(e0[j])*np.outer(al[j], al[j]) for j in range(h))
    def rhs(x, y):
        ep = psi(chi, x)[0]
        M2 = sum(ep[j]*np.outer(al[j], al[j]) for j in range(h))
        return np.concatenate([y[h:], (M2 - Om**2*np.eye(h)) @ y[:h]])
    w = np.exp(2j*np.pi/h)
    rho = {c: w**(np.arange(h)*c) for c in range(1, h)}
    sols = []
    for c in range(1, h):                      # basis of BC-satisfying data inside the sum-zero subspace
        v = rho[c].astype(complex)
        y0 = np.concatenate([v, -H @ v])
        s = solve_ivp(rhs, (0, -X), y0, rtol=1e-12, atol=1e-14, method='DOP853')
        sols.append(s.y[:, -1])
    # decompose at x = -X into modes rho_c e^{+-ikx}
    modes, labels = [], []
    for c in range(1, h):
        if mass(c) < Om:
            p = np.sqrt(Om**2 - mass(c)**2)
            for sg, lab in ((1, 'in'), (-1, 'out')):
                modes.append(np.concatenate([rho[c]*np.exp(1j*sg*p*(-X)), 1j*sg*p*rho[c]*np.exp(1j*sg*p*(-X))])); labels.append((lab, c, p))
        else:
            k = np.sqrt(mass(c)**2 - Om**2)
            for sg, lab in ((1, 'dec'), (-1, 'grow')):
                modes.append(np.concatenate([rho[c]*np.exp(sg*k*(-X)), sg*k*rho[c]*np.exp(sg*k*(-X))])); labels.append((lab, c, k))
    Mm = np.array(modes).T
    coef = np.array([np.linalg.lstsq(Mm, s, rcond=None)[0] for s in sols]).T   # rows: modes, cols: basis solutions
    gi = [i for i, l in enumerate(labels) if l[0] == 'grow']
    ii = [i for i, l in enumerate(labels) if l[0] == 'in']
    oi = [i for i, l in enumerate(labels) if l[0] == 'out']
    # combos with no growing part: null space of coef[gi]
    if gi:
        _, _, Vh = np.linalg.svd(coef[gi]); N = Vh[len(gi):].conj().T
    else:
        N = np.eye(coef.shape[1])
    IN, OUT = coef[ii] @ N, coef[oi] @ N
    K = OUT @ np.linalg.inv(IN)
    return K, [labels[i][1] for i in ii]

for chi in (0.3, 1.0, 2.0, 2.9, 3.6, 5.0):
    pats = patterns(chi)
    if pats[0][0] > 1e-8 or pats[1][0] < 1e-3:
        print(f'chi={chi}: no unique pattern {pats[0][0]:.1e} {pats[1][0]:.1e}'); continue
    C = pats[0][1]
    th = 0.6; Om = mass(1)*np.cosh(th)
    K, chans = kmat(chi, C, Om)
    # class b of the pattern: compare with DG for b = 1, 2
    dg = {b: complex(PR['K_uniform'](1, n)(th, 0)*PR['M_sol'](1, b, n)(th, 0)) for b in (1, 2)}
    off = abs(K[0, 1]) if K.shape[0] > 1 else 0
    print(f'chi={chi} C={C} channels {chans}: K11={K[0,0]:.10f} off={off:.1e}  DG b=1 {dg[1]:.10f}  b=2 {dg[2]:.10f}')
