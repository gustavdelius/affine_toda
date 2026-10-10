QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
"""c_3^(1): fused amplitudes from S_33 = F33 R33.  R-only products give exact eigenvalue moduli for symmetric fusions;
F33 (with Richardson-extrapolated truncation) is included where needed (asymmetric fusions, 3*)."""
import numpy as np, os
from scipy.special import loggamma
from common import Spinor33, Ws, n
omega = float(os.environ.get('OMEGA', '2.37')); q = QS*np.exp(-1j*np.pi*omega); T = 3*omega + 1; mu = 2*T
S = Spinor33(q)
A = [omega, 2*omega + 0.5, 3*omega]
def c_of(t): return np.prod([np.sin(np.pi*(t - a)) for a in A], axis=0)
def lf1(z):
    z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A) - len(A)*np.log(np.pi)
def logf_J(t, J):
    js = np.arange(1, J+1); Aj = 2*T*(js-1); d = lambda z, z0: lf1(z) - lf1(z0)
    return np.sum(d(t + Aj, Aj) - d(t + Aj + T, Aj + T) + d(-t + Aj + T, Aj + T) - d(-t + Aj + 2*T, Aj + 2*T))
JJ = int(os.environ.get('JTRUNC', '2000'))
def logF(theta, J=None, rich=True):
    t = mu*theta/(2j*np.pi)
    if rich:
        l1, l2 = logf_J(t, JJ), logf_J(t, 2*JJ); lf = 2*l2 - l1
    else: lf = logf_J(t, J or JJ)
    return np.log(c_of(t) + 0j) + lf - np.log(c_of(0.0) + 0j)
def R33(theta): return S.R(np.exp(mu*theta))
def S33(theta): return np.exp(logF(theta))*R33(theta)
u1 = np.pi*(2*omega + 0.5)/T; u2 = np.pi*omega/T
P64 = np.zeros((64, 64))
for i in range(8):
    for j in range(8): P64[j*8+i, i*8+j] = 1
def orthimg(M, r):
    U_, sv, _ = np.linalg.svd(M); return U_[:, :r], sv
# multiplets of 1 and 2: image of residue of R33 = range of the projectors with the pole
P1img = orthimg(S.Pk[2] + S.Pk[3], 8)[0]                      # 3+3 -> 1 at x = -q^-4 : channels 7 + 1
P2img = orthimg(S.Pk[1] + S.Pk[2] + S.Pk[3], 29)[0]           # 3+3 -> 2 at x = q^-2 : channels 21 + 7 + 1
part = {3: dict(off=[0.0], B=np.eye(8, dtype=complex)),
        1: dict(off=[-u1/2, u1/2], B=P1img),
        2: dict(off=[-u2/2, u2/2], B=P2img)}
def dims(a): return part[a]['B'].shape[1]
def swaps(a, b):
    sites = [('a', o) for o in part[a]['off']] + [('b', o) for o in part[b]['off']]
    seq = []; na = len(part[a]['off'])
    for j in range(len(part[b]['off'])):
        pos = na + j
        for p in range(pos - 1, j - 1, -1):
            L, R = sites[p], sites[p+1]; seq.append((p, L[1] - R[1])); sites[p], sites[p+1] = R, L
    return seq
def apply(a, b, theta, X, Rf):
    Ba, Bb = part[a]['B'], part[b]['B']; Na, Nb = len(part[a]['off']), len(part[b]['off']); N = Na + Nb; m = X.shape[1]
    V = np.einsum('ia,jb,abm->ijm', Ba, Bb, X.reshape(dims(a), dims(b), m)).reshape((8,)*N + (m,))
    for p, off in swaps(a, b):
        M = Rf(theta + 1j*off)
        V = np.moveaxis(V, [p, p+1], [0, 1]); sh = V.shape
        V = np.moveaxis((M @ V.reshape(64, -1)).reshape(sh), [0, 1], [p, p+1])
    V = V.reshape(8**Nb, 8**Na, m)
    Y = np.einsum('ib,ja,ijm->bam', Bb.conj(), Ba.conj(), V)
    resid = np.linalg.norm(V - np.einsum('ib,ja,bam->ijm', Bb, Ba, Y))/max(np.linalg.norm(V), 1e-300)
    return Y.reshape(dims(b)*dims(a), m), resid
def flip(Da, Db):
    P = np.zeros((Da*Db, Da*Db))
    for i in range(Da):
        for j in range(Db): P[j*Da + i, i*Db + j] = 1
    return P
def Spl(a, b, theta, Rf=R33):
    """particle-labelled amplitude on W_a (x) W_b"""
    Y, res = apply(a, b, theta, np.eye(dims(a)*dims(b), dtype=complex), Rf)
    return flip(dims(b), dims(a)) @ Y, res
def scalar_logF(a, b, theta):
    return sum(logF(theta + 1j*off) for _, off in swaps(a, b))
def add_bound(name, a, b, t0, rank, ubar=None, Rf=S33):
    """bound state in a+b at t0 (pole of S_ab); constituents at offsets -ubar_a (a... ) following allS conventions:
       image of residue lives in W_b(left) (x) W_a(right) [output of S_ab]; offsets [o - ub for o in b] + [o + ua for o in a]"""
    u = np.pi*t0/T
    if ubar is None: ua = ub = u/2
    else: ua, ub = ubar
    d = 1e-7
    Y, _ = apply(a, b, 1j*(u + d), np.eye(dims(a)*dims(b), dtype=complex), Rf)
    img, sv = orthimg(d*Y, rank)
    vec = np.kron(part[b]['B'], part[a]['B']) @ img
    part[name] = dict(off=[o - ub for o in part[b]['off']] + [o + ua for o in part[a]['off']], B=vec)
    return sv[:rank+2]/sv[0]
