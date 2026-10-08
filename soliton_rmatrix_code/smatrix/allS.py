import numpy as np, itertools, time
import os
src = open('s33.py').read().replace('J = 6000', 'J = ' + os.environ.get('JTRUNC', '1500'))
exec(src[:src.index("# ---------------- pole scan")].split("th = 0.31+0.17j\nlhs = P @ Sfull")[0].split("# --- unitarity and crossing")[0])
S0 = Sbraid; F33 = F                                   # S_33 = c f R_33 ; F33 = c f
u1 = np.pi*(2*omega + 0.5)/T; u2 = np.pi*omega/T       # 3+3 -> 1 and 3+3 -> 2
wts64 = (Ws[:, None, :] + Ws[None, :, :]).reshape(64, 3)
def multiplet(u):
    d = 1e-7; U_, sv, _ = np.linalg.svd(d*S0(1j*(u + d))); r = int(np.sum(sv > 1e-6*sv[0])); return U_[:, :r]
def hw_coeff(B, w):
    mask = np.all(np.isclose(wts64, w), axis=1); Q = B[~mask, :]
    c = np.linalg.svd(Q)[2][-1].conj(); return c/np.linalg.norm(c)
B1, B2 = multiplet(u1), multiplet(u2)
part = {3: dict(off=[0.0], B=np.eye(8, dtype=complex), hw=np.eye(8)[0].astype(complex)),
        1: dict(off=[-u1/2, u1/2], B=B1, hw=hw_coeff(B1, [1, 0, 0])),
        2: dict(off=[-u2/2, u2/2], B=B2, hw=hw_coeff(B2, [1, 1, 0]))}
dims = {a: part[a]['B'].shape[1] for a in part}
def swaps(a, b):
    """list of (position, imaginary rapidity offset) for moving b's constituents left through a's"""
    sites = [('a', o) for o in part[a]['off']] + [('b', o) for o in part[b]['off']]
    seq = []; na = len(part[a]['off'])
    for j in range(len(part[b]['off'])):
        pos = na + j
        for p in range(pos - 1, j - 1, -1):
            L, R = sites[p], sites[p+1]                      # L is an a-site, R the moving b-site
            seq.append((p, L[1] - R[1])); sites[p], sites[p+1] = R, L
    return seq
def apply_S(a, b, theta, X):
    Ba, Bb = part[a]['B'], part[b]['B']; Na, Nb = len(part[a]['off']), len(part[b]['off']); N = Na + Nb
    m = X.shape[1]
    V = np.einsum('ia,jb,abm->ijm', Ba, Bb, X.reshape(dims[a], dims[b], m)).reshape((8,)*N + (m,))
    for p, off in swaps(a, b):
        M = S0(theta + 1j*off)
        V = np.moveaxis(V, [p, p+1], [0, 1]); sh = V.shape
        V = np.moveaxis((M @ V.reshape(64, -1)).reshape(sh), [0, 1], [p, p+1])
    V = V.reshape(8**Nb, 8**Na, m)
    Y = np.einsum('ib,ja,ijm->bam', Bb.conj(), Ba.conj(), V)
    resid = np.linalg.norm(V - np.einsum('ib,ja,bam->ijm', Bb, Ba, Y))/max(np.linalg.norm(V), 1e-300)
    return Y.reshape(dims[b]*dims[a], m), resid
def K_top(a, b, theta):
    X = np.kron(part[a]['hw'], part[b]['hw'])[:, None]; Y, _ = apply_S(a, b, theta, X)
    tgt = np.kron(part[b]['hw'], part[a]['hw'])
    K = np.vdot(tgt, Y[:, 0]); return K, np.linalg.norm(Y[:, 0] - K*tgt)/max(abs(K), 1e-300)
def prodF(a, b, theta):
    return np.prod([F33(theta + 1j*off) for _, off in swaps(a, b)])
if __name__ == "__main__":
    print("multiplet dimensions:", dims)
    th = 0.23 + 0.13j
    for a, b in [(3, 1), (3, 2), (1, 2), (1, 1), (2, 2)]:
        X = np.eye(dims[a]*dims[b], dtype=complex)
        Yab, r1 = apply_S(a, b, th, X); Yba, r2 = apply_S(b, a, -th, np.eye(dims[a]*dims[b], dtype=complex))
        U = Yba @ Yab
        K, rk = K_top(a, b, th)
        print(f"S_{a}{b}: closes on multiplets ({r1:.0e}); unitarity S_{b}{a}(-th)S_{a}{b}(th) = 1 to {np.abs(U - np.eye(len(U))).max():.1e};"
              f" highest-weight state maps to highest weight ({rk:.0e}); swaps {[(p, round(o/np.pi*T, 3)) for p, o in swaps(a, b)]} (offsets in t)")
