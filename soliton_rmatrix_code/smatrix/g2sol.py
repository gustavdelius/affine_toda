import numpy as np, os, time, itertools
from scipy.special import loggamma
exec(open('g2alg.py').read().split('if __name__ == "__main__":')[0])
Bc = list(np.load(f"{ALG}_comm_w{omega}.npy")); I = np.eye(n); N2 = n*n; t0 = time.time()
K0 = K(0); K0i = np.linalg.inv(K0)
EA, EB = np.kron(E[0], I), np.kron(K0, E[0]); FA, FB = np.kron(F[0], K0i), np.kron(I, F[0])
def Rcoef(x):
    cols = []
    for Bk in Bc:
        cols.append(np.concatenate([(Bk @ (x*EA + EB) - (EA + x*EB) @ Bk).ravel(), (Bk @ (FA/x + FB) - (FA + FB/x) @ Bk).ravel()]))
    Mx = np.array(cols).T; _, _, Vh = np.linalg.svd(Mx, full_matrices=False); c = Vh[-1].conj()
    R_ = sum(ck*Bk for ck, Bk in zip(c, Bc)); return R_/R_[0, 0]
def label(z):
    s = np.angle(z)/np.pi; best = None
    for L in range(-30, 31):
        for sig, nm in ((0, '+'), (1, '-'), (0.5, '+i'), (-0.5, '-i')):
            d_ = (s - (sig - L*omega)) % 2; d_ = min(d_, 2 - d_)
            if best is None or d_ < best[0]: best = (d_, L, nm, sig)
    return best
poles = []
for eps, nm, sig in ((1, '+', 0), (-1, '-', 1), (1j, '+i', 0.5), (-1j, '-i', -0.5)):
    for l in range(1, 31):
        xs = eps*q**(-l); n1 = np.linalg.norm(Rcoef(xs*(1 + 1e-5))); n2 = np.linalg.norm(Rcoef(xs*(1 + 1e-7)))
        if n2/n1 > 30:
            sv = np.linalg.svd(1e-7*Rcoef(xs*(1 + 1e-7)), compute_uv=False); poles.append((nm, l, sig, int(np.sum(sv > 1e-6*sv[0]))))
print(f"{ALG}: R-matrix poles {[(p[0] + 'q^-' + str(p[1]), 'rank', p[3]) for p in poles]}  [{time.time()-t0:.0f}s]")
xa = 0.83+0.29j; tc = None
def mats(x): return [((x if i == 0 else 1)*E[i], F[i]/(x if i == 0 else 1), K(i)) for i in range(3)]
for eps, nm, sig in ((1, '+', 0), (-1, '-', 1)):
    for l in range(-30, 31):
        tt = eps*q**l; G = 0
        for (Ei, Fi, Kk), (Ej, Fj, Kj) in zip(mats(xa), mats(xa*tt)):
            Kii = np.linalg.inv(Kk)
            for A_, B_ in (((-Kii@Ei).T, Ej), ((-Fi@Kk).T, Fj), (Kii.T, Kj)):
                r_ = np.kron(np.eye(n), A_.T) - np.kron(B_, np.eye(n)); G = G + r_.conj().T @ r_
        ev_ = np.linalg.eigvalsh(G)
        if ev_[0] < 1e-16*ev_[-1]: tc = (nm, l, sig); print(f"crossing point: V*(x) = V(x * {nm}q^{l})")
np.save(f"{ALG}_poles_w{omega}.npy", np.array([poles, tc], dtype=object), allow_pickle=True)
