"""Independent numerical check on a periodic Fourier (Galerkin) grid:
  (i) spectrum of the folded fluctuation operator A = -d^2 + W^T M(x) W below threshold, counted with algebraic
      multiplicity (eigenvalue clusters), versus the net order of prod_b X_b;
  (ii) heat trace Tr(e^{-tA} - e^{-tA0}) versus Proposition 8.5 part 2 with the product-rule X_b;
  (iii) the one-loop energy (8.1) evaluated directly from the eigenvalues, versus the closed-form mass formula.
The soliton is built from the Hirota tau functions (composite: tau-recursion with summed delta)."""
import numpy as np, sys
from scipy.integrate import quad
from toda import soliton, degree, perm_matrix
from fold import load, sigma_action, orbits, folded
from particles import parent_roots
pi = np.pi

def build(alg, r, perm, cons_index):
    D = load(alg, r); C, n, h = D['C'], D['n'], D['h']
    masses, vecs = D['masses'], D['vecs']
    roots, _ = parent_roots(alg, r)
    N = r + 1
    if perm is None:
        W = np.eye(r); G = [np.eye(N)]
    else:
        P = perm_matrix(perm)
        G = [np.eye(N)]; g = P.copy()
        while not np.allclose(g, np.eye(N)): G.append(g); g = g @ P
        # invariant subspace of root space
        A_ = roots[1:]; B_ = np.array([roots[perm[j]] for j in range(1, N)])
        S = B_.T @ np.linalg.inv(A_.T)
        Gs = [np.eye(r)]; gs = S.copy()
        while not np.allclose(gs, np.eye(r)): Gs.append(gs); gs = gs @ S
        Pi = sum(Gs)/len(Gs); u, s, vt = np.linalg.svd(Pi); W = u[:, s > 0.5]
    a0 = cons_index
    mu = masses[a0]
    delta = sum(gm @ vecs[a0] for gm in G)
    delta = delta/np.max(np.abs(delta))
    c = soliton(C, n, mu, delta, Kmax=int(4*max(n)) + 4)
    deg = degree(c)
    c = c[:, :max(deg)+1]
    return dict(C=C, n=n, h=h, roots=roots, W=W, c=c, mu=mu, deg=deg, masses=masses)

def Mfield(B, x, xi):
    C, n, roots, W, c, mu = B['C'], B['n'], B['roots'], B['W'], B['c'], B['mu']
    E = np.exp(mu*x + xi)
    deg = B['deg']
    logt = np.zeros((len(n), len(x)), complex); tau = np.zeros((len(n), len(x)), complex)
    big = np.abs(E) > 1
    for j in range(len(n)):
        cj = c[j][:deg[j]+1]
        small = np.polyval(cj[::-1], E[~big])
        # for |E|>1: tau_j = E^D (c_D + c_{D-1}/E + ...)
        Ei = 1/E[big]
        large = np.polyval(cj, Ei)            # sum_k c_k E^{k-D} = sum_k c_k Ei^{D-k}
        logt[j, ~big] = np.log(small); logt[j, big] = deg[j]*np.log(E[big]) + np.log(large)
        tau[j, ~big] = small; tau[j, big] = large
    
    ex = np.exp(-(C @ logt))                                              # e^{i beta alpha_k.phi} = prod tau^{-C}
    Ar = roots @ W                                                        # projected roots (N, d)
    M = np.einsum('j,ja,jb,jx->xab', n, Ar, Ar, ex)
    M0 = np.einsum('j,ja,jb->ab', n, Ar, Ar)
    return M, M0, np.min(np.abs(tau))

def spectrum(B, L=40.0, Npts=256, xi=None):
    x = (np.arange(Npts) - Npts/2)*L/Npts
    if xi is None:
        args = []
        for j in range(len(B['n'])):
            cj = B['c'][j][:B['deg'][j]+1]
            if len(cj) > 1: args += list(np.angle(np.roots(cj[::-1])))
        grid = np.linspace(-pi, pi, 2001)
        dist = [min(abs(np.angle(np.exp(1j*(g - a)))) for a in args) for g in grid]
        xi = 1j*grid[int(np.argmax(dist))]
    M, M0, mn = Mfield(B, x, xi)
    d = M0.shape[0]
    k = 2*pi*np.fft.fftfreq(Npts, d=L/Npts)
    F = np.fft.fft(np.eye(Npts), axis=0)/np.sqrt(Npts)        # unitary DFT
    Fi = F.conj().T
    A = np.zeros((Npts*d, Npts*d), complex)
    for a in range(d):
        for b in range(d):
            Vab = Fi @ np.diag(M[:, a, b]) @ F if False else None
    # build in Fourier basis: kinetic diag, potential convolution
    Mk = np.fft.fft(M, axis=0)/Npts                           # Fourier coefficients of M(x)
    idx = (np.arange(Npts)[:, None] - np.arange(Npts)[None, :]) % Npts
    for a in range(d):
        for b in range(d):
            blk = Mk[idx, a, b]
            if a == b: blk = blk + np.diag(k**2)
            A[a*Npts:(a+1)*Npts, b*Npts:(b+1)*Npts] = blk
    lam = np.linalg.eigvals(A)
    m2 = np.linalg.eigvalsh(M0)
    lam0 = np.concatenate([k**2 + m for m in m2])
    trV = np.einsum('xaa->x', M - M0[None]).mean()            # (1/N) sum_x tr(M - M0) -> for counterterm
    # counterterm Tr[(A-A0) A0^{-1/2}] in the channel basis of M0
    w0, U = np.linalg.eigh(M0)
    Vch = np.einsum('ia,xab,bj->xij', U.T, M - M0[None], U).mean(axis=0)   # (1/N) sum_x (M-M0) in channel basis
    ct = sum(Vch[i, i].real*np.sum(1/np.sqrt(k**2 + w0[i])) for i in range(d)) + 1j*sum(Vch[i, i].imag*np.sum(1/np.sqrt(k**2 + w0[i])) for i in range(d))
    dM = 0.5*np.sum(np.sqrt(lam + 0j) - np.sqrt(lam0 + 0j)) - 0.25*ct
    return lam, lam0, dM, xi, mn

def heat_formula(chans, h, t):
    tot = 0.0
    for co, mb, roots in chans:
        for x, o in roots:
            w = np.cos(pi*x/h)
            if x < h/2 - 1e-9: tot += o*np.exp(-t*mb*mb*(1 - w*w))
            g = lambda kk: (o/(pi*1j)/(kk + 1j*mb*w)*np.exp(-t*(mb*mb + kk*kk)))
            tot += quad(lambda kk: g(kk).real, 0, np.inf, limit=200)[0]
    return tot

if __name__ == "__main__":
    from analysis import dspin, D4_Z3, E6_Z2
    cases = [("g_2 triple composite", 'd', 4, D4_Z3, 0), ("b_3 spinor composite", 'd', 4, dspin(4), 1),
             ("b_4 spinor composite", 'd', 5, dspin(5), 1), ("f_4 composite {1,6}", 'e', 6, E6_Z2, 0)]
    sel = sys.argv[1:] and [cases[int(i)] for i in sys.argv[1:]] or cases
    for name, alg, r, perm, a0 in sel:
        B = build(alg, r, perm, a0)
        F = folded(alg, r, perm, name, verbose=False)
        s = [s_ for s_ in F['sols'] if a0 in s_['cons']][0]
        print(f"\n=== {name}: constituents {s['cons']}, tau degrees {B['deg']}, closed-form Delta M = {s['dM']:.6f}")
        thr = min(mb**2 for co, mb, roots in s['chans'])
        # predicted below-threshold eigenvalues with net order
        from collections import defaultdict
        pred = defaultdict(int)
        for co, mb, roots in s['chans']:
            for x, o in roots:
                if x < F['h']/2 - 1e-9:
                    lam_ = mb*mb*(1 - np.cos(pi*x/F['h'])**2)
                    pred[round(lam_, 8)] += o
        pred = {k_: v for k_, v in pred.items() if v != 0 and k_ < thr - 1e-9}
        print(f"  predicted isolated eigenvalues (lambda: net order): {pred}   lowest threshold {thr:.6f}")
        for L, Np in [(40.0, 192), (48.0, 256), (56.0, 320)]:
            lam, lam0, dM, xi, mn = spectrum(B, L=L, Npts=Np)
            below = np.sort_complex(lam[lam.real < thr - 0.02])
            clusters = {}
            for lv in pred:
                clusters[lv] = int(np.sum(np.abs(lam - lv) < 0.02))
            print(f"  L={L}, N={Np}, Im xi={xi.imag:.3f}, min|tau|={mn:.2e}: eigenvalues below threshold: {np.round(below, 5)}")
            print(f"      cluster counts near predicted eigenvalues: {clusters};  lattice (8.1) Delta M = {dM.real:.6f} {dM.imag:+.1e}i")
        ts = [0.3, 0.6, 1.0, 1.5]
        ht_num = [np.sum(np.exp(-t*lam)) - np.sum(np.exp(-t*lam0)) for t in ts]
        ht_th = [heat_formula(s['chans'], F['h'], t) for t in ts]
        print("  heat trace numeric :", np.round(np.real(ht_num), 5))
        print("  heat trace formula :", np.round(ht_th, 5))
