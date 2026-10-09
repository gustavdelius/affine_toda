"""Numerical checks of Theorem G and Proposition 8.5(1) for a list of sectors.
 (1) Hirota Wronskian matrix W(psi_b^-, psi_d^+) at complex lambda: x-independent, = -2 kappa_b X_b u_b^T u_d.
 (2) Fourier-grid determinant ratio det(A_N - lam)/det(A_0N - lam), cutoff-corrected, vs prod_b X_b(lam).
 (3) isolated eigenvalues of A_N (with multiplicities) vs net count of prod_b X_b.
Usage: python3 theoremG_check.py 'a4:1' 'd4:1' 'b3(1)<d4:3' ...  [env N, LP]"""
import sys, os, time
import numpy as np, mpmath as mp
from sectors import parent_sector, fold_sector, Xclosed, net_isolated

N = int(os.environ.get('N', '384'))
LPF = float(os.environ.get('LP', '40'))


def check(spec):
    nm, a = spec.rsplit(':', 1)
    a = int(a)
    S = fold_sector(nm, a) if '<' in nm else parent_sector(nm, a)
    mu = float(S.mu)
    Xf = {c[0]: Xclosed(S, c[0]) for c in S.channels}
    print(f"### {S.label}: mu={mu:.4f} xi={S.xi:.4f} min|tau|={S.bg.min_tau:.3f} channels m^2={ {k: round(v,4) for k,v in S.mb2.items()} }", flush=True)
    # (1) Wronskians
    xs = np.array([-2.0, 0.0, 1.5]) / mu
    for lam in [0.3 + 0.7j, 0.5 * min(S.mb2.values()) - 0.4j, 1.1 * max(S.mb2.values()) + 0.9j]:
        names, W, Xs = S.wronskian_matrix(lam, xs)
        const = np.max(np.abs(W - W[:, :, :1])) / np.max(np.abs(W))
        conj = {}
        for b in names:
            for d in names:
                if abs(S.mb2[b] - S.mb2[d]) < 1e-9 and abs(S.raw[b] @ S.raw[d]) > 1e-8:
                    conj[b] = d
        pred = np.zeros(W.shape[:2], complex)
        for i, b in enumerate(names):
            j = names.index(conj[b])
            pred[i, j] = -2 * S.kappa(b, lam) * Xs[b] * (S.raw[b] @ S.raw[conj[b]])
        err = np.max(np.abs(W[:, :, 0] - pred)) / np.max(np.abs(pred))
        xcl = max(abs(Xs[b] / Xf[b](S.kappa(b, lam) / np.sqrt(S.mb2[b])) - 1) for b in names)
        print(f"  W at lam={lam:.3f}: x-variation {const:.1e}; |W - diag(-2 kappa X u^T u)|/|W| = {err:.1e}; Hirota X vs closed form {xcl:.1e}", flush=True)
    # (2),(3) grid
    Lp = LPF / min(mu, 1.0) if LPF > 0 else 40
    H, H0, trV, kmax, x = S.grid(N, Lp)
    t0 = time.time()
    for lam in [0.3 + 0.7j, -0.5 + 0.2j, 1.7 - 1.1j]:
        s1, l1 = np.linalg.slogdet(H - lam * np.eye(len(H)))
        s0, l0 = np.linalg.slogdet(H0 - lam * np.eye(len(H)))
        ratio = s1 / s0 * np.exp(l1 - l0)
        corr = np.exp(trV / (np.pi * kmax))
        px = S.Xprod(lam, Xf)
        print(f"  det ratio at lam={lam}: raw/prodX={ratio/px:.6f}  corrected/prodX={ratio*corr/px:.8f}  (prod X={px:.5f})", flush=True)
    ev = np.linalg.eigvals(H)
    iso, emb, mmin2 = net_isolated(S)
    sel = np.sort_complex(ev[ev.real < mmin2 - 0.02])
    print(f"  N={N} Lp={Lp:.1f}: eigenvalues below lowest threshold {mmin2:.4f}: {np.round(sel, 5).tolist()}")
    print("  net count prediction:", [(round(l, 6), n) for l, n, _ in iso if n != 0], " cancelled:", [(round(l, 6), [(e[2], e[3], e[1]) for e in g]) for l, n, g in iso if n == 0])
    # compare
    ok = True
    used = np.zeros(len(sel), bool)
    for l, n, g in iso:
        near = np.where(np.abs(sel - l) < float(os.environ.get("CLUSTER", "0.05")))[0]
        used[near] = True
        if len(near) != n:
            ok = False
            print(f"   MISMATCH at {l:.6f}: predicted {n}, found {len(near)}")
        elif n > 0:
            print(f"   lam={l:.6f}: predicted {n}, found {len(near)}, max dev {np.max(np.abs(sel[near]-l)):.1e}")
    if not used.all():
        ok = False
        print("   UNPREDICTED eigenvalues:", sel[~used])
    print("  SPECTRUM", "OK" if ok else "FAIL", f"({time.time()-t0:.0f}s)", flush=True)


if __name__ == '__main__':
    for spec in sys.argv[1:]:
        check(spec)
