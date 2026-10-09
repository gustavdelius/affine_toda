"""Jordan structure of isolated fluctuation eigenvalues from the Fourier-grid operator.
For each target lambda0: cluster of grid eigenvalues near lambda0 (count = algebraic multiplicity),
ordered Schur form, singular values of (B - lambda_bar)^k on the invariant subspace (k = 1..m).
For Jordan type (p_1,...,p_g): rank (B - lambda_bar)^k = sum_i max(p_i - k, 0).
Usage: python3 grid_jordan.py 'g2(1)<d4:1@1.5' 'b3(1)<d4:3@1.5' ...   [env N (list), LP]"""
import sys, os
import numpy as np
import scipy.linalg as sla
from sectors import parent_sector, fold_sector, net_isolated

Ns = [int(x) for x in os.environ.get('N', '256,384').split(',')]
LPF = float(os.environ.get('LP', '40'))


def run(spec):
    s, lam0 = spec.split('@')
    lam0 = float(lam0)
    nm, a = s.rsplit(':', 1)
    a = int(a)
    S = fold_sector(nm, a) if '<' in nm else parent_sector(nm, a)
    iso, emb, mmin2 = net_isolated(S)
    pred = [n for l, n, g in iso + emb if abs(l - lam0) < 1e-9]
    pred = pred[0] if pred else 0
    print(f"### {S.label}, lambda0={lam0}: net count {pred}; lowest threshold {mmin2:.4f}", flush=True)
    for N in Ns:
        Lp = LPF / min(float(S.mu), 1.0)
        H, H0, trV, kmax, x = S.grid(N, Lp)
        R = float(os.environ.get("R", "0.1"))
        T, Z, sdim = sla.schur(H, output='complex', sort=lambda z: abs(z - lam0) < R)
        B = T[:sdim, :sdim]
        ev = np.diag(B)
        lb = ev.mean()
        sv = []
        Mk = np.eye(sdim, dtype=complex)
        for k in range(1, sdim + 1):
            Mk = Mk @ (B - lb * np.eye(sdim))
            sv.append(np.linalg.svd(Mk, compute_uv=False))
        spread = np.max(np.abs(ev - lb)) if sdim else 0
        print(f"  N={N} Lp={Lp:.0f}: {sdim} eigenvalues within {R} of lambda0; mean-lambda0={lb-lam0:.2e}, spread={spread:.2e}")
        print(f"     eigenvalues - lambda0: {np.round(ev - lam0, 8).tolist()}")
        for k, s_ in enumerate(sv, 1):
            print(f"     singular values of (B-lb)^{k}: {np.array2string(s_, precision=2)}")
        sys.stdout.flush()


if __name__ == '__main__':
    for spec in sys.argv[1:]:
        run(spec)
