QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
import numpy as np, sys, math, time
args = sys.argv[1:]; sys.argv = ['x', args[0], args[1]]
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
t0 = time.time()
for om in [float(v) for v in args[2].split(',')]:
    q = QS*np.exp(-1j*np.pi*om); Pk = setup(q); Rpl = lambda x: P @ Rmat(x, q, Pk)
    lc = None
    for lx in np.linspace(0.001, 12, 4801):
        if np.abs(np.abs(np.linalg.eigvals(Rpl(np.exp(lx)))) - 1).max() > 1e-7: lc = lx; break
    rng = np.random.default_rng(1); worst_in = (0, None); worst_all = (0, None); kin = kall = 0
    for k in range(150):
        # restricted: all pairwise |l_ij| < 0.95 lc  (theta_1, theta_2, theta_3 within a window)
        th = rng.uniform(0, 0.95*lc, 3); ls = [th[0] - th[1], th[0] - th[2]]
        e = np.linalg.eigvals(transferN(Rpl, D, ls)); d = np.abs(np.abs(e) - 1).max()
        if d > worst_in[0]: worst_in = (d, np.round(ls, 3)); kin = krein_err(e)
        ls2 = rng.uniform(-6, 6, 2)
        e = np.linalg.eigvals(transferN(Rpl, D, ls2)); d = np.abs(np.abs(e) - 1).max()
        if d > worst_all[0]: worst_all = (d, np.round(ls2, 3)); kall = krein_err(e)
    print(f"{alg}_{N} spinor, omega={om}: 2-body unimodular for |log x| < {lc:.4f}; 3-body (all pairs inside window) max dev {worst_in[0]:.1e} at {worst_in[1]} (Krein {kin:.0e}); "
          f"3-body unrestricted max dev {worst_all[0]:.1e} (Krein {kall:.0e})  [{time.time()-t0:.0f}s]", flush=True)
