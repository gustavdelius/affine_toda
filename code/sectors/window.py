QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
import numpy as np, sys
args = sys.argv[1:]; sys.argv = ['x', args[0], args[1]]
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
out = []
for om in (2.37, 1.61, 0.3013):   # 0.3 is a root of unity (q^20 = 1): rank-3 projectors degenerate
    q = QS*np.exp(-1j*np.pi*om); Pk = setup(q)
    lc = None
    for lx in np.linspace(0.001, 12, 2401):
        if np.abs(np.abs(np.linalg.eigvals(P @ Rmat(np.exp(lx), q, Pk))) - 1).max() > 1e-7: lc = lx; break
    T = {'c': N*om + (N-1)/2, 'a': (2*N-1)*om + N - 1}[alg]
    out.append(f"omega={om}: l_c={lc:.3f} (theta_c = l_c/2T = {lc/(2*T):.4f})")
print(f"{alg}_{N} spinor:", "; ".join(out))
