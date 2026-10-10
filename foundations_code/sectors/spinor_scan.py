QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
import numpy as np, sys, math, time
from lib import *
N = int(sys.argv[1]); alg = sys.argv[2]
D = 2**N; Wn = spin_weights(N); P = flip(D, D)
roots = {i: np.eye(N)[i-1] - np.eye(N)[i] for i in range(1, N)}; roots[N] = np.eye(N)[N-1]
def rhos(x, q):
    if alg == 'c':   # Lambda^{N-k}: prod_{j<=k} <2j>_{(-1)^{j+1}}
        out = [1.0]
        for j in range(1, N+1): out.append(out[-1]*br(x, 2*j, (-1)**(j+1), q))
    else:            # Lambda^{N-k}: prod_{i=1}^{ceil(k/2)} <4k-8i+6>_+
        out = []
        for k in range(N+1):
            v = 1.0
            for i in range(1, (k+1)//2 + 1): v = v*br(x, 4*k - 8*i + 6, 1, q)
            out.append(v)
    return out
def setup(q):
    e, f, k, qi = spinor_rep(N, q, 1.0, alg)
    comps = components(e, f, k, Wn, range(1, N+1), roots)
    byd = {d: Pm for d, Pm, _ in comps}
    Pk = [byd[math.comb(2*N+1, N-kk)] for kk in range(N+1)]
    return check_projectors(Pk, Wn)
def Rmat(x, q, Pk): return sum(r*Pm for r, Pm in zip(rhos(x, q), Pk))
if __name__ == "__main__":
    # verify against direct intertwiner for small N
    if N <= 3 and len(sys.argv) > 3:
        from rsolve import intertwiner
        q = QS*np.exp(-1j*np.pi*0.731); Pk = setup(q); x0 = 0.77+0.3j
        Rc = intertwiner(spinor_rep(N, q, x0, alg), spinor_rep(N, q, 1.0, alg), range(N+1), Wn, Wn)[0]; Rc /= Rc[0, 0]
        print("closed form vs intertwiner:", np.abs(Rmat(x0, q, Pk) - Rc).max())
    lxs = np.concatenate([np.linspace(0.005, 8, 400), np.linspace(8, 40, 81)])
    oms = [float(v) for v in __import__("os").environ["OMS"].split(",")] if "OMS" in __import__("os").environ else np.arange(0.0125, 2.0, 0.025)
    t0 = time.time()
    for om in oms:
        q = QS*np.exp(-1j*np.pi*om); Pk = setup(q)
        dp = maxdev_scan(lambda x: P @ Rmat(x, q, Pk), lxs, (1,))
        dm = maxdev_scan(lambda x: P @ Rmat(x, q, Pk), lxs, (-1,))
        print(f"omega {om:.4f}: x>0 maxdev {dp[0]:.2e} (logx {dp[2]})   x<0 maxdev {dm[0]:.2e} (logx {dm[2]})  [{time.time()-t0:.0f}s]", flush=True)
