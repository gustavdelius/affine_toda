QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
import numpy as np, sys, itertools
N_, alg_ = sys.argv[1], sys.argv[2]; sys.argv = ['x', N_, alg_]
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
m = np.array([[int(b) for b in format(i, f'0{N}b')] for i in range(D)])     # 1 = down spin (index bit)
lx = np.concatenate([np.linspace(0.05, 6, 25), [8, 11]])
oms = (2.37, 1.61, 0.3)
setups = {om: setup(QS*np.exp(-1j*np.pi*om)) for om in oms}
Rs = {om: [Rmat(np.exp(l), QS*np.exp(-1j*np.pi*om), setups[om]) for l in lx] for om in oms}
good = []
for bits in itertools.product([0, 1], repeat=N*N):
    M = np.array(bits).reshape(N, N)
    ph = np.array([(-1)**(int(m[i] @ M @ m[j]) % 2) for i in range(D) for j in range(D)])   # index i*D+j
    Q = np.diag(ph)
    dev = 0
    for om in oms:
        for R in Rs[om]:
            dev = max(dev, np.abs(np.abs(np.linalg.eigvals(P @ R @ Q)) - 1).max())
            if dev > 1e-6: break
        if dev > 1e-6: break
    if dev < 1e-6: good.append((M.tolist(), dev))
print(f"{alg}_{N}: {len(good)} of {2**(N*N)} Z2 bicharacter twists give unimodular 2-body spectra at sampled points:")
for g in good[:20]: print("  ", g)
