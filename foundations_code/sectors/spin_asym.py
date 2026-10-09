import numpy as np, sys, math
sys.argv = ['x', sys.argv[1], sys.argv[2]]
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
for om in [0.05, 0.15, 0.25, 0.35, 0.45, 0.49, 0.55, 0.75]:
    q = np.exp(-1j*np.pi*om); Pk = setup(q)
    ds = []
    for lx in (8, 12, 16, 20):
        e = np.linalg.eigvals(P @ Rmat(np.exp(lx), q, Pk)); ds.append(np.abs(np.abs(e) - 1).max())
    p = [-(np.log(ds[i+1]) - np.log(ds[i]))/4 for i in range(3)]
    print(f"omega {om}: dev at logx 8,12,16,20: {np.array(ds)}  ; local exponent p (dev ~ x^-p): {np.round(p, 3)}")
