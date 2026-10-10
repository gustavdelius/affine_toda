QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
import numpy as np, sys
N_ = sys.argv[1]; sys.argv = ['x', N_, 'c']
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
g = np.array([bin(i).count('1') % 2 for i in range(D)])
Pg = np.zeros((D*D, D*D))
for i in range(D):
    for j in range(D): Pg[j*D+i, i*D+j] = (-1)**(g[i]*g[j])
for om in (2.37, 1.61, 0.3):
    q = QS*np.exp(-1j*np.pi*om); Pk = setup(q)
    for nm, PP in (("plain flip", P), ("graded flip", Pg)):
        dev = max(np.abs(np.abs(np.linalg.eigvals(PP @ Rmat(np.exp(l), q, Pk))) - 1).max() for l in np.linspace(0.01, 12, 200))
        print(f"n={N} omega={om}: 2-body with {nm}: maxdev {dev:.2e}")
    # does the plain-flip 2-body spectrum change sign structure? check eigenvalue sets of P R and Pg R at one point
