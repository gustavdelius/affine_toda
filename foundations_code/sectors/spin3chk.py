QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
import numpy as np, sys
sys.argv = ['x', '2', 'c']
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
om = 2.37; q = QS*np.exp(-1j*np.pi*om); Pk = setup(q); Rpl = lambda x: P @ Rmat(x, q, Pk)
for l in (1.0, 2.0, 3.0, 5.0):
    print("2-body logx", l, "maxdev", np.abs(np.abs(np.linalg.eigvals(Rpl(np.exp(l)))) - 1).max())
for ls in ([3.0, 1.0], [3.0, 5.0], [3.0, 3.1], [3.0, 0.0001], [-3, 2], [4, -4], [2.5, 6]):
    T1 = transferN(Rpl, D, ls); e = np.linalg.eigvals(T1)
    print("3-body", ls, "maxdev T1", np.abs(np.abs(e) - 1).max())
# 4-body
rng = np.random.default_rng(3); w = 0
for k in range(40):
    ls = rng.uniform(-6, 6, 3); e = np.linalg.eigvals(transferN(Rpl, D, ls)); w = max(w, np.abs(np.abs(e) - 1).max())
print("4-body random max dev", w)
# 2-body is P R ; what about R alone (positional)?
print("positional Rcheck eigen maxdev at logx 3:", np.abs(np.abs(np.linalg.eigvals(Rmat(np.exp(3.0), q, Pk))) - 1).max())
