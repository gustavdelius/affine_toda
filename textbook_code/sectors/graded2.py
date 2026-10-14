QS = float(__import__("os").environ.get("QSIGN", "1"))   # QSIGN=-1: the physical q = -e^{-i pi omega}
import numpy as np, sys
from lib import apply_pair, krein_err
sys.argv = ['x', '2', 'c']
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
g = np.array([bin(i).count('1') % 2 for i in range(D)])
Pg = np.zeros((D*D, D*D))
for i in range(D):
    for j in range(D): Pg[j*D+i, i*D+j] = (-1)**(g[i]*g[j])
def transfer_graded(Rpl, ls):
    Nn = len(ls) + 1; T = np.eye(D**Nn, dtype=complex)
    for j, l in enumerate(ls, start=1):
        # R_{1,j+1} = (Pg chain) R_{12} (Pg chain)^{-1}, built by conjugating with graded swaps of factors (2..j+1)
        def op(v):
            for k in range(j, 1, -1): v = apply_pair(Pg, v, k-1, k, D, Nn)        # move factor j+1 -> position 2 (0-based 1)
            v = apply_pair(Rpl(np.exp(l)), v, 0, 1, D, Nn)
            for k in range(2, j+1): v = apply_pair(Pg, v, k-1, k, D, Nn)
            return v
        T = op(T)
    return T
for om in (2.37, 1.61, 0.3):
    q = QS*np.exp(-1j*np.pi*om); Pk = setup(q); Rpl = lambda x: Pg @ Rmat(x, q, Pk)
    rng = np.random.default_rng(5)
    for Nn in (2, 3, 4):
        w = 0; kr = 0
        for k in range({2: 200, 3: 150, 4: 60}[Nn]):
            ls = rng.uniform(-8, 8, Nn - 1); e = np.linalg.eigvals(transfer_graded(Rpl, ls)); d = np.abs(np.abs(e) - 1).max()
            if d > w: w = d; kr = krein_err(e)
        print(f"c_2 spinor, graded (fermionic for odd # down spins), omega={om}: N={Nn} max dev {w:.1e} (Krein {kr:.0e})", flush=True)
