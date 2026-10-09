import numpy as np, sys, time
args = sys.argv[1:]; sys.argv = ['x', args[0], args[1]]
exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
from lib import transferN, krein_err
om = float(args[2]); lc = float(args[3]); q = np.exp(-1j*np.pi*om); Pk = setup(q); Rpl = lambda x: P @ Rmat(x, q, Pk)
rng = np.random.default_rng(1); t0 = time.time(); wi = wa = 0; ki = ka = 0
for k in range(int(args[4])):
    th = rng.uniform(0, 0.95*lc, 3); e = np.linalg.eigvals(transferN(Rpl, D, [th[0]-th[1], th[0]-th[2]])); d = np.abs(np.abs(e)-1).max()
    if d > wi: wi, ki = d, krein_err(e)
    e = np.linalg.eigvals(transferN(Rpl, D, rng.uniform(-6, 6, 2))); d = np.abs(np.abs(e)-1).max()
    if d > wa: wa, ka = d, krein_err(e)
print(f"{alg}_{N} spinor omega={om}: 3-body inside 2-body window (|l|<{lc}): max dev {wi:.1e} (Krein {ki:.0e}); unrestricted: {wa:.1e} (Krein {ka:.0e})  [{time.time()-t0:.0f}s]")
