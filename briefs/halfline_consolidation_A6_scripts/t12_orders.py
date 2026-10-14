"""Orders of the zeros/poles of c_sigma(mu) at candidate real points (winding number of c on a small circle)."""
import sys
import numpy as np
from t10_fitc import cval, n

om = float(sys.argv[2])
name = sys.argv[3]
pts = [float(x) for x in sys.argv[4].split(',')]
s = 1 if name == 's+' else -1
sig = s*(-1)**n
r = 0.006
K = 12
for p in pts:
    zs = p + r*np.exp(2j*np.pi*np.arange(K + 1)/K) + 1e-4j*0
    cs = np.array([cval(name, om, z, sig) for z in zs])
    ph = np.unwrap(np.angle(cs))
    wnd = (ph[-1] - ph[0])/(2*np.pi)
    print(f'{name} om={om} mu0={p}: winding {wnd:.3f}', flush=True)
