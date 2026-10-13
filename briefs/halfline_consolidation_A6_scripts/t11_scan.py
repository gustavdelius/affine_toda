"""Locate zeros/poles of c_sigma(mu) on the real mu axis (all candidate positions are real for real omega), and check
periodicity."""
import sys
import numpy as np
from t10_fitc import cval, n

om = float(sys.argv[2]) if len(sys.argv) > 2 else 0.53
T = (n - 1)*om
for name, s in (('s+', 1), ('s-', -1)):
    sig = s*(-1)**n
    m0 = 0.123 + 0.071j
    print(f'{name}: c(mu+1)/c(mu) = {cval(name, om, m0 + 1, sig)/cval(name, om, m0, sig):.8f}, '
          f'c(mu+2)/c(mu) = {cval(name, om, m0 + 2, sig)/cval(name, om, m0, sig):.8f}, '
          f'c(-mu)c(mu) = {cval(name, om, -m0, sig)*cval(name, om, m0, sig):.8f}')
    xs = np.linspace(-1.0, 1.0, 801) + 1e-3j
    vals = np.array([cval(name, om, x, sig) for x in xs])
    la = np.log10(np.abs(vals))
    ext = []
    for i in range(1, len(xs) - 1):
        if la[i] > la[i - 1] and la[i] > la[i + 1] and la[i] > np.median(la) + 1.0:
            ext.append(('pole', round(xs[i].real, 4)))
        if la[i] < la[i - 1] and la[i] < la[i + 1] and la[i] < np.median(la) - 1.0:
            ext.append(('zero', round(xs[i].real, 4)))
    print(f'  om={om} T={T}: {ext}')
