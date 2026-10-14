import sys
import numpy as np
from t10_fitc import cval, n
om = float(sys.argv[2]); name = sys.argv[3]
s = 1 if name == 's+' else -1
sig = s*(-1)**n
xs = np.linspace(0.0, 1.0, 241)[:-1] + 1e-3j
vals = np.array([cval(name, om, x, sig) for x in xs])
la = np.log10(np.abs(vals))
med = np.median(la)
ext = []
L = len(xs)
for i in range(L):
    a, b, c = la[(i-1) % L], la[i], la[(i+1) % L]
    if b > a and b > c and b > med + 0.7: ext.append(('pole', round(xs[i].real, 4)))
    if b < a and b < c and b < med - 0.7: ext.append(('zero', round(xs[i].real, 4)))
print(name, om, ext)
