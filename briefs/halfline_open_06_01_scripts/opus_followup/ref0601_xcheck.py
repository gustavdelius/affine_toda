"""Cross-check the tau-method K (ref0601_a5_tauK.py) against the worker's ODE solver linrefl2.refl_matrix2 at delta=0.1
for the regular pattern and the target pattern ++---- (theta=0.6)."""
import sys, numpy as np, mpmath as mp
W = sys.argv[1]
sys.path.insert(0, W); import ref0601_a5_tauK as T
sys.path.insert(0, '/home/gustav/Git/affine_toda/briefs/halfline_open_06_01_scripts'); from linrefl2 import refl_matrix2
chi = 5*np.pi/3 + 0.1; Om = np.cosh(0.6); d = 4-2*np.sqrt(3)
def tf(x):
    j = np.arange(6); e = np.exp(np.sqrt(3)*x)
    return 1 + 2*d*np.cos(chi + 2*np.pi*j/3)*e + d*d/4*e*e
for C in [(1,-1,-1,1,-1,-1), (1,1,-1,-1,-1,-1), (1,1,-1,1,1,-1)]:
    Kt, op, _ = T.Kmatrix(5*mp.pi/3 + mp.mpf('0.1'), C, mp.cosh(mp.mpf('0.6')))
    Kt = np.array(Kt.tolist(), dtype=complex)
    Kw, lab = refl_matrix2(5, tf, C, Om, X=16.0)
    print(C, 'tau-method K=', np.round(Kt, 6).tolist(), '\n   linrefl2 K=', np.round(Kw, 6).tolist(),
          '\n   |K| entries tau:', np.round(np.abs(Kt), 6).tolist(), ' linrefl2:', np.round(np.abs(Kw), 6).tolist(), flush=True)
