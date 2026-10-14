"""Is the K_3 pole of d_5 at Omega = 1 (theta = 3 i pi/8) a bound state?  Bound-state determinant of the
full linearized problem (all channels closed) on a fine grid around Omega^2 = 1, and over (0, m_min^2)."""
import numpy as np
import b2_dn_classical as b2
n = 5
psi, p0 = b2.vacuum(n)
for lo, hi, N in ((0.0, 1.1716 - 1e-4, 60), (0.95, 1.05, 41)):
    g = np.linspace(lo, hi, N)
    d = np.array([b2.bound_det(n, psi, p0, x) for x in g])
    ch = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1])))
    i = np.argmin(np.abs(d))
    print(f"Omega^2 in [{lo:.3f},{hi:.3f}]: {ch} sign changes; min |det| = {np.abs(d).min():.3e} at {g[i]:.4f}; det(1.0) ~ {b2.bound_det(n, psi, p0, 1.0):.3e}", flush=True)
