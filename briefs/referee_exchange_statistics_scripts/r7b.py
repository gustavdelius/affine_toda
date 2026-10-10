import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from spin_common import *
for n, alg in ((3, 'c'),):
    D, Wn, setup = make(n, alg); P = flip(D, D)
    for om in (0.3, 0.7, 1.2, 1.61, 0.45):
        q = np.exp(-1j*np.pi*om); Rp, Rm = setup(q), setup(-q)
        x = 1.7+0.4j; Cp, Cm = P@Rp(x), P@Rm(x)
        nzp, nzm = np.abs(Cp) > 1e-10, np.abs(Cm) > 1e-10
        r = Cm[nzp & nzm]/Cp[nzp & nzm]
        # eigenvalues comparison (Rcheck spectra are x-dependent phases; compare)
        ep = np.sort_complex(np.round(np.linalg.eigvals(Cp), 8)); em = np.sort_complex(np.round(np.linalg.eigvals(Cm), 8))
        print(f"om={om}: same sparsity {np.array_equal(nzp, nzm)}, ratio |r| in [{np.abs(r).min():.3f},{np.abs(r).max():.3f}], "
              f"max|Im r| {np.abs(r.imag).max():.2e}, spectra equal {np.allclose(ep, em, atol=1e-6)}, "
              f"Herm check |Cp Cp^dag - 1| {np.abs(Cp@Cp.conj().T - np.eye(D*D)).max():.1e}")
