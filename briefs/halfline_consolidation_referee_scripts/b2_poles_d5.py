"""Direct check of the physical-strip poles of prod_b S_ab(theta, B=1) for d_5 (and d_4): evaluate
|prod S| at theta = i pi z/h + eps for eps -> 0; a simple pole grows like 1/eps."""
import numpy as np
import b2_dn_classical as b2
sm = b2.load_smatrix()
for n in (4, 5):
    h, bl, pfm, _ = sm.setup('d', n)
    al, mk = b2.dn_roots(n); masses = np.sqrt(np.linalg.eigvalsh((al.T*mk) @ al)); mmin = masses.min()
    for a in range(n):
        m = pfm[a]*mmin
        for z in range(1, (h + 1)//2):
            v1, v2 = [abs(np.prod([sm.S(bl[(a, b)], h, 1.0, 1j*np.pi*z/h + e) for b in range(n)])) for e in (1e-3, 1e-5)]
            if v2 > 1e3:
                print(f"d_{n} node {a} mass {m:.4f}: pole at theta = i pi {z}/{h} (|prodS| {v1:.1e} -> {v2:.1e}); "
                      f"2|p| = {2*m*np.sin(np.pi*z/h):.6f}, Omega^2 = {m**2*np.cos(np.pi*z/h)**2:.6f}")
    print(f"d_{n} masses {np.round(sorted(set(np.round(masses, 6))), 6)}")
