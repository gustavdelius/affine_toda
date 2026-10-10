"""Quantum dimensions at roots of unity (chapter 27, exr-quantum-dims).

For q = e^{i pi r/p}, gcd(r,p)=1, the U_q(sl2) quantum dimensions [k]_q = sin(pi r k/p)/sin(pi r/p),
1 <= k <= p-1, are all positive iff r = +-1 mod 2p, alternate in sign iff r = p+-1 mod 2p, and the
ratios [k+1][k-1]/[k]^2 (2 <= k <= p-2) are all positive iff r = +-1 mod p.
Also checks the parameter identities used in the exercises of chapters 21-25.
"""
import numpy as np


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


ok = True
for p in range(5, 60):
    for r in range(1, 2*p):
        if np.gcd(r, p) != 1:
            continue
        qd = lambda k: np.sin(np.pi*r*k/p)/np.sin(np.pi*r/p)
        allpos = all(qd(k) > 0 for k in range(1, p))
        alt = all(qd(k)*(-1)**(k+1) > 0 for k in range(1, p))
        ratio = all(qd(k+1)*qd(k-1)/qd(k)**2 > 0 for k in range(2, p-1))
        ok &= allpos == (r % (2*p) in (1, 2*p-1)) and alt == (r % (2*p) in (p-1, p+1)) and ratio == (r % p in (1, p-1))
report('quantum-dimension sign criteria for 5 <= p < 60', ok)

# parameter identities (chapters 22-25), checked at random couplings
rng = np.random.default_rng(0); ok = True
for _ in range(20):
    w = rng.uniform(0.3, 5); b2 = 4*np.pi/(w+1)                     # omega = 4 pi/beta^2 - 1
    ok &= abs((3*w+1) - (12*np.pi/b2 - 2)) < 1e-12                  # g2: 2T = 4 pi h/beta^2 - h_vee
    ok &= abs((6*w+3) - (0.5*(4*np.pi*4/(b2/3) - 6))) < 1e-12        # d4^(3): beta_lit^2 = beta^2/3
    for n in range(2, 7):
        T = n*w + (n-1)/2; B = -1/(w+0.5)
        ok &= abs(2*T/(w+0.5) - (2*n + B)) < 1e-12                   # c_n: H = 2n + B
        ok &= abs(B + 2*b2/(8*np.pi - b2)) < 1e-12                   # B = -2 beta^2/(8 pi - beta^2)
        T2 = (2*n-1)*w + n - 1
        ok &= abs(T2/(w+0.5) - (2*n - 1 - b2/(8*np.pi - b2))) < 1e-12  # a_{2n-1}^(2) = lattice formula
    H = rng.uniform(12.5, 30); Hp = 1/(1/6 - 1/H)
    m1, m2, m3 = np.sin(np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(np.pi/Hp), np.sin(2*np.pi/H)*np.sin(2*np.pi/Hp)
    ok &= abs(m3/m1 - 2*np.cos(np.pi/H)) < 1e-12 and abs(m2/m1 - 2*np.cos(np.pi/H + np.pi/6)) < 1e-12
report('T, H and mass identities of chapters 22-25 (g2, d4^(3), c_n, a_{2n-1}^(2), CDS e6^(2))', ok)
