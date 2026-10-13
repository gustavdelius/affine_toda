"""Items 1 and 3 (ch. 41, prp-an-transmission (3), (5)), numerically from the chapter's formulas (numpy only):
 - sector determinants of eq-dt-reversed, q^r (1 + (-1)^r q X^h), so the type-I_+ matrix for V (its inverse on the
   sector -Q, offset -s, X -> -X_+) has det q^-r/(1 - q X_+^h): a pole, not a zero -> 'defined and invertible';
 - for that inverse, the entries of a column have moduli proportional to 1, t, ..., t^(h-1), t = |X_+| (in some
   cyclic order), with the common denominator |1 - q X_+^h|: at theta = theta_d (t = 1) all h^2 entries have the
   same modulus 1/|1 - q X_+^h|, not 1; the relative size of entry d is min(t^d, t^(d-h+1)) as for M_Vbar."""
import numpy as np
rng = np.random.default_rng(7)
ok = True
for h in (3, 4, 5):
    r = h - 1
    for om in (0.6, 1.37, 3.1):
        q = -np.exp(-1j*np.pi*om)
        s = rng.normal(size=h)*0.4; s -= s.mean()
        Q = np.zeros(h); Q[0] = 1.0                   # a V sector of type I_+: total weight Q = e_0, sum Q = 1
        for (th, thd) in ((0.0, 0.0), (0.3, 0.3), (0.7, 0.2), (-0.5, 0.4)):
            Xp = -np.exp(om*(thd - th - 1j*np.pi/2))
            Xr = -Xp                                 # eq-dt-reversed at X -> -X_+
            sm, Qm = -s, -Q                          # sector -Q, offset -s
            idx = {i: k for k, i in enumerate(range(h))}
            M = np.zeros((h, h), complex)
            for i in range(h):
                N = Qm + np.eye(h)[i]; n = N + sm        # Vbar input vb_i (x) |N>, total weight -e_i + N = -Q, sum N = 0
                M[i, i] += q**n[i]
                M[(i + 1) % h, i] += Xr*q**(n[i] - n[(i + 1) % h])
            det = np.linalg.det(M)
            ok &= abs(det - q**r*(1 + (-1)**r*q*Xr**h)) < 1e-10*max(1, abs(det))
            Minv = np.linalg.inv(M)
            ok &= abs(np.linalg.det(Minv) - q**(-r)/(1 - q*Xp**h)) < 1e-9*max(1, abs(np.linalg.det(Minv)))
            t = abs(Xp)
            mods = np.abs(Minv)
            if abs(th - thd) < 1e-12:
                ok &= np.allclose(mods, 1/abs(1 - q*Xp**h), rtol=1e-10)
                if h == 3 and om == 1.37 and th == 0.0:
                    print(f'h = 3, omega = 1.37, theta = theta_d: all |entries| = {mods[0,0]:.6f} = 1/|1 - q X_+^h| = {1/abs(1 - q*Xp**h):.6f}')
            for col in range(h):
                c = np.sort(mods[:, col]/mods[:, col].max())[::-1]
                want = np.sort([min(t**d, t**(d - h + 1)) for d in range(h)])[::-1]
                ok &= np.allclose(c, want, rtol=1e-9)
print('PASS' if ok else 'FAIL', ': items 1, 3 - type-I_+ V matrix: det q^-r/(1 - q X_+^h) (poles), equal moduli 1/|1 - q X_+^h| at theta = theta_d, relative sizes min(t^d, t^(d-h+1))')
