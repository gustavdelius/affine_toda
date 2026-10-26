"""Integration agent C: does DG's soliton-bootstrap condition sigma(mu) = prod_{k=1}^n sigma(mu + T - k omega)
(DG hep-th/9904002 eq. (ext); book ch. 34, before eq-an-cdd) hold for the CDD factors P_x = (T/2-x)_I (T/2+x)_I,
x = omega/2 (the book's choice in A^(+)) and x = (omega+1)/2 (the alternative of ch. 34), T = h omega/2?
Control: DG's sigma_c = -(c omega/2)_I (T - c omega/2)_I must pass for every n.
Result (printed): both P_x pass for odd n (h even) and fail for even n (h odd), and so does their ratio."""
import mpmath as mp
mp.mp.dps = 40


def blkI(a, mu, h, om):
    return mp.sin(mp.pi * (mu + a) / (h * om)) / mp.sin(mp.pi * (mu - a) / (h * om))


def Px(mu, x, h, om):
    T = h * om / 2
    return blkI(T / 2 - x, mu, h, om) * blkI(T / 2 + x, mu, h, om)


def sig(mu, c, h, om):
    return -blkI(c * om / 2, mu, h, om) * blkI(h * om / 2 - c * om / 2, mu, h, om)


def defect(f, mu, h, om):
    T = h * om / 2
    return f(mu) / mp.fprod(f(mu + T - k * om) for k in range(1, h))


omg = mp.sqrt(2) + mp.mpf(1) / 7
ok = True
for h in (3, 4, 5, 6):
    n = h - 1
    for mu in (mp.mpc(0.23, -0.17), mp.mpc(-0.41, 0.3)):
        dg = max(abs(defect(lambda m: sig(m, c, h, omg), mu, h, omg) - 1) for c in range(1, n + 1))
        d1 = defect(lambda m: Px(m, omg / 2, h, omg), mu, h, omg)
        d2 = defect(lambda m: Px(m, (omg + 1) / 2, h, omg), mu, h, omg)
        expect = (n % 2 == 1)
        ok &= dg < 1e-25 and ((abs(d1 - 1) < 1e-25) == expect) and ((abs(d2 - 1) < 1e-25) == expect) \
            and ((abs(d1 / d2 - 1) < 1e-25) == expect)
        print(f"n={n} mu={mp.nstr(mu, 3)}: DG sigma_c defect-1 {mp.nstr(dg, 2)}; P_(om/2) defect {mp.nstr(d1, 6)}; "
              f"P_((om+1)/2) defect {mp.nstr(d2, 6)}")
print(("PASS" if ok else "FAIL") + " soliton-bootstrap condition: DG sigma_c pass for all n; P_(om/2), P_((om+1)/2) and "
      "their ratio pass exactly for odd n (n=3,5) and fail for even n (n=2,4)")
