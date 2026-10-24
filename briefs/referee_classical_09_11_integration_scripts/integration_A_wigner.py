"""Item 12 (integration A, Part III): sign of the leading semiclassical phase of the
a_n^(1) soliton S-matrix against ch. 10's d delta/d theta = (2h/beta^2) ln A_11(theta).

Sources, read from the PDFs:
- Book eq-sg-s0 (ch. 3) and eq-sg-smatrix: S_T = sinh(lam th)/sinh(lam(i pi - th)) S0, lam = 4 pi/beta^2 - 1.
- Book sec-scalar-factor (ch. 21): F = c f, which for h = 2 (T = a = lam) must give S0.
- Hollowood hep-th/9203076 eq. (4.5) (PDF p. 23), his n = book h, his lambda = 4 pi/beta^2 + O(1)
  (his (5.6), (5.10); his beta^2 is the book's). S_{ik->ki} is the transmission amplitude.
Checks:
1. book c f at h = 2 reproduces -i dln S0/dtheta (calibration of the ch. 21 formula);
2. Hollowood (4.5) at n = 2 reproduces the book's S_T exactly (so no inversion or theta -> -theta
   between his conventions and the book's);
3. beta^2/(2h) (-i dln S_{ik->ki}/dtheta) -> ln A_11(theta) < 0 as lambda -> inf, for h = 2, 3, 4,
   the error shrinking like 1/lambda (book sign; Hollowood's own section 5 has the opposite sign).
"""
import mpmath as mp
mp.mp.dps = 20
half = mp.mpf(1)/2

def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)

def sg_S0(theta, lam):                     # -i d ln S0/d theta from eq-sg-s0, xi = pi/lam
    xi = mp.pi/lam
    g = lambda s: mp.sinh((mp.pi-xi)*s/2)*mp.cos(theta*s)/(mp.sinh(xi*s/2)*mp.cosh(mp.pi*s/2))
    return -mp.quadosc(g, [0, mp.inf], omega=theta)

def book_cf_h2(theta, lam):                 # -i d ln F/d theta, F = c f of ch. 21 with h = 2, T = a = lam
    T = a = lam; t = T*theta/(1j*mp.pi)
    p1 = lambda s: mp.digamma(s-a) + mp.digamma(1+s+a)
    term = lambda m: p1(t+2*T*(m-1)) - p1(-t+2*T*(m-half)) - p1(t+2*T*(m-half)) + p1(-t+2*T*m)
    d = mp.pi*mp.cot(mp.pi*(t-a)) + mp.nsum(term, [1, mp.inf])
    return -1j*d*T/(1j*mp.pi)

def holl_ik(theta, n, lam):                 # -i d ln S_{ik->ki}/d theta, Hollowood (4.5); real part = phase
    s = 1j*theta/(2*mp.pi); y = lambda j: j/(n*lam); r = mp.mpf(1)/n
    term = lambda j: (mp.digamma(s+1+y(j)) + mp.digamma(s+1+y(j-1)) + mp.digamma(-s+y(j)) + mp.digamma(-s+y(j-1))
                      - mp.digamma(-s+r+y(j)) - mp.digamma(-s-r+y(j-1)) - mp.digamma(s+1+r+y(j)) - mp.digamma(s+1-r+y(j-1)))
    return -1j*1j/(2*mp.pi)*mp.nsum(term, [1, mp.inf])   # complex: S_T alone is not unitary

th0 = mp.mpf('0.83'); lams = [mp.sqrt(3), mp.e]
d1 = max(abs(book_cf_h2(th0, l) - sg_S0(th0, l)) for l in lams)
report(f'ch. 21 c f at h = 2 equals S0 of eq-sg-s0 (phase derivative, lam = sqrt3, e)  [diff {mp.nstr(d1, 3)}]', d1 < 1e-12)
sT = lambda th, l: sg_S0(th, l) - 1j*mp.diff(lambda x: mp.log(mp.sinh(l*x)/mp.sinh(l*(1j*mp.pi-x))), th)
d2 = max(abs(holl_ik(th0, 2, l) - sT(th0, l)) for l in lams)
report(f'Hollowood (4.5) at n = 2 equals the book S_T (phase derivative, lam = sqrt3, e)  [diff {mp.nstr(d2, 3)}]', d2 < 1e-12)

ok = True
for h in (2, 3, 4):
    for th in (mp.mpf('0.7'), mp.sqrt(2)):
        lnA = mp.log((mp.cosh(th)-1)/(mp.cosh(th)-mp.cos(2*mp.pi/h)))
        v = [4*mp.pi/(l+1)/(2*h)*mp.re(holl_ik(th, h, l)) for l in (20*mp.sqrt(2), 40*mp.sqrt(2))]
        e = [abs(x - lnA) for x in v]
        print(f'   h={h} theta={mp.nstr(th, 4)}: beta^2/(2h) * phase derivative = {mp.nstr(v[0], 8)}, {mp.nstr(v[1], 8)} '
              f'(lam = 28.3, 56.6); ln A_11 = {mp.nstr(lnA, 8)}')
        ok &= lnA < 0 and v[1] < 0 and e[1] < 0.6*e[0] and e[1] < 0.03*abs(lnA)
report('Hollowood (4.5): -i dlnS/dtheta -> +(2h/beta^2) ln A_11 < 0 for h = 2, 3, 4 (the sign of ch. 10)', ok)
