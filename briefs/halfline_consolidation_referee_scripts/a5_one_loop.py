"""Referee check of the one-loop attributions of the consolidation round (A5 a-c, B1 f).

Book blocks: (x) = sinh(th/2 + i pi x/2h)/sinh(th/2 - i pi x/2h), {x} = (x-1)(x+1)/((x-1+B)(x+1-B)).
Kim's block [x] = (x-1/2)(x+1/2)/((x-1/2+B/2)(x+1/2-B/2))  (Kim hep-th/9506031, Kim-Cho hep-th/9505138).
B = (beta^2/2pi)/(1+beta^2/4pi), so beta^2 = 2 pi B + O(B^2).

d_n amplitudes of the book (ch. 35):
  K_1   = (h-B/2)(h/2-1)(h/2+1-B/2) / [(h-1)(1-B/2)(h/2)(h/2-B/2)]
  K^N_1 = K_1|_{B->2-B},  K^(-)_1 = K_1/K_1^cl,  K_1^cl = -(h/2-1)(h/2+1)/[(1)(h-1)(h/2)^2]
Checks
  1. d_4: O(B) of K^N_1 and of K^N_2 = K^N_1(th+i pi/6) K^N_1(th-i pi/6) S_11(2th), S_11={1}{5},
     against Kim-Cho (29), (40) with i beta^2/24 = i pi B/12;
  2. K^(-)_1(th) and K^(-)_1(th+i pi) against Kim-Cho (29) at O(B);
  3. Kim-Cho's exact K_1 = [1/2][3/2][5/2]^2[7/2][9/2] and Kim's d_n vector rule (fig. 7):
     agreement with K^N_1 at O(B), disagreement at O(B^2) (n = 4..7);
  4. Kim's d_n vector rule = CAF R~_1(th + i pi) (their (dn1)), as functions, n = 4..7.
All at complex rapidities chosen here.
"""
import mpmath as mp

mp.mp.dps = 40


def blk(x, th, h):
    return mp.sinh(th / 2 + 1j * mp.pi * x / (2 * h)) / mp.sinh(th / 2 - 1j * mp.pi * x / (2 * h))


def curly(x, th, h, B):
    return blk(x - 1, th, h) * blk(x + 1, th, h) / (blk(x - 1 + B, th, h) * blk(x + 1 - B, th, h))


def kimblk(x, th, h, B):
    return blk(x - mp.mpf(1) / 2, th, h) * blk(x + mp.mpf(1) / 2, th, h) / (
        blk(x - mp.mpf(1) / 2 + B / 2, th, h) * blk(x + mp.mpf(1) / 2 - B / 2, th, h))


def K1(th, h, B):
    b = lambda x: blk(x, th, h)
    return b(h - B / 2) * b(h / 2 - 1) * b(h / 2 + 1 - B / 2) / (b(h - 1) * b(1 - B / 2) * b(h / 2) * b(h / 2 - B / 2))


def K1cl(th, h):
    b = lambda x: blk(x, th, h)
    return -b(h / 2 - 1) * b(h / 2 + 1) / (b(1) * b(h - 1) * b(h / 2) ** 2)


def KN1(th, h, B):
    return K1(th, h, 2 - B)


def Km1(th, h, B):
    return K1(th, h, B) / K1cl(th, h)


def kim_vector(th, h, B):
    p = kimblk(mp.mpf(h - 1) / 2, th, h, B)
    for k in range(1, 2 * h - 2, 2):
        p *= kimblk(mp.mpf(k) / 2, th, h, B)
    return p


def caf_norm(x, th, h, B):
    # CAF ||x|| for simply laced (h = H): <x>_th/<x>_{-th}
    def ang(t):
        s = lambda a: mp.sinh((t + 1j * mp.pi * a / (2 * h)) / 2)
        return s(x - 1) * s(x + 1) / (s(x - 1 + B) * s(x + 1 - B))
    return ang(th) / ang(-th)


def caf_R1_shift(th, h, B):
    """CAF (dn1) with j=1: R~_1(th+i pi) = ||h-1|| prod_{k=1}^{h-1} ||2k-1||"""
    p = caf_norm(h - 1, th, h, B)
    for k in range(1, h):
        p *= caf_norm(2 * k - 1, th, h, B)
    return p


def kimcho_K1(th, b2):
    sh, ch = mp.sinh(th), mp.cosh(th)
    return 1 + 1j * b2 / 24 * (sh / ch - sh / (ch - 1) + 2 * sh / (2 * ch + mp.sqrt(3)) - 2 * sh / (2 * ch - 1))


def kimcho_K2(th, b2):
    sh, ch = mp.sinh(th), mp.cosh(th)
    s3 = mp.sqrt(3)
    return 1 + 1j * b2 / 24 * (sh / ch - sh / (ch - 1) - 4 * sh / (2 * ch - s3) + 2 * sh / (2 * ch + s3) + 2 * sh / (2 * ch + 1))


def coeffs(f, order=2):
    """Taylor coefficients in B of log f(B) at B=0 (f(0) = 1 assumed for B-coefficient extraction)"""
    eps = mp.mpf('1e-12')
    g = lambda B: mp.log(f(B))
    c1 = (g(eps) - g(-eps)) / (2 * eps)
    c2 = (g(eps) - 2 * g(0) + g(-eps)) / (2 * eps ** 2)
    return g(0), c1, c2


if __name__ == '__main__':
    ok = True
    pts = [mp.mpc('0.37', '0.21'), mp.mpc('1.3', '-0.4'), mp.mpc('-0.8', '0.95')]
    h = 6
    print("1. d_4 one loop: O(B) coefficient of log K^N vs Kim-Cho (29), (40)")
    for th in pts:
        g0, c1, c2 = coeffs(lambda B: KN1(th, h, B))
        kc = (kimcho_K1(th, 2 * mp.pi * mp.mpf('1e-12')) - 1) / mp.mpf('1e-12')     # d/dB at 0
        KN2 = lambda B: KN1(th + 1j * mp.pi / 6, h, B) * KN1(th - 1j * mp.pi / 6, h, B) * curly(1, 2 * th, h, B) * curly(5, 2 * th, h, B)
        g02, c12, _ = coeffs(KN2)
        kc2 = (kimcho_K2(th, 2 * mp.pi * mp.mpf('1e-12')) - 1) / mp.mpf('1e-12')
        print(f"   th={mp.nstr(th, 4)}: K^N_1(B=0)={mp.nstr(mp.exp(g0), 6)}, dlogK^N_1={mp.nstr(c1, 12)}, Kim-Cho {mp.nstr(kc, 12)};"
              f" K^N_2: {mp.nstr(c12, 12)} vs {mp.nstr(kc2, 12)}")
        ok &= abs(mp.exp(g0) - 1) < 1e-20 and abs(c1 - kc) < 1e-8 and abs(c12 - kc2) < 1e-8
    print("2. K^(-)_1 vs Kim-Cho at O(B): at th and at th + i pi")
    for th in pts:
        g0, c1, _ = coeffs(lambda B: Km1(th, h, B))
        g0s, c1s, _ = coeffs(lambda B: Km1(th + 1j * mp.pi, h, B))
        kc = (kimcho_K1(th, 2 * mp.pi * mp.mpf('1e-12')) - 1) / mp.mpf('1e-12')
        print(f"   th={mp.nstr(th, 4)}: K^(-)_1(B=0)={mp.nstr(mp.exp(g0), 6)} d={mp.nstr(c1, 10)};"
              f" shifted: K(B=0)={mp.nstr(mp.exp(g0s), 6)} d={mp.nstr(c1s, 10)}; Kim-Cho d={mp.nstr(kc, 10)}")
        ok &= abs(c1s - kc) < 1e-8 and abs(c1 - kc) > 1e-3
    print("3. Kim's self-dual vector amplitudes vs K^N_1: O(B) and O(B^2) of log")
    for n in (4, 5, 6, 7):
        hh = 2 * n - 2
        for th in pts[:2]:
            a = coeffs(lambda B: KN1(th, hh, B))
            b = coeffs(lambda B: kim_vector(th, hh, B))
            print(f"   d_{n} th={mp.nstr(th, 4)}: B^0 {mp.nstr(mp.exp(a[0]), 5)} {mp.nstr(mp.exp(b[0]), 5)};"
                  f" B^1 {mp.nstr(a[1], 10)} {mp.nstr(b[1], 10)}; B^2 {mp.nstr(a[2], 8)} {mp.nstr(b[2], 8)}")
            ok &= abs(a[1] - b[1]) < 1e-8 and abs(a[2] - b[2]) > 1e-3 and abs(mp.exp(a[0]) - mp.exp(b[0])) < 1e-15
    print("4. Kim's d_n vector rule vs CAF R~_1(th + i pi) at B = 0.37")
    for n in (4, 5, 6, 7):
        hh = 2 * n - 2
        d = max(abs(kim_vector(th, hh, mp.mpf('0.37')) - caf_R1_shift(th, hh, mp.mpf('0.37'))) for th in pts)
        print(f"   d_{n}: max diff {mp.nstr(d, 3)}")
        ok &= d < 1e-25
    print("PASS" if ok else "FAIL", "one-loop attributions (Kim-Cho d_4, Kim d_n vector, CAF)")
