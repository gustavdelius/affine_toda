"""Task 5: continued SP particle amplitudes L_a = K_{a -> h-a} (real coupling), with particle_reflection.py's blocks
and S-matrix (exec of its definitions).
 L_1(theta) = K0(mu) E(theta) [sigma_a(mu_a) sigma_b(mu_b)] at omega = -2/B, mu = -i T theta/pi (direct substitution),
 E = (sinh th - i sin phi)/(sinh th + i sin(phi + 2 pi l/h - pi)) with phi = pi A/T kept as a free parameter,
 L_n(theta) = K0 Ebar [sigma_b(mu_a) sigma_a(mu_b)] (the soliton-side Bbar_1 -> B_1 amplitude; the Z_h phase dropped).
 Higher L_a from the bootstrap L_{a+1}(th) = L_a(th + i pi/h) S_{a,h-1}(2th + i pi(1-a)/h) L_1(th - i pi a/h).
Checks: continuation of K0 to the book blocks; L_n from the bootstrap = soliton-side L_n; unitarity
L_a(th) L_{h-a}(-th) = 1; crossing L_a(i pi/2 - th) = S_aa(2 th) L_a(i pi/2 + th); bootstrap at every fusing
(a+b < h and a+b > h); B -> 0 limit."""
import sys
import mpmath as mp

src = open('/home/gustav/Git/affine_toda/textbook_code/particle_reflection.py').read()
ns = {}
exec(src[:src.index('THS = [')], ns)
blk, S = ns['blk'], ns['S']
mp.mp.dps = 30
I, PI = mp.mpc(0, 1), mp.pi
h_ = mp.mpf(1)/2


def bI(x, mu, T):
    return mp.sin(PI*(mu + x)/(2*T))/mp.sin(PI*(mu - x)/(2*T))


def K0(mu, T, om):
    d = (T - om)/2
    return (bI(T/2, mu, T)*bI(T/2 - h_, mu, T)*bI(d, mu, T)*bI(h_ - d, mu, T)*bI(h_ + d, mu, T)
            * bI(T - h_ - d, mu, T))


CDD = {
    'none': (lambda m, T, om: 1, lambda m, T, om: 1),
    'sigma_A': (lambda m, T, om: bI(-T/2, m, T)**2, lambda m, T, om: bI(-T/2, m, T)**2),
    'sigma_B': (lambda m, T, om: -1/(bI(T/2, m, T)**2*bI(T/2 - om/2, m, T)*bI(T/2 + om/2, m, T)),
                lambda m, T, om: -1/(bI(T/2, m, T)**2*bI(T/2 - om/2, m, T)*bI(T/2 + om/2, m, T))),
    '(om)_I pair': (lambda m, T, om: bI(om, m, T), lambda m, T, om: -bI(T - om, m, T)),
    '(-om/2)^2 pair': (lambda m, T, om: bI(-om/2, m, T)**2, lambda m, T, om: bI(T + om/2, m, T)**2),
}


def make_L(n, B, l, phi, cdd):
    h = n + 1
    om = -2/mp.mpf(B)
    T = h*om/2
    sa, sb = CDD[cdd]
    s2 = mp.sin(phi + 2*PI*l/h - PI)

    def L1(th, bar=False):
        mu = -I*T*th/PI
        ma, mb = mu - (T - 1)/2, mu + (T - 1)/2
        if not bar:
            E = (mp.sinh(th) - I*mp.sin(phi))/(mp.sinh(th) + I*s2)
            c = sa(ma, T, om)*sb(mb, T, om)
        else:
            E = (mp.sinh(th) - I*s2)/(mp.sinh(th) + I*mp.sin(phi))
            c = sb(ma, T, om)*sa(mb, T, om)
        return K0(mu, T, om)*E*c

    cache = {}

    def L(a):
        if a == 1:
            return lambda th: L1(th)
        if a not in cache:
            prev = L(a - 1)
            cache[a] = lambda th, prev=prev, a=a: (prev(th + I*PI/h)*S(a - 1, h - 1, n)(2*th + I*PI*(2 - a)/h, B)
                                                   * L1(th - I*PI*(a - 1)/h))
        return cache[a]
    return L, L1


THS = [mp.mpc('0.37', '0.11'), mp.mpc('-0.83', '0.29'), mp.mpc('1.21', '-0.47')]
if __name__ == '__main__':
    # 1. continuation of K0 to the book blocks
    err = 0
    for n in range(1, 8):
        h = n + 1
        for B in (mp.mpf('0.317'), mp.mpf('1.43')):
            om = -2/B
            T = h*om/2
            for th in THS:
                mu = -I*T*th/PI
                book = (blk(mp.mpf(h)/2, th, h)*blk(mp.mpf(h)/2 + B/2, th, h)*blk(mp.mpf(h)/2 - 1, th, h)
                        * blk(mp.mpf(h)/2 - 1 - B/2, th, h)*blk(mp.mpf(h)/2 + 1 + B/2, th, h)
                        / blk(mp.mpf(h)/2 - 1 + B/2, th, h))
                err = max(err, abs(K0(mu, T, om)/book - 1))
    print(f'1. K0 at omega = -2/B = (h/2)(h/2+B/2)(h/2-1)(h/2-1-B/2)(h/2+1+B/2)/(h/2-1+B/2), n = 1..7: {mp.nstr(err, 2)}')
    # 2-5. axioms
    for cdd in CDD:
        for n in (1, 2, 3, 4, 5):
            h = n + 1
            worst = {'Ln=bootstrap': 0, 'unitarity': 0, 'crossing': 0, 'bootstrap': 0}
            for l in sorted({0, 1, n}):
                for B in (mp.mpf('0.317'), mp.mpf('1.43')):
                    phi = mp.mpc('0.7', '0.2')
                    L, L1 = make_L(n, B, l, phi, cdd)
                    for th in THS:
                        worst['Ln=bootstrap'] = max(worst['Ln=bootstrap'], abs(L(n)(th)/L1(th, bar=True) - 1))
                        for a in range(1, h):
                            worst['unitarity'] = max(worst['unitarity'], abs(L(a)(th)*L(h - a)(-th) - 1))
                            worst['crossing'] = max(worst['crossing'], abs(L(a)(I*PI/2 - th)
                                                    / (S(a, a, n)(2*th, B)*L(a)(I*PI/2 + th)) - 1))
                            for b in range(1, h):
                                if a + b < h:
                                    c, ua, ub = a + b, PI*b/h, PI*a/h
                                elif a + b > h:
                                    c, ua, ub = a + b - h, PI*(h - b)/h, PI*(h - a)/h
                                else:
                                    continue
                                rhs = L(a)(th + I*ua)*S(a, h - b, n)(2*th + I*ua - I*ub, B)*L(b)(th - I*ub)
                                worst['bootstrap'] = max(worst['bootstrap'], abs(L(c)(th)/rhs - 1))
            print(f'cdd={cdd:15s} n={n}: ' + '  '.join(f'{k} {mp.nstr(v, 2)}' for k, v in worst.items()), flush=True)
