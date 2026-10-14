"""Independent check (referee, 2026-10-12) of the pole orders at nu_0 and nu_{-1} of the Dirichlet
soliton reflection factors on MD's excited states, coded from Mattsson-Dorey (hep-th/0008071, bsgdht.tex)
directly: P^{+-}(u) Gamma-product (TeX l.288-300), R_0 (l.296-305), a^1_n product form (l.975-980),
a^0_n = a^1_n(xi -> xi_bar), P_{|c;n1..nk>} = P_(c) a^1_{n1} a^0_{n2} ... (l.1213-1218), P_(1)^{+-} = Pbar^{-+} (l.958).
Products are truncated at L terms: the tail is analytic and non-vanishing near the strip points, so it does
not change pole orders.  Order = -d log|f| / d log(eps) measured between eps=1e-7 and 1e-9."""
import mpmath as mp
mp.mp.dps = 40
PI = mp.pi
L = 40

def lgP(sign, xi, u, lam):
    # log P^{sign}(u) without R_0 (MD l.288-292)
    s = 0
    xp = xi/PI
    for l in range(1, L+1):
        def t(uu):
            v = lam*uu/PI
            return (mp.loggamma(0.5 + 2*l*lam + sign*xp + v) + mp.loggamma(0.5 + (2*l-2)*lam - sign*xp + v)
                    - mp.loggamma(0.5 + (2*l-1)*lam + xp + v) - mp.loggamma(0.5 + (2*l-1)*lam - xp + v))
        s += t(u) - t(-u)
    return s

def lgR0(u, lam):
    s = 0
    for k in range(1, L+1):
        def t(uu):
            v = 2*lam*uu/PI
            return (mp.loggamma(1 + lam*(4*k-4) - v) + mp.loggamma(4*lam*k - v)
                    - mp.loggamma(lam*(4*k-3) - v) - mp.loggamma(1 + lam*(4*k-1) - v))
        s += t(u) - t(-u)
    return s

def blk(x, u):   # MD's (x) = sin(u/2 + pi x/2)/sin(u/2 - pi x/2)
    return mp.sin(u/2 + PI*x/2)/mp.sin(u/2 - PI*x/2)

def a1(n, xi, u, lam):   # MD l.975-980
    r = mp.mpf(1)
    for x in range(1, n+1):
        A = xi/(lam*PI) + 1/(2*lam) - mp.mpf(x)/lam
        B = xi/(lam*PI) - 1/(2*lam) - mp.mpf(x)/lam
        r *= blk(A, u)*blk(B, u)/(blk(A+1, u)*blk(B+1, u))
    return r

def P_state(sign, labels, xi, u, lam):
    xib = PI*(lam+1) - xi
    k = len(labels)
    if k % 2 == 0:
        val = mp.exp(lgP(sign, xi, u, lam) + lgR0(u, lam))
    else:   # P^{+-}_{|1;..>} = Pbar^{-+}
        val = mp.exp(lgP(-sign, xib, u, lam) + lgR0(u, lam))
    for i, n in enumerate(labels, start=1):
        val *= a1(n, xi if i % 2 == 1 else xib, u, lam)
    return val

def order(f, u0):
    e1, e2 = mp.mpf('1e-7'), mp.mpf('1e-9')
    o = []
    for d in (1, -1):
        r = mp.log(abs(f(u0 + d*e2))/abs(f(u0 + d*e1)))/mp.log(e1/e2)
        o.append(r)
    return o

def nu(n, xi, lam): return xi/lam - (2*n+1)*PI/(2*lam)
def w(n, xi, lam): return PI - xi/lam - (2*n-1)*PI/(2*lam)

def predict(labels, sign, pos):
    """analytic count: prefactor + a-factor orders (my derivation)"""
    k = len(labels)
    if pos == 0:
        pre = 1 if k % 2 == 1 else (1 if sign == +1 else 0)
        return pre + sum(1 for n in labels if n >= 1)
    if pos == -1:
        pre = 1 if k % 2 == 1 else (0 if sign == +1 else -1)
        ev = labels[1::2]
        return pre + sum(2 if n > 1 else (1 if n == 1 else 0) for n in ev)

if __name__ == '__main__':
    for lam, xi, states in [(6.3, 5.0, [(), (0,), (1,), (0, 5), (0, 5, 1)]),
                            (12.7, 15.0, [(), (0,), (1,), (2,), (0, 6), (1, 6), (1, 6, 2), (0, 6, 2), (1, 6, 2, 9)])]:
        xi = mp.mpf(xi); lam = mp.mpf(lam)
        print(f'lambda={lam} xi=eta={xi}  MD range 0<xi<{float(PI*(lam+1)/2):.4f};  eta_N-pi={float(PI*(lam-1)/2):.4f}')
        print('  nu_0=%.5f nu_-1=%.5f (pi/2=%.5f)' % (nu(0, xi, lam), nu(-1, xi, lam), PI/2),
              ' nu_n:', [round(float(nu(n, xi, lam)), 4) for n in range(0, 5)],
              ' w_n:', [round(float(w(n, xi, lam)), 4) for n in range(0, 11)])
        for S in states:
            # validity: nu_{n1} > w_{n2} > nu_{n3} > ... all in (0, pi/2)
            pos = [nu(n, xi, lam) if i % 2 == 0 else w(n, xi, lam) for i, n in enumerate(S)]
            ok = all(0 < p < PI/2 for p in pos) and all(pos[i] > pos[i+1] for i in range(len(pos)-1))
            for sign in (+1, -1):
                for p in (0, -1):
                    u0 = nu(p, xi, lam)
                    o = order(lambda u: P_state(sign, S, xi, u, lam), u0)
                    print(f'  state c={len(S)%2} {S!s:14} valid={ok!s:5} P{"+" if sign>0 else "-"} at nu_{p}: '
                          f'measured {float(o[0]):+.4f}/{float(o[1]):+.4f}  predicted {predict(S, sign, p):+d}')
