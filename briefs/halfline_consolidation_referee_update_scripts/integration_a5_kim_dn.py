"""A5(b,c): Kim's d_n vector rule (hep-th/9506031, (20) and Fig. 7: exponents 1,...,1,2,1,...,1 over k = 1,3,..,2h-3,
middle number 2) vs K^N_1 (eq-dn-neumann-amplitude) at O(B), n = 4..7; and vs CAF R~_1(theta + i pi) = ||h-1|| prod_{k=1}^{h-1} ||2k-1||
(hep-th/0308004 (4.19)) with ||2x|| = Kim's [x]. Also: equal at B=0, differ at B = 0.6."""
import mpmath as mp
mp.mp.dps = 30
I = mp.mpc(0, 1)
blk = lambda x, th, h: mp.sinh((th + I*mp.pi*x/h)/2)/mp.sinh((th - I*mp.pi*x/h)/2)
def kimb(x, th, h, B):
    half = mp.mpf(1)/2
    return blk(x-half, th, h)*blk(x+half, th, h)/(blk(x-half+B/2, th, h)*blk(x+half-B/2, th, h))
def kim_vec(th, n, B):
    h = 2*n - 2
    ks = list(range(1, 2*h - 2, 2))          # 1, 3, ..., 2h-3
    mid = len(ks)//2
    r = mp.mpc(1)
    for i, k in enumerate(ks):
        r *= kimb(mp.mpf(k)/2, th, h, B)**(2 if i == mid else 1)
    return r, ks[mid]
def KN1(th, h, B):
    b = lambda x: blk(x, th, h)
    return b(h/2-1)*b(h/2+B/2)*b(h-1+B/2)/(b(h-1)*b(h/2)*b(h/2-1+B/2)*b(B/2))
for n in (4, 5, 6, 7):
    h = 2*n - 2
    worst, fin = 0, []
    for th in (mp.mpc('0.37', '0.21'), mp.mpc('-0.8', '0.5')):
        f = lambda B: kim_vec(th, n, B)[0]
        g = lambda B: KN1(th, h, B)
        worst = max(worst, abs(f(0) - 1), abs(g(0) - 1), abs(mp.diff(f, 0) - mp.diff(g, 0)))
        fin.append(abs(f(mp.mpf('0.6'))/g(mp.mpf('0.6')) - 1))
    print(f'd_{n}: middle k = {kim_vec(0.3, n, 0)[1]} (= h-1 = {h-1}); Kim vector vs K^N_1: value and B-derivative at 0 {mp.nstr(worst, 3)}; '
          f'|ratio-1| at B=0.6: {[mp.nstr(x, 3) for x in fin]}')
