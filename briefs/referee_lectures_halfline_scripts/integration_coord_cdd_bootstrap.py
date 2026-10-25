# Test the soliton-bootstrap condition sigma(mu) = prod_{k=1}^n sigma(mu + T - k*omega), T = h*omega/2,
# for CDD factors P_x(mu) = (T/2 - x)_I (T/2 + x)_I, with (a)_I(mu) = sin(pi(mu+a)/(h omega)) / sin(pi(mu-a)/(h omega)).
import mpmath as mp
mp.mp.dps = 30
def blk(a, mu, h, w): return mp.sin(mp.pi*(mu+a)/(h*w))/mp.sin(mp.pi*(mu-a)/(h*w))
def P(x, mu, h, w):
    T = h*w/2
    return blk(T/2-x, mu, h, w)*blk(T/2+x, mu, h, w)
w = mp.sqrt(2); mu = mp.mpf('0.31') + mp.mpf('0.17')*1j
for n in range(1, 6):
    h = n+1; T = h*w/2
    for name, x in [("P_{w/2}", w/2), ("P_{(w+1)/2}", (w+1)/2)]:
        lhs = P(x, mu, h, w)
        rhs = mp.fprod(P(x, mu+T-k*w, h, w) for k in range(1, n+1))
        print(f"n={n} {name:12s} |lhs-rhs|={mp.nstr(abs(lhs-rhs),3):>10s}  lhs/rhs={mp.nstr(lhs/rhs,6)}")
