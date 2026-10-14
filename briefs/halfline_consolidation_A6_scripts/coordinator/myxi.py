import mpmath as mp
mp.mp.dps = 30
blk = lambda x, mu, T: mp.sin(mp.pi*(mu + x)/(2*T))/mp.sin(mp.pi*(mu - x)/(2*T))
H = mp.mpf(1)/2
for n, om in ((4, mp.mpf('0.37')), (5, mp.mpf('0.53'))):
    T = (n - 1)*om
    Xi = lambda m: blk(T/2, m, T)**2*blk(T/2 - om/2, m, T)*blk(T/2 + om/2, m, T)*blk(om/2, m, T)*blk(T - om/2, m, T)
    mu = mp.mpc('0.13', '0.21'); b = lambda x: blk(x, mu, T)
    KBp = (b(om/2 - T)*b(H - T)**2*b(-om/2 - H)*b(T/2 - om/2)*b(T/2 + om/2 + H)*b(T/2 - H)
           / (b(T/2)*b(H)*b(H - om/2)*b(T + om/2 - H)*b(T/2 - H + om/2)*b(T/2 + H - om/2)))
    KN = b(om/2 - T)*b(H)*b(T - om/2 - H)*b(T/2 - om/2)*b(T/2 - H)/(b(T/2)*b(T/2 - om/2 - H))
    print(n, 'Xi unitary', mp.nstr(abs(Xi(mu)*Xi(-mu) - 1), 3), ' Xi crossing', mp.nstr(abs(Xi(T/2 - mu)/Xi(T/2 + mu) - 1), 3),
          ' dressed breather', mp.nstr(abs(Xi(mu - T/2 + H)*Xi(mu + T/2 - H)*KBp/KN - 1), 3))
