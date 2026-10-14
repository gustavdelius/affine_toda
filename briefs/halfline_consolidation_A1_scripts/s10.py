import sys
sys.path.insert(0, '/home/gustav/Git/affine_toda/textbook_code')
src = open('/home/gustav/Git/affine_toda/textbook_code/particle_reflection.py').read()
exec(src.split('# ---------------------------------------------------------------- 1. blocks')[0])

def blkI(a, mu, T):
    return mp.sin(PI*(mu + a)/(2*T))/mp.sin(PI*(mu - a)/(2*T))

def Kminus(a, n):
    h = n + 1
    return Prod(h, sign=-1).mul((a, 0)).mul((h - a, 0))*K_uniform(a, n)

ok_cont = True
for n in range(1, 8):
    h = n + 1
    for B in BS:
        w = -2/B; T = h*w/2
        for th in THS:
            mu = -I*T*th/PI
            v = -blkI(w/2, mu, T)*blkI(mp.mpf(1)/2 - T, mu, T)*blkI(-w/2 - mp.mpf(1)/2, mu, T)
            ref = -blk(1, th, h)*blk(h - B/2, th, h)*blk(-1 + B/2, th, h)
            ok_cont &= abs(v - ref) < TOL and abs(v - Kminus(1, n)(th, B)) < TOL
print('continuation', ok_cont)
res = [axioms(lambda a, n=n: Kminus(a, n), n, 'minus') for n in range(1, 8)]
print('axioms', res)
print('B->0', all(abs(Kminus(a, n)(th, 0) - 1) < TOL for n in range(1, 8) for a in range(1, n + 1) for th in THS))
print('dual', all(close(Kminus(a, n).dual(), Prod(n + 1, sign=-1).mul((a, 0)).mul((n + 1 - a, 0))*K_neumann(a, n))
                  for n in range(1, 8) for a in range(1, n + 1)))
G64 = Prod(3, {(1, 0): 1, (-1, Fr(1, 2)): 1, (3, Fr(-1, 2)): 1}, sign=-1)
print('G98 (6.4)', close(Kminus(1, 2), G64), close(Kminus(2, 2), G64))
# n = 1: sinh-Gordon; what is K^(-)_1 there?  -(1)(1)(2-B/2)(-1+B/2)
print('n=1 value at sample', [mp.nstr(Kminus(1, 1)(th, mp.mpf('0.6')), 8) for th in THS[:2]])
