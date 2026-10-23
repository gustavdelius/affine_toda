"""Conjectured closed form of the dressed source of VJS (37):
   s_phys_a = cosh(w c_a)/cosh(w c),  a<n;   s_phys_n = 1/(2 cosh(w c)),
   c = ((2n-1)pi - 2n g)/4,  c_a = c - a (pi-g)/2.  Test at random real and complex w."""
import mpmath as mp, random
exec(open('bethe_masses.py').read().split('def check_n2')[0])
random.seed(1); worst = 0
for n in range(2, 9):
    for trial in range(6):
        g = mp.mpf(random.uniform(0.05, 1.5)); w = mp.mpc(random.uniform(-6, 6), random.uniform(-0.3, 0.3))
        v = sphys(n, g, w); c = ((2*n-1)*mp.pi - 2*n*g)/4
        pred = [mp.cosh(w*(c - a*(mp.pi - g)/2))/mp.cosh(w*c) for a in range(1, n)] + [1/(2*mp.cosh(w*c))]
        worst = max(worst, max(abs(v[j] - pred[j]) for j in range(n)))
print("max deviation over n=2..8, random gamma in (0.05,1.5), random w:", mp.nstr(worst, 3))
