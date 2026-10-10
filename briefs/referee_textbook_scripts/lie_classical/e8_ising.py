"""e_8^(1) Toda S_11 = {1}{11}{19}{29} (h=30) factorizes as S_min * CDD(B); S_min = Zamolodchikov's Ising-in-field S_11."""
import numpy as np
h = 30
blk = lambda x, th: np.sinh(th/2 + 1j*np.pi*x/(2*h))/np.sinh(th/2 - 1j*np.pi*x/(2*h))
curly = lambda x, th, B: blk(x-1, th)*blk(x+1, th)/(blk(x-1+B, th)*blk(x+1-B, th))
f = lambda a, th: np.tanh((th + 1j*np.pi*a)/2)/np.tanh((th - 1j*np.pi*a)/2)
th = np.linspace(-3, 3, 13) + 0.1j
Smin = np.prod([blk(x-1, th)*blk(x+1, th) for x in (1, 11, 19, 29)], axis=0)
Szam = f(2/3, th)*f(2/5, th)*f(1/15, th)
print('S_min == Zamolodchikov f_{2/3} f_{2/5} f_{1/15}:', np.allclose(Smin, Szam))
for B in (0.0, 0.3, 1.0):
    St = np.prod([curly(x, th, B) for x in (1, 11, 19, 29)], axis=0)
    print(f'B={B}: S_Toda/S_min = 1/prod (x-1+B)(x+1-B); S_Toda==1? {np.allclose(St,1)}; S_Toda==S_min? {np.allclose(St, Smin)}')
