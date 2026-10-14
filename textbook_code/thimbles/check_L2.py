"""a_2, L=2 (n=4): identify the dominant contributing saddle from the small-hbar asymptotics of the exact Z."""
import numpy as np, sys
from lie import An
from model import Chain
from exactZ import Z_grid
from flows2 import realmin_ReS
g = An(2); lam = g.fundamental_weights_rep(1)[(1,)]
ch = Chain(g, 2, 1.0, lam)
cps = np.load('cps_a2L2_k1.00.npy')
S = np.array([ch.action(u) for u in cps])
print("min_R Re S =", realmin_ReS(ch, nstart=30), " Imax bound on R = %.3f" % (2*1.0*3*np.sqrt(3)/2))
for i in range(7):
    u = cps[i]; th = 2*np.pi*lam - np.conj(u[::-1])
    j = np.argmin([np.max(np.abs(c - th)) for c in cps])
    print("S=%.5f%+.5fi  u/2pi (eps basis) site1=%s site2=%s   Theta'-partner index %d" % (S[i].real, S[i].imag,
          np.round(g.eps_coords(u[0])/(2*np.pi), 3), np.round(g.eps_coords(u[1])/(2*np.pi), 3), j))
U0 = cps[0]; mu, V = ch.tangent(U0); c1 = ch.two_loop_c1(U0)
for hbar in (0.4, 0.3, 0.2, 0.15):
    Z = Z_grid(ch, hbar, h=0.025)
    Z1 = ch.one_loop(U0, hbar, V, mu)
    print("hbar=%.2f Z_exact=%.8e%+.1ei  Z/Z1loop(S=8.728)=%.6f%+.1ei   1+hbar c1=%.6f" % (hbar, Z.real, Z.imag, (Z/Z1).real, (Z/Z1).imag, (1+hbar*c1).real)); sys.stdout.flush()
