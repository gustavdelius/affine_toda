"""A4: d_n uniform vacuum (prp-dn-uniform-vacuum, ch. 35).  Units m = beta = 1.
V = sum_j n_j (exp(alpha_j.phi) - 1)                       (eq-imaginary-potential, real coupling)
B = 2 sum_j C_j sqrt(2 n_j/alpha_j^2) (exp(alpha_j.phi/2) - 1)   (eq-boundary-potential), alpha^2 = 2, C = +1
Field equation phi'' = grad V; boundary condition phi'(0) = -grad B(phi(0)).
1. Diagram-invariant subspace and reduced equations for d_4 (phi = f(e2-e3)) and d_5 (phi = f(e2-e4)).
2. f(0) from the first integral and the boundary condition: exp f0 = 9/8 and 4 - 2 sqrt 2.
3. Closed forms: d_4  u' = (u-1) sqrt(2u(u+2)); d_5 tanh(f/4) = exp(2(x-x0)); check they solve the full
   n-component field equation and the boundary condition.
4. Energies by quadrature E = int_{-inf}^0 (phi'^2/2 + V) dx + B = int_0^{f0} 2 sqrt(W) df + B, compared with
   (16 sqrt2 - 8 - 6 sqrt6) and (32(1+sqrt2) sin(pi/8) - 24 - 4 sqrt2).
5. Coercivity constant: min over unit u of max_j alpha_j.u > 0 (d_4, d_5, d_6), and V >= c |phi|^2 sampled.
"""
import mpmath as mp
import numpy as np
mp.mp.dps = 30


def dn_roots(n):
    e = np.eye(n)
    roots = [-e[0]-e[1]] + [e[i]-e[i+1] for i in range(n-1)] + [e[n-2]+e[n-1]]
    marks = [1, 1] + [2]*(n-3) + [1, 1]
    return np.array(roots), np.array(marks)


for n, direction in ((4, (0, 1, -1, 0)), (5, (0, 1, 0, -1, 0))):
    R, nm = dn_roots(n)
    d = np.array(direction, float)
    print(f'd_{n}: alpha_j . v =', (R@d).astype(int), ' marks', nm, ' |v|^2 =', d@d)
    # invariance: every diagram symmetry permutes nodes with equal marks; the projections are equal on
    # each orbit iff v is invariant.  Here: end nodes all -1, middle nodes 2 (d4) / 1,1 (d5).

# ---- d_4: W(f) = V on the invariant line; kinetic phi'^2/2 = f'^2 (|v|^2 = 2)
W4 = lambda f: 4*(mp.e**(-f) - 1) + 2*(mp.e**(2*f) - 1)
B4 = lambda f: 2*(4*1*(mp.e**(-f/2) - 1) + 1*mp.sqrt(2)*(mp.e**f - 1))      # end nodes sqrt(1), centre sqrt(2)
W5 = lambda f: 4*(mp.e**(-f) - 1) + 2*2*(mp.e**f - 1)
B5 = lambda f: 2*(4*(mp.e**(-f/2) - 1) + 2*mp.sqrt(2)*(mp.e**(f/2) - 1))
for name, W, B, guess, closed_f0, closed_E in (
        ('d4', W4, B4, 0.1, mp.log(mp.mpf(9)/8), 16*mp.sqrt(2) - 8 - 6*mp.sqrt(6)),
        ('d5', W5, B5, 0.15, mp.log(4 - 2*mp.sqrt(2)), 32*(1 + mp.sqrt(2))*mp.sin(mp.pi/8) - 24 - 4*mp.sqrt(2))):
    # E = int (f'^2 + W) + B; boundary condition 2 f'(0) = -B'(f0); first integral f'^2 = W, f' > 0
    g = lambda f: 2*mp.sqrt(W(f)) + mp.diff(B, f)
    f0 = mp.findroot(g, guess)
    E = mp.quad(lambda f: 2*mp.sqrt(W(f)), [0, f0]) + B(f0)
    print(f'{name}: f0 = {mp.nstr(f0, 20)}, exp f0 = {mp.nstr(mp.e**f0, 20)} (closed {mp.nstr(mp.e**closed_f0, 20)});'
          f' E = {mp.nstr(E, 18)}  closed form {mp.nstr(closed_E, 18)}')

# ---- closed-form profiles against the full field equation and boundary condition
def full_check(n, direction, prof, dprof, d2prof):
    R, nm = dn_roots(n)
    v = np.array(direction, float)
    worst = 0
    for x in (-3.0, -1.2, -0.3, 0.0):
        f, fp, fpp = prof(x), dprof(x), d2prof(x)
        phi = f*v
        gradV = (R.T*(nm*np.exp(R@phi)))@np.ones(len(nm))
        worst = max(worst, np.abs(fpp*v - gradV).max())
    phi0 = prof(0.0)*v
    gradB = (R.T*(2*np.sqrt(nm)*np.exp(R@phi0/2)/2))@np.ones(len(nm))
    worst = max(worst, np.abs(dprof(0.0)*v + gradB).max())
    return worst


# d_4: u = e^f, u' = (u-1) sqrt(2u(u+2)); integrate backwards from u(0) = 9/8 numerically (mpmath odefun)
sol = mp.odefun(lambda x, u: -(u - 1)*mp.sqrt(2*u*(u + 2)), 0, mp.mpf(9)/8)   # in s = -x
prof4 = lambda x: float(mp.log(sol(-x)))
dprof4 = lambda x: float((sol(-x) - 1)*mp.sqrt(2*sol(-x)*(sol(-x) + 2))/sol(-x))
def d2prof4(x):
    f = prof4(x)
    return 2*(np.exp(2*f) - np.exp(-f))
print('d4 closed form: full field eq + boundary residual', full_check(4, (0, 1, -1, 0), prof4, dprof4, d2prof4))
# d_5: tanh(f/4) = exp(2(x - x0)), exp f0 = 4 - 2 sqrt 2
f0 = np.log(4 - 2*np.sqrt(2)); c = np.tanh(f0/4)
prof5 = lambda x: 4*np.arctanh(c*np.exp(2*x))
dprof5 = lambda x: 4*np.sinh(prof5(x)/2)
d2prof5 = lambda x: 4*np.sinh(prof5(x))
print('d5 closed form: full field eq + boundary residual', full_check(5, (0, 1, 0, -1, 0), prof5, dprof5, d2prof5))

# ---- coercivity constant
rng = np.random.default_rng(0)
for n in (4, 5, 6):
    R, nm = dn_roots(n)
    u = rng.normal(size=(200000, n)); u /= np.linalg.norm(u, axis=1)[:, None]
    m = (u@R.T).max(axis=1).min()
    # V/|phi|^2 on random points of radius 0.01..20
    rad = np.exp(rng.uniform(np.log(0.01), np.log(20), 200000))
    phi = u*rad[:, None]
    V = (nm*(np.exp(phi@R.T) - 1)).sum(axis=1)
    print(f'd_{n}: min_u max_j alpha_j.u ~ {m:.4f} (> 0); min V/|phi|^2 sampled = {np.min(V/rad**2):.4f}')
