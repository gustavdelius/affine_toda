"""Referee check for item 5 (ch43 prp-dpt-doubled), written by the second referee.

No book code imported.  Three checks:
 (A) For h=2 the a_r matrix (ch41 eq-dt-tv) times q^{s/2}, conjugated by the
     diagonal gauge g_n = tau^n q^{-n(n-1)/2 - n s/2}, tau = -i e^{i pi w/2},
     is the ch37 matrix (eq-sg-transmission-matrix).
 (B) The same gauge maps the h=2 case of the explicit a_r second module
        rho~(e_j)|N> = C q^{-(N+s)_j} |N+alpha_j>,  rho~(k_j)|N> = |N+alpha_j-alpha_{j-1}>,
        C = e^{w(theta_d+i pi/2)}/(q-q^{-1}),
     to ch42 eq-dbd-second-module.
 (C) For h=3,4,5 the free constant c_j q^{+(N+s)_j} of the e-type solution violates
     rho~(k_j) Y_j rho~(k_j)^{-1} = q^2 Y_j unless c_j = 0.
"""
import numpy as np

w, s, td = np.sqrt(3) - 0.2, 0.37, 0.29
q = -np.exp(-1j*np.pi*w)
tau = -1j*np.exp(1j*np.pi*w/2)
assert abs(tau**2 - 1/q) < 1e-15
ns = range(-6, 7)
g = {n: tau**n * q**(-n*(n-1)/2 - n*s/2) for n in range(-9, 10)}

# (A) entries as (coefficient, label shift); A: v+ -> v+, C: v+ -> v-, B: v- -> v+, D: v- -> v-
worst_A = 0
for theta in (-0.4, 0.3 + 0.2j):
    Y = np.exp(w*(theta - td)); X = np.exp(w*(theta - td - 1j*np.pi/2))
    for n in ns:
        a41 = {'A': (q**(-n - s/2), 0), 'C': (X, 1), 'B': (X, -1), 'D': (q**(n + s/2), 0)}
        a37 = {'A': q**(-n), 'C': -1j*Y*q**(-n), 'B': -1j*Y*q**(n + s), 'D': q**(n + s)}
        for key, (c, k) in a41.items():
            new = q**(s/2) * c * g[n + k] / g[n]
            worst_A = max(worst_A, abs(new - a37[key]) / abs(a37[key]))
print(f'(A) gauge maps h=2 eq-dt-tv x q^(s/2) to eq-sg-transmission-matrix: max rel. error {worst_A:.2e}')

# (B) the second module, h=2: N=(n,-n), s=(s/2,-s/2), alpha_1 = e1-e2 (n -> n+1), alpha_0 = -alpha_1
yd = np.exp(w*td)
C = np.exp(w*(td + 1j*np.pi/2)) / (q - 1/q)
a = 1j*yd*q**(-s)/(q**2 - 1); b = 1j*yd*q**s/(q**2 - 1)
worst_B = 0
for n in ns:
    ar = {'e1': (C*q**(-(n + s/2)), 1), 'k1': (1, 2), 'e0': (C*q**(n + s/2), -1), 'k0': (1, -2)}
    ch42 = {'e1': a*q**(-2*n), 'k1': q**(-2*n - 2 - s), 'e0': b*q**(2*n), 'k0': q**(2*n - 2 + s)}
    for key, (c, k) in ar.items():
        new = c * g[n + k] / g[n]
        worst_B = max(worst_B, abs(new - ch42[key]) / abs(ch42[key]))
print(f'(B) gauge maps h=2 second module to eq-dbd-second-module: max rel. error {worst_B:.2e}')

# (C) Cartan relation selects c_j = 0 for h >= 3
for h in (3, 4, 5):
    offs = np.linspace(-0.21, 0.21, h)
    alpha = [np.roll(np.eye(h, dtype=int)[0] - np.eye(h, dtype=int)[1], j - 1) for j in range(h)]  # alpha[j] = e_j - e_{j+1}, 1-based j, alpha[0] = e_h - e_1
    shift_k = [alpha[j] - alpha[(j - 1) % h] for j in range(h)]
    N = np.array([1, -2] + [1] + [0]*(h - 3))
    for c_free in (0.0, 0.3):
        worst = 0
        for j in range(h):
            comp = (j - 1) % h  # 0-based index of the component e_j
            Yj = lambda M: C*q**(-(M[comp] + offs[comp])) + c_free*q**(M[comp] + offs[comp])
            # k_j Y_j k_j^{-1} |N> = Y_j(N - shift) |N + alpha_j>
            lhs = Yj(N - shift_k[j]); rhs = q**2 * Yj(N)
            worst = max(worst, abs(lhs - rhs) / abs(rhs))
        print(f'(C) h={h} free constant c={c_free}: max violation of k_j Y_j k_j^-1 = q^2 Y_j: {worst:.2e}')
