"""Integration check, ch. 43 item 5: the explicit second module rho~ of prp-dpt-doubled.

Own code (exact propagation on dicts of labels, no window), from eq-dt-weyl, eq-dt-tv, eq-dt-tvbar and the e-type equation
   T [pi(e_j) (x) 1 + pi(k_j) (x) Y] = [Y (x) 1 + kappa (x) pi(e_{psi(j)})] T,   psi(e_j) = e_{j-1},
with the candidate (0-based indices: alpha_j = e_{j-1} - e_j, first end a = j-1, second end b = j)
   rho~(e_j)|N> = (t/(q - 1/q)) q^{-(N+s)_a} |N + alpha_j>  (+ c_j q^{+(N+s)_a} |N + alpha_j>),  t = e^{omega(theta_d + i pi/2)},
   rho~(k_j)|N> = |N + alpha_j - alpha_{j-1}>.
Checks: V and Vbar, every node, h = 3, 4, 5, several theta; the kappa equation T(pi(k_j) (x) kappa) = (kappa (x) pi(k_{j-1})) T;
c_j != 0 still solves the V equations but breaks k_j rho~(e_j) k_j^-1 = q^2 rho~(e_j); the U_q(b^+) relations of rho~;
exact y-power obstructions for the other twists on V; h = 2 against eq-sg-transmission-matrix and eq-dbd-second-module."""
import numpy as np
om, TD = np.sqrt(2), 0.27
q = -np.exp(-1j*np.pi*om)
t = np.exp(om*(TD + 1j*np.pi/2))

def unit(h, i):
    v = [0]*h; v[i % h] = 1; return tuple(v)
def add(N, *vs, sign=1):
    N = list(N)
    for v in vs:
        for k in range(len(N)): N[k] += sign*v[k]
    return tuple(N)
def alpha(h, j):
    return tuple(np.array(unit(h, j - 1)) - np.array(unit(h, j)))
def acc(d, k, v):
    d[k] = d.get(k, 0) + v

# V(y): e_j v_j = y v_{j-1}; Vbar(y): e_j vb_{j-1} = y vb_j; k_j = q^{+-H_j}, H_j = E_{j-1,j-1} - E_{jj}
def rep_e(h, conj, y, j, i):
    a, b = (j - 1) % h, j % h
    if not conj:
        return [(a, y)] if i == b else []
    return [(b, y)] if i == a else []
def rep_k(h, conj, j, i):
    a, b = (j - 1) % h, j % h
    H = (i == a) - (i == b)
    return q**(-H if conj else H)

def T_apply(h, s, X, conj, state):
    """state: dict (i, N) on V (x) W -> dict (N, i) on W (x) V, by eq-dt-tv or eq-dt-tvbar."""
    out = {}
    for (i, N), c in state.items():
        Ns = np.array(N) + s
        if not conj:
            acc(out, (N, i), c*q**(-Ns[i]))
            acc(out, (add(N, unit(h, i), unit(h, i - 1), sign=1) if False else tuple(np.array(N) + np.array(unit(h, i)) - np.array(unit(h, i - 1))), (i - 1) % h), c*X)
        else:
            kk = np.zeros(h)
            for jj in range(h):
                b = (i + jj) % h
                kk[b] += 1
                acc(out, (tuple(np.array(N) - np.array(unit(h, i)) + np.array(unit(h, b))), b), c*X**jj*q**(Ns @ kk))
    return out

def Yop(h, s, j, c=0.0):
    a = (j - 1) % h
    def f(N):
        n = N[a] + s[a]
        return [(tuple(np.array(N) + np.array(alpha(h, j))), t/(q - 1/q)*q**(-n) + c*q**n)]
    return f
def Kop(h, j):
    w = np.array(alpha(h, j)) - np.array(alpha(h, j - 1))
    return lambda N: [(tuple(np.array(N) + w), 1.0)]

def lhs_e(h, s, y, conj, j, Y, inp):
    """[pi(e_j) (x) 1 + pi(k_j) (x) Y] on V (x) W, input dict (i, N)."""
    out = {}
    for (i, N), c in inp.items():
        for (i2, v) in rep_e(h, conj, y, j, i):
            acc(out, (i2, N), c*v)
        for (N2, v) in Y(N):
            acc(out, (i, N2), c*rep_k(h, conj, j, i)*v)
    return out
def rhs_e(h, s, y, conj, jpsi, Y, K, inp):
    """[Y (x) 1 + kappa (x) pi(e_psi(j))] on W (x) V, input dict (N, i)."""
    out = {}
    for (N, i), c in inp.items():
        for (N2, v) in Y(N):
            acc(out, (N2, i), c*v)
        for (i2, v) in rep_e(h, conj, y, jpsi, i):
            for (N2, w) in K(N):
                acc(out, (N2, i2), c*v*w)
    return out
def diff(d1, d2):
    keys = set(d1) | set(d2)
    num = max([abs(d1.get(k, 0) - d2.get(k, 0)) for k in keys] + [0])
    den = max([abs(v) for v in list(d1.values()) + list(d2.values())] + [1e-300])
    return num/den

def e_type_residual(h, s, th, conj, j, jpsi, c=0.0, Nlabels=None):
    y = np.exp(om*th); X = np.exp(om*(th - TD - 1j*np.pi/2))
    Y, K = Yop(h, s, j, c), Kop(h, j)
    worst = 0
    for N in Nlabels:
        for i in range(h):
            inp = {(i, N): 1.0}
            L = T_apply(h, s, X, conj, lhs_e(h, s, y, conj, j, Y, inp))
            R = rhs_e(h, s, y, conj, jpsi, Y, K, T_apply(h, s, X, conj, inp))
            worst = max(worst, diff(L, R))
    return worst

def kappa_residual(h, s, th, conj, j, Nlabels):
    """T (pi(k_j) (x) kappa) = (kappa (x) pi(k_{j-1})) T."""
    X = np.exp(om*(th - TD - 1j*np.pi/2)); K = Kop(h, j)
    worst = 0
    for N in Nlabels:
        for i in range(h):
            L = T_apply(h, s, X, conj, {(i, K(N)[0][0]): rep_k(h, conj, j, i)})
            R0 = T_apply(h, s, X, conj, {(i, N): 1.0})
            R = {}
            for (N2, i2), v in R0.items():
                acc(R, (K(N2)[0][0], i2), v*rep_k(h, conj, j - 1, i2))
            worst = max(worst, diff(L, R))
    return worst

rng = np.random.default_rng(1)
ths = (-0.32, 0.41, 0.82 + 0.11j)
for h in (3, 4, 5):
    s = rng.normal(size=h)*0.3; s -= s.mean()
    Ns = [tuple([0]*h), tuple(np.array(unit(h, 0)) - np.array(unit(h, 2))), tuple(2*np.array(unit(h, 1)) - 2*np.array(unit(h, h - 1)))]
    rV = max(e_type_residual(h, s, th, False, j, (j - 1) % h, 0.0, Ns) for th in ths for j in range(h))
    rVb = max(e_type_residual(h, s, th, True, j, (j - 1) % h, 0.0, Ns) for th in ths for j in range(h))
    rK = max(kappa_residual(h, s, th, cj, j, Ns) for th in ths for j in range(h) for cj in (False, True))
    rVc = max(e_type_residual(h, s, th, False, j, (j - 1) % h, 0.3, Ns) for th in ths for j in range(h))
    rVbc = max(e_type_residual(h, s, th, True, j, (j - 1) % h, 0.3, Ns) for th in ths for j in range(h))
    # other twists fail on V (y-power obstruction): e-type psi = id, pi^-1
    rfail = min(e_type_residual(h, s, 0.41, False, j, jp, 0.0, Ns) for j in range(h) for jp in (j, (j + 1) % h))
    print(f'h={h}: e-type psi=pi, c=0: V {rV:.1e}, Vbar {rVb:.1e}; kappa eq {rK:.1e}; with c_j=0.3: V {rVc:.1e}, Vbar {rVbc:.1e}; '
          f'with rho~ the e-type psi=id, pi^-1 fail by >= {rfail:.2f}')

# relations of rho~ (as operators on labels): Cartan conjugation, quadratic relations, product of kappas
def compose(*ops):
    def f(N):
        cur = [(N, 1.0)]
        for op in reversed(ops):
            nxt = []
            for (M, c) in cur:
                for (M2, v) in op(M):
                    nxt.append((M2, c*v))
            cur = nxt
        return cur
    return f
def opdiff(A, B, Ns):
    w = 0
    for N in Ns:
        a, b = {}, {}
        for k, v in A(N): acc(a, k, v)
        for k, v in B(N): acc(b, k, v)
        w = max(w, diff(a, b))
    return w
def scal(c, op):
    return lambda N: [(M, c*v) for (M, v) in op(N)]
def Kinv(h, j):
    w = np.array(alpha(h, j)) - np.array(alpha(h, j - 1))
    return lambda N: [(tuple(np.array(N) - w), 1.0)]
for h in (3, 4, 5):
    s = rng.normal(size=h)*0.3; s -= s.mean()
    Ns = [tuple([0]*h), tuple(np.array(unit(h, 0)) - np.array(unit(h, 2)))]
    cart = 0; cartc = 0; quad = 0; serre = 0
    for i in range(h):
        for j in range(h):
            aij = 2 if i == j else (-1 if (i - j) % h in (1, h - 1) else 0)
            cart = max(cart, opdiff(compose(Kop(h, i), Yop(h, s, j), Kinv(h, i)), scal(q**aij, Yop(h, s, j)), Ns))
            cartc = max(cartc, opdiff(compose(Kop(h, j), Yop(h, s, j, 0.3), Kinv(h, j)), scal(q**2, Yop(h, s, j, 0.3)), Ns)) if i == j else cartc
            if (j - i) % h == 1:   # adjacent: e_i e_j = q^{-1} e_j e_i  with j = i+1
                quad = max(quad, opdiff(compose(Yop(h, s, i), Yop(h, s, j)), scal(1/q, compose(Yop(h, s, j), Yop(h, s, i))), Ns))
                Yi, Yj = Yop(h, s, i), Yop(h, s, j)
                for (A, B) in ((Yi, Yj), (Yj, Yi)):
                    t1 = compose(A, A, B); t2 = scal(-(q + 1/q), compose(A, B, A)); t3 = compose(B, A, A)
                    for N in Ns:
                        tot = {}
                        for k, v in t1(N) + t2(N) + t3(N): acc(tot, k, v)
                        scale = max(abs(v) for k, v in t1(N))
                        serre = max(serre, max(abs(v) for v in tot.values())/scale)
            elif i != j and (i - j) % h != 1:
                quad = max(quad, opdiff(compose(Yop(h, s, i), Yop(h, s, j)), compose(Yop(h, s, j), Yop(h, s, i)), Ns))
    prodK = opdiff(compose(*[Kop(h, j) for j in range(h)]), lambda N: [(N, 1.0)], Ns)
    print(f'h={h}: Cartan conj {cart:.1e}; c_j=0.3 breaks it by {cartc:.2f}; quadratic relations {quad:.1e}; q-Serre {serre:.1e}; prod kappa = 1 {prodK:.0e}')

# h = 2: basis d_n = tau^-n q^{n(n-1)/2 + n s/2}, tau = -i e^{i pi om/2}, N = (n, -n), s = (s/2, -s/2), T times q^{s/2}
sS = 0.37; s2 = np.array([sS/2, -sS/2]); tau = -1j*np.exp(1j*np.pi*om/2); yd = np.exp(om*TD)
d = lambda n: tau**(-n)*q**(n*(n - 1)/2 + n*sS/2)
devT = devR = 0
for th in ths:
    y = np.exp(om*th); X = np.exp(om*(th - TD - 1j*np.pi/2))
    for n in range(-3, 4):
        N = (n, -n)
        o1 = T_apply(2, s2, X, False, {(0, N): 1.0}); o2 = T_apply(2, s2, X, False, {(1, N): 1.0})
        A = q**(sS/2)*o1[(N, 0)]; C = q**(sS/2)*o1[((n + 1, -n - 1), 1)]*d(n)/d(n + 1)
        D = q**(sS/2)*o2[(N, 1)]; B = q**(sS/2)*o2[((n - 1, -n + 1), 0)]*d(n)/d(n - 1)
        ch37 = (q**-n, -1j*(y/yd)*q**(n + sS), -1j*(y/yd)*q**-n, q**(n + sS))     # A, B, C, D of eq-sg-transmission-matrix
        devT = max(devT, max(abs(u - v)/abs(v) for u, v in zip((A, B, C, D), ch37)))
for n in range(-3, 4):
    N = (n, -n)
    e1 = Yop(2, s2, 1)(N)[0]; k1 = Kop(2, 1)(N)[0]; e0 = Yop(2, s2, 0)(N)[0]; k0 = Kop(2, 0)(N)[0]
    a = 1j*yd*q**(-sS)/(q**2 - 1); b = 1j*yd*q**sS/(q**2 - 1)
    got = (e1[1]*d(n)/d(e1[0][0]), k1[1]*d(n)/d(k1[0][0]), e0[1]*d(n)/d(e0[0][0]), k0[1]*d(n)/d(k0[0][0]))
    want = (a*q**(-2*n), q**(-2*n - 2 - sS), b*q**(2*n), q**(2*n - 2 + sS))
    shifts = (e1[0][0] - n, k1[0][0] - n, e0[0][0] - n, k0[0][0] - n)
    assert shifts == (1, 2, -1, -2), shifts
    devR = max(devR, max(abs(u - v)/abs(v) for u, v in zip(got, want)))
print(f'h=2: q^(s/2) M_V in the basis d_n|n> = eq-sg-transmission-matrix [{devT:.1e}]; rho~ = eq-dbd-second-module, shifts +1, +2, -1, -2 [{devR:.1e}]')
