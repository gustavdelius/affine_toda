"""Integration check, ch. 43 item 6: residue ranks of T_V (species-1 pole) and T_Vbar (species-r pole) on a charge sector;
the Lee-Yang poles are in Tbar_b; the imaginary-theta_d particle factor is a phase only for v in pi Z.
Own code from eq-dt-tv, eq-dt-tvbar (sector blocks built directly), eq-ly-transmission, eq-dpt-particle; the scalar
factor gbar from defect_common (book helper) only to confirm that its pole at the species-r position is not cancelled."""
import sys, numpy as np
sys.path.insert(0, 'code')
from defect_common import an_scalar_gbar, an_scalar_g, an_scalar_constant

def blocks(h, om, X, Q, s):
    """Sector of total charge Q (sum Q_i = 1 for V, -1 for Vbar).  V: basis v_i (x) |Q - e_i>; Vbar: vb_i (x) |Q + e_i>."""
    q = -np.exp(-1j*np.pi*om)
    I = np.eye(h)
    MV = np.zeros((h, h), complex)
    for i in range(h):                       # input v_i (x) |N>, N = Q - e_i
        N = Q - I[i]
        MV[i, i] += q**(-(N + s)[i])         # output |N> (x) v_i, label index i
        MV[(i - 1) % h, i] += X              # output |N + e_i - e_{i-1}> (x) v_{i-1}
    Qb = Q - 2*I[0] + I[0]                   # Vbar sector: total charge Qb with sum -1 (use Q - e_0 shape)
    Qb = Q - np.ones(h)*0 - I[0] - I[0]
    MVb = np.zeros((h, h), complex)
    for i in range(h):                       # input vb_i (x) |N>, N = Qb + e_i
        N = Qb + I[i]
        kk = np.zeros(h)
        for jj in range(h):
            b = (i + jj) % h
            kk[b] += 1
            MVb[b, i] += X**jj*q**((N + s) @ kk)   # output |N - e_i + e_b> = |Qb + e_b> (x) vb_b
    return MV, MVb, q

rng = np.random.default_rng(3)
for h in (3, 4, 5, 6):
    r = h - 1
    for om in (np.sqrt(2), 1.37, 3.1):
        s = rng.normal(size=h)*0.3; s -= s.mean()
        Q = np.zeros(h); Q[1] = 1                       # sum 1 for V
        q = -np.exp(-1j*np.pi*om)
        # species-1 pole (eq-dt-pole): X = e^{omega(-i pi r/h - i pi/(h omega))}; species r: X = e^{-i pi omega/h - i pi/h}
        X1 = np.exp(om*(-1j*np.pi*r/h) - 1j*np.pi/h)
        Xr = np.exp(-1j*np.pi*om/h - 1j*np.pi/h)
        MV, _, _ = blocks(h, om, X1, Q, s)
        _, MVb, _ = blocks(h, om, Xr, Q, s)
        sv1 = np.linalg.svd(MV, compute_uv=False); svr = np.linalg.svd(MVb, compute_uv=False)
        print(f'h={h} om={om:.3f}: (-1/q)^r X1^h = {((-1/q)**r*X1**h):.3f}, Xr^h/q = {Xr**h/q:.3f}; '
              f'M_V(theta_1) sv/sv0 = {np.round(sv1/sv1[0], 3)}; M_Vbar(theta_r) sv/sv0 = {np.round(svr/svr[0], 3)}')
# pole of gbar at the species-r position is simple and not cancelled (book helper)
for h, om in ((3, 1.37), (4, 1.37), (5, 1.37)):
    r = h - 1; q = -np.exp(-1j*np.pi*om); c = an_scalar_constant(r, om, q)
    Dr = 1j*np.pi/2 - 1j*(np.pi/h + np.pi/(h*om))
    vals = [eps*an_scalar_gbar(Dr + eps, r, om, 400, c) for eps in (1e-5, 5e-6)]
    print(f'h={h}: eps*gbar(theta_r + eps) = {vals[0]:.4f}, {vals[1]:.4f}')

# Lee-Yang: T_b = f(theta - i pi(3-b)/6), Tbar_b = f(theta + i pi(3-b)/6), f = f_{2/3}
f = lambda th: (np.sinh(th) + 1j*np.sin(2*np.pi/3))/(np.sinh(th) - 1j*np.sin(2*np.pi/3))
for b in (-0.5, 0.5, 1.5):
    for y0, name in (((b + 1)/6, '|1+>'), ((b - 1)/6, '|2+>')):
        th = 1j*np.pi*y0; eps = 1e-7
        print(f'LY b={b} {name} y={y0:.3f}: eps*T_b={abs(eps*f(th + eps - 1j*np.pi*(3 - b)/6)):.1e}, eps*Tbar_b={abs(eps*f(th + eps + 1j*np.pi*(3 - b)/6)):.2e}')
    # T_b poles in the strip 0 < y < 1/2?
    ys = np.linspace(0.001, 0.499, 4000)
    print(f'   max |T_b| on the strip axis: {max(abs(f(1j*np.pi*y - 1j*np.pi*(3 - b)/6)) for y in ys):.2f}')

# imaginary theta_d = i v: |T_a| = 1 on the real line iff v in pi Z; v = pi inside the window iff a < h/2
def Ta(th, thd, a, h):
    g = np.pi*a/h
    return (np.exp(-(th - thd)) + np.exp(1j*g))/(np.exp(-(th - thd)) + np.exp(-1j*g))
ths = np.linspace(-3, 3, 61)
for h in (3, 4, 5):
    for a in range(1, h):
        g = np.pi*a/h
        dev_pi = max(abs(abs(Ta(t, 1j*np.pi, a, h)) - 1) for t in ths)
        dev_gen = max(abs(abs(Ta(t, 1j*(np.pi - 0.3), a, h)) - 1) for t in ths)
        inwin = np.pi - g < np.pi < 1.5*np.pi - g
        # pole position at v = pi: theta = i gamma_a
        res = abs(1/Ta(1j*g + 1e-9, 1j*np.pi, a, h))
        print(f'h={h} a={a}: v=pi phase dev {dev_pi:.0e}, v=pi-0.3 dev {dev_gen:.2f}; v=pi in window: {inwin} (a<h/2: {a < h/2}); 1/T at i gamma_a: {res:.0e}')
# sine-Gordon row: window pi/(2 lam) < v < pi/2 + pi/(2 lam) contains v = pi iff 1/2 < lam < 1
for lam in (0.4, 0.6, 0.9, 1.2):
    print(f'SG lam={lam}: v=pi in window: {np.pi/(2*lam) < np.pi < np.pi/2 + np.pi/(2*lam)}')
