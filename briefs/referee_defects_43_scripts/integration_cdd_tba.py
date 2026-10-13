"""Integration check, ch. 43 item 7: the soliton CDD condition and the TBA fugacities.
A soliton CDD factor f(theta - theta_d) multiplies g; by eq-dt-bootstrap gbar gets F(D) = prod_{m<r} f(D + i pi(2m-r+1)/h).
Crossing (ch. 41): g(D) gbar(D - i pi) c_I = 1 and gbar(D) g(D - i pi) c_II = 1 then need
   (I)  f(D) F(D - i pi) = 1,    (II)  F(D) f(D - i pi) = 1.
Claim: both are  prod_{k<h} f(D - i pi + 2 pi i k/h) = 1 (up to a shift)."""
import numpy as np
def Tcl(D, a, h):
    g = np.pi*a/h
    return (np.exp(-D) + np.exp(1j*g))/(np.exp(-D) + np.exp(-1j*g))
def F(f, D, h):
    r = h - 1
    return np.prod([f(D + 1j*np.pi*(2*m - r + 1)/h) for m in range(r)], axis=0)
def cond(f, D, h):
    return np.prod([f(D - 1j*np.pi + 2j*np.pi*k/h) for k in range(h)], axis=0)
Ds = np.array([0.41, -0.7 + 0.2j, 1.3 - 0.1j])
# the two crossing conditions are the string condition, at D and at a shifted D
for h in (2, 3, 4, 5):
    for name, f in (('Tcl_1', lambda D: Tcl(D, 1, h)), ('sinh', lambda D: (np.sinh(D) + 0.3j)/(np.sinh(D) - 0.3j)),
                    ('e^{ic sinh}', lambda D: np.exp(0.19j*np.sinh(D)))):
        I_ = f(Ds)*F(f, Ds - 1j*np.pi, h); II = F(f, Ds, h)*f(Ds - 1j*np.pi)
        shift = -1j*np.pi + 2j*np.pi/h
        print(f'h={h} {name:12s}: |(II) - string(D)| = {abs(II - cond(f, Ds, h)).max():.1e}, |(I) - string(D + i pi - 2 pi i/h)| = '
              f'{abs(I_ - cond(f, Ds - shift, h)).max():.1e}; string - 1 = {abs(cond(f, Ds, h) - 1).max():.2e}; '
              f'f(D) f(D + i pi) - 1 = {abs(f(Ds)*f(Ds + 1j*np.pi) - 1).max():.2e}')
# breather multiplier: product of f over the h constituent rapidities of B^(a)_1 (sec-breather-transmission)
def constituents(h, a, om, D):
    ub = np.pi/2 - np.pi/(h*om)
    return [D - 1j*ub + 1j*np.pi*(2*m - h + a + 1)/h for m in range(h - a)] + [D + 1j*ub + 1j*np.pi*(2*m - a + 1)/h for m in range(a)]
for h, a in ((4, 2), (3, 1), (4, 1)):
    for name, f in (('Tcl_1', lambda D: Tcl(D, 1, h)), ('e^{ic sinh}', lambda D: np.exp(0.19j*np.sinh(D)))):
        vals = [np.prod([f(c) for c in constituents(h, a, om, 0.41)]) for om in (np.sqrt(2), np.sqrt(5))]
        print(f'h={h} a={a} {name:12s}: breather multiplier at omega=sqrt2: {vals[0]:.5f}, sqrt5: {vals[1]:.5f}')
    # closed form for e^{ic sinh}: exp(2 i c sin(pi a/h)/sin(pi/h) sin(pi/(h om)) sinh D)
    for om in (np.sqrt(2), np.sqrt(5)):
        cf = np.exp(2j*0.19*np.sin(np.pi*a/h)/np.sin(np.pi/h)*np.sin(np.pi/(h*om))*np.sinh(0.41))
        num = np.prod([np.exp(0.19j*np.sinh(c)) for c in constituents(h, a, om, 0.41)])
        print(f'    e^(ic sinh) closed form check at omega={om:.3f}: {abs(cf - num):.1e}')
# T^cl_1 satisfies the string condition for every h (hand: prod (u - x e^{-2 pi i k/h})/(u^-1 - x e^{-2 pi i k/h}) = (u^h - x^h)/(u^-h - x^h) = 1)
print('T^cl_1 string condition h = 2..12:', max(abs(cond(lambda D: Tcl(D, 1, h), Ds, h) - 1).max() for h in range(2, 13)))
print('T^cl_1 f(D) f(D+i pi) - 1, h = 3..6:', [round(abs(Tcl(0.41, 1, h)*Tcl(0.41 + 1j*np.pi, 1, h) - 1), 3) for h in range(3, 7)])

# TBA fugacities F_a(theta) = T_abar(i pi/2 + theta) for real theta, theta_d
thd = 0.27
for h in (3, 4, 5, 6):
    th = np.linspace(-3, 3, 13)
    dev = 0; mx = 0
    for a in range(1, h):
        Fa = Tcl(1j*np.pi/2 + th - thd, h - a, h); Fab = Tcl(1j*np.pi/2 + th - thd, a, h)
        dev = max(dev, abs(Fab - np.conj(Fa)).max()); mx = max(mx, abs(Fa).max())
        x = np.exp(thd - th); g = np.pi*a/h
        dev = max(dev, abs(abs(Fa)**2 - (1 + x**2 - 2*x*np.sin(g))/(1 + x**2 + 2*x*np.sin(g))).max())
    line = f'h={h}: F_abar = conj F_a and |F_a|^2 formula [{dev:.1e}], max |F_a| = {mx:.3f}'
    if h % 2 == 0:
        Fh = Tcl(1j*np.pi/2 + th - thd, h//2, h)
        line += f'; F_h/2 + tanh((theta-theta_d)/2) = {abs(Fh + np.tanh((th - thd)/2)).max():.1e}, imag {abs(Fh.imag).max():.1e}, sign change: {Fh.min() < 0 < Fh.max()}'
    print(line)
