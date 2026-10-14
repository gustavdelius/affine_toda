"""tau(lam) conservation along REGULAR three-soliton image solutions (phases from the scan in e2b_singular.log)."""
import numpy as np
src = open('e2_solitons.py').read()
ns = {'__name__': 'e2c'}
exec(compile(src[:src.index('# zero-curvature sanity check')], 'e2_solitons.py[defs]', 'exec'), ns)
Sol, interaction, transport, compat, class_of, fields = ns['Sol'], ns['interaction'], ns['transport'], ns['compat'], ns['class_of'], ns['fields']


def three(n, a, b, th, xi1, y3):
    h = n + 1
    Aab = interaction(a, b, th, 0.0, h).real; Aaa = interaction(a, a, th, -th, h).real
    return Sol(n, [a, a, b], [th, -th, 0.0], [xi1, -np.conj(xi1) - np.log(Aaa) - np.log(Aab), -np.log(Aab) + 1j * y3])


xs_full = np.linspace(-25, 25, 10001); i0 = 5000
cases = [(1, 1, 1, 0.5, 1.0 + 0.7j, 1.1), (2, 1, 1, 0.5, 1.0 + 2.094j, 2.094), (2, 2, 1, 0.6, 1.0 + 0.0j, 4.189),
         (3, 1, 2, 0.5, 1.0 + 5.76j, 1.571), (3, 2, 1, 0.45, 1.0 + 1.571j, 5.76), (3, 1, 1, 0.5, 1.0 + 2.618j, 2.618),
         (3, 2, 2, 0.5, 1.0 + 1.571j, 1.571)]
for (n, a, b, th, xi1, y3) in cases:
    s = three(n, a, b, th, xi1, y3)
    img = 0; vs = []; dIm = 0
    for t in np.linspace(-8, 8, 17):
        ph, dph, T = s.phi(xs_full, t)
        c = ph + np.conj(ph[:, ::-1]); img = max(img, np.max(np.abs(c - c[:, [i0]])))
        vs.append(ph[:, i0].real); dIm = max(dIm, np.max(np.abs(dph[:, i0].imag)))
    vs = np.array(vs); nu, off, triv = class_of(vs[0], n)
    D = np.diag(np.exp(1j * np.pi * nu)); E = np.eye(n + 1)
    cD = max(compat(s, t, lam, D) for t in (-4.0, -1.0, 0.5, 3.0) for lam in (0.6, 1.5))
    c1 = max(compat(s, t, lam, E) for t in (-4.0, -1.0, 0.5, 3.0) for lam in (0.6, 1.5))
    line = (f'a_{n} (a={a}, b={b}, th={th}): image res {img:.1e}, Re phi(0,t) spread {np.ptp(vs, axis=0).max():.1e}, |d_x Im phi(0)| {dIm:.1e}, '
            f'nu={np.round(nu, 3)} in 2Lambda: {triv}; compat K=D {cD:.1e}, K=1 {c1:.1e}')
    print(line)
    for lam in (0.6, 1.5):
        tD, t1 = [], []
        for t in (-6.0, -3.0, 0.0, 3.0, 6.0):
            T1 = transport(s, t, lam, 16.0, 1600); T2 = transport(s, t, 1 / lam, 16.0, 1600)
            tD.append(np.trace(D @ T1 @ T2.conj().T)); t1.append(np.trace(T1 @ T2.conj().T))
        rel = lambda z: np.max(np.abs(np.array(z) - z[0])) / abs(z[0])
        print(f'    lam={lam}: max_t |tau(t) - tau(-6)|/|tau| with K_+=D: {rel(tD):.1e};  with K_+=1: {rel(t1):.1e}')
