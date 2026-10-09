"""Topological charges of continuum Hirota solutions of a_n^(1) (imaginary coupling)."""
import numpy as np, itertools, collections
from lie import An
from hirota import tau, field

def charge_label(g, species, thetas, xis, t=0.0, X=None, N=20001):
    m = 2*np.sin(np.pi*np.array(species)/g.h)
    if X is None:
        X = 25.0
    x = np.linspace(-X, X, N)
    T = tau(g, species, thetas, xis, x, t=t)
    if np.min(np.abs(T)) < 1e-6 * np.max(np.abs(T), axis=1).min()**0 and np.min(np.abs(T))<1e-6:
        return None
    # check sampling fine enough for unwrap: phase jumps
    ph = np.angle(T)
    d = np.diff(np.unwrap(ph, axis=1), axis=1)
    if np.max(np.abs(d)) > 1.0:
        return None
    u = field(g, T)
    q = (u[:, -1] - u[:, 0]) / (2*np.pi)
    # remove the non-quantized part: u(+inf) should be 2 pi * weight (real)
    lab = g.label(q.real, tol=1e-3)
    if lab is None or np.max(np.abs(q.imag)) > 1e-3:
        return ('bad', tuple(np.round(q, 3)))
    return lab

def rep_name(g, lab):
    if lab is None or lab[0]=='bad': return str(lab)
    return ''.join(str(i+1) for i, c in enumerate(lab) for _ in range(c)) or '0'

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    for n in (2, 3, 4, 5):
        g = An(n)
        print("=== a_%d^(1), h=%d" % (n, g.h))
        for a in range(1, n+1):
            labs = collections.Counter()
            for th in np.linspace(0, 2*np.pi, 721)[:-1] + 1e-3:
                lab = charge_label(g, [a], [0.0], [1j*th])
                labs[rep_name(g, lab)] += 1
            print(" species %d: charges (eps-subset labels, count of Im xi samples):" % a, dict(labs),
                  " #distinct =", len([k for k in labs if k not in ('None',)]), " dim Lambda^a =", len(list(itertools.combinations(range(n+1), a))))
