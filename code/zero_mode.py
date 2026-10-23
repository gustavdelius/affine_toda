"""Zero-mode model of sec-solvable-check (chapter 45): H = |p|^2/2 - sign*g sum_j n_j exp(i alpha_j.x) on the field-space
torus, in the plane-wave basis p = c_1 alpha_1 + c_2 alpha_2 with |c_i| <= P.  The potential shifts p by the simple
affine roots, so H is real and not symmetric.

  a_2: |alpha|^2 = 4, shifts alpha_1, alpha_2, alpha_0 = -alpha_1 - alpha_2, marks (1,1,1).
  c_2: alpha_1 = e1-e2, alpha_2 = 2 e2, alpha_0 = -2 e1 (e_i^2 = 1; long roots have squared length 4),
       marks (alpha_1, alpha_2, alpha_0) = (2,1,1).

Prints the table of the chapter: lowest level, number of complex levels among the 40 with smallest real part,
Z(T) = Tr exp(-T H) for T = 0.5, 2, 8; then the opposite sign for a_2 at g = 1, and the convergence in P.
Exits non-zero if a claim of the chapter fails (lowest level real, simple and of minimal real part with a gap of at
least 1.3 to the next real part, Z > 0, P-independence,
complex levels present, a_2 opposite sign negative Z(2) with complex lowest pair, c_2 sign equivalence).
Runtime: 6 min measured on a heavily loaded machine (the P = 20 diagonalizations dominate)."""
import itertools
import sys

import numpy as np

MODELS = {
    "a2": (2 * np.array([[2.0, -1], [-1, 2]]), [(1, 0), (0, 1), (-1, -1)], (1, 1, 1)),
    "c2": (np.array([[2.0, -2], [-2, 4]]), [(1, 0), (0, 1), (-2, -1)], (2, 1, 1)),
}


def spectrum(model, g, P, sign=1):
    G, shifts, marks = MODELS[model]
    keys = list(itertools.product(range(-P, P + 1), repeat=2))
    idx = {c: i for i, c in enumerate(keys)}
    H = np.zeros((len(keys), len(keys)))
    for c, i in idx.items():
        v = np.array(c)
        H[i, i] = 0.5 * v @ G @ v
        for s, n in zip(shifts, marks):
            c2 = (c[0] + s[0], c[1] + s[1])
            if c2 in idx:
                H[idx[c2], i] -= sign * g * n
    ev = np.linalg.eigvals(H)
    return ev[np.argsort(ev.real, kind="stable")]


ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


Ts = (0.5, 2, 8)
print("| algebra | g | lowest level | complex among lowest 40 | Z(T), T=0.5,2,8 |")
for model in ("a2", "c2"):
    for g in (0.5, 1.0):
        ev = spectrum(model, g, 16)
        nc = int((np.abs(ev[:40].imag) > 1e-6).sum())
        Z = [np.exp(-T * ev).sum().real for T in Ts]
        print(f"| {model} | {g} | {ev[0].real:.4f} | {nc} | " + ", ".join(f"{z:.3g}" for z in Z) + " |")
        check(f"{model} g={g}: lowest level real, minimal real part", abs(ev[0].imag) < 1e-9 and np.all(ev.real >= ev[0].real - 1e-12))
        check(f"{model} g={g}: vacuum simple, every other level has real part larger by >= 1.3 "
              f"(gap {ev[1].real - ev[0].real:.3f})", ev[1].real - ev[0].real >= 1.3)
        check(f"{model} g={g}: Z(T) > 0", all(z > 0 for z in Z))
        check(f"{model} g={g}: complex levels among the lowest 40", nc > 0)
        e8, e20 = spectrum(model, g, 8)[:10], spectrum(model, g, 20)[:10]
        check(f"{model} g={g}: lowest 10 levels unchanged for P = 8 and 20", np.abs(e8 - e20).max() < 1e-6)

evm = spectrum("a2", 1.0, 16, sign=-1)
Zm = np.exp(-2 * evm).sum().real
print(f"a2, g=1, opposite sign: Z(2) = {Zm:.3f}, lowest levels {evm[0]:.3f}, {evm[1]:.3f}")
check("a2 opposite sign: Z(2) < 0 and lowest level a complex pair", Zm < 0 and abs(evm[0].imag) > 1e-3 and abs(evm[0] - np.conj(evm[1])) < 1e-6)
e_p, e_m = spectrum("c2", 1.0, 12), spectrum("c2", 1.0, 12, sign=-1)
check("c2: g and -g give the same spectrum (sum n_j even)", np.abs(e_p[:20] - e_m[:20]).max() < 1e-8)
sys.exit(0 if ok else 1)
