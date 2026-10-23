"""Regression checks for the independently re-derived PR #14 referee findings.

Run from the repository root: python3 referee_scratch/halfline_updated/checks.py. The charge enumeration uses rational
phase intervals, not a phase grid. CDD ambiguity is tested after deriving its
breather image directly in real-coupling blocks, without Gamma products.
The d5 pole is counted from Dorey's blocks and tested against the decaying
Jost subspace; the Hessian also checks the proposed spinor-channel explanation.
See briefs/halfline_consolidation_integration.md for the derivations and scope.
"""
from fractions import Fraction as Q
from itertools import product
from math import gcd
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp


def check(name, condition, detail=''):
    print(f'{"PASS" if condition else "FAIL"} {name}' + (f': {detail}' if detail else ''))
    assert condition, name


def charge_classes(n):
    h = n + 1
    reached = set()
    for b in range(1, h):
        # tau_j = 1 + exp(m_b x + 2 pi i (phase + jb/h)).
        # A charge changes only when a tau ray passes through the negative axis.
        cuts = sorted({(Q(1, 2) - Q(j*b, h)) % 1 for j in range(h)})
        charges = set()
        for left, right in zip(cuts, cuts[1:] + [cuts[0] + 1]):
            phase = (left + right)/2
            args = [(phase + Q(j*b, h) + Q(1, 2)) % 1 - Q(1, 2) for j in range(h)]
            # -sum alpha_j arg(tau_j)/(2 pi), alpha_0=e_h-e_1.
            weight = [args[j] - args[(j+1) % h] for j in range(h)]
            charges.add(tuple(weight))
            dynkin = [weight[j] - weight[j+1] for j in range(n)]
            assert all(x.denominator == 1 for x in dynkin)
            reached.add(tuple(int(x) % 2 for x in dynkin))
        assert len(charges) == h // gcd(b, h)
    all_classes = set(product((0, 1), repeat=n)) - {(0,)*n}
    return all_classes - reached


def block(x, theta, h):
    return mp.sinh((theta + 1j*mp.pi*x/h)/2) / mp.sinh((theta - 1j*mp.pi*x/h)/2)


def pair_image(q, theta, h, B):
    """Breather image of P_x, q=2x/omega, omega=-2/B.

    b_x(theta+iv)b_x(theta-iv)=b_(x+hv/pi)b_(x-hv/pi)(theta).
    The fusion shifts are v=pi(1/2+B/(2h)).
    """
    return mp.fprod(block(x, theta, h) for x in
                    (h+q+B/2, q-B/2, h-q+B/2, -q-B/2))


def cdd_ratio(theta, h, B, replace_zero, replace_one):
    result = mp.mpc(1)
    if replace_zero:
        result *= pair_image(-B/2, theta, h, B)/pair_image(0, theta, h, B)
    if replace_one:
        result *= pair_image(1-B/2, theta, h, B)/pair_image(1, theta, h, B)
    return result


def dn_uniform(theta, h, B):
    b = lambda x: block(x, theta, h)
    return b(h-B/2)*b(h//2-1)*b(h//2+1-B/2) / (
        b(h-1)*b(1-B/2)*b(h//2)*b(h//2-B/2))


def dn_s(a, b, theta, h, B):
    def curly(x):
        return block(x-1, theta, h)*block(x+1, theta, h) / (
            block(x-1+B, theta, h)*block(x+1-B, theta, h))
    return mp.fprod(curly(abs(a-b)+p)*curly(h-abs(a-b)-p)
                    for p in range(1, 2*min(a, b), 2))


def bootstrap(first, n, B):
    h = 2*n-2
    family = {1: first}
    for a in range(1, n-2):
        old = family[a]
        family[a+1] = lambda t, a=a, old=old: (
            first(t+1j*mp.pi*a/h)*old(t-1j*mp.pi/h)*dn_s(1, a, 2*t+1j*mp.pi*(a-1)/h, h, B))
    return family


def check_cdd():
    worst = loop = mp.mpf(0)
    departures = []
    points = [mp.mpc('.47', '.19'), mp.mpc('-.93', '.27')]
    for n in (4, 5, 6):
        h = 2*n-2
        for flags in ((True, False), (False, True), (True, True)):
            B = mp.mpf('.6')
            first = lambda t: dn_uniform(t, h, 2-B)*cdd_ratio(t, h, B, *flags)
            family = bootstrap(first, n, B)
            for t in points:
                ratio = lambda coupling: cdd_ratio(t, h, coupling, *flags)
                loop = max(loop, abs(ratio(0)-1), abs(mp.diff(ratio, 0)))
                departures.append(abs(ratio(B)-1))
                for a, Ka in family.items():
                    worst = max(worst, abs(Ka(t)*Ka(-t)-1),
                                abs(Ka(t)*Ka(t+1j*mp.pi)/dn_s(a, a, 2*t, h, B)-1))
                for a in range(4, n-1):
                    alternative = (family[2](t+1j*mp.pi*(a-2)/h)*family[a-2](t-2j*mp.pi/h)
                                   *dn_s(2, a-2, 2*t+1j*mp.pi*(a-4)/h, h, B))
                    worst = max(worst, abs(alternative/family[a](t)-1))
    check('three alternative d_n CDD choices preserve one loop and vector-particle axioms',
          loop < mp.mpf('1e-35') and worst < mp.mpf('1e-35') and min(departures) > mp.mpf('.001'),
          f'loop {mp.nstr(loop, 3)}, axioms {mp.nstr(worst, 3)}, finite-B differences '
          f'{mp.nstr(min(departures), 3)}–{mp.nstr(max(departures), 3)}')
    # For a3 the P_(omega/2) replacement gives the same pair-image ratio.
    h = 4
    bI = lambda x, mu, T: mp.sin(mp.pi*(mu+x)/(2*T))/mp.sin(mp.pi*(mu-x)/(2*T))
    om = mp.mpf('.61'); T = 2*om; mu = mp.mpc('.19', '.13')
    sig = lambda m: bI(T/2-(om+1)/2, m, T)*bI(T/2+(om+1)/2, m, T)
    dev = abs(sig(mu)/mp.fprod(sig(mu+T-k*om) for k in range(1, 4))-1)
    for t in points:
        dev = max(dev, abs(mp.diff(lambda B: cdd_ratio(t, h, B, False, True), 0)))
        # Kim (33), coefficient in B (beta^2=2 pi B+O(B^2)).
        kim_loop = lambda z: 1j*mp.pi/8*(mp.sinh(z)/(mp.cosh(z)+1/mp.sqrt(2))
                                              -mp.sinh(z)/(mp.cosh(z)-1))
        minus = lambda B: -block(1, t, h)*block(h-B/2, t, h)*block(-1+B/2, t, h)
        dev = max(dev, abs(mp.diff(minus, 0)-kim_loop(t+1j*mp.pi)))
        def kim_exact(B):
            z = t+1j*mp.pi
            return mp.fprod(block(x-mp.mpf('.5'), z, h)*block(x+mp.mpf('.5'), z, h) /
                            (block(x-mp.mpf('.5')+B/2, z, h)*block(x+mp.mpf('.5')-B/2, z, h))
                            for x in map(mp.mpf, ('.5', '1.5', '2.5')))
        assert abs(minus(mp.mpf('.6'))/kim_exact(mp.mpf('.6'))-1) > mp.mpf('.01')
    check('a3 CDD ambiguity survives soliton bootstrap; shifted Kim equality holds only at one loop', dev < mp.mpf('1e-35'))


def check_dn_poles():
    # Load only root/block definitions; the source's unrelated E8 bootstrap
    # suite runs on import and is deliberately not part of this regression.
    ns = {}
    source = (Path(__file__).parent / 'snapshot/code/smatrix_simply_laced.py').read_text().split('\nB = 0.37;')[0]
    exec(compile(source, 'smatrix_simply_laced.py', 'exec'), ns)
    h, blocks, _, _ = ns['setup']('d', 5)
    powers = {}
    for b in range(5):
        for x, power in blocks[2, b].items():
            for y, k in ((x-1, power), (x+1, power), (x, -2*power)):
                powers[y] = powers.get(y, 0) + k
    poles = {x: powers.get(x, 0)-powers.get(2*h-x, 0) for x in range(1, h//2)}
    check('d5 particle 3 has simple poles at i pi/8 and 3i pi/8', poles == {1: 1, 2: -2, 3: 1}, str(poles))
    t0 = 3j*mp.pi/8
    amplitude = lambda t: mp.fprod(block(x, t, h)**p for x, p in powers.items() if p)
    residues = [e*amplitude(t0+e) for e in (mp.mpf('1e-8'), mp.mpf('1e-10'))]
    check('extra d5 pole has nonzero small residue', abs(residues[1]) > mp.mpf('.008')
          and abs(residues[1]-residues[0]) < mp.mpf('1e-8'), mp.nstr(residues[1], 9))

    eye = np.eye(5)
    roots = np.array([-eye[0]-eye[1]] + [eye[i]-eye[i+1] for i in range(4)] + [eye[3]+eye[4]])
    marks = np.array([1, 1, 2, 2, 1, 1])
    f0 = np.log(4-2*np.sqrt(2)); c = np.tanh(f0/4)
    field = lambda x: 4*np.arctanh(c*np.exp(2*x))*np.array([0, 1, 0, -1, 0])
    mass = lambda x: (roots.T*(marks*np.exp(roots@field(x))))@roots
    R = (roots.T*(np.sqrt(marks)*np.exp(roots@field(0)/2)))@roots/2
    spinors, middle = [0, 4], [1, 2, 3]
    coupling = max(np.linalg.norm(mass(x)[np.ix_(spinors, middle)]) for x in (-4, -1, 0))
    check('d5 spinor channels decouple in both bulk and boundary Hessians',
          coupling == 0 and np.linalg.norm(R[np.ix_(spinors, middle)]) == 0)
    eig, U = np.linalg.eigh(mass(-30))
    initial = (U*np.sqrt(eig-1))@U.T
    # Y=epsilon' epsilon^-1 for the decaying Jost subspace, Omega^2=1.
    def riccati(x, y):
        Y = y.reshape(5, 5)
        return (mass(x)-np.eye(5)-Y@Y).ravel()
    sol = solve_ivp(riccati, (-24, 0), initial.ravel(), method='DOP853',
                    first_step=.01, max_step=.2, rtol=2e-11, atol=2e-13)
    J = sol.y[:, -1].reshape(5, 5)+R
    smallest = np.linalg.eigvalsh((J+J.T)/2)[0]
    check('d5 extra pole at Omega=m is not a decaying bound state', sol.success and smallest > .1,
          f'min eigenvalue of Y(0)+R = {smallest:.9f}')


if __name__ == '__main__':
    mp.mp.dps = 45
    counts = {n: len(charge_classes(n)) for n in (2, 3, 4, 5)}
    check('exact missing single-soliton class counts', counts == {2: 0, 3: 2, 4: 5, 5: 21}, str(counts))
    check_cdd()
    check_dn_poles()
