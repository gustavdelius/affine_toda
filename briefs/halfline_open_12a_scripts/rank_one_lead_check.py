"""Check of the AKR lead (half-line open problem 12, part a).

Conventions
-----------
AKR (hep-th/0110218, eq. 1):  (1/8pi)(d phi)^2 + mu sum_i e^{b e_i.phi} + mu_B sum_i eta_i e^{b e_i.phi/2}
For a_1: e_1=sqrt2, e_0=-sqrt2.  With phi = sqrt2 phi' this is the FZZ/BRZ normalisation
    (1/4pi)(d phi')^2 + 2 mu cosh(2 b phi') + mu_B (e^{b phi'} + e^{-b phi'}),
i.e. BRZ (0710.4789) eq. (AB) with M^+ = M^- = mu_B.
FZZ / BRZ:  M = M0 cosh(pi b s),  M0 = sqrt(mu / sin(pi b^2)).
p = b^2/(1+b^2) = B_book/2 (book's B ranges over (0,2)).
Sine-Gordon continuation: b^2 -> -beta_SG^2/(8 pi) = -1/(lambda+1).

Run:  python3 rank_one_lead_check.py
"""
import mpmath as mp

mp.mp.dps = 30
pi = mp.pi


# ---------------------------------------------------------------- 1. convention chain
def chain(b2):
    b = mp.sqrt(b2)
    p = b2 / (1 + b2)
    muB2_over_mu = mp.cot(pi * b2 / 2) / 2                       # AKR (10)
    cosh2 = muB2_over_mu * mp.sin(pi * b2)                        # FZZ: cosh^2(pi b s)
    s = mp.acosh(mp.sqrt(cosh2)) / (pi * b)                       # principal branch
    C = mp.sqrt(muB2_over_mu * mp.sin(pi * b2))                   # M/M0
    # CT / continued UV-IR:  C = cos(pi a),  E = 2a(1-B/2) = 2a(1-p)
    a = mp.acos(C) / pi
    E_CT = 2 * a * (1 - p)
    # BRZ: cos(pi E/2) = cosh(p eta),  b eta = pi s
    E_BRZ = 2 / pi * mp.acos(mp.cosh(p * pi * s / b))
    return dict(b2=b2, s=s, s_expected=1j * b / 2, C=C, C_expected=mp.cos(pi * b2 / 2),
                E_CT=E_CT, E_BRZ=E_BRZ, p=p)


def chain_sg(lam):
    """Sine-Gordon side: continue AKR (10) with b^2=-1/(lam+1), mu_SG=-2mu, mutilde=-2mu_B."""
    b2 = -1 / (lam + 1)
    mu = mp.mpf(-1)                     # continued bulk coupling is negative (vacuum at 0)
    muB2 = mu / 2 * mp.cot(pi * b2 / 2)   # AKR (10), still positive
    mu_SG, mut2 = -2 * mu, 4 * muB2
    beta2_over_8 = pi / (lam + 1)
    mutcrit2 = 2 * mu_SG / mp.sin(beta2_over_8)                   # BPT (3.3)
    ratio = mp.sqrt(mut2 / mutcrit2)
    eta = (lam + 1) * mp.acos(ratio)                              # BPT (3.2) at phi0=0, theta=0
    eps1 = -mp.cos(eta) / mp.sin(pi * lam)                        # book (eq-bsg-eta-theta), vartheta=0
    return dict(lam=lam, muB2=muB2, ratio=ratio, target=mp.cos(pi / (2 * (lam + 1))),
                eta=eta, eps_hat=eps1)


# ---------------------------------------------------------------- 2. boundary energies
def E_bdy_BPT(E, p, m=1):
    """BPT (hep-th/0108157) eq. (sinh_benergy) = BRZ (Eb), symmetric boundary, as a function
    of the book's E (cos(pi E/2) = cos(eta_BPT/lambda) = cosh(p eta_BRZ))."""
    return m / (4 * mp.sin(pi * p)) * (2 * mp.cos(pi * E / 2) + 1 - mp.cos(pi * p / 2) - mp.sin(pi * p / 2))


def E_bdy_Fateev(p, h=2, m_particle=1):
    """Fateev hep-th/0103014 eq. (eQ) for the (+) boundary, A_1: Sigma_m = particle mass."""
    Sigma = m_particle
    return mp.sin(pi / (2 * h)) * Sigma / (8 * mp.sin(pi * p / (2 * h)) * mp.cos(pi * (1 - p) / (2 * h)))


def E_bdy_AKR_printed(p, h=2, m_particle=1):
    """AKR's printed Delta c_boun converted with E = pi Delta c/(24 R), mbar = m_particle/2."""
    mbar = m_particle / 2
    return pi / 24 * 6 * mbar * mp.cos(pi / (2 * h)) / (mp.sin(pi * p / (2 * h)) * mp.sin(pi * (1 - p) / (2 * h)))


# ---------------------------------------------------------------- 3. AKR delta_1
def delta1_AKR(b):
    """AKR eq. (deltas), delta_1 for the (++) strip, s1=s2=ib/2."""
    g = mp.euler
    f = lambda t: (2 * t * mp.cosh(b * t / 2) ** 2 / (mp.sinh(b * t) * mp.sinh(t / b))
                   - t * (mp.exp(-b * t) + mp.exp(-t / b)) / ((1 - mp.exp(-b * t)) * (1 - mp.exp(-t / b)))) / t
    def fs(t):
        with mp.workdps(60):
            return f(t)
    t0 = mp.mpf('1e-8')
    I = mp.quad(fs, [t0, mp.mpf('1e-3'), 1, 10, mp.inf]) + t0 * fs(t0)
    return -2 * g * (b + 1 / b) - 2 * (b - 1 / b) * mp.log(b) + 2 * I


def ddelta1(b, s):
    """delta_1(s1=s2=ib/2) - delta_1(s1=s2=s): from BRZ (Dint), d/dP at P=0."""
    f = lambda t: 2 * (mp.cosh(b * t / 2) ** 2 - mp.cos(s * t) ** 2) / (mp.sinh(b * t) * mp.sinh(t / b))
    t0 = mp.mpf('1e-8')
    return 2 * (mp.quad(f, [t0, 1, 10, mp.inf]) + t0 * f(t0))


if __name__ == "__main__":
    print("1. Convention chain, sinh-Gordon side")
    for b2 in [mp.mpf('0.1'), mp.mpf('0.4286'), mp.mpf('0.8086')]:
        d = chain(b2)
        print(f"  b^2={float(b2):.4f}  s={mp.nstr(d['s'],12)} (ib/2={mp.nstr(d['s_expected'],12)})"
              f"  M/M0={mp.nstr(d['C'],12)} cos(pi b^2/2)={mp.nstr(d['C_expected'],12)}"
              f"  E_CT={mp.nstr(d['E_CT'],12)} E_BRZ={mp.nstr(mp.re(d['E_BRZ']),12)} p={mp.nstr(d['p'],12)}")
    print("   sine-Gordon side (book's cnj-bsg-uv-ir at phi0=0)")
    for lam in [mp.mpf('1.7'), mp.mpf('5.3'), mp.mpf('12.3')]:
        d = chain_sg(lam)
        print(f"  lambda={float(lam)}  mu_B^2={mp.nstr(d['muB2'],8)}>0  mutilde/mutilde_crit={mp.nstr(d['ratio'],14)}"
              f"  cos(pi/(2(lam+1)))={mp.nstr(d['target'],14)}  eta={mp.nstr(d['eta'],14)} (pi/2={mp.nstr(pi/2,14)})"
              f"  eps_hat={mp.nstr(d['eps_hat'],3)}")

    print("\n2. Boundary energy: Fateev (+) at h=2 vs BPT/BRZ at E=p (eps_hat=0) and E=0 (M=M0)")
    for p in [mp.mpf('0.05'), mp.mpf('0.2'), mp.mpf('0.3'), mp.mpf('0.4471'), mp.mpf('0.7')]:
        F = E_bdy_Fateev(p)
        print(f"  p={float(p):.4f}  Fateev={mp.nstr(F,15)}  BPT(E=p)={mp.nstr(E_bdy_BPT(p,p),15)}"
              f"  BPT(E=0)={mp.nstr(E_bdy_BPT(0,p),10)}  AKR-printed={mp.nstr(E_bdy_AKR_printed(p),10)}"
              f"  dual: Fateev(1-p)={mp.nstr(E_bdy_Fateev(1-p),12)} BPT(E=1-p)={mp.nstr(E_bdy_BPT(1-p,p),12)}")
    print("   classical limit m/(2 pi p) at p=0.05:", mp.nstr(1 / (2 * pi * mp.mpf('0.05')), 10))

    print("\n3. AKR delta_1: s=ib/2 (AKR table, TBA agreed to 1e-5) versus s=0 (M=M0)")
    table = {0.25: -8.12340, 0.30: -6.15569, 0.35: -4.73422, 0.40: -3.68012, 0.45: -2.89175,
             0.50: -2.30886, 0.55: -1.89627, 0.60: -1.63628, 0.65: -1.52627, 0.70: -1.58079,
             0.75: -1.84022, 0.80: -2.39472}
    for B, d1tab in table.items():
        b = mp.sqrt(mp.mpf(B) / (1 - mp.mpf(B)))
        d1 = delta1_AKR(b)
        d0 = d1 - ddelta1(b, 0)
        print(f"  B={B:.2f}  delta_1(ib/2)={mp.nstr(d1,7)} (AKR table {d1tab})   delta_1(s=0)={mp.nstr(d0,7)}")
