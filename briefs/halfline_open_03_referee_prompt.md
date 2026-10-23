# Referee prompt: half-line open problem 03
Tier: B
Brief: briefs/halfline_open_03_neumann_soliton_route.md
Report: briefs/halfline_open_03_report.md
Scripts: briefs/halfline_open_03_scripts/ (the worker's own files only)
- **Setup.** The scripts import unchanged copies of the book modules, so copy `code/*.py` into the same scratch directory first.
- **Main entry points.** `check_all.py` runs checks 1–9 in about 4 minutes; `python3 check_all.py full` adds the $n=5$ bootstrap. `closed_forms.py` gives $A^{(\pm)}$, and `particle_checks.py` does the continuation and the axioms.
- **Other files.** The `t*.py` and `fit*.py` files are exploration.

Scope: the full brief, including the $s=-1$ bonus. This is the most consequential report of the batch. If it stands, it moves `cnj-neumann-coideal` and `cnj-neumann-amplitude` from conjecture to [Numerical] for $n\le5$. Tier B because it claims a proof (claim 3) and corrects premises of the brief and the book (claims 1, 2).

## Claims to referee

1. [Correction of the brief and book] **Boundary unitarity for the non-diagonal $K$.** It is $A(\theta)A(-\theta)\,K(-\theta)^{\rm T}K(\theta)=1$, not $K(\theta)K(-\theta)\propto\mathbf 1$. The antisoliton is reflected by $JKJ=K^{\rm T}$ ($J$ = basis reversal = charge conjugation, from the diagram automorphism $j\leftrightarrow n{+}1{-}j$, which fixes $\hat\epsilon$). $K(\theta)K(-\theta)$ is not scalar for $n\ge2$.
2. [Correction of book code] **The branch.** `k_closed` in `code/boundary_coideal.py` uses `np.sqrt(-q)` and a positive $\epsilon_*$. For odd $N$ and $\cos\frac{\pi\omega}2<0$ that is the wrong sign class. The book checks only at $\omega=0.37$, so no statement in the book is wrong.
3. [Theorem] **The unitarity function.** $K(-\theta)^{\rm T}K(\theta)=U\mathbf 1$, $U=(s\cosh2T\theta-\cos\pi T)/(2\sin^2\frac{\pi\omega}2)$, $T=h\omega/2$. The proof goes through $K$ as a polynomial in the twisted shift $Z$ ($Z^N=s$) and the eigenvalue $\kappa(\zeta)$. The same $\kappa$ reproduces the $\det K$ formula of ch. 34.
4. [Numerical] **The crossing function.** $\Phi_+=\sin\pi(\frac\omega2-\mu)/\sin\pi(\frac\omega2+\mu)$ and $\Phi_-=\cos\pi(\mu-\frac\omega2)/\cos\pi(\mu+\frac\omega2)$, $\mu=-iT\theta/\pi$, for $n=1..5$ (1.3e-14). For $n=2$ the resulting functional equations are Gandenberger's (4.4), (4.6).
5. [Numerical] **The singlet projection.** $(1\otimes K)\check R(1\otimes K^{\rm T})$ maps the $B_1$ singlet to $c_\pm(\mu)$ times itself (closed forms in report §2), for $n=2..5$ (1.7e-14).
6. [Numerical] **The scalar factors $A^{(\pm)}$** (report §3), solving $A(\mu)A(-\mu)=1/U$ and $A(\frac T2-\mu)/A(\frac T2+\mu)=F_{11}(2\mu)/\Phi(\mu)$ to 3.5e-8 for $n=1..5$. Checks:
   - $n=1$: they equal GZ/BPTT at the Neumann point and at $(\frac\pi2,\frac{i\pi\omega}2)$.
   - $n=2$: $A^{(+)}$ is Gandenberger's (4.14) times $\sigma^+$ of (5.3).
   - $A^{(-)}$ is said to differ from Gandenberger's (4.15) by a constant $-i$; the worker blames garbled PDF text and did not pursue it.
7. [Numerical] **The soliton bootstrap** $V^{\otimes n}\to\bar V$, measured relative to $K=1$, holds for $A^{(+)}$ at $n=2..5$ and for $A^{(-)}$ at $n=2..4$.
8. [Claim about CDD factors] **What fixes the CDD factor.**
   - The zero-momentum pole at $\theta=i\pi/2$ requires $(\frac T2)_I^2$.
   - The factor $(\frac T2-\frac\omega2)_I(\frac T2+\frac\omega2)_I$ is forced by the bootstrap for even $n$, and is allowed but not forced for odd $n$, where it is DG's $\sigma_{(n-1)/2}$. Without it $K_1$ changes by DG's $k_{h/2-1}$ (3.25).
   - Kim's one-loop $a_3$ result requires the factor at $n=3$.
   - An unexplained pole at $\theta=i\pi(\frac12-\frac1h)$ lies in the physical strip.
9. [Numerical] **The main result.**
   - $A(\mu_a)F_{11}(2\mu)A(\mu_b)c_+(\mu)=(-\frac{n\omega}2)_I(\frac{n\omega}2-\frac12)_I(\frac12)_I$ for $n=1..5$ (9.5e-8).
   - Under $\omega\to-2/B$ this becomes $(-n)(n+\frac B2)(-\frac B2)=K^N_1$, a block identity for $n\le7$.
   - $K^N_a$ satisfies the axioms of `prp-uniform-amplitude` and tends to 1.
10. [Numerical] **The $s=-1$ bonus.**
    - The breather amplitude is DG(3.11)$\times\sigma_1$, which continues to $K^{(-)}_a=-(a)(h-a)K_a^{\rm unif}$.
    - At $n=2$ this equals Gandenberger's (6.4), $-(1)(-1+\frac B2)(3-\frac B2)$.
    - The axioms hold for $n\le7$, and the amplitude tends to 1.
11. [Numerical] **Semiclassics.** $\frac{\beta^2}{4h}\frac d{d\theta}(-i\ln A^{(+)})$ tends to the Neumann time advance $\frac12\ln(1-\sin^2\frac\pi h/\cosh^2\theta)$, with errors falling like $1/\omega$.
12. [Text] **The proposed replacement text and status labels** (report §7), including the new `prp-an-nondiagonal-scalar` and the change of `prp-an-nondiagonal` (iii) to "$n\le5$".

## Coordinator notes (suspected weak points)

- **Claim 8 is where a wrong answer would hide.** The CDD factor is chosen by the worker. Check the three steps separately.
  - Is the rank-one-$K$ argument for the $\theta=i\pi/2$ pole sound, and is the pole present in GZ at $n=1$?
  - Does $(\frac T2)_I^2$ alone really fail the bootstrap for even $n$? Test it yourself.
  - The unexplained physical-strip pole at $\theta=i\pi(\frac12-\frac1h)$ is a genuine open issue. A pole in the physical strip needs a boundary bound state or a Coleman–Thun explanation. Check whether Gandenberger explains it for $n=2$.
- **Claim 9's precision is only 1e-7**, because of truncated infinite products. Evaluate the breather identity independently, for example from an integral representation of the blocks or with mpmath and a much longer product, at one $(n,\omega,\theta)$ for $n=3$ and $n=4$. Test the block identity after $\omega\to-2/B$ symbolically.
- **Claim 7, normalization.** The bootstrap is "measured relative to $K=1$, so that bulk factors and embedding normalizations cancel". Check that this ratio method cannot hide an error, for example a missing overall scalar that the ratio cancels.
- **Claim 1.** Find what the book currently says about boundary unitarity for non-diagonal $K$ (`sec-boundary-unitarity-crossing`, `sec-an-nondiagonal-k`), and decide whether it is wrong or only incomplete. Check the claim that $K^{\rm T}$ is the intertwiner $\bar V(\theta)\to V(-\theta)$ in the book's conventions, including the sign of $q$.
- **Claims 6 and 10 against Gandenberger** (hep-th/9806003). Read his TeX source rather than PDF text. Settle the $-i$ in (4.15), and check (4.14), (5.3), (5.16)–(5.18) and (6.4) directly.
- **Claim 8, Kim.** Check against `kim1995` eqs. (33) and (39) at $O(\beta^2)$.
- **Prior art.** Has anyone published the scalar factor of the non-diagonal $a_n^{(1)}$ K-matrix (Gandenberger 1999b left it open), or the Neumann amplitudes from the soliton side? Search for later work, for example by Doikou, Nepomechie, Avan–Doikou, Regelskis–Vlaar, Appel–Vlaar, and Delius–MacKay.
