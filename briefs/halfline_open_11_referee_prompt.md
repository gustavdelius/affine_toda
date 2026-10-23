# Referee prompt: half-line open problem 11
Tier: B
Brief: briefs/halfline_open_11_soliton_preserving_scalar.md
Report: briefs/halfline_open_11_report.md
Scripts: briefs/halfline_open_11_scripts/
- **Setup.** The scripts import unchanged copies of `code/*.py` and of `gz.py`, `nk_common.py`, `rhonp.py` from `briefs/halfline_open_03_scripts/`; copy all of these alongside.
- **Core modules:** `sp_core.py` (the matrices), `scal.py` (the scalar factors), `bproj.py` and `breather.py` (the breathers). The `e*.py` and `t*.py` files are the individual checks listed in the report's table.

Scope: the full brief. Its step 5 was meant to test the tension "soliton-preserving at real coupling looks like Neumann classically, but continues to Dirichlet for $n=1$". The worker was given problem 3's (Neumann = soliton-conjugating $s=+1$) report as a method template. Tier B because of two claims proved by hand (claims 3, 5) and two corrections to the book (claims 2, 4).

## Claims to referee

1. [Numerical] **The SP matrices.** In book conventions (principal gradation) the diagonal solutions of the braid reflection equation form one family: a cyclic block $S$ of length $l$ and a parameter $\eta$,
   $K_V=\mathrm{diag}(e^{a_j\omega\theta})\,\mathrm{diag}(\sinh(\eta\pm T\theta))$ (report §1).
   - The braid reflection equation holds to 1e-15 for $n=2,3$.
   - The gradation change maps it exactly to de Vega–González-Ruiz.
   - The coideal is $X=$ all nodes except $0,l$, with $\tau:j\to l-j$ (by hand, unverified).
2. [Book correction] **The DVGR statement in `sec-preserving-conjugating`, marked "unchecked in book conventions", is correct.** Gives the replacement formula.
3. [Theorem, by hand] **The $\mathbb Z_N$ obstruction.** For $n\ge2$ no invertible diagonal solution of the braid reflection equation is $\mathbb Z_N$-covariant: a Möbius function of $y^N$ would have to equal $y^{2-N}$. Since Delius's boundary at a vacuum is $\mathbb Z_N$-invariant and charge conservation forces diagonal $K$, that boundary has no $\mathbb Z_N$-symmetric quantum $K$ in $V$.
4. [Book correction, by hand] **`sec-soliton-preserving` is wrong to exclude coweights outside $2\Lambda_W^\vee$.** With $\sigma=\dagger$ and constant diagonal $K_+$, $\operatorname{Re}\phi(0)\in\frac\pi\beta\Lambda_W^\vee$ is allowed with $K_+=e^{i\pi\mu\cdot H}$. Coweights outside $2\Lambda_W^\vee$ fail only with $K_+=1$.
5. [Theorem, by hand] **Real-coupling Dirichlet is not classically integrable for $n\ge2$.** $\phi(0)=\phi_0$ is not integrable with a field-independent $K$, any linear $\sigma$, in the fundamental representation. The argument: $\sigma=-\mathrm{Ad}_K$ on $\mathfrak h$, so $\psi_j=\lambda^2$ at every node, and $\lambda^{2N}=1$ fails. Sine-Gordon escapes because $E_{\alpha_0}=E_{-\alpha_1}$.
6. [Numerical] **The antisoliton matrix** is fixed uniquely by the mixed reflection equation. It is $K_{\bar V}$ with the shifted parameter $\bar\eta=\eta+i\pi\omega(l-N/2)$. Unitarity is $K(\theta)K(-\theta)=U\mathbf 1$.
7. [Numerical] **Crossing.** The function is $\Phi_{1,2}=-(-q)^{\pm(1-l)}\sin\pi(2\mu-T)/\sin\pi(2\mu+T)$, independent of $\eta$ (fitted, $n\le4$, two values of $\omega$), and the transmission entries of $\check R_{\bar VV}$ equal 1.
8. [Numerical] **The minimal scalar factors $a,b$** (report §3, $\hat E_P$ and $\Lambda_P$ products) solve the four functional equations to about 1e-7 for six $(n,l)$. The identity $F_{11}(x-2T)=F_{11}(x)\cdots$ is stated.
   - For $n=1$, GZ Dirichlet is the minimal solution times $(-T/2)_I^2$.
   - The CDD factor for $n\ge2$ is not fixed; the candidates are $\sigma_A$ and $\sigma_B$.
9. [Numerical] **The breathers.**
   - The singlet in $V(\theta_a)\otimes\bar V(\theta_b)$ maps exactly onto the singlet in $\bar V(-\theta_b)\otimes V(-\theta_a)$, i.e. $B_1\to\bar B_1$.
   - The amplitude is $K_{1\to n}=-K_0E_\eta\times$CDD, with $K_0$ a product of six blocks fitted for $n\le3$ (residual ≤5e-6) and $E_\eta$ derived analytically.
   - For $n=1$ it equals Ghoshal's Dirichlet $K_{B_1}$ (1e-6), and at $\eta=i\pi/2$ it continues to $K_D$.
10. [Claim] **Step 5 verdict.**
    - At $\eta=i\pi/2$ the $B\to0$ limit of $K_{1\to n}$ is $\theta$-dependent for $h\ge3$ for every $l$ and $\eta$, so it is neither Dirichlet ($-1$) nor Neumann ($+1$), and particles reflect $a\to h-a$.
    - The $n=1$ Dirichlet result is a sine-Gordon accident.
    - There is no conflict with problem 3, and no second soliton route to the Neumann amplitudes.
11. [Text] **The proposed book text**: `sec-preserving-conjugating`, a new ch. 35 paragraph, open problem 11 marked partially solved, and the `sec-soliton-preserving` correction.

## Coordinator notes (suspected weak points)

- **Claim 9, $B_1\to\bar B_1$, is the physical crux.** Check that the singlet in $V\otimes\bar V$ and the singlet in $\bar V\otimes V$ really are the two different particles 1 and $n$ in the book's conventions (`sec-breather-reflection`), and not the same particle $B_1$ written in the other tensor order. The order of the rapidities is reversed by reflection. If the identification is a labeling artefact, claim 10 changes substantially. Compare with how the book projects onto $B_1$ for the soliton-conjugating matrices, where $V\to\bar V$ under reflection.
- **Claim 3.** Re-derive the obstruction. Check whether $\mathbb Z_N$ covariance needs to hold exactly, or only up to a diagonal gauge, a twist or a rescaling of $K$ by a scalar function of $\theta$. A gauge-covariant version might not be obstructed. The book's $Z$ in problem 3 is a *twisted* shift for $s=-1$.
- **Claim 4.** Recheck with `code/boundary_classical.py`, or by hand from the book's derivation in `sec-soliton-preserving` (ch. 30), whether $K_+=e^{i\pi\mu\cdot H}$ with $\operatorname{Re}\phi(0)\in\frac\pi\beta\Lambda_W^\vee\setminus2\Lambda_W^\vee$ satisfies the classical boundary condition. Check also whether such a $\phi(0)$ is a vacuum, as the brief's trap requires.
- **Claim 5 needs a prior-art check.** Bowcock–Corrigan–Dorey–Rietdijk (hep-th/9501098) and later work by Corrigan and collaborators classified the integrable boundary conditions of affine Toda theory; Dirichlet conditions may already be treated there. Check whether "real Dirichlet is not integrable for $n\ge2$" is known, and whether the restriction "field-independent $K$, linear $\sigma$, fundamental representation" leaves a loophole.
- **Claims 8–9.** $K_0$ is a *fit* (residual 5e-6). Evaluate the breather projection independently at higher precision at one point each for $n=2,3$, and confirm the block content. Recompute the $n=1$ comparison with Ghoshal's $K_{B_1}$ yourself.
- **Claim 6.** Verify $\bar\eta=\eta+i\pi\omega(l-N/2)$ independently by solving the mixed reflection equation at one $(n,l)$.
- **Prior art.** Doikou's and Doikou–Nepomechie's soliton-preserving work (`doikou2000`), and Arnaudon–Avan–Crampé–Doikou–Frappat–Ragoucy (math-ph/0406021), may already contain the SP crossing function, the antisoliton shift $\bar\eta$, or the scalar factors. Check.
