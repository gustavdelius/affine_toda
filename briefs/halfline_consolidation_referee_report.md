# Updated referee report: half-line consolidation (PR #14)

**Date:** 2026-10-14; outstanding numerical checks completed later the same day.
**Referee:** Codex (GPT-6).
**Book commit reviewed:** `3a75d316d3679c9dae6f9026c69a3679ecbafbcd`, on
`half-line-consolidation`. **Report re-examined:** the Opus report of
2026-10-13, preserved in commit `86ffbda`. The original item numbers A1–A5,
B1–B6 and O.1–O.5 are retained. They are not the consolidation workers' A1–A6.

This report supersedes that report's verdicts where they differ. It records
completed work from this review; it does not claim to review subsequent book
commits. The current checkout is `main`; the historical sources and results
were recovered read-only. No book edit, commit, push or PR update accompanies
this report, as requested by the editor.

**Recommendation:** retain the formulas with the qualifications below. The
main corrections concern the scope of the Dirichlet obstruction, missing
single-soliton charges, CDD nonuniqueness, an extra classical pole, and the
rapidity at which the mixed-multiplet instability appears. Reject the original
referee's proposed spinor-coupling explanation of the extra pole. Do not count
its inaccurate singular-endpoint calculation as independent confirmation of
a limiting value.

**Evidence policy.** The derivations below were reconstructed from the Lax
condition, energy functional, root Hessians and block definitions. The new
charge/CDD/pole checks, matrix Riccati scan, and high-precision Hirota
fluctuation calculation were written during this review. Re-executions of the
original referee's and book's scripts are supplementary checks, not independent
re-derivations. Their failures and incomplete runs are recorded explicitly.
Primary sources were read in TeX or extracted PDF text. Numerical scans and
identities between conjectured amplitudes do not establish their physical
boundary identifications or arbitrary-rank validity.

**Labels.** A1, A2 and A4 admit [Theorem] only with the hypotheses and proofs
below. B1–B6 retain [Numerical] for the tested identities or spectra; continuation
and boundary identifications retain their conjectural assumptions. The CDD
counterexamples disprove uniqueness from one loop, not the displayed amplitudes.

## Summary and proposed changes

| Item | Verdict and re-derived conclusion | Proposed change |
|---|---|---|
| A1 | **CONFIRMED WITH CHANGES.** The evaluation-representation obstruction is correct. The lead-in overstates its scope; individual non-simply-laced boundary parameters can be unbounded. The proof needs the invariant-image ideal argument. | Ch. 30: qualify the scope, say “cannot all”, specify evaluation representations, supply the cycle and ideal arguments. Ch. 35 problem 7: same scope. |
| A2 | **CONFIRMED WITH CHANGES.** Half-coweight compatibility and the reflection-equation statement survive. The inference from a fundamental weight to a single-soliton charge does not. Missing-class counts are 0, 2, 5, 21 for ranks 2–5. | Qualify the half-soliton and three-soliton constructions in chs. 30–31 and `tbl-boundary-types`; explicitly record the missing classes. |
| A3 | **CONFIRMED WITH CHANGES.** Delius (1.4) gives a half-coweight value; (1.5) follows only for even coweights. A vacuum shift changes the image coweight by twice a coweight. | Retain the correction in Appendix D, with A2's qualification. |
| A4 | **CONFIRMED WITH CHANGES.** Strict convexity gives uniqueness. Coercivity and the direct method also give existence for every rank when all C_j are nonnegative. The d4/d5 kinks and energies agree. | Add the function space and existence argument; problem 2 should ask for closed forms of higher-rank uniform vacua. |
| A5(a–d) | **CONFIRMED WITH CHANGES.** Kim–Cho is the d4 one-loop source; Kim's vector rule and CAF's shifted self-dual solution coincide; Gandenberger's minus breather and CDRS's conjectures have the stated conventions. Kim's figure ends at d6. | Retain the attributions; specify that d7 uses Kim's stated rule. |
| B1(a–f) | **CONFIRMED at the tested points.** Matrix unitarity, crossing, singlet coefficients, scalar equations, breather identities, continuation, and one-loop comparisons survive at the tested ranks and points. | Keep the formulas and numerical labels. |
| B1(g) | **CONFIRMED WITH CHANGES.** One loop does not determine Xi. The three alternatives preserve the tested scalar and vector-particle constraints and differ at finite coupling. The rho-plus comparison also needs its extra CDD factor. | Correct chs. 33 and 35 and problems 2–3. Make boundary identifications conditional on continuation and the chosen CDD factor; give the precise rho ratios. |
| B2 | **CONFIRMED WITH CHANGES; proposed mechanism WRONG.** The diagonal classical reflection calculation survives the reruns. An independent stable scan now excludes bound-state zeros at 110 subthreshold frequencies each for d4 and d5, with domain and tolerance checks. Particle 3 of d5 has the extra simple pole at 3i pi/8. **The referee's proposed inheritance from a coupled spinor channel is false:** both Hessians decouple the spinors from particle 3. | Correct the pole list, state the coincident frequency without asserting that mechanism, and keep the non-bound-state check. Keep the absence-of-bound-states statement numerical and limited to d4/d5; d6 was not scanned. |
| B3 | **CONFIRMED at the tested points.** The SP crossing functions, scalar equations and phased breather coefficient survive. | Keep the proposition and amplitude. |
| B4 | **CONFIRMED WITH CHANGES.** A new independent Hirota calculation verifies all 56 regular backgrounds, the middle-soliton limits, distinct stationary-pair limits, and the exceptional a4 limit along the specified paths. The original ODE evidence was insufficient; the new convergent computations supply the missing verification. | Keep the withdrawal of “BP Table 2 corrected”. Specify the approach and target Hessian, and limit path-independence claims to the paths tested. These calculations do not implement BP's prescription at the singular wall. |
| B5 | **CONFIRMED WITH CHANGES.** The mirror convention, cyclic sectors and numerical thresholds survive. The slow-rapidity-zero description is wrong for V Vbar. The unmirrored convention coincides with the mirror convention for a uniform left wall. | Correct the onset and mirror descriptions. Round the finer V Vbar estimates to 0.667 and 0.500. Add a regression exhibiting breaking at finite slow rapidity while the near-zero point remains unbroken. |
| B6 | **CONFIRMED at the tested points.** The minus scalar, breather, continuation, bootstrap and CDRS dual comparison survive. | Keep the formulas. |
| O.1 | **CONFIRMED.** CAF explicitly gave the B=0 product and associated it with fixed boundary conditions. | Credit CAF for that formula; distinguish the numerical identification with the uniform kink. |
| O.2 | **CONFIRMED.** The shifted a3 Kim agreement is perturbative, not an exact equality. | Say “agrees at one loop with Kim's (33)”. |
| O.3 | **CONFIRMED.** The a3 alternative P_(omega+1)/2 survives the soliton bootstrap and one loop. | Correct the CDD discussion and its dependent ch. 35 statements. |
| O.4 | **CONFIRMED WITH CHANGES.** The report did not check semiclassics. An additional Gamma-product calculation supports equality of the **leading phase derivatives**, not equality of the finite-coupling phases. | State the observable actually checked: leading semiclassical phase derivatives. |
| O.5 | **historical build CONFIRMED.** Cross-references and citations require a fresh render after these edits. | The earlier HTML render completed; this report makes no changes requiring a new book render. |

## Independent re-derivations

### A1: full Dirichlet

For fixed phi(0), its normal derivative is free. The Cartan coefficient of
compatibility therefore vanishes independently of the step coefficients.
After absorbing conjugation and replacing transpose by the dual representation,
K preserves weights between rho and rho-prime. For rank at least two the
2(r+1) root shifts are distinct. Each nonzero node coefficient cancels, leaving

`K rho(E_±j) = lambda^(∓2) rho'(E_±j) K`.

The finite negative simple generators construct E_-theta in h-1 brackets.
Their factor lambda^(2h-2), compared with lambda^-2 from the affine node,
gives `(lambda^(2h)-1) rho'(E_-theta) K = 0`. The image is invariant; its
annihilator in simple g is an ideal containing E_-theta, hence all of g.
Finite-dimensional complete reducibility excludes a nonzero trivial image.
This is an evaluation-representation argument, not the same proof for arbitrary
tensor products at different evaluation points. At rank one the root shifts
overlap and the coefficient cancellation fails, explaining the control case.

For antilinear sigma, independent real and imaginary normal derivatives make
both Cartan products vanish. K has only zero-weight support; the remaining
step terms map to different weight spaces, forcing the root generators to
annihilate its image as well. The spin-2 flux gives an independent obstruction:
the velocity term of Q2+c Q-2 has nonzero curl unless c=1; at c=1 a nonzero
cubic polynomial in the free normal derivative remains. This is not a theorem
excluding every conceivable definition of integrability.

### A2–A3: half-coweights and missing charges

For K=exp(i pi nu.H), conjugation of E_j has sign
`(-1)^(alpha_j.nu)`. The boundary value gives exactly the same sign to
`conjugate(w_j)/w_j`; the imaginary normal derivative vanishes. In the defining
a_n representation the free real derivative forces K diagonal. A nonzero entry
propagates around the affine cycle, and the step equations then force these
signs and determine K up to scale. Diagonal K reduces the reflection equation
term by term to its identity-matrix case. This establishes compatibility and
the reflection equation; the pre-existing caveat about involution remains.

The static one-soliton tau functions are
`tau_j=1+exp(m_b x+2 pi i (r+jb/h))`. A topological charge changes only when a
ray crosses the negative real axis. Enumerating the *exact rational intervals*
between `r=1/2-jb/h mod 1`, rather than scanning arbitrary phases, gives h/gcd(b,h)
charges. Reducing their Dynkin coordinates modulo 2 and taking the union over b
gives the missing-class counts in the table. In a3 the missing parity vectors
are (0,1,0) and (1,0,1), the classes of epsilon1+epsilon2 and epsilon1+epsilon4.
A complementary pair of weights need not occur among these charges.

Self-mirroring makes the energy density at x the conjugate of that at -x,
so the real half-line energy is half the soliton mass whenever the solution is
regular. The three-soliton image needs both interaction shifts, including
`Re xi3=-log A_ab(theta)`; compatibility alone does not supply a soliton in
each class. At x=0 Delius's image condition gives Im phi_D=-pi lambda, and
`phi_D -> phi_D+2 pi i mu` changes lambda to lambda-2mu. This independently
recovers the factor-of-two correction and the limitation of its half-soliton
interpretation.

### A4: variational vacuum

The Hessian of V is `sum n_j exp(alpha_j.phi) alpha_j alpha_j^T`, positive
definite; the boundary Hessian is positive semidefinite for C_j>=0. The positive
affine relation among spanning roots implies a positive lower bound for
`max_j alpha_j.u` over unit u. Thus V grows exponentially in every direction
and, together with its positive Hessian at its unique minimum, satisfies
`V(phi)>=c|phi|^2` globally. The boundary potential is bounded below. The energy
is coercive on H1 of the half-line. A weakly convergent minimizing subsequence
has a locally convergent subsequence and a convergent trace; lower
semicontinuity and Fatou give a minimizer. H1 embeds in L-infinity in one
dimension and its elements tend to zero at the infinite end, legitimizing
the field variation. Strict convexity proves uniqueness and global minimality.
Convexity alone, without coercivity, would not prove existence.

Diagram invariance reduces d4 to phi=f(e2-e3), d5 to phi=f(e2-e4).
The equations are f''=2(exp(2f)-exp(-f)) and f''=4 sinh(f). The first integral
and boundary condition give exp(f0)=9/8 and 4-2 sqrt(2). Quadrature gives
energies -0.0695214587295478 and -0.0927092091312040 in m/beta^2 units.

### B1 and O.2–O.3: what one loop can fix

Let P_x=(T/2-x)_I(T/2+x)_I and q=2x/omega. The breather constituents have
opposite shifts `v=pi(1/2+B/(2h))` after continuation. The block identity

`(z)(theta+iv) (z)(theta-iv) = (z+hv/pi)(theta) (z-hv/pi)(theta)`

therefore maps P_x to

`D(q;B)=(h+q+B/2)(q-B/2)(h-q+B/2)(-q-B/2)`.

Here parentheses are real-coupling blocks. At B=0, D(q;0)=1 for every q.
Replacing q=0 by -B/2 or q=1 by 1-B/2 consequently changes D only at O(B^2).
Their ratios have value 1 and derivative zero at B=0, so they cannot be
distinguished at one loop. Direct 45-digit checks at new complex rapidities
give residuals below 4e-45 for the continued axioms, while finite-B ratios
depart from 1 by 0.004–0.23. The a3 replacement obeys the soliton-bootstrap
product too. Exact duality specifies the displayed *breather amplitude*;
this argument does not prove uniqueness of the underlying soliton CDD factor
against factors invisible to that breather.

The prefactor comparison is literal: dividing the a_n minus amplitude by its
diagonal amplitude gives rho-minus. For plus, dividing out P_(omega/2) leaves
one inverse (-T/2)_I, whereas rho-plus contains one positive such factor. Their
ratio is (-T/2)_I^2, which is homogeneous crossing-symmetric and unitary.

Kim (33) is the coefficient of beta^2. Differentiating the continued minus
amplitude at B=0 agrees after theta -> theta+i pi; evaluating the full
self-dual product at B=.6 does not. Only that perturbative agreement is justified. These distinctions also apply to conditional boundary identifications.

### B2: the extra d5 pole and the incorrect spinor explanation

At B=1, expanding every curly block as `(x-1)(x+1)/(x)^2` gives orders +1,
-2, +1 for particle 3 at theta=i pi z/8, z=1,2,3. The z=3 pole has residue
approximately `-0.00814729956 i`; a magnitude cutoff can easily miss it.
Since m3=2 sqrt(2) sin(3 pi/8) m, its momentum is i(1+sqrt(2))m and its
frequency is m. This kinematic correction is verified.

However, in the stated orthonormal root convention its fluctuation Hessian is

```
M^2(f) = diag_spinors(2 exp(-f), 2 exp(-f)) plus the middle block
         [[4 cosh(f), -2 exp(f), 0],
          [-2 exp(f), 4 exp(f), -2 exp(f)],
          [0, -2 exp(f), 4 cosh(f)]].
```

The spinors are e1 and e5. Their off-block entries vanish identically, as do
those of the boundary Hessian. Particle 3 belongs to the middle block, which
can mix the vector-like masses; it cannot inherit this pole by coupling to a
spinor. Record the equal frequency, but leave the detailed pole mechanism unassigned. Independently evolving the Riccati matrix of decaying
solutions at Omega^2=1 gives minimum eigenvalue 0.897061476 for Y(0)+R,
excluding a bound-state zero there. This also avoids exponentially ill-scaled
raw determinants. The completed independent scan described below extends this
check to 110 subthreshold frequencies each for d4 and d5, including points
within 1e-10 of the continuum threshold. It is numerical evidence, not an
interval-certified proof. No independent d6 bound-state scan is available.

### B3–B4 and B6: identities and limiting backgrounds

Crossing for an SP matrix pairs V with Vbar, not V with V. Contracting its
diagonal entries with the transmission-normalized mixed R-matrix gives the
two factors `-(-q)^(1-l) s(mu)` and `-(-q)^(l-1) s(mu)`. The lowest-breather
singlet changes order to Vbar tensor V, retaining the convention-dependent
phase exp(-2 pi i l/h). The scalar equations and breather identity are checked
with independently summed Gamma products.

For a real continued Hirota background, adding one infinitesimal soliton
gives `epsilon=-sum alpha_j delta(tau_j)/tau_j`. Differentiating the field
equation verifies this is an exact linearized solution. Imposing
`epsilon'(0)+Hess(B) epsilon(0)=0` gives its reflection coefficients; direct
ODE integration provides an independent comparison on regular backgrounds.
At a singular endpoint both this ratio and the boundary Hessian are singular,
so a path must be specified. Finite-distance ODE values that differ from a
limiting value do not by themselves prove different limits. The high-precision
Hirota calculation now independently reproduces the chi-path number, both
sides of the endpoint, the distinct d-path limit, and the exceptional a4
limit. The original referee's lost-accuracy ODE remains insufficient evidence;
the explicit new construction and convergence data are given below.

For the minus a_n amplitude, continuing each imaginary block gives
`-(1)(h-B/2)(-1+B/2)=-(1)(h-1)K1`. The bootstrap generates the factor
`-(a)(h-a)` for every particle. Replacing B by 2-B yields the CDRS conjectured
C=-1 family, not a proof of the classical identification of the original minus
coideal.

### B5 and O.4: numerical scope

The mirror rule follows by reflecting the entire double row, including the
particle label at the left wall. At new rapidities the independently built
matrices commute, preserve the Krein form, and commute with Z tensor ... tensor Z
to 1.1e-13. The wrong convention fails only when the left wall is non-diagonal;
for U it is the same identity matrix. Bisection with theta=(12,.001) gives VV
thresholds .698194 (UN) and .536227 (NN). Testing only this rapidity limit in
V Vbar is insufficient; the additional regression explicitly separates the
finite slow-rapidity instability from the near-zero control.

For O.4 the independent Gamma-product sums were differentiated at theta=.83,
omega=40.37,80.37,160.37. The derivatives divided by omega for plus, minus and
half the bulk phase derivative are respectively

| omega | plus | minus | half bulk |
|---|---|---|---|
| 40.37 | -0.311522015 | -0.363190154 | -0.388859198 |
| 80.37 | -0.346000148 | -0.371911889 | -0.384826128 |
| 160.37 | -0.363337560 | -0.376312913 | -0.382790135 |

First-order extrapolation leaves differences 7.9e-5 and 4.0e-5. This supports
the leading semiclassical time delay, not an exact identity of phases. The
finite CDD blocks have O(1) phase derivatives and do not change the leading
O(omega) term.

## Primary-source checks

- Delius writes “The boundary condition corresponding to (1.4) is” before (1.5). [Delius, hep-th/9809140](https://arxiv.org/abs/hep-th/9809140): (1.4), (1.5), (5.1), (5.2); the factor of two and the allowed diagonal K were checked directly.
- Corrigan writes “Dirichlet boundary conditions are not generally permitted”. [Corrigan, hep-th/9612138](https://arxiv.org/abs/hep-th/9612138): section 4 limits the Dirichlet assertion to the integrability framework being discussed. BCDR/CDRS use boundary potentials, not a separate general Dirichlet exclusion theorem.
- [Kim–Cho, hep-th/9505138](https://arxiv.org/abs/hep-th/9505138): (29), (40), and the self-dual exact products; beta^2=2 pi B+O(B^2).
- [Kim, hep-th/9506031](https://arxiv.org/abs/hep-th/9506031): the a3 appendix (33), (39), and the d_n vector exponent rule below Figure 7. Higher-rank use is a conjectural rule, not additional perturbation theory.
- CAF say the amplitude “tends to a theory with fixed boundary conditions”. [Castro-Alvaredo–Fring, hep-th/0308004](https://arxiv.org/abs/hep-th/0308004): (5.1)–(5.3) explicitly give the classical product and discuss fixed versus Neumann limits; the shifted vector self-dual solution agrees with Kim's rule.
- [Gandenberger, hep-th/9806003](https://arxiv.org/abs/hep-th/9806003): (5.17) and (6.4), translated block by block.
- [CDRS, hep-th/9404108](https://arxiv.org/abs/hep-th/9404108): (5.7), (5.16), and the tentative status of the bootstrap picture. No new claim of exhaustive originality is made here.

## A5: source attributions and explicit comparisons

**(a) Correct.** Kim–Cho (hep-th/9505138), equations (29) and (40), are the
$d_4$ one-loop sources. Kim (hep-th/9506031) refers back to that calculation;
his appendix calculates $a_3$. With $\beta^2=2\pi B+O(B^2)$, the derivative of
the proposed Neumann amplitude agrees with those coefficients. For example,
at $\theta=.37+.21i$, the supplementary rerun gives
$\partial_B\log K_1=-.595711473453-1.12736942855i$ and
$\partial_B\log K_2=-.562297025118-2.01707295533i$ from either expression.

**(b) Correct.** Translate Kim's block $[x]$ to CAF's $\|2x\|$. His vector
rule becomes
$\|h-1\|\prod_{k=1}^{h-1}\|2k-1\|=\widetilde R_1(\theta+i\pi)$.
This is the same shift stated by CAF's equation labelled (KimCF).
Kim's phrase “only the middle number is two” identifies the extra central
factor. Figure 7 ends at $d_6$; the $d_7$ comparison applies the stated rule.

**(c) Correct, with that rank qualification.** The self-dual products and
$K^N$ have equal constant and linear terms but unequal quadratic terms.
For $d_4$ at the preceding rapidity, the coefficients of $B^2$ in $\log K$
are $-.0045981527+.02655555i$ for $K^N$ and
$.055469185+.5401141i$ for the self-dual product. Thus one-loop agreement
does not imply exact agreement. Kim–Cho explicitly write: “The exact boundary
reflection matrix is determined uniquely if we assume the strong-weak coupling
'duality'.” The correction concerns that uniqueness claim.

**(d) Correct.** Translate Gandenberger's sine block in (5.17) with the same
$T$ and $\mu$ conventions; its three factors give the $a_2$ minus breather.
Continuing $\omega=-2/B$ and bootstrapping gives
$K_a^{(-)}=-(a)(h-a)K_a$. Its dual is the general $a_n$ family in CDRS (5.16);
CDRS call the earlier $a_2$ choice a “minimal hypothesis” and their final
picture “tentative”. Describe this as a conjectured identification, not a
proof of the boundary condition.

## Numerical qualifications and replacement wording

The file locations below refer to the reviewed commit; named statement IDs
are included because subsequent line numbers have changed.

| Location | Replacement or qualification |
|---|---|
| `part7-half-line/30-boundary-conditions.qmd:185`, `prp-no-dirichlet` and its lead-in; ch. 35 open problem 7 | “For rank at least two, full Dirichlet is excluded in this Sklyanin construction with a field-independent K in finite-dimensional evaluation representations without trivial summands.” Say that the boundary parameters cannot **all** be sent to infinity. Do not claim exclusion of every notion of integrability. |
| Ch. 30 line 272, “Half-solitons at the wall”; ch. 31 line 91, `sec-images-preserving`; `tbl-boundary-types`; Appendix D's Delius entry | “When its class contains the topological charge of a regular single soliton, the boundary admits the corresponding static half-soliton and three-soliton image construction.” This is not a classification of all static solutions. |
| Ch. 35 line 304, `prp-dn-uniform-vacuum` and open problem 2 | State the H1 function space, coercivity and the direct-method existence argument before uniqueness. The remaining higher-rank problem is finding closed forms. |
| Ch. 33 line 257, `prp-dn-nondiagonal-scalar`; ch. 35, `eq-dn-neumann-amplitude` | “One loop requires a CDD correction but does not determine Xi uniquely. The displayed choice gives exact duality of the continued vector breather amplitude.” The alternatives were checked against scalar and vector-particle constraints; do not say that every spinor or full soliton bootstrap constraint was checked. |
| Ch. 35 line 320 onwards, discussion following `prp-dn-uniform-classical` | “Particle 3 of d5 also has a simple pole at theta=3i pi/8, with momentum i(1+sqrt(2))m and frequency m. The spinor frequency coincides, but the spinors decouple in both fluctuation Hessians; this coincidence does not establish an inheritance mechanism.” |
| Ch. 35 line 339, `sec-boundary-krein` | “For VV the onset is seen with the slow rapidity tending to zero. For V Vbar the lower threshold occurs at finite slow rapidity.” Round the finer mixed estimates to 0.667 (UN) and 0.500 (NN). State the failure of the alternative mirror convention only for a non-diagonal left wall. |
| Ch. 35, shifted Kim comparison | “The a3 minus amplitude agrees at one loop with Kim's (33) after theta is shifted by i pi.” |
| Ch. 35, d4 semiclassical comparison | “The calculations support equality of the leading semiclassical phase derivatives.” Do not assert equality of the full phases. |

**B1.** The original matrix checks cover $d_4,d_5$; scalar checks cover
$d_4,d_5,d_6$, with continued vector-particle closure through $d_7$.
The rerun of the Gamma-product sums gives worst residual $8.76\times10^{-30}$;
continued identities give residuals of order $10^{-29}$. These are
supplementary results. The new block calculation independently establishes
the one-loop ambiguity: at 45 digits its derivative residual is
$4.56\times10^{-55}$, vector-particle residual $3.12\times10^{-45}$, and
finite-$B$ differences range from .00369 to .226 at the selected points.
The comparison of $\rho_+$ requires the extra $(-T/2)_I^2$ factor derived above.

**B2.** Supplementary multichannel scattering reruns reproduce the classical
product to $2.3\times10^{-11}$ ($d_4$), $2.7\times10^{-11}$ ($d_5$), and
$1.8\times10^{-11}$ ($d_6$, separate run). The earlier raw-column bound-state
runner failed to finish. This has now been replaced by an independently
written matrix Riccati scan, which completes and passes for $d_4,d_5$.
The new calculation tests the eigenvalues of the boundary mismatch, rather
than reproducing a normalization-dependent determinant. Its domain, precision,
frequency coverage and output are recorded below.

**B3.** The matrix contraction and singlet-order derivation above give the
two crossing factors and the phase $e^{-2\pi i l/h}$. Supplementary checks
at $(n,l)=(2,1),(3,2),(3,1)$ give scalar residuals below
$3\times10^{-25}$; the matrix/breather errors are about $4\times10^{-15}$.
This is not a check of the full soliton bootstrap or the rank-one GZ
identification.

**B4.** The original report overstates what its ODE proves. Its historical rerun reports
**FAIL**: diagonal amplitudes agree to about $5\times10^{-12}$, but
spurious off-diagonal entries reach $1.2\times10^{-8}$ for $a_4$ and
$4.4\times10^{-8}$ for $a_5$, above its $10^{-8}$ tolerance. That rerun supported
the regular-background formula at lower precision but did not meet the
requested precision for every matrix entry. The independent calculation below
now resolves this issue and reproduces the formula beyond the claimed 27 digits.

The middle-soliton endpoint errors decrease from $.005$ to $.0005$ as
$1-d$ decreases from $10^{-2}$ to $10^{-3}$, then encounter an ODE floor.
For the stationary-pair endpoint, the d-path values of $|K_{15}|$ are
$.927785,.936485,.936746$. The chi-path integration loses accuracy: at
$\delta=10^{-3}$ one value is $5.062077$, and at $10^{-4}$ it is $72.627245$.
Comparing these with finite-$\delta$ values of another path cannot prove
unequal limits. The new independent high-precision calculation now establishes
the missing numerical convergence: the chi-path value is $.719193205484$,
the d-path value tends to approximately $.937437212$, and the exceptional
$a_4$ value is approximately $.147136511$. The book's check script is no
longer the sole source of these results. Keep the withdrawal of the proposed
correction to Bowcock–Perkins Table 2, whose singular-wall prescription is
different from the limits computed here.

**B5.** The additional weight-preserving, row-equilibrated intertwiner solve
reproduces the mirror, Krein and cyclic-sector identities with worst residual
$1.16\times10^{-13}$. VV bisections give .698194 (UN) and .536227 (NN).
The supplementary mixed scan gives deviations .049 at $(20,1.5)$,
$\omega=.672$ (UN), and .062 at $(20,2)$, $\omega=.505$ (NN), while the
near-zero slow-rapidity controls remain unbroken. These are finite-grid
threshold estimates, not analytic transition values. Three-soliton breaking
is supported by the 45-triple supplementary scan; it does not mean every
rapidity triple is broken. An initial dense solve falsely suggested a
$2.26\times10^{-8}$ departure on a uniform-wall grid; the restricted solve
reduces the worst mixed-grid deviation to $1.03\times10^{-11}$.

**B6.** The block continuation and bootstrap derivation above confirms the
minus-family factor. The supplementary scalar/breather checks through $n=4$
and particle bootstrap through $n=7$ pass at about $10^{-30}$. The physical
identification after duality remains conditional on the continuation conjecture.

## Completion of the outstanding numerical checks

### B2: stable subthreshold spectral scan

The fluctuation equations and Robin boundary condition are
$F''=(M^2(x)-zI)F$, $F'(0)+RF(0)=0$, where $z=\Omega^2$ and
$R=\operatorname{Hess}\mathcal B$. For the matrix of solutions decaying at
$-\infty$, set $Y=F'F^{-1}$. Direct differentiation gives

$$
Y'=M^2(x)-zI-Y^2,\qquad
Y(-L)=U\operatorname{diag}\sqrt{m_a^2-z}\,U^T.
$$

Here $U$ diagonalizes the asymptotic mass matrix. A bound state requires a
zero eigenvalue of $Y(0)+R$, provided this chart of the decaying subspace remains
regular. The integration succeeds without a chart singularity at every sampled
frequency; symmetry is checked before taking the eigenvalues. This removes the
exponential column scaling that defeated the earlier determinant calculation.

`b2_completion.py` reconstructs both Hessians directly from the roots and the
verified closed-form $d_4,d_5$ vacua. Each scan has 101 evenly spaced values of
$z$ from $-4$ to $m_{\min}^2-.001$, nine threshold-focused values
$m_{\min}^2-10^{-k}$ ($k=2,\ldots,10$), and $z=1$, with duplicates removed:
110 distinct points. All eigenvalues of $Y(0)+R$ remain positive.

| Algebra | $m_{\min}^2/m^2$ | Smallest sampled eigenvalue of $(Y+R)/m$ | Value at $\Omega^2=m^2$ |
|---|---:|---:|---:|
| $d_4$ | 2 | 0.848538937418 | 1.891805812446 |
| $d_5$ | 1.171572875254 | 0.463269881452 | 0.897061476243 |

The smallest values occur at the point $10^{-10}m^2$ below threshold.
Changing $L$ from 20 to 30 and relative tolerance from $2\times10^{-10}$ to
$2\times10^{-12}$ changes the eigenvalues by at most $3.52\times10^{-11}$
($d_4$) and $2.59\times10^{-11}$ ($d_5$), at the pole frequency and three
near-threshold points. The positive quadratic fluctuation form excludes
nonpositive eigenvalues analytically. The positive-frequency result remains a
numerical scan, not a theorem covering the continuous interval. No $d_6$
absence claim is added.

### B4: independent high-precision fluctuation construction

`b4_completion.py` imports neither the book's nor the original referee's
calculation. Starting from `eq-n-soliton` and `eq-interaction-coefficient`,
write each term of a static tau function as $c_{j,S}e^{r_Sx}$, with
$r_S=\sum_{a\in S}m_a$. Adding one infinitesimal soliton of particle species
$c$ and spatial exponent $k$ gives

$$
\delta\tau_j(x)=e^{kx+2\pi ijc/h}
\sum_S c_{j,S}e^{r_Sx}
\prod_{a\in S}
\frac{k/m_c-\cos(\pi(a-c)/h)}{k/m_c-\cos(\pi(a+c)/h)},
\qquad k^2=m_c^2-\Omega^2.
$$

The fluctuation is $\epsilon=-\sum_j\alpha_j\delta\tau_j/\tau_j$.
Derivatives are calculated from the finite exponential sums and the quotient
rule, without numerical differentiation. Normalizing its asymptotic Fourier
vector fixes the particle amplitudes. Incoming/outgoing open channels use
$k=\pm i\sqrt{\Omega^2-m_c^2}$; closed channels use the positive real root.
The matrix of their boundary residuals $\epsilon'+R\epsilon$ determines the
outgoing coefficients and the closed-channel coefficients by a linear solve.

At 80-digit arithmetic, two interior phases in **every** phase interval of
all pairs in $a_2$–$a_5$, plus the self-conjugate backgrounds, give 56 regular
backgrounds and three frequencies each. The worst residuals are
$1.52\times10^{-79}$ for the static equation,
$2.57\times10^{-79}$ for the tested fluctuation equations,
$1.60\times10^{-78}$ for the boundary condition, and
$9.33\times10^{-79}$ for agreement of the full open-channel reflection matrix
with the diagonal DG amplitude. Thus the original ODE tolerance failure is
resolved by an independent calculation.

For the $a_5$ target $(+,+,-,-,-,-)$ at $\theta=.6$, the limiting calculations
are:

| Distance $\delta$ | $|K_{15}|$, $\chi=5\pi/3+\delta$ | $|K_{15}|$, $\chi=5\pi/3-\delta$ | $|K_{15}|$, $d=(1-\delta)/(1+s)$ |
|---|---:|---:|---:|
| $10^{-3}$ | 0.719192521391 | 0.719192521391 | 0.936483685827 |
| $10^{-5}$ | 0.719193205416 | 0.719193205416 | 0.937427690600 |
| $10^{-7}$ | 0.719193205484 | 0.719193205484 | 0.937437116694 |
| $10^{-9}$ | 0.719193205484 | 0.719193205484 | 0.937437210954 |

The zeroth-order boundary defect tends to zero along all three approaches:
its final values are $1.15\times10^{-9}$, $2.31\times10^{-9}$ and
$4.81\times10^{-5}$. The d-path has the slower square-root convergence.
For the alternative target $(+,+,-,-,+,-)$, the two phase sides instead give
$.091302814886$ and $.272069513902$ at $\delta=10^{-9}$. Repeating the first
target's phase and d-path calculations at 120 digits changes its amplitudes
by at most $8.87\times10^{-26}$ relative to 80 digits. Fluctuation equations
for all species also pass inside the wall layer. These convergent sequences,
not the original finite-distance ODE observations, verify the differing limits.

For middle-soliton endpoints, both the single-soliton family and a family with
an added species-1/conjugate pair of amplitude $(1-d)^2$ approach the diagonal
DG value. At $1-d=10^{-9}$, the maximum of the tested diagonal and off-diagonal
errors is $4.97\times10^{-10}$ ($a_3$) and $1.42\times10^{-9}$ ($a_5$),
scaling linearly with $1-d$. This verifies the **tested** paths; it is not
proof of independence under arbitrary perturbations of the background.

For exceptional $a_4$, take the static species $(1,4,2,3)$ with phases
$(9\pi/5,-9\pi/5,18\pi/5,-18\pi/5)$ and pair amplitudes $d_1,d_2$.
The symmetry makes $\tau_2(0)=\tau_4(0)$. Solve $\tau_2(0)=0$ and
$\tau'_0(0)/\tau_0(0)=\tau'_3(0)/\tau_3(0)$, the equality required by the
boundary equation on the nonvanishing nodes. This independently gives

$$
d_1=.257789025012387020138024812130,\qquad
d_2=.563669884984393415822712186658.
$$

The boundary defect evaluated at $x=-10^{-30}$ is $9.46\times10^{-16}$,
consistent with its square-root approach. Reducing $d_1$ or reducing $d_2$
by the relative fraction $\delta$ gives:

| $\delta$ | $|K_{14}|$, reduce $d_1$ | $|K_{14}|$, reduce $d_2$ |
|---|---:|---:|
| $10^{-4}$ | 0.147126544275 | 0.146877326806 |
| $10^{-6}$ | 0.147136412477 | 0.147133965847 |
| $10^{-8}$ | 0.147136509883 | 0.147136485462 |
| $10^{-10}$ | 0.147136510856 | 0.147136510612 |

The same non-diagonal limit is therefore independently reproduced for these
two paths. These are approaches through static bulk configurations; they need
not satisfy the target boundary condition at finite distance from the endpoint.
The calculation specifies a limiting Hessian prescription. It does not show
that a singular background has a unique reflection problem, and it does not
test Bowcock–Perkins's distinct prescription of setting the perturbation of
a vanishing tau function to zero at the wall.

**Additional replacement wording** for ch. 31, the bullets in
`sec-boundary-missing-charges`: replace an unrestricted “the answer does not
depend on the approach” by “the approaches tested here have the same limit”.
For the exceptional background retain the existing qualification “the limits
we tried agree”. The stationary-pair path dependence and the quoted numerical
values can now be retained with independent verification.

## Integration

No integration into the book or PR is requested with this updated report.
The table above is the editor's disposition record: **retain** the qualified
A1–A5 and B1/B3/B5/B6 findings; **modify** B2 to keep the extra pole while
rejecting its proposed spinor mechanism; **confirm with numerical scope**
the completed B2 bound-state scans and B4 limits. Retain the criticism of the
original B4 ODE argument, but the previously missing independent computations
are now supplied. Qualify claims of approach independence to the paths tested. Retain
O.1–O.3, narrow O.4 to leading phase derivatives, and treat O.5 as a historical
build result. Subsequent commits on the PR branch are outside this review's
claimed coverage.

## Completed checks and reproducibility

Evidence is preserved under `referee_scratch/halfline_updated/`:

- `b2_completion.py`: independent matrix Riccati scans and domain/tolerance refinement.
- `b4_completion.py`: independent Hirota subset sums, differential-equation checks,
  regular reflection, three kinds of singular limits, and precision refinement.
- `checks.py`: this review's exact charge enumeration, independent CDD block
  calculation, d5 pole count/residue, decoupling test and Riccati check.
- `halfline_consolidation_integration_checks.py`: the additional Krein and
  semiclassical checks; it uses the original referee's representation and
  Gamma-product helpers, with the independent stable intertwiner solver.
- `snapshot/`: read-only historical copies of the helper sources, avoiding
  dependence on the current main checkout's renamed files.
- `logs/`: the completed numerical outputs, including failures;
  `book-halfline_consolidation_checks.py.log`, `integration-krein.log` and
  `integration-semiclassical.log` contain the new checks' results.

The 20-script supplementary referee runner completed 17 scripts successfully,
one script with a printed FAIL (B4), one killed process, and one unsuccessful
helper invocation. Nine book check scripts passed, including the new regression;
this is regression evidence, not nine independent proofs. The earlier HTML
render completed with an existing TeX-to-MathML warning. No render was run to
prepare this updated report, and no claim is made that all referee checks passed.

The older test logs predate the checkout's move to main. The new completion
logs are `b2_completion.log`, `b4_regular.log`, `b4_pairs.log`, `b4_middle.log`,
`b4_exceptional.log`, and `b4_precision.log`. All six completion runs pass.
The earlier raw-script failures are retained as provenance; their outstanding
scientific checks are resolved by the new independent implementations.

The saved checks can be run from the repository root (the broad, slow suites
are not needed to reproduce the new counterexamples):

```bash
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/b2_completion.py
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/b4_completion.py regular
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/b4_completion.py pairs
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/b4_completion.py middle
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/b4_completion.py exceptional
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/b4_completion.py precision
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/checks.py
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/halfline_consolidation_integration_checks.py krein
OPENBLAS_NUM_THREADS=1 python3 referee_scratch/halfline_updated/halfline_consolidation_integration_checks.py semiclassical
```
