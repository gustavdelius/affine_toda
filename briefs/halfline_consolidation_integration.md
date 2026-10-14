> **Superseded (2026-10-14).** This file records edits that Codex made in commit 904f924 and that Claude had not verified; they were verified, and partly rewritten, in the `## Integration (Claude, 2026-10-14)` section of [halfline_consolidation_referee_report.md](halfline_consolidation_referee_report.md), which supersedes this record.

# PR #14: re-derivation and integration of the consolidation referee report

2026-10-14, Codex. Source: `halfline_consolidation_referee_report.md` (Opus,
2026-10-13). The referee's item numbers below are distinct from the earlier
consolidation workers' A1–A6 numbers. The report is retained as a historical
record; this file records which conclusions survive re-derivation.

The proof arguments were reconstructed from the Lax compatibility condition,
the energy functional and the block definitions. Numerical verification uses
the referee's independent implementations as well as the book implementations,
with additional independent checks in
`textbook_code/halfline_consolidation_checks.py` and
`briefs/halfline_consolidation_integration_checks.py`. Primary papers were read
in arXiv TeX or extracted PDF text. Numerical scans remain numerical evidence,
not proofs at arbitrary rank or coupling.

## Item-by-item disposition

| Referee item | Re-derived conclusion | Integration |
|---|---|---|
| A1 | The evaluation-representation obstruction is correct. The lead-in overstates its scope; individual non-simply-laced boundary parameters can be unbounded. The proof needs the invariant-image ideal argument. | Ch. 30: qualify the scope, say “cannot all”, specify evaluation representations, supply the cycle and ideal arguments. Ch. 35 problem 7: same scope. |
| A2 | Half-coweight compatibility and the reflection-equation statement survive. The inference from a fundamental weight to a single-soliton charge does not. Missing-class counts are 0, 2, 5, 21 for ranks 2–5. | Qualify the half-soliton and three-soliton constructions in chs. 30–31 and `tbl-boundary-types`; explicitly record the missing classes. |
| A3 | Delius (1.4) gives a half-coweight value; (1.5) follows only for even coweights. A vacuum shift changes the image coweight by twice a coweight. | Retain the correction in Appendix D, with A2's qualification. |
| A4 | Strict convexity gives uniqueness. Coercivity and the direct method also give existence for every rank when all C_j are nonnegative. The d4/d5 kinks and energies agree. | Add the function space and existence argument; problem 2 now asks for closed forms of higher-rank uniform vacua. |
| A5(a–d) | Kim–Cho is the d4 one-loop source; Kim's vector rule and CAF's shifted self-dual solution coincide; Gandenberger's minus breather and CDRS's conjectures have the stated conventions. Kim's figure ends at d6. | Retain the attributions; specify that d7 uses Kim's stated rule. |
| B1(a–f) | Matrix unitarity, crossing, singlet coefficients, scalar equations, breather identities, continuation, and one-loop comparisons survive at the tested ranks and points. | Keep the formulas and numerical labels. |
| B1(g) | One loop does not determine Xi. The three alternatives preserve the tested scalar and vector-particle constraints and differ at finite coupling. The rho-plus comparison also needs its extra CDD factor. | Correct chs. 33 and 35 and problems 2–3. Make boundary identifications conditional on continuation and the chosen CDD factor; give the precise rho ratios. |
| B2 | The diagonal classical reflection calculation survives. Particle 3 of d5 has the extra simple pole at 3i pi/8. **The referee's proposed inheritance from a coupled spinor channel is false:** both Hessians decouple the spinors from particle 3. | Correct the pole list, state the coincident frequency without asserting that mechanism, and keep the non-bound-state check. Limit the bound-state scan claim to d4/d5; neither original implementation scanned d6. |
| B3 | The SP crossing functions, scalar equations and phased breather coefficient survive. | Keep the proposition and amplitude. |
| B4 | Regular-background reflection and the middle-soliton limit survive. Singular pair limits require the path and target boundary Hessian to be specified. The referee's ODE does not independently resolve the chi-path limit. | Keep the qualified numerical statements and the withdrawal of “BP Table 2 corrected”. Do not promote finite-delta ODE data to a proof of distinct limiting values. |
| B5 | The mirror convention, cyclic sectors and numerical thresholds survive. The slow-rapidity-zero description is wrong for V Vbar. The unmirrored convention coincides with the mirror convention for a uniform left wall. | Correct the onset and mirror descriptions. Round the finer V Vbar estimates to 0.667 and 0.500. Add a regression exhibiting breaking at finite slow rapidity while the near-zero point remains unbroken. |
| B6 | The minus scalar, breather, continuation, bootstrap and CDRS dual comparison survive. | Keep the formulas. |
| O.1 | CAF explicitly gave the B=0 product and associated it with fixed boundary conditions. | Credit CAF for that formula; distinguish the numerical identification with the uniform kink. |
| O.2 | The shifted a3 Kim agreement is perturbative, not an exact equality. | Say “agrees at one loop with Kim's (33)”. |
| O.3 | The a3 alternative P_(omega+1)/2 survives the soliton bootstrap and one loop. | Correct the CDD discussion and its dependent ch. 35 statements. |
| O.4 | The report did not check semiclassics. An additional Gamma-product calculation supports equality of the **leading phase derivatives**, not equality of the finite-coupling phases. | State the observable actually checked: leading semiclassical phase derivatives. |
| O.5 | Cross-references and citations require a fresh render after these edits. | Full HTML render and source/reference checks recorded below. |

## Re-derivations

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
self-dual product at B=.6 does not. The book now states only that perturbative
agreement. These distinctions also apply to conditional boundary identifications.

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
spinor. The book records the equal frequency but leaves the detailed pole
mechanism unassigned. Independently evolving the Riccati matrix of decaying
solutions at Omega^2=1 gives minimum eigenvalue 0.897061476 for Y(0)+R,
excluding a bound-state zero there. This also avoids exponentially ill-scaled
raw determinants. Absence on an entire spectral interval remains a numerical
scan, with d6 explicitly outside its verified scope.

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
Hirota calculation is the evidence for the reported chi-path number; the
referee's lost-accuracy ODE is not additional confirmation of that number.

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
V Vbar is insufficient; the new textbook regression explicitly separates the
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

- [Delius, hep-th/9809140](https://arxiv.org/abs/hep-th/9809140): (1.4), (1.5), (5.1), (5.2); the factor of two and the allowed diagonal K were checked directly.
- [Corrigan, hep-th/9612138](https://arxiv.org/abs/hep-th/9612138): section 4 limits the Dirichlet assertion to the integrability framework being discussed. BCDR/CDRS use boundary potentials, not a separate general Dirichlet exclusion theorem.
- [Kim–Cho, hep-th/9505138](https://arxiv.org/abs/hep-th/9505138): (29), (40), and the self-dual exact products; beta^2=2 pi B+O(B^2).
- [Kim, hep-th/9506031](https://arxiv.org/abs/hep-th/9506031): the a3 appendix (33), (39), and the d_n vector exponent rule below Figure 7. Higher-rank use is a conjectural rule, not additional perturbation theory.
- [Castro-Alvaredo–Fring, hep-th/0308004](https://arxiv.org/abs/hep-th/0308004): (5.1)–(5.3) explicitly give the classical product and discuss fixed versus Neumann limits; the shifted vector self-dual solution agrees with Kim's rule.
- [Gandenberger, hep-th/9806003](https://arxiv.org/abs/hep-th/9806003): (5.17) and (6.4), translated block by block.
- [CDRS, hep-th/9404108](https://arxiv.org/abs/hep-th/9404108): (5.7), (5.16), and the tentative status of the bootstrap picture. No new claim of exhaustive originality is made here.

## Validation record

Validation is in progress. Final commands, results and numerical limitations
will replace this paragraph before the changes are committed to PR #14.
