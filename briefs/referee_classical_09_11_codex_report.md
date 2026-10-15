# Codex referee report: classical chapters 9–11, phase 2

Date: 2026-10-15. Model: GPT-6 (Codex; no more specific identifier exposed). Reviewed commit: `bcd82b9ee9667975eed02109f08a257d50ba38b1`, worktree `.claude/worktrees/referee-briefs-part-iii`. Line numbers throughout refer to this version, not the later main checkout.

This continues [the frozen Codex phase-1 report](referee_classical_09_11_codex_phase1_report.md), abbreviated **P1** below. Its SHA-256 is `26c44d6352239cde46bb1ae60c9f66567a41fe7b79b752a4993b68bdc4faa24b`. The coordinator list's SHA-256 is `7e46e33bf78567a291bb5350625a8a67b9139c75665e02d0a9de3252de1f41b6`. The other referee's two reports were not read in any checkout. Neither earlier Claude reports nor author scripts supply evidence. No book file, git index or branch was changed; no rendering or subagents were used. All new filenames contain `codex`.

**Overall:** retain the established classical constructions with the corrections below. Do not retain the Appendix E theorem/originality claim or the assertion that H6 is closed. The new direct integrations support spin-2 and spin-4 additivity on the regular A2 breathers tested; they do not establish unrestricted additivity. The requested centralizer and A4/E6 coupling checks pass. The comparison with the general quantum soliton phase remains **unclear in the book's convention**, with an explicit sign discrepancy documented below, rather than an asserted correction to Hollowood.

Evidence scripts and outputs are in `referee_scratch/classical_09_11/codex/` and archived in [the Codex evidence directory](referee_classical_09_11_codex_scripts/). P1 sections already supply complete computations for many subitems; those are referenced rather than repeated. A finite numerical check is never presented as a universal proof.

## 1. Theta selection rule, Appendix E, H6 and prior art

**Overlap:** P1 ranks/entries 3 and 18, both checked, with 1, 2 and 4 supplying regularity counterexamples. **Verdict:** the universal rule and H6 closure are **wrong**; the restricted implication is **correct**; exact novelty of the reduced-data wording is **unclear**.

P1 §2.3 gives the conditional proof and a topological-charge counterexample. The coordinator's requested proof of the two missing premises cannot be supplied as stated: they are not consequences of additivity and one of them is false.

Let `D={(a_i,theta_i)}`. For a *specified* hierarchy charge assume (i) a convergent charge on the solution domain, (ii) `Q[phi]=F_Q(D)`, independent of omitted phases/positions, and (iii) `Q[Theta phi]=epsilon_Q conjugate(Q[phi])`. Then

\[
D=\Theta D\ \Longrightarrow\ Q[\phi]=F_Q(D)=F_Q(\Theta D)
=Q[\Theta\phi]=\epsilon_Q\overline{Q[\phi]}.
\]

Additivity gives premise (ii) only if the one-constituent charges are themselves phase independent. On a regular single soliton, differentiation with respect to its complex translation parameter is `(1/k) partial_x` on any local density, where `k=m_a cosh theta`. Its integrated derivative vanishes when the density decays and no pole crosses the integration line. This proves phase independence within such a regular domain, not across singular phase boundaries. For scattering solutions with distinct real velocities, separation and locality give additivity if mixed-density contributions tend to zero integrably. Extending this to bound states requires analytic continuation in a domain with no poles and uniform control of the integrals. These hypotheses and that extension are not proved in the chapter.

**Parity tests.** The printed sinh-/sine-Gordon current has even Theta parity: each derivative changes sign and is squared, and `cosh(beta phi)` or `cos(beta phi)` is even. For A2, put `U=i beta phi=(u,v,-u-v)`, `m=beta=1`, `D=partial_t-partial_x`, `u_k=D^k u`, `v_k=D^k v`. Theta acts as `U -> conjugate(U)`. Define

\[
z_0=e^{u-v},\quad z_1=e^{u+2v},\quad z_2=e^{-2u-v},\qquad
\partial_+D(u,v)=(-z_0+z_2,z_0-z_1).
\]

`codex_phase2_currents_v3.py` solves the polynomial conservation equations from these field equations, rather than importing a current. For spin 2 it finds

\[
T=-u_1^2v_1-u_1v_1^2-u_2v_1-2v_1v_2+v_3,\qquad X=3z_0v_1,
\quad \partial_+T+DX=0.
\]

For spin 4 its explicit 18-term `T` is recorded in `codex_phase2_currents_v3.json`, with

\[
X=3z_0(u_1v_2+v_1^3-v_1v_2),\qquad \partial_+T+DX=0.
\]

The exact rational nullspace dimensions are 6 and 21; quotienting by the constructed total derivatives leaves one primitive current at each spin (derivative ranks 5 and 20). The script verifies both differential identities exactly. Both representatives have real coefficients in `U` and are Theta-even. Their species-1 rest charges are respectively `3 sqrt(3) i` and `-9 sqrt(3) i/2`. Dividing by these imaginary constants makes the one-soliton eigenvalues real, and changes both charges to Theta-odd. Thus the book's chosen spin-2 normalization is consistent.

But `(1+i)` times either raw current is still a local conserved current and has neither parity `+1` nor `-1`: the proposed parity would require `1+i=epsilon(1-i)`. An antilinear involution permits a real basis; it does not make every vector of the complex current space a parity eigenvector. P1's real sine-Gordon kink also disproves the phrase “every conserved charge” when phase-dependent topological charges are included.

**Khastgir–Sasaki example.** Their [hep-th/9507001](https://arxiv.org/abs/hep-th/9507001), eqs. (3.10)–(3.12), uses `sigma x-lambda t+x0`. Choose the book's A2 species 1, `theta=.3+.2i`, `xi=.4i`; then `sigma=sqrt(3)cosh(theta)`, `lambda=sqrt(3)sinh(theta)`, `x0=xi`, and `sigma^2-lambda^2=3`. These are their c-soliton data (with the field-sign convention translated). The slice `t=0` is regular: the modulus condition forces `x=0`, where none of the three phases is an odd multiple of pi. Independent canonical integration at `t=.1` gives

`E=10.646930242183003067385 + 0.628722418343127100433 i`.

The one-element reduced multiset maps to `{(2,.3-.2i)}`, so it is not Theta invariant. Solving the two real equations for a tau zero gives the first positive singular time `t=4.3963606744302846`, at `x=1.2807153133477245`, node 2; residual `1.23e-16` (`codex_phase2_extra.py`). The source's §5 says “all c-soliton solutions develop singularities after a finite time.” This example verifies the claim about that family, but calling it excluded from *the physical sector* is an interpretation, not a classical theorem.

**H6.** At `part9-foundations/44-twelve-gaps.qmd:33`, H6 asks which solutions should count as physical, not just which energies are real. A sufficient condition on selected integrated charges supplies neither necessity, stability, regularity nor a quantum-state selection. The “closed” table entry at line 81 and chapter 11:49 overclaim.

**Prior art read and checked.** The relevant constructions, energy arguments and conclusions of all four requested papers were examined. None proves the full general affine-Toda higher-charge statement used by the book.

- **Cen–Fring, [1602.05465](https://arxiv.org/abs/1602.05465), §§2–3.** Complex KdV, mKdV and sine-Gordon solutions; real integrated energies are obtained by integration and by spacetime shifts restoring antilinear symmetry. Equation (2.60) is a numerical two-sine-Gordon-soliton additivity check. The conclusion describes the general N-soliton extension as a conjecture, using the words “we conjecture that one might be able”. Their dispersion parameters are real. This already supplies symmetry-plus-phase-insensitivity reasoning, but not a proof for complex-rapidity Toda breathers. Their PT includes spacetime reflection, unlike the field-only Theta here.
- **Cen–Correa–Fring, [1608.01691](https://arxiv.org/abs/1608.01691), §3, eqs. (3.4), (3.10)–(3.18).** Gardner recursion gives KdV local charges. The stated conclusion is “Thus PT-symmetry and integrability guarantee the reality of all charges.” The argument uses even density parity and asymptotic separation into localized one-solitons, yielding `I_n=2 sum_k alpha_k^(2n-1)/(2n-1)`. The support-disjointness statements must be read as asymptotic limits, not exact compact support at finite time. This is prior art for the logic of the book's hierarchy-charge criterion, in a different model and parameter domain.
- **Fring–Taira, [2007.15425](https://arxiv.org/abs/2007.15425), §2, eqs. (2.7)–(2.9), and the three model applications.** They derive BPS equations, evaluate boundary energies, and use conjugation combined with spatial reflection. They explicitly allow symmetry to exchange equal-energy self-dual and anti-self-dual solutions: “the energy is the same for either choice.” Independent reasoning is `E[Theta phi]=conjugate(E[phi])` plus equality of the two energies, exactly the logical step relevant here. It does not require the individual field to be fixed. The examples include coupled/extended sine-Gordon models; no general affine-Toda hierarchy result follows.
- **Correa–Fring–Taira, [2110.06825](https://arxiv.org/abs/2110.06825), §2 and §§3.1–3.2.** The general criterion allows a map to “itself or possibly to a second degenerate solution”. Their Bullough–Dodd equation is related to rank-one twisted Toda by field/scaling conventions. Type-I energies are evaluated as `-6|k|`. Type-II energies require a complex vacuum shift and corresponding constant subtraction in the Hamiltonian; they are not finite energies relative to the unshifted vacuum. The stability analysis distinguishes the two families despite their real energies. This is particularly direct prior art for the distinction between fixed fields and symmetry-related degenerate solutions, and evidence that energy reality alone does not settle physical acceptability. The source's spacetime PT and energy normalization must not be identified without conversion with the book's Theta and imaginary-coupling energy.

The additional search found Cen–Fring's [Multicomplex solitons](https://www.atlantis-press.com/journals/jnmp/125950459/view), which explicitly uses antilinear symmetries as a selection principle for real conserved quantities in multicomplex KdV systems. This strengthens the prior-art objection to a broad “selection rule from symmetry” originality claim. It does not establish priority for the exact Toda species/rapidity formulation. No novelty inference is drawn from a negative search.

**Paste-ready replacements.** At `part3-classical/11-reality.qmd:45`, replace the selection-rule bullet by:

```markdown
::: {#prp-theta-reduced-data}
**[Established; reformulated, conditional]** Let $Q$ be a local hierarchy charge whose convergent value is determined by the reduced Hirota data $D=\{(a_i,\theta_i)\}$ alone. Choose its normalization so that $Q[\Theta\phi]=\epsilon_Q\overline{Q[\phi]}$, with $\epsilon_Q=\pm1$. If $\Theta D=D$, then $Q=\epsilon_Q\bar Q$. This applies to the energy and momentum in the regular decay sector of @prp-hirota-energy. For higher charges, additivity and independence of the internal phases are additional assumptions. Topological charges are not included.
:::

::: {.proof}
Write $Q[\phi]=F_Q(D)$. Invariance of $D$ gives $Q[\phi]=F_Q(D)=F_Q(\Theta D)=Q[\Theta\phi]=\epsilon_Q\overline{Q[\phi]}$.
:::
```

This label records a conditional reformulation, not a new theorem. The label should match an explicit statement id, not a section id. At 11:49 use P1 §2.3's replacement: a sufficient reality criterion, with physical selection and singularities still open. At 11:47 replace the last sentence by “Their complex energies are consistent with the failure of this sufficient reality condition; no physical-sector classification follows from this observation alone.” At 11:5 use “Derive a sufficient reality or imaginary-value criterion for charges of definite Theta parity.” At 11:6 use “Distinguish energy reality from regularity and physical-state selection.”

At `part9-foundations/44-twelve-gaps.qmd:81` replace the status cell by:

> partly addressed: sufficient reality conditions for regular classical charges; a general physical-solution criterion and singularities remain open

At `part3-classical/11-reality.qmd:96`, replace the novelty sentence by the following (bibliography keys shown are proposed new keys, not claims that entries already exist):

```markdown
Antilinear symmetries have also been used to explain real integrated energies of complex solitons [@cen_fring2016], and, together with asymptotic separation, real higher charges [@cen_correa_fring2017]. A symmetry may relate distinct solutions of equal energy rather than fix each field [@fring_taira2020; @correa_fring_taira2022]. Here these ideas are organized using the field-only symmetry $\Theta$ and reduced Hirota species/rapidity data. The higher-charge statement remains conditional on additivity, phase independence and regularity; we make no originality claim for the general symmetry argument.
```

At `appendices/E-original-results.qmd:8`, the preferred replacement is **deletion of this item from the list of original results**, retaining the reformulated statement in chapter 11. If the editor deliberately retains an organizational contribution in this appendix, replace it with:

```markdown
- **[Established; reformulated, conditional]** For a convergent local hierarchy charge of definite $\Theta$ parity whose value depends only on the reduced species/rapidity data, invariance of those data implies $Q=\epsilon_Q\bar Q$. This organizes earlier antilinear reality arguments; it does not establish higher-charge additivity or select a physical sector (@prp-theta-reduced-data).
```

P1 §§2.3–2.4 already supply the field/data distinction and breather-domain replacements at 11:26,28,54.

## 2. Energy, momentum and the stress tensor

**Overlap:** P1 entries/ranks 2 and 13, checked; entry 1 excludes singular integrations. **Verdict:** the formula is **correct on the stated regular decay sector**, **wrong without that domain**. The normalization of the tested tensor is **correct**. Its precise formula's attribution to the unavailable `olive1993a` remains **unclear**.

P1 §§2.2 and 2.13 give the regular antikink counterexample to unrestricted complex rapidity, boundary derivation, static first integral and trace. They fully cover those subitems, including the recommended restricted **[Established; reformulated]** label.

**Exact cited source.** [Takács–Watts, hep-th/9810006](https://arxiv.org/abs/hep-th/9810006), p.11, eq. (4.6), displays

\[
E=M\sum_i\cosh\theta_i,\qquad P=M\sum_i\sinh\theta_i.
\]

This is their A2 expression; their `M=8m/beta^2` uses their conventions. It is not a source for all algebras with unrestricted complex data. [Olive–Turok–Underwood, hep-th/9305160](https://arxiv.org/abs/hep-th/9305160), §8.3, p.27, says that continuation works “providing that its real part is positive”, referring to `exp(eta_i)`, and states independence from the constants `Q_i`. This is positive-real spatial decay, not arbitrary complex rapidity. The original [1993 energy–momentum paper](https://doi.org/10.1016/0550-3213(93)90318-J), whose abstract was accessible, establishes the total derivative and additive real-soliton contributions; Takács–Watts later uses the complex continuation. The full original source formula was not available, so I do not manufacture a quotation of its prefactor.

**Direct canonical integrals.** `codex_phase2_integrals_v2.py` evaluates derivatives of the finite tau polynomial analytically, constructs `phi_t`, `phi_x` and the exponential potential, then integrates

\[
T_{00}=\tfrac12(\phi_t^2+\phi_x^2)+V,\qquad
P=-\int\phi_t\cdot\phi_x\,dx.
\]

There is no complex conjugation in the kinetic terms and no replacement of the integrand by boundary data. Units are `m=beta=1`, so `M_1=M_2=6sqrt(3)`.

| Data, integration slice | Canonical E | Canonical P |
|---|---|---|
| species 1, `theta=0`, `xi=.3i`, `t=0` | 10.392304845413263761165 | 0 |
| species 1, `theta=.3+.2i`, `xi=.4i`, `t=.1` | 10.646930242183003067385 + 0.628722418343127100433 i | 3.101585063440180077815 + 2.158239607451275474235 i |
| species 1,2, `theta=.23 +/- i*pi/3`, both `xi=i*pi/2`, `t=.17` | 10.668395194383174381312 | 2.411359787244120294605 |
| same breather, `xi_1=.2+i*pi/2`, `xi_2=-.2+i*pi/2`, `t=.43` | same to shown digits | same to shown digits |

The unshifted breather is regular globally: in its rest frame `A=-2` and `Re tau_j=1+2r^2>0`. For the shifted data set `s=exp(.2)`. Directly, `Re tau_j >= 1+2r^2-(s-s^(-1))r`, whose minimum is `1-(s-s^(-1))^2/8>0`. Thus both tested breather families are globally regular; this does not assert regularity for arbitrary unequal phases.

At 45 digits on `[-90,90]`, the maximum energy/momentum residual against the additive prediction is `1.5e-45`. All three covariant stress components are checked independently at five points per case; maximum residual is below `2.4e-45`:

\[
T_{00}=2\sum_j\partial_x^2\log\tau_j,\quad
T_{01}=2\sum_j\partial_x\partial_t\log\tau_j,\quad
T_{11}=2\sum_j\partial_t^2\log\tau_j.
\]

These are precisely `(eta_munu Box-partial_mu partial_nu)C` for `C=-2 sum log tau` and signature `(+,-)`. Also `T^mu_mu=2V=Box C`; the trace equality follows immediately by multiplying the bilinear log equations by `-4/(beta^2 alpha_j^2)` and summing. For static solutions the trace and first integral suffice; for general moving multi-solitons the trace alone is not a proof of the full tensor identity.

**Locality and tau gauge.** The field determines coroot combinations of the logs, not the remaining redundant combination. In A type, `log tau_j -> log tau_j+f` leaves the field unchanged but changes `C` by `-2hf/beta^2`. If `Box f=0`, it also preserves the log-tau field equations. Thus `C` is not a gauge-independent local function of `phi`. The tensor identity must use the standard normalized soliton tau gauge; it cannot be asserted for arbitrary redundant tau representatives. For example `f=x^2+t^2` is harmonic and leaves the field unchanged but changes the proposed stress improvement. This does not contradict the identity on the normalized finite Hirota polynomials used in the book.

[Zhu–Caldi, hep-th/9404097](https://arxiv.org/abs/hep-th/9404097), eqs. (1.4) and (1.7), gives the redundant `ln tau_j-l_j ln tau_0` field parametrization and the static identity `1/2 (partial_x varphi)^2=U(varphi)`. I did not locate the full displayed tensor formula there. Aratyn et al., [hep-th/9212086](https://arxiv.org/abs/hep-th/9212086), eq. (2.22), independently displays the conformal-to-affine tensor improvement and derives the boundary mass. This supports the construction, but is not substituted for the exact OTU quote requested and unavailable.

**Replacement at `10-solitons.qmd:174`:** keep the printed sign/prefactor, add “in the standard normalized soliton tau gauge”, replace “and the energies are real” by “and the integrated charges follow from spatial boundary terms when the integrals converge; regular real-rapidity solitons have real energies.” Add “The trace check alone does not establish the full tensor identity.” Restrict `prp-hirota-energy` using P1 §2.2's paste-ready text.

## 3. Higher-charge additivity tested directly

**Overlap:** P1 entries 3 and 15 flagged the assumptions and Cartan conventions, but did not construct and integrate the A2 currents. **Verdict:** the normalized spin-2 example is **correct in the cases tested**; universal complex-rapidity additivity remains **unclear**.

The exact construction in §1 and canonical-field evaluation in §2 are independent evidence. With `Q_s=integral(T+X) dx`, direct integration gives

| case | raw Q2 | raw Q4 |
|---|---|---|
| species 1 at rest | `+3sqrt(3)i` | `-9sqrt(3)i/2` |
| species 2 at rest | `-3sqrt(3)i` | `+9sqrt(3)i/2` |
| complex single soliton | `-3.687015598347877987820 + 8.720611956246317320633i` | `18.563561859336835032209 - 18.029202317866202495819i` |
| breather, both choices of xi in §2 | `-14.256665864950335973763` | `-33.875420264140018561644` |

Normalize each charge by its species-1 rest value `q_1`. On the breather, the integrated results are

\[
Q_2/q_1=2.7436966249585417700284\,i
=2i e^{.46}\sin(2\pi/3),
\]
\[
Q_4/q_1=-4.3462184463139873857525\,i
=2i e^{.92}\sin(4\pi/3).
\]

Largest spin-2 error: `7.3e-34`; spin-4 error: `5.7e-45`. The complex single-soliton integrals also equal `q_a exp(s theta)` for both spins. This directly tests additivity at a bound-state rapidity without substituting it into the integrand.

**Convergence audit.** The first quadrature used 32 digits and `[-40,40]`; it printed two spin-2 FAILs at the overambitious `1e-16` tolerance (errors about `1e-14`). The breather's slow tail has decay rate `sqrt(3)cosh(.23)/2`, so the finite cutoff matters. The independently extended range `[-90,90]` and 45-digit calculation reduces those errors to `7e-34`, with unchanged real parts. Both outputs are retained. This is a resolved tail-truncation effect, not suppressed failed evidence.

[Freeman, hep-th/9408092](https://arxiv.org/abs/hep-th/9408092), introduction p.2, explicitly limits the calculation to “all conserved charges for the single soliton solutions”; the eigenvalue is `2(1-cos(pi N/h))`. For A2 and `N=2`, this is 3, with eigenvector `(1,-1)`. Freeman's Cartan convention has the root-length denominator on the other index from the book, so left/right terminology outside simply-laced type must be translated. His construction does not establish breather additivity, and the chapter correctly says it is assumed.

No new theorem label follows from these finite tests. A check script for the charge and an explicitly conditional statement would be appropriate if the editor integrates the example.

## 4. Time delays and the semiclassical phase

**Overlap:** P1 entry/rank 8 checked trajectories, boosted counterexample and sine-Gordon phase. General A-series quantum amplitudes were not checked there. **Verdict:** trajectory formula **correct**, unrestricted time-advance inequality **wrong**; Wigner chain rule and sine-Gordon comparison **correct**; general source-to-book phase comparison **unclear**.

P1 §2.8 contains the frame counterexample and SG transform derivation. For completeness, following both trajectories for `theta_1>theta_2` gives

\[
\Delta x_1=-\log A/k_1,\qquad \Delta x_2=+\log A/k_2,
\quad k_i=m_{a_i}\cosh\theta_i.
\]

Consequently `E_1 Delta x_1=-(2h/beta^2)log A`, `E_2 Delta x_2=+(2h/beta^2)log A`, and `Delta t_i=-Delta x_i/v_i` when `v_i` is nonzero. In the caption's frame (`theta_1=.6`, `theta_2=-.6`) both are forward along their own motion and both time shifts are negative. The caption's positive expression is the *magnitude* of each displacement, not the signed x displacement of particle 2. Outside such a frame the inequalities do not persist. At identical species and zero rapidity difference, `A=0`, so `0<A<1` needs the distinct-velocity qualification.

In the rest frame of particle 2, `E_1=M_1 cosh(theta)`, so `dE_1/dtheta=p_1`. Given the convention `S=exp(i delta)` and Wigner's `Delta t=d delta/dE_1`, the displayed chain rule is exactly `delta'=p_1 Delta t`. It predicts `(2h/beta^2)log A`. P1 derives the specified chapter-3 integral's leading SG value `(16/beta_SG^2)log tanh(theta/2)`.

**Independent quantum-side expansion and unresolved sign.** I read [Hollowood, hep-th/9203076](https://arxiv.org/abs/hep-th/9203076), §§3–5, including the Gamma product and its weak-coupling expansion. His `n` equals this book's `h`; his WKB convention is `S=exp(2i delta_H/hbar)`, eq. (5.2). Hence a factor of two must be included before comparing phases. Telescoping the Euler–Maclaurin integral of the Gamma factors in (5.7), with `lambda_0=4pi`, gives for elementary transmission/identical-particle processes

\[
\log S_H(\theta)=\text{constant}+
\frac{2ih}{\beta^2}\int_0^\theta
\log\!\left[
\frac{\sin(\theta'/2i+\pi/h)\sin(\theta'/2i-\pi/h)}{\sin^2(\theta'/2i)}
\right]d\theta'+O(1).
\]

The bracket is `1/A_11(theta')`, not `A_11(theta')`. This reproduces the source's eq. (5.8), and gives the opposite derivative from the book **before a convention conversion**. The independent script `codex_phase2_phase.py` integrates the derivative of the finite sine-log integral resulting from (5.7); it does not merely compare the two printed formulas.

For bound-state fusion into species `a,b`, the constituent shifts are `i pi(a+1-2r)/h` and `i pi(b+1-2q)/h`. Multiplying their elementary interaction factors telescopes to

\[
\prod_{r=1}^a\prod_{q=1}^b
A_{11}\!\left(\theta+\frac{i\pi}{h}(a-b-2r+2q)\right)=A_{ab}(\theta).
\]

The same shifts applied to the source leading logarithmic derivative give `-log A_ab`, on a continuous branch at positive real theta. Direct quadrature and product checks at `theta=sqrt(2)` give the following coefficients after dividing `-i d log S_H/dtheta` by `2h/beta^2`:

| theory | `(a,b)` | source coefficient | book's `log A_ab` |
|---|---|---|---|
| A2 | (1,1) | +0.821164893701769 | -0.821164893701769 |
| A2 | (1,2) | +0.638597832801906 | -0.638597832801906 |
| A2 | (2,2) | +0.821164893701769 | -0.821164893701769 |
| A3 | (1,1) | +0.614517404916466 | -0.614517404916466 |
| A3 | (1,2) | +0.673630898983532 | -0.673630898983532 |
| A3 | (2,2) | +0.992335931335434 | -0.992335931335434 |

Errors in these integral/product identities are below `3e-41`. This verifies the required same- and different-species dependence on the quantum side, but not the final orientation identification. Hollowood's time-delay sign in (5.3) is also opposite to the book's passage-time convention, while his text says the rank-one amplitude is the sine-Gordon amplitude. An inverse amplitude or reversed in/out convention would resolve the derivative sign; an additive phase constant or the WKB factor of two would not. I have not established the required in/out conversion from the source's operator convention to the book's. The fusion computation concerns the leading scalar contribution and does not audit every non-diagonal channel or nonuniform limit near poles.

Thus I do **not** report a correction to the published S-matrix, or a verified all-A-series match. At `10-solitons.qmd:122`, the safe replacement for the quantum comparison is:

> With `S=exp(i delta)` and the passage-time convention used here, Wigner's relation predicts `d delta_ab/dtheta=(2h/beta^2)log A_ab(theta)` at leading order. The sine-Gordon scalar factor printed in this book has this limit. A comparison with higher-rank soliton amplitudes must first match their in/out and rapidity conventions.

Use P1 §2.8's replacement for the frame inequality, and say “displacement magnitude” in the caption at 10:124. [Fring–Johnson–Kneipp–Olive, hep-th/9405034](https://arxiv.org/abs/hep-th/9405034), eqs. (4.10), (4.14), explicitly gives the opposite signs for the two spatial shifts, agreeing with the independent asymptotics.

## 5. Lax pair, boosts, currents and allowed spins

**Overlap:** P1 ranks 5, 9 and 15. The first two were checked; rank 15's twisted centralizers were previously unclear. **Verdict:** curvature and principal-boost identity **correct**; twisted-generator description and spin sign **wrong**; all requested centralizer examples **correct**. Full quantum scope/source coverage remains partly **unclear**.

**Twisted matrix check.** `codex_phase2_extra.py` uses

\[
\alpha_1=a,\quad\alpha_0=-2a,\quad(n_0,n_1)=(1,2),\quad
H=a\,\operatorname{diag}(1,0,-1),
\]
\[
e_1=\sqrt2(E_{12}+E_{23}),\quad e_0=E_{31},\quad f_j=e_j^T,
\quad(c_0,c_1)=(\sqrt2 a,a).
\]

These obey exactly `[H,e_j]=alpha_j e_j` and `[e_j,f_k]=delta_jk (2alpha_j/alpha_j^2)H`. Substitution of the printed Lax pair gives Cartan curvature

\[
-\beta\Box\phi\,H-m^2\big(2a e^{\beta a\phi}-2a e^{-2\beta a\phi}\big)H.
\]

All off-diagonal entries cancel. Substituting the twisted field equation gives the zero 3-by-3 matrix. Keep `a` arbitrary to avoid confusing finite-root normalization (`a=sqrt(2)`) with affine-long-root normalization (`a=1/sqrt(2)`). The affine generator is in the twisted loop algebra; it is not the negative highest root vector of the fixed A1 finite algebra. P1 §2.5 supplies the needed wording replacement at 09:57.

**Boosts/gradation.** P1's derivative calculation proves the printed principal-parameter identity with no compensating gauge. In this twisted representation take `R=diag(1,0,-1)`. Since `[R,e_1]=e_1`, `[R,e_0]=-2e_0`,

\[
\lambda^R(e_1+\lambda^3e_0)\lambda^{-R}=\lambda(e_1+e_0).
\]

The homogeneous parameter is therefore `z=lambda^3`. In that parameter the boost includes `z -> exp(3vartheta)z` and constant Cartan conjugation by `G=exp(vartheta R)`, together with the light-cone prefactor `exp(∓vartheta)`. This explicitly verifies the remark about another gradation.

**Spin-3 current.** P1 §2.9 derives the identity, its conserved sum and kink integral. The identity is `partial_-T+partial_+X=0`, so `partial_t(T+X)=partial_x(T-X)`. To match the earlier generic convention `partial_-T=partial_+X_generic`, set `X_generic=-X`; otherwise the generic density is `T-X_generic`, not its sum. This is a notation inconsistency, not a failure of conservation.

A free-field one-particle check fixes the rapidity sign independently of the kink. With standard mode normalization `1/sqrt(2E_p)`, the normal-ordered quadratic part of the displayed charge has eigenvalue

\[
\frac{(E-p)^4+m^2(E-p)^2}{E}=2(E-p)^3=2m^3e^{-3\theta}.
\]

The quartic normal-ordered term does not contribute to this one-particle matrix element at tree level. Quantum renormalization may change the overall coefficient, not the Lorentz weight. This agrees with the independent kink value `-16m_0^3 exp(-3theta)/(3 beta_SG^2)`. At 09:96 replace “spin 3” by “spin -3” and `exp(3theta)` by `exp(-3theta)`; say that exchanging `partial_+` and `partial_-` gives the spin +3 partner. At 09:41 use the same conservation-law sign convention as the example or explicitly distinguish `X_generic`.

**Centralizers.** `codex_phase2_centralizers.py` builds loop matrices, graded bases and twisting automorphisms independently, then solves `[E_1,X_s]=0` with the twisting constraints. It uses A-type defining matrices, the symplectic fixed subalgebra for C2, and a Chevalley-bracket construction of so(8) with triality for G2/D4. It does not enter the expected exponent lists into the linear solve.

| algebra | tested residues with centralizer dimension 1 | period |
|---|---|---|
| A2, A3, A4 untwisted | all residues except 0 | 3, 4, 5 respectively |
| A4 twisted | 1,3,7,9 | 10 |
| C2 untwisted | 1,3 | 4 |
| D3 twisted (= A3 twisted) | 1,3,5 | 6 |
| G2 untwisted | 1,5 | 6 |
| D4 triality twisted | 1,5,7,11 | 12 |

Other tested graded dimensions are zero. Maximum residual is `3.7e-15`; the smallest nonzero singular value is at least `.77`, so the numerical nullities are well separated. The two dual pairs yield the same infinite spin sets despite different convenient periods. For general A type, a cyclic regular matrix `E` has centralizer polynomials in `E`, with `E^h=zI`; tracelessness removes degrees divisible by h. The twisted anti-transpose involution on A_(2n) removes even degrees as well, leaving precisely odd degrees not divisible by `2n+1`. These arguments support the finite checks. E6 twisted and arbitrary twisted multiplicities were not independently reconstructed.

**Drinfeld–Sokolov mechanism.** Split each graded loop-algebra component into `ker ad E_1` and `im ad E_1`. After the leading potential is brought to constant semisimple form, the image part at each successive inverse-spectral order is removed by solving the algebraic `ad E_1` equation; the source at that order is a differential polynomial in previous coefficients. Thus no spatial integration is required and the remaining coefficients are local in the fields and their derivatives. Zero curvature gives their conservation laws. In the centrally extended algebra the Heisenberg commutator is central; it is abelian after passing to the loop quotient used by the matrix Lax pair. The chapter's description is correct at this level, with ordinary monodromy boundary assumptions understood. Freeman §2.3 independently describes the two chiral abelianizations. The original [Drinfeld–Sokolov paper](https://math.berkeley.edu/~frenkel/Drinfeld-Sokolov.pdf), pp.2011–2013, Propositions 6.2 and 6.6, was also located and checked: the former constructs the formal gauge and differential-polynomial coefficients, while the latter identifies the residual coefficients as densities modulo total derivatives. Its exact statement includes: ‘these densities do not depend on the arbitrariness in the choice of U’.

**Quantum references.** Feigin–Frenkel [hep-th/9310022](https://arxiv.org/abs/hep-th/9310022), §4.7, Theorem 4.7.3, treats quantum deformation of local integrals. Its supplied proof has the explicit restriction “all the exponents of g are odd and the Coxeter number is even”; it says the general proof will appear separately. The construction is at generic coupling, not a nonperturbative Hilbert-space theorem at every special coupling. Delius–Grisaru–Zanon [hep-th/9201067](https://arxiv.org/abs/hep-th/9201067), §2, explicitly verifies a quantum spin-3 current and says “Assuming that higher-spin quantum-conserved currents” before using the larger hierarchy. Niedermaier [hep-th/9401078](https://arxiv.org/abs/hep-th/9401078) constructs charge eigenvalues in the factorized-scattering description and compares them with perturbation theory; it is not an independent all-model construction of local operators from the interacting Lagrangian.

Therefore the broad physical assertion may be established elsewhere, but these citations cannot all be described as proofs of the full unrestricted sentence. At 09:32 a more precise replacement is:

> Quantum conserved charges are supported by screening-operator constructions, explicit current calculations and the factorized-scattering description [@feigin1993; @dgz1992; @niedermaier1994]. These arguments have different scopes: generic-coupling algebraic constructions, particular quantum currents, and scattering-state eigenvalues, respectively.

Use P1 §2.15's simply-laced qualification for the Cartan eigenvector sentence at 09:98.

## 6. Masses, cubic couplings and Dorey's rule

**Overlap:** P1 entries 7, 14 and 20; all checked for the stated examples, with entry 11 covering C_n masses. **Verdict:** mass formulas, cubic factor, representatives and D4 figure **correct**; the unqualified area identity and unit wording **wrong**; Dorey's **[Established]** label appropriate. Some original-source attributions remain **unclear**, as listed in §9.

P1 §§2.11, 2.20 independently construct the C_n, D_n and E8 mass matrices and check the table and trace sum rule. P1 §2.7 proves the missing nonzero-coupling qualification and the need to use the actual masses in the area, or include an extra `m^2` when the area is dimensionless. Taylor expansion gives `-(i beta)^3/beta^2=i beta`, so the imaginary-coupling cubic tensor is exactly i times the real one in the same basis.

**New A4 and E6 checks.** `codex_phase2_couplings.py` constructs the potential's cubic tensor directly from roots. For A4 it uses the four normalized discrete Fourier modes of the five cyclic differences. The Fourier sum vanishes unless `a+b+c=0 mod 5`. Of all 20 unordered triples, four are nonzero: `(1,1,3)`, `(1,2,2)`, `(2,4,4)`, `(3,3,4)`; maximum area error is `1.78e-15`.

For E6, roots are obtained by Cholesky factorization of the explicitly specified Cartan matrix with chain `1-3-4-5-6` and node 2 attached to 4. The highest-root marks are `(1,2,2,3,2,1)`. The affine centre rotation cycles `0,1,6` and `2,3,5`, fixing 4. Its orthogonal action commutes with the independently constructed mass matrix. Joint mass/centre eigenvectors, with characters assigned using fundamental weights, give the charge basis. An independently generated Coxeter orbit test predicts precisely the 18 nonzero couplings among all 56 unordered triples, with no mismatches and maximum area error `6.11e-15`.

The computed E6 masses in that node order at `m=1` are

`(1.1260325, 1.59245043, 2.17532775, 3.0763780, 2.17532775, 1.1260325)`.

Real orthonormal mass eigenvectors give counterexamples even among *nonzero* couplings: in A4 a representative ratio `|C|/(4Delta/sqrt(h))` is `1/sqrt(2)`, and in the E6 real basis it is `.4809652855442468` for the first cubic triple. These are basis effects, not violations of the charge-basis rule. No full all-algebra claim is inferred solely from these finite checks.

For `w=s_B s_W`, `w^-1=s_W s_B`. If a is white, `s_B lambda_a=lambda_a` and `(1-w^-1)lambda_a=alpha_a`; if black, `s_B lambda_a=lambda_a-alpha_a` and `s_W lambda_a=lambda_a`, giving `s_W alpha_a`. This verifies the representatives with the book's order. P1 §2.14 fully checks D4's roots and figure: the `(1,1,central)` triangle has sides `(1,1,sqrt(3))`, interior angles `(pi/6,pi/6,2pi/3)`, and hence the particle fusing angle for the equal outer pair is `u=pi/3`, complementary to the interior angle opposite the central mass. The printed caption correctly gives the *interior* angles, not an assertion that `u=2pi/3`.

**Additional original source inspected.** The scanned [Fring–Liao–Olive preprint](https://openaccess.city.ac.uk/id/eprint/963/2/Smat.pdf) was read directly, including pp.5–7. Equation (34) gives the root-sum selection, followed by “This is the confirmation of Dorey's rule”; equation (44) gives `C_ijbar(k)=4i beta epsilon(i,j,q) Delta_ijk/sqrt(h)`, with the simply-laced structure-constant modulus 1. This verifies the chapter's attribution of a Lie-algebra derivation, independently of citation lists. Dorey's sequel was checked in P1 §2.14. No new theorem label is warranted for either established result.

At `09-lagrangian.qmd:130–136` use P1 §2.7's exact replacement: “every **nonzero** cubic coupling” and “triangle of the **actual** masses”, removing “with m=1”.

## 7. Hirota solitons, charges, degrees and Lie duality

**Overlap:** P1 entries 10–13 and 16, with 6 for fusion and 4 for breathers. **Verdict:** bilinear factors, examples, charge count and tested mass duality **correct**; vertex formula without spacetime factors **wrong**; universal exact degrees remain **unclear/conditional**.

P1 §2.12 completely derives the non-simply-laced and twisted bilinear equations, including the transpose, the factor `n_j alpha_j^2/2` and the A1 double neighbour. P1 §2.16 supplies the missing `W_k(x,t)` in the displayed vertex matrix element and verifies the printed C2 polynomial with residuals `[0,0,0,0]`.

**Embedded real sine-Gordon kinks.** For `h=2n`, species `a=n`, let `E=exp(2mx)` be positive and set `Im xi=pi/2`. Then `tau_even=1+iE`, `tau_odd=1-iE`. With cyclic roots and `v=sum_even alpha_j`,

\[
\sum_{\rm odd}\alpha_j=-v,\quad v^2=2n,\quad
\alpha_j\cdot v=2(-1)^j,\qquad
\phi=-\frac{2v}{\beta}\arctan E.
\]

This field is real. Writing `phi=(v/|v|)psi`, the reduced potential is

\[
V=\frac{2nm^2}{\beta^2}(1-\cos(b\psi)),\qquad b=2\beta/|v|,
\]

so its sine-Gordon mass parameter is `m_0=2m`. The canonical kink mass is `8m_0/b^2=8nm/beta^2`, equal to `2h m_a/beta^2`. `codex_phase2_extra.py` constructs the roots for n=2,3, obtains `v^2=4,6`, alternating root products `+2,-2`, and energies 16,24 at unit parameters. This verifies both requested examples analytically as well as by finite root arithmetic.

**G2 orbit-degree check.** Fold D4 by triality. Denote the affine tau by F, the common outer tau by G, and the central tau by H. The static folded equations at m=1 are

\[
FF''-(F')^2=F^2-H,\quad GG''-(G')^2=G^2-H,
\quad HH''-(H')^2=2(H^2-FG^3).
\]

For the invariant central species, with `E=exp(sqrt(6)x)`,

\[
F=G=1+E,\qquad H=1-4E+E^2,
\]

giving degrees `(1,1,2)=n^vee`. For the three-element outer orbit, coefficient recursion from initial coefficients `(F_1,G_1,H_1)=(3,-1,0)`, `E=exp(sqrt(2)x)`, gives

\[
F=1+3E+E^2+E^3/27,\quad
G=1-E-E^2/3+E^3/27,
\]
\[
H=1+3E^2-16E^3/27+E^4/3+E^6/729.
\]

Substitution in all three equations gives exact polynomial residuals `[0,0,0]`. The degrees are `(3,3,6)=3n^vee`. The C2 fixed species has degree one at every node and the orbit pair degree two, as P1 already checked. These examples confirm both constructions; they do not prove nonvanishing top coefficients for every folded orbit. Nilpotency alone gives only an upper bound.

The assertion that Part III relies on Hypothesis D only for the degree claims is reasonable **if the general mass/energy results are independently quoted from the established soliton literature**. A derivation of the general mass formula solely by counting tau degrees would inherit D. No extra proof of the universal exact-degree assertion is provided by the finite examples; keep the conditional/numerical wording at 10:77 and do not upgrade it.

**Topological charges and masses.** P1 §§2.10 and 2.13 fully check the endpoint formula, phase chambers, count, A3 figure and mass-exercise first integral. [McGhee, hep-th/9307035](https://arxiv.org/abs/hep-th/9307035), p.2, displays `h_tilde_a=h/gcd(a,h)` and explicitly says it “turns out to be equal to the number of topological charges”. Thus the exact count, not only an upper bound, is indeed in the cited source.

**Lie duality source and independent check.** P1 §2.11 derives the C_n/B_n mass comparison. The two G2 families above have masses proportional to `sqrt(6)` and `3sqrt(2)`, giving ratio `sqrt(3)`, as the dual rank-two particle system requires. [Kneipp–Olive, hep-th/9404030](https://arxiv.org/abs/hep-th/9404030) does not discuss only twisted algebras: §6, pp.20–21, explicitly treats untwisted non-simply-laced folding and derives

\[
M_{\langle i\rangle}=|\langle i\rangle|M_i.
\]

It explains that a single folded soliton comes from coincident parent constituents, while the invariant configuration has unchanged energy and momentum. The attribution in 10:104 is therefore supported; it is not a mistaken citation merely because the title emphasizes twisted algebras. General all-type Lie duality has not been independently rederived in this pass beyond the stated root/folding checks.

**Part VIII pointers.** In `prp-stored-charge`, transmission changes only `xi -> xi+log z`. Substitution into the independently derived endpoint expression can only choose another phase chamber of the *same* `h/gcd(a,h)` cycle; it cannot add a weight absent from that family. This verifies the pointer for regular one-soliton transmission in that proposition's scope, not for arbitrary defects or absorption. In `prp-bianchi`, `a_i=-exp(-theta_i)`, and direct division gives `(a_1-a_2)/(a_1+a_2)=-tanh((theta_1-theta_2)/2)`. Squaring yields `A_11`; hence it is exactly the square of the inverse of the displayed Bianchi coefficient.

## 8. Other chapter-11 statements

**Overlap:** P1 ranks 1–4, 17 and 19. **Verdict:** `i phi`, the data conjugation and corrected regular-breather statements are **correct**; ordinary-charge continuation through poles and physical-sector interpretations stated as facts are **wrong**. The scoped A3 static obstruction is **correct**.

For `psi=i phi`, multiplying the imaginary equation by i gives

\[
\Box\psi=-\frac{m^2}{\beta}\sum_j n_j\alpha_j e^{\beta\alpha_j\cdot\psi},
\]

which is the real-coupling equation. Thus **+i**, as printed, is right. For `psi=-i phi` the right side has the opposite sign with exponent `-beta alpha_j.psi`; this is not the same theory except where an additional field inversion symmetry exists. A Theta-fixed field makes `psi` real, and its vacuum asymptotics have no topological difference (with a common vacuum shift removed).

The interaction identity follows from `cos(pi((h-a)-(h-b))/h)=cos(pi(a-b)/h)` and `cos(pi((h-a)+(h-b))/h)=cos(pi(a+b)/h)`, with real coefficients in cosh. Combined with `conjugate(omega^{ja})=omega^{j(h-a)}`, it proves `eq-theta-data` term by term in the finite tau sum. Log branches must be continued consistently; a vacuum shift accounts for branch changes.

P1 §2.1 completely refutes the zero-residue inference with a double pole and divergent symmetric cutoff integral. The **tau polynomials** are entire in x,t and their parameters away from parameter poles of the interaction coefficient, so their algebraic continuation through a spacetime zero is harmless. Their logarithms, the field and the real-line charge integrals are not thereby regular. “The Hirota data continue analytically” should be replaced by this precise distinction, using P1's replacement at 11:60 and 91.

Takács–Watts explicitly relates complex time delays to non-unitary quantum scattering, but does not prove an identification of every singular final solution with a Theta-broken sector in the book's later Krein framework. At `11-reality.qmd:62`, replace the last two sentences by:

> The relation to a Theta-broken quantum sector is an interpretation motivated by the later scattering analysis. The classical calculation establishes complex time delays and possible loss of regularity; it does not by itself identify a physical quantum sector.

The A3 missing-charge statement is covered by P1 §2.19's independent unstable-manifold/Hirota completeness argument and subset-energy enumeration. Its legitimate scope is regular static solutions approaching vacua at both ends. Replace “no static saddle” at 11:72 by that phrase. The subsequent bound-state proposal is a conjecture, and the chapter correctly calls it one; no intersection-number or quantum-completeness claim has been verified here.

For the sine-Gordon breather take conjugate rapidities and `xi_2=conjugate(xi_1)+i*pi`. Then `E_2=-conjugate(E_1)`, `A=-tan^2 u`, and `tau_1=conjugate(tau_0)` with positive real part `1+tan^2(u)|E_1|^2`. This proves real, regular fields and verifies the printed sign. Without the extra i*pi, the tau functions are real; a rest-frame choice of time makes the two linear terms negative and produces roots. P1 §2.4 gives explicit tau functions and the general neutral-breather range. Calling the nonzero real sine-Gordon breather a Theta-fixed *field* is nevertheless wrong: it is the reduced species/rapidity multiset that is fixed.

## 9. Citation audit and all 17 exercises

**Overlap:** distributed across P1 entries 1–20. The coordinator specifically added systematic citation coverage, one-particle/current constructions and exercise coverage. **Verdict:** mixed as itemized below; inaccessible or only abstract-level sources are **unclear**, not endorsed from reputation.

### Citation-by-citation audit

The locations below cover every bibliography-bearing sentence in chapters 9–11. Repeated uses are grouped, with their separate locations shown. “Correct” applies to the specified claim and scope, not every assertion of the cited paper. Exact supporting quotations/equations already given in this report or P1 are not repeated.

| Citation and book locations | Verdict and direct evidence/scope |
|---|---|
| `olive1985`: 09:32,82,188 | **Correct at abstract level; full exponent derivation unclear.** The [publisher's abstract](https://doi.org/10.1016/0550-3213(85)90347-5) says “All the corresponding hamiltonians have vanishing mutual Poisson brackets.” It describes abelianizing the gauge potentials to obtain local densities. The full paper's precise spin classification was unavailable; independent centralizers are in §5. |
| `feigin1993`: 09:32 | **Correct with qualifications; unrestricted scope unclear.** §4.7 and its stated restriction/generic coupling are discussed in §5 above. |
| `dgz1992`: 09:32 | **Correct for explicit quantum-current evidence; not an all-spin proof.** §2 distinguishes the spin-3 calculation from the assumed higher hierarchy. |
| `niedermaier1994`: 09:32 | **Correct for quantum eigenvalue construction and perturbative checks.** Abstract and introduction specify a free-field realization of the ZF algebra. No unrestricted nonperturbative local-operator theorem was verified. |
| `drinfeld_sokolov1985`: 09:82 | **Correct.** The original paper, Propositions 6.2 and 6.6 (pp.2011–2013), supplies the formal gauge and local densities modulo total derivatives; §5 gives the independent recursive argument. |
| `kac1990`: 09:88 | **Unclear source check:** the actual relevant pages of the third edition were not obtained. The requested examples and general A argument are independently checked; this does not certify every twisted exponent or the source's h convention. |
| `klassen1990`: 09:98 | **Unclear:** full original not read. The book needs the simply-laced/left-right/classical-quantum qualifications already identified in P1 §2.15. |
| `dorey1991`: 09:98,140,151,154,188 | **Unclear full original; attribution supported by its accessible abstract and the verified sequel/Fring–Liao–Olive derivation.** The [original abstract](https://doi.org/10.1016/0550-3213(91)90428-Z) describes a proposed root fusing rule. Do not claim to have checked its complete proof. |
| `dorey1992`: 09:154 | **Correct.** [hep-th/9110058](https://arxiv.org/abs/hep-th/9110058), eq. (2.8), orbit construction and Coxeter-plane discussion; P1 §2.14 and §6 here independently verify representatives and examples. |
| `braden1990`: 09:110,130,188 | **Unclear original-source audit:** publication record found, full source unavailable. Masses and corrected area rule independently verified in P1 and §6; the unqualified area formula in the citing paragraph is wrong. |
| `freeman1991`: 09:110,188 | **Unclear original-source audit:** full paper not obtained. Its PF mass result is independently supported by the inspected Fring–Liao–Olive derivation; that is not a substitute quotation from Freeman 1991. |
| `mikhailov1981`: 09:188 | **Correct for the historical construction claim.** The actual paper text, [Two-dimensional generalized Toda lattice](https://doi.org/10.1007/BF01209308), abstract and §§2–4, constructs the root potentials, zero curvature and mass matrix. Exact abstract excerpt: “The reduction group and conservation laws are found and the mass spectrum is calculated.” This is not a proof of historical priority over every earlier special case. |
| `olive1983`: 09:188 | **Unclear:** full paper was not read; bibliographic attribution alone is insufficient evidence. |
| `fring1991`: 09:188 | **Correct.** Original scanned preprint pp.5–7, eqs. (24), (34), (44), checked directly; see §6. |
| `hollowood1992`: 10:22,65,114,172 | **Correct for A-series solitons and the Hirota method; unclear as the sole support for every D/E polynomial claim.** [hep-th/9110010](https://arxiv.org/abs/hep-th/9110010) develops simply-laced solitons; its explicit general N-soliton construction is A series. The book's all-type statement relies on the collection of later papers, not this citation alone. The A interaction coefficient follows by substituting the dispersion relation into its characteristic-polynomial ratio, as independently verified in P1. |
| `olive1993a`: 10:22,69,96,172,174; 11:33,96 | **Abstract supports construction/total derivative/additivity; exact full-source check unclear.** Original full text unavailable, especially the C prefactor. Do not cite its abstract as an exact quotation of that equation. |
| `olive1993b`: same locations, plus 10:77 | **Correct for vertex construction, nilpotency bound and additive boundary momenta in its decay domain.** §8.3 explicitly supplies the complex-rapidity restriction. P1 §2.16 records the exact nilpotency qualification. The displayed book matrix element is still missing the spacetime factors. |
| `kneipp1996`: 10:22,104,172 | **Correct for folding and both twisted and untwisted non-simply-laced mass rules.** §6, especially pp.20–21, supports the coincident-orbit formula; see §7. |
| `bernard1991`, `felder1992`: 10:27 | **Unclear exact source-to-book convention check:** the full original texts were not obtained/read. [Bernard–LeClair's original abstract](https://doi.org/10.1007/BF02099173) confirms non-local affine quantum-group symmetry and resulting formal S-matrices, but not the exact sign conversion. The specified physical q and opposite coproduct require a separate quantum-group audit; nothing in this classical pass silently certifies their sign for every algebra. |
| `mcghee1994`: 10:88,90,178; 11:66,96 | **Correct.** p.2 explicitly gives the equality `h/gcd(a,h)`, the charge cycle and missing weights; own derivation and enumeration in P1 §2.10. |
| `fring1994`: 10:114,122,172 | **Correct for interaction-factor interpretation and time shifts.** Eqs. (4.10), (4.14) match both spatial signs; the negative time delay is a centre-of-momentum statement. The book's unrestricted lab-frame inequality needs correction. |
| `harder1995`: 10:136,172 | **Correct attribution, but the book omits essential restrictions.** [hep-th/9409035](https://arxiv.org/abs/hep-th/9409035), §3.2, excludes positive interaction coefficient for neutral type-B breathers and also treats charged families. P1 §2.4 independently derives the missing band. |
| `tw1999`: 10:138,181; 11:39,62,96 | **Correct for A2 additive formula and problematic scattering, with domain/convention caveats.** Eq. (4.6) quoted in §2; §4 gives complex delays/singular scattering, and the discussion distinguishes quantum excited solitons without established classical counterparts. The new Theta-sector identification at 11:62 is an interpretation, not a result proved there. |
| `mackay1993`: 10:172 | **Correct contribution to the multi-source existence claim; insufficient alone for every species.** [hep-th/9208057](https://arxiv.org/abs/hep-th/9208057) constructs solitons by diagram automorphisms. Kneipp–Olive §6 explicitly notes the earlier omission of coincident-orbit species. The book's collective citation is therefore better than assigning completeness to MacKay–McGhee alone. |
| `aratyn1993`: 10:172 | **Correct for recursive Hirota construction and mass analysis.** [hep-th/9212086](https://arxiv.org/abs/hep-th/9212086), §§2–3, derives the conformal tensor improvement and recursive scheme; the explicit folded examples here independently test its mechanism. No independent all-algebra completeness proof is claimed in this audit. |
| `zhu1995`: 10:172,174; 11:96 | **Correct for reality/regularity analysis and pairing conditions, not unrestricted real energies.** Eq. (2.5) concerns paired tau functions and pointwise-real density; §1 gives the static first integral. I did not locate the full C tensor formula there. P1 §2.18 explains why its condition is not identical to the reduced-data criterion. |
| `freeman1995`: 10:174; 11:39,46,96 | **Correct for single solitons, after normalization/convention qualification.** p.2 explicitly says single solitons, left Cartan eigenvectors and a possible overall phase; §3 computes boundary charges. It does not prove complex-breather additivity. |
| `khastgir1996`: 10:180; 11:96 and the named claim 11:47 | **Correct for c-soliton blow-up and negative-energy initial data.** §5, eqs. (5.5)–(5.7), gives the tau-zero linear system; own example is in §1. Eq. (5.13) chooses `phi(x,0)=0`, `phi_t(x,0)=i f(x)` with real square-integrable f, yielding `E=-1/2 integral f^2<0` immediately from the canonical density. The paper does not prove global regularity of these negative-energy solutions; neither should the book. |

The last negative-energy computation also separates *negative energy* from *complex energy*: the former follows from an imaginary initial velocity without using any soliton energy formula. All citations used by the Notes sentences are included in the table; none is endorsed merely because it appears in an older report.

### Exercise outlines

The following covers all 17 named exercises. “Solvable” means with the chapter formulas and the indicated earlier background, after the corrections explicitly listed. All calculations are independent derivations here or in the cited P1 section.

| Exercise | Verdict, outline and solvability |
|---|---|
| 09 `exr-lax` | **Correct; solvable.** Insert the supplied 2-by-2 raising/lowering matrices, `[e,f]=sqrt(2)H`, `alpha_0=-alpha_1`. Root terms cancel and the Cartan coefficient yields `Box phi=-(2sqrt(2)m^2/beta)sinh(sqrt(2)beta phi)`. P1 §2.5 gives the exact zero-matrix check. |
| 09 `exr-spin3` | **Identity and task correct; surrounding spin sign wrong.** Differentiate T and X using `partial_-partial_+phi=-m^2 sinh(beta phi)/beta`; terms cancel pairwise. Replace sinh/cosh by sin/cos and the quartic coefficient by its negative. With `phi_z=2 sech(z)/beta_SG`, integrate sech squared and fourth powers to get `Q=-16m_0^3 exp(-3theta)/(3 beta_SG^2)`. Solvable after specifying the SG parameter dictionary; P1 §2.9. |
| 09 `exr-dn-masses` | **Correct; solvable.** Let N=n-1. `sum_(a=1)^(N-1) sin^2(pi a/(2N))=(N-1)/2`; the sine masses contribute `4(n-2)m^2`, the two spinor masses `4m^2`, totaling `4(n-1)m^2=2hm^2`. D4 eigenvalues are `(2,2,2,6)m^2`. |
| 09 `exr-cn-masses` | **Correct; solvable using folding background.** Reflection exchanges Fourier modes a and 2n-a, so each pair has one invariant combination; mode n is itself invariant. The restricted eigenvalues are `4m^2 sin^2(pi a/(2n))`, a=1,…,n. P1 §2.11 independently constructs the matrix. |
| 09 `exr-shg-cubic` | **Correct after the nonzero-coupling area qualification.** The ±root contributions cancel in rank one. In A2 the Fourier sum vanishes except `a+b+c=0 mod3`, so only 111 and 222 survive, each with modulus `3m^2 beta`. Equal masses `sqrt(3)m` give `Delta=3sqrt(3)m^2/4`; `m_2^2=2m_1^2(1+cos u)` yields `u=2pi/3`. P1 §2.7. |
| 09 `exr-dorey-a2` | **Correct; solvable.** Choose node 1 black, 2 white. Then `gamma_1=alpha_1+alpha_2`, `gamma_2=alpha_2`; the first orbit is `{alpha_1+alpha_2,-alpha_1,-alpha_2}`, the other its negatives. Each orbit's three vectors sum to zero; no mixed closing triple exists. This matches precisely 111 and 222. |
| 10 `exr-bilinear` | **Correct; solvable.** `exp(i beta alpha_k.phi)=product_j tau_j^(-a_jk)`. Multiply each log equation by `2alpha_j/alpha_j^2`; the added constant vanishes because `sum n_j alpha_j=0`. Finally `1/2(D_x^2-D_t^2)tau.tau=-tau^2 Box log tau`. P1 §2.12. |
| 10 `exr-moving-soliton` | **Correct; solvable.** Differentiating a single exponential gives `k^2-nu^2=m_a^2(cosh^2theta-sinh^2theta)=m_a^2`; the right side gives `4m^2 sin^2(pi a/h)omega^(ja)E`. No additional factor of two occurs. |
| 10 `exr-soliton-mass` | **Correct on the regular slice specified; solvable.** Dot the static field equation with phi' to get the first integral with zero vacuum constant. Then `E=2 integral V=(2/beta^2)sum_j(2/alpha_j^2)[partial_x log tau_j]`. For A2 each boundary jump is `m_a`; hence `E=6m_a/beta^2`, independent of xi within a regular phase sector. P1 §2.13 and direct integration §2 here. |
| 10 `exr-charges-a3` | **Correct; solvable.** Continue each ray's phase into `(-pi,pi)`. For a=1 the four jumps give the four defining weights; for a=2, tau_even and tau_odd alternate, giving only `±(1/2,-1/2,1/2,-1/2)` in the standard cyclic root coordinates (up to relabelling). Four other two-form weights are absent. P1 §2.10. |
| 10 `exr-lie-duality` | **Correct; solvable with D-to-B folding background.** Identifying the two D_(n+1) spinor nodes retains one spinor mass and all chain masses. Ratios to the spinor are `2sin(pi a/(2n))` plus 1, exactly the orbit-pair C_n soliton ratios. P1 §2.11. |
| 10 `exr-fusion` | **Wrong without species/range restriction; corrected version solvable.** Euler's formula gives `sin A exp(iB)+sin B exp(-iA)=sin(A+B)`. For A2 specify `1+1 -> 2`: take rapidities `theta +/- i*pi/3`, send the two linear amplitudes to zero and hold their A-weighted product fixed. The surviving tau is species 2. For `2+2 -> 1` use the upper-channel prescription in P1 §2.6; `1+2` is not fusion to species 3. Paste “for a+b<h” and specify the channel at 10:163. |
| 10 `exr-breather-sg` | **Correct; solvable after choosing regular phases.** Use `A=-tan^2u`, `xi_1=-log(tan u)`, `xi_2=xi_1+i*pi`. With `y=m_0 x cos u`, `z=m_0 t sin u`, `tau_0=1+e^(2y)-2i cot(u)e^y sin z`, `tau_1=conjugate(tau_0)`. Their log ratio gives the usual real breather, up to time origin/sign and the stated complementary-angle convention. P1 §2.4. |
| 11 `exr-theta-solutions` | **Correct; solvable.** Conjugate the imaginary equation and multiply by -1; its right side becomes the equation evaluated at `-conjugate(phi)`. The real equation is preserved by plain conjugation. Under `psi=i phi`, Theta becomes `psi -> conjugate(psi)`. A2's affine root set is not inversion invariant, disproving plain field parity there. The A1 exception is correct. |
| 11 `exr-theta-interaction` | **Correct; solvable.** Use cos-evenness and its 2pi period as in §8, then conjugate each subset term in the N-soliton polynomial. |
| 11 `exr-breather-real` | **Correct for the stated nondegenerate bound-state setting, with regularity/domain distinctions.** Conjugate species, rapidities and xi make tau real. In the rest frame, a time with destructive interference gives `1-2r+Ar^2`, which is negative somewhere since A<1; its zeros sweep curves. In the regular neutral domain A<0 the large-r sign already proves this. Adding i*pi in SG gives `E_2=-conjugate(E_1)` and conjugate taus of positive real part. It is reduced-data invariance, not invariance of the real field. P1 §2.4 and §8 here. |
| 11 `exr-singular-points` | **Wrong as asked; corrected version solvable.** Rank-two real Jacobian of `(Re tau,Im tau)` gives an isolated spacetime zero. A simple zero on a slice gives `partial_x^2 log tau=-1/(x-x_*)^2+O(1)`, no residue. The integral diverges nonetheless: P1's SG example has `E_epsilon=8(1-coth(2epsilon))=-4/epsilon+8+O(epsilon)`. Replace the requested conclusion by explaining the distinction between a divergent ordinary integral and an optional finite-part prescription (P1 §2.1). |

## 10. Overlap table

Rows C1–C9 cover each coordinator item; rows P1–P20 cover each entry of the frozen phase-1 risk list separately. P1 entry numbers are also their ranks. Composite items have explicit mixed verdicts rather than a misleading single “correct”.

| Claim | Coordinator list | Phase-1 list | Verdict |
|---|---|---|---|
| C1: selection rule, H6, parity and novelty | 1 | 3,18; domains 1,2,4 | wrong unrestricted; conditional implication correct; exact novelty unclear |
| C2: energy proposition and tensor | 2 | 2,13,1 | wrong unrestricted; regular-domain formula/tensor tests correct |
| C3: A2 spin-2 additivity | 3 | 3,15 (flagged, not integrated) | correct in new direct tests; general extension unclear |
| C4: shifts and semiclassical phase | 4 | 8 | shifts/SG correct; frame inequality wrong; higher-rank convention match unclear |
| C5: Lax/boost/current/spins | 5 | 5,9,15 | curvature/boost/tested spins correct; generator wording/spin sign wrong |
| C6: masses/couplings/Dorey | 6 | 7,14,20,11 | masses/tested rule correct; unconditional area identity wrong |
| C7: solitons/charges/degrees/duality | 7 | 10,11,12,13,16 | examples/count/factors correct; vertex formula wrong; universal degrees unclear |
| C8: other reality claims | 8 | 1,3,4,19 | +i and conjugation correct; pole-integral/sector overclaims wrong |
| C9: citations and exercises | 9 | distributed; 1,4,6,9,10,12,13,18 especially | mixed, detailed in §9 |
| P1: singular energy integrals | 2,8,9 | 1 | wrong zero-residue inference |
| P2: arbitrary-complex-rapidity energy formula | 2 | 2 | wrong domain |
| P3: universal selection rule and label | 1,3 | 3 | wrong scope; conditional proof correct |
| P4: breather existence band | 8,9 indirectly; not explicitly requested as a band | 4 | wrong unrestricted band |
| P5: twisted Lax generator wording | 5 | 5 | wrong description; curvature correct |
| P6: fusion range and endpoint | 9 (exercise); not a separate fusion-range prompt | 6 | wrong unrestricted species formula |
| P7: area-rule selection/unit qualification | 6 | 7 | wrong unconditional formula |
| P8: time delay/frame and SG phase | 4 | 8 | frame inequality wrong; SG result correct |
| P9: spin-3 sign/current | 5,9 | 9 | rapidity sign wrong; identity correct |
| P10: topological-charge cycle | 7,9 | 10 | correct |
| P11: C_n/B_n folded masses | 6,7,9 | 11 | correct tested family; broader full derivation incomplete |
| P12: bilinear transpose and factors | 7,9 | 12 | correct |
| P13: regular soliton mass/tensor | 2,7,9 | 13 | correct with convergence/gauge scope; full source formula unclear |
| P14: Dorey/degenerate basis/D4 figure | 6 | 14 | correct tested cases; additional A4/E6 checks pass |
| P15: twisted spins and Cartan convention | 5; 3 for A2 | 15 | requested centralizers now correct; general quantum eigenvector wording needs qualification |
| P16: missing vertex spacetime factor/degrees | 7 degree question; no explicit missing-factor question | 16 | displayed vertex formula wrong; exact universal degrees conditional |
| P17: vacua and complex “minimum” | no | 17 | vacua correct; unqualified minimum wording wrong |
| P18: novelty/source support | 1,9 | 18 | broad originality unsupported; exact novelty unclear |
| P19: missing static saddles | 8 | 19 | correct for regular vacuum-to-vacuum A3 solutions |
| P20: exceptional mass table | 6 | 20 | correct |

## 11. What each list found

“Found” below distinguishes an explicitly targeted defect from mere occurrence in a broad section. Section 10 also records indirect coverage, so the comparison does not depend on silently changing the overlap definition.

**Errors explicitly identified through both lists:** singular-charge convergence; the missing domain of the energy proposition; selection-rule scope/label and H6 overclaim; twisted affine-generator wording; spin-3 rapidity sign; unqualified time-advance sign; missing area-rule selection/scale qualification. The coordinator's questions supplied substantially deeper stress tests, especially parity normalization and the higher-charge integrals, for concerns already present in the blind risk list.

**Errors first targeted by the phase-1 list:** the regular neutral-breather mass band, the fusion species/range/endpoint problem, the missing spacetime factors in the vertex formula, and the use of “minimum” for a complex potential. The coordinator touches the first two via broad exercise/citation checks and the third via the degree paragraph, so under a broad “any overlap” definition only the minimum wording is strictly phase-1-only. The precise distinction is visible in the table. P1 supplies paste-ready corrections for all four.

**Additional errors/overclaims exposed by the coordinator's targeted follow-up:** the mismatch between the generic conservation-law sign and the displayed current's X convention; treating the Takács–Watts-to-Theta-sector interpretation as an established identification; and conflating the scopes of the quantum-conservation references when presenting them as an unrestricted proof. The last is a citation-scope overclaim, not a disproof of quantum integrability. Phase 1 flagged related broad physical/quantum language but did not perform these source/current checks. The explicit higher-rank S-matrix discrepancy is **not** counted as an error because the convention match remains unclear.

The blind list was better at domain and degenerate-limit failures of compact formulas. The coordinator list was better at demanding the missing constructions, auditing attribution, and distinguishing a computationally successful example from a theorem. Their combination substantially strengthens the evidence without turning unproved extensions into established results.

## 12. Claims not fully checked, and reasons

- **Full higher-rank quantum phase match in the book's convention:** the leading elementary/fused scalar calculation is done, including A2/A3 equal and unequal species, but its sign is opposite before converting source in/out conventions. No proved conversion or exhaustive non-diagonal-channel audit was obtained. Do not claim either a confirmed match or a published-paper correction.
- **A general higher-charge additivity/phase-independence theorem for complex multi-solitons:** exact A2 currents and converged integrations strongly support the stated examples. They do not establish a global analytic continuation across singular parameter loci or all algebras.
- **Full stress-tensor formula for every exceptional/folded solution:** trace/static identities are derived, all components checked for the requested A2 examples, and the normalized tau gauge clarified. The original OTU prefactor passage was inaccessible. Zhu–Caldi's paper was read, but that full formula was not located there.
- **Universal exact polynomial degrees and coincident-orbit products:** C2/G2 examples and nilpotency bounds are checked; the full nonvanishing statement remains Hypothesis D. No upgrade is recommended.
- **All twisted exponents and all quantum charge eigenvectors:** the requested centralizers and general A argument pass. E6 twisted was not newly constructed, and special-coupling/operator-domain issues are not settled by generic screening algebra or numerical centralizers.
- **Original sources not fully available/read:** `braden1990`, `freeman1991`, `dorey1991` (abstract only), `olive1993a` (abstract only), `olive1985` (abstract only), `olive1983`, `kac1990`, `klassen1990`, `bernard1991`, `felder1992`. Searches and available primary texts are documented in §9; missing exact source checks are marked unclear rather than reconstructed from later citations.
- **Exact priority of the reduced Toda-data formulation:** the general symmetry mechanism has clear prior art, including a rank-one affine-Toda model. No exhaustive priority theorem for the exact formulation can be inferred from this search.
- **Quantum physical-sector selection and the later missing-charge conjecture:** only the regular static A3 obstruction is independently checked. Thimble intersections, the full bound-state charge enumeration, asymptotic completeness and state selection are outside the established evidence here.

## 13. Other observations and evidence handoff

**Unverified:** the position sentence at 10:59 should probably specify the rest-frame centre; its printed exponent suggests lab intercept `-Re xi/(m_a cosh theta)`. P1 already recorded this minor issue; it is not counted as a verified correction here. Likewise, the broad statement that convex real-coupling Toda has precisely the mass-matrix particles “and their bound states” is not established by convexity alone; its full quantum meaning was not audited.

**Evidence retained:** final new scripts are `codex_phase2_currents_v3.py` (with JSON), `codex_phase2_integrals_v2.py`, `codex_phase2_centralizers.py`, `codex_phase2_couplings.py`, `codex_phase2_extra.py`, and `codex_phase2_phase.py`, each with its output. The first integration and output are retained to document the resolved tail-cutoff FAILs. The old slow current solver was stopped after the exact v3 calculation had completed; it is not evidence for a verdict. Failed preliminary versions are retained only in scratch. Final evidence is archived under `briefs/referee_classical_09_11_codex_scripts/`; source PDFs and their text extracts remain in the Codex scratch directory.

All nine coordinator items have been addressed, every frozen risk-list entry is mapped, and all 17 exercises have outlines. Unresolved claims are explicitly marked above. Phase 1 remains frozen. Phase 2 complete as a referee report; no integration into the book has been attempted.
