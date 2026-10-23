<!-- Worker report, 2026-10-12. Model: Claude Opus 5.5, effort max. Scope: brief halfline_open_08_masked_bound_states.md with coordinator priorities (steps 1-2, then step 3; steps 4-5 if time) and a warning about the unexplained nu_{-j} poles of boundary sine-Gordon (12b). Interrupted once by the usage limit and resumed. Saved verbatim by the coordinator. Scripts copied to briefs/halfline_open_08_scripts/ (self-contained; mpmath only). -->

# Problem 8 report: boundary Coleman–Thun residues and the a₂ "unexplained" pole

Work done in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/p08/`. No repository file was touched.

## Verdict

- **The a₂ pole is not a bound state, so DG's spectrum stands.** The pole is p = ½+(2n−1)B/2 of K₁ on b₋ₙ,ₙ, in the window 1/(2n−1) < B < 2/(2n−1). It is a boundary Coleman–Thun pole: the double-triangle diagram D1, with its middle reflection at a zero. Three independent pieces of evidence, all checked numerically:
  1. **Sign.** −i·Res is negative throughout the window, so in a unitary theory the pole cannot be a pure bound state. The values are listed in section 2.
  2. **Threshold.** At B = 1/(2n−1) the state b₋ₙ₊₁,ₙ reaches the threshold of b₋ₙ₊₁,ₙ₋₁ plus a particle 2 at rest, and becomes a virtual state.
  3. **Diagram.** D1 can be drawn with positive lengths (is realizable) exactly in the window. Its middle reflection sits on the zero the dissolved state leaves behind.
- **The residue has an exact closed form.** Res K₁^{b₋ₙ,ₙ}(iu₁) = ½ g²Γ² K₂^{α}′(θ_R). Here α = b₋ₙ₊₁,ₙ₋₁, (i/2)g² is the residue that creates b₋ₙ,ₙ, iΓ² is the residue of S₁₁ at 2πi/3, and θ_R = i((2n−1)πB/6 − π/6) is the zero. It follows from the boundary and bulk bootstraps, and holds to working precision for every case tested.
- **The naive "K′ at the zero" rule is incomplete.** My first-principles evaluation of the asymptotic part of D1 gives −½g²Γ²K′, the opposite sign. When a zero lowers the order, configurations where a sub-triangle collapses onto the boundary contribute at the same order. Their coefficient is a local boundary amplitude, not a product of on-shell couplings, so the rule "the residue picks up K′ at the zero" does not by itself fix the residue.
- **No masked bound state found.** One formal loophole remains (section 5).

## 1. Residue formula: derivation and calibration

**Formalism.** I Fourier-transform in time and keep the position x ≤ 0, with exact large-distance propagators (sign conventions: K ≈ (i/2)g²/(θ−iv), S ≈ iΓ²/(θ−iu)).
- **Bulk line.** A line of energy E is e^{iκ|x−y|}/2κ, with κ = √(E²−m²) and Im κ > 0.
- **Reflected line.** It is K(θ_E)·e^{−iκ(x+y)}/2κ. The reflection factor depends only on the line's energy, because energy is conserved at the wall, so its off-shell continuation is unambiguous.
- **Bulk vertex.** −iC, with C² = 8m_a²m_b² sin²u·Γ².
- **Boundary emission/absorption vertex.** γ = κg. A factor ½ appears when the boundary state is unchanged; see the Lee–Yang test below.
- **Boundary propagator.** i/(ω−e_γ).
- **Amplitude.** δK = 𝒜/(2κ₀).
- **Where poles come from.** Poles arise from the x-integrals ∫₋∞⁰ e^{ifx}dx = 1/(if), where f = 0 is momentum conservation at a vertex, pinched against boundary propagators in the loop-energy integral.

**Test 1: trees reproduce the bootstrap signs.** The s-channel tree gives (i/2)g². The u-channel (crossed) tree gives −(i/2)g²·S(2iv). This matches the bootstrap value Res K^β(iv) = S(0)·S(2iv)·Res K^α.

**Test 2: boundary Lee–Yang (`p08_calib_ly.py`).**
- At θ = iπ/2 the s- and u-channel trees add, because 1 − S(iπ) = 2, so the s-channel-only coupling satisfies g_s² = g²/2.
- With that factor, the one-point pole at iπ/6 comes out as |Res| = |Γg|/2 = 1.26794919, which is Ghoshal–Zamolodchikov's (3.41).
- Without the doubling it would be |Γg|/√2 = 1.793. This also explains Ghoshal–Zamolodchikov's normalization of g.
- The test fixes magnitudes only; the relative sign of g and Γ is a convention.

**Double-triangle D1, general formula.** The boundary β emits particle c at rapidity −iv and drops to α. The incoming particle a at iu fuses with c into d. Particle d reflects off α at θ_d, splits back into a and c, and c is reabsorbed.
- **Without a zero:** a double pole with coefficient −½ g²Γ² K_d^α(θ_d)·κ_d/κ_a.
- **With a zero:** the momentum factors cancel and the naive residue is −½ g²Γ² K_d′(θ_d).
- **Why a zero is always there:** whenever β = α + c(iv), K_d^α(θ) contains the factor 1/K_c^α(iπ/3 − θ), which vanishes at θ_d. So diagrams of this type always come with a zero.

**Caveat.** The overall sign of this loop evaluation has not been cross-checked against an independently known loop residue. Only the trees and the Lee–Yang magnitude are calibrated. The structural statement does not depend on that sign: near-wall (collapsed) regions add (regular)×1/f terms, which become simple poles after the pinch, so they are of the same order as the zero-reduced term.

**Sine-Gordon calibration and the ν₋ⱼ poles.** Not done. The analogous hypothesis is untested: a bound-state pole that crosses θ = 0 at threshold leaves a zero, and a zero-reflection diagram uses it. The residue-sign test used here could decide it.

## 2. The a₂ test, with numbers (B = 9/20 unless stated)

**Sign of −i·Res K at each strip pole** (`p08_signs.py`; poles at θ = iπp/3):
- **Every creation pole of a state in DG's spectrum is positive:** 11/40 and 51/40 on b₀,₀; 33/40 and 29/40 on b₋₁,₁; 3/8 on b₋₂,₂; 11/40 on b₁,₀ (for K₂) and on b₀,₁.
- **The pole 47/40 of K₁ on b₋₂,₂ is −1.18749.**
- **The other non-creation poles have either sign:**
  - crossed-channel poles: −2.603 on b₋₁,₁ and +0.302 on b₋₂,₂;
  - the DG fig. 12c pole on b₃,₋₃: +0.059.

**Scan of the pole p₁ across B** (`p08_scan.py`). The residue is positive below the window and changes sign exactly at B = 1/(2n−1). There p₁ collides with the crossed pole at 3/2 − (2n−1)B/2 and the pole is double; at the grid point B = 1/5 (n = 3) the scan prints order 2. Inside the window:

| n | window | −i·Res inside the window |
|---|---|---|
| 2 | 1/3 < B < 2/3 | −11.59 at B = 7/20, rising to −0.119 at B = 13/20 |
| 3 | 1/5 < B < 2/5 | −4.71 to −0.243 |
| 4 | 1/7 < B < 2/7 | −12.9 to −0.235 |

**Threshold mechanism** (`p08_threshold.py`).
- The creation pole of b₋ₙ₊₁,ₙ in K₂^{b₋ₙ₊₁,ₙ₋₁} sits at x = ½ − (2n−1)B/2.
- Its residue goes to zero linearly as B approaches 1/(2n−1) from below. For n = 2, −i·Res is 0.283, 0.076 and 0.0156 at B = 7/30, 37/120 and 197/600.
- Above B = 1/(2n−1) the strip contains a simple zero at the mirror point, which is θ_R.
- At B = 9/20 the would-be state b₋₁,₂ has energy 1.2167 m₁. That is below the threshold e(b₋₁,₁) + m₁ = 1.2334 m₁, but only on the unphysical sheet, so it is a virtual state.

**Diagram D1 geometry.** The particle-2 angle is π/3 − v = (2n−1)πB/6 − π/6. It is positive, so the triangles close with positive lengths, if and only if B > 1/(2n−1). The incoming angle u₁ < π/2 gives the upper edge B < 2/(2n−1). No internal lines cross. Counting gives order 2, lowered to 1 by the zero.

**Exact-residue identity** (`p08_identity.py`). Exact residue / (½g²Γ²K′) = 1 to working precision for:
- n = 2 at B = 9/20, 2/5, 3/5;
- n = 3 at B = 3/10 and 3/20;
- n = 4 at B = 1/5.

At B = 9/20: g² = 3.36623, Γ² = 0.319032, Res = −1.187489724 i.

## 3. Coincident diagrams (partial)

At n = 2, B = 1/4 = 1/(2k), the pole p₁ = 7/8 coincides with the pole that creates b₋₃,₃.
- The pole stays simple.
- Its residue is twice the D1 structure (ratio 2.000000), with −i·Res = +5.86. This is a sum over a degenerate pair of states, consistent with the level crossing.

The a₄ case n₂ = −n₃ and a₃, b = 2 were not started.

## 4. Masked bound states

None found. The remaining loophole: a state at e(b₋ₙ,ₙ) + m₁cos u₁ that couples to b₋ₙ,ₙ + particle 1 but not to b₋ₙ₊₁,ₙ₋₁ + particle 2. Such a coupling must vanish because K₂ has a zero there. That state would have to appear at B = 1/(2n−1), hidden under a coupling-dependent CT contribution of opposite sign. Residues alone cannot exclude it, but continuity and the conserved charges, which match the virtual state's, make it implausible.

## 5. Proposed text

**Replacement for "This range leaves one pole unexplained…":**

> In the window $\frac1{2n-1}<B<\frac2{2n-1}$ the pole $p=\frac12+(2n-1)\frac B2$ of $K_1$ on $b_{-n,n}$ is not a bound state. At $B=\frac1{2n-1}$ the state $b_{-n+1,n}$ reaches the threshold of $b_{-n+1,n-1}$ plus a particle 2 at rest. Its creation pole in $K_2^{b_{-n+1,n-1}}$, at $p=\frac12-(2n-1)\frac B2$, passes through $\theta=0$ with vanishing residue, and above the threshold the strip contains instead a zero at $\theta_R=\frac{i\pi}3\big((2n-1)\frac B2-\frac12\big)$, a virtual state. Exactly in the window there is a boundary Coleman–Thun diagram. $b_{-n,n}$ emits a particle 1 and drops to $b_{-n+1,n-1}$. The incoming particle 1 fuses with it into a particle 2, which reflects at $\theta_R$, the zero, and splits back, and the emitted particle is reabsorbed. The diagram has two loops and would give a double pole; the zero makes it simple. The residue agrees: $-i\operatorname{Res}$ is positive below the window, changes sign at $B=\frac1{2n-1}$, where the pole meets the pole at $\frac32-(2n-1)\frac B2$, and is negative in the window. At $B=\frac9{20}$ the pole $\frac{47}{40}$ on $b_{-2,2}$ has $-i\operatorname{Res}K_1=-1.187$, impossible for a bound state in a unitary theory. The spectrum of Delius and Gandenberger therefore stands.

**Replacement for open problem 8:**

> **Masked bound states.** In the bulk of non-simply-laced theories Coleman–Thun poles can mask genuine bound states [DG sec. 7]. On the half-line the leading coefficient of an on-shell diagram is a product of couplings $g$, $\Gamma$, reflection and scattering factors and a geometric factor; this reproduces Ghoshal and Zamolodchikov's one-point residue in boundary Lee–Yang. For the diagrams that are simple only because a reflection happens at a zero, which dominate @sec-boundary-spectrum-a2, the naive rule ($K\to K'$ at the zero) fails: configurations in which a sub-diagram collapses onto the boundary contribute at the same order, through a local boundary amplitude. The exact residues have the bootstrap form $\frac12g^2\Gamma^2K'$, but residues alone cannot exclude a masked state. For $a_2^{(1)}$, residue signs and threshold behaviour exclude one. A general residue rule, the $a_3^{(1)}$ ($b=2$) and $a_4^{(1)}$ ($n_2=-n_3$) cases, and the unexplained $\nu_{-j}$ poles of boundary sine-Gordon remain open.

## Status

- **Done:** step 3 (the a₂ verdict); the leading-order formula; the Lee–Yang calibration; the exact-residue identity.
- **Partial:**
  - step 2, the rule at a zero: shown to be incomplete, with the local term not computed;
  - step 4: only the B = 1/4 coincidence.
- **Not done:**
  - the sine-Gordon calibration and the ν₋ⱼ poles;
  - a₄ and a₃ (b = 2);
  - Riva's Eₙ cases;
  - an independent cross-check of the sign of the D1 loop evaluation.

## Scripts

All in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/p08/`. Each runs in seconds.
- `p08_common.py`: block products with exact Laurent coefficients, DG (6.7), energies.
- `p08_signs.py`: residue-sign table at B = 9/20.
- `p08_scan.py`: sign scan of the pole p₁ across B, for n = 1–4.
- `p08_identity.py`: the identity Res = ½g²Γ²K₂′(θ_R) and the order of the zero.
- `p08_calib_ly.py`: Lee–Yang normalization test.
- `p08_threshold.py`: threshold crossing of b₋ₙ₊₁,ₙ.

## Integration (branch half-line-open-problems-integration)

Basis: the worker report, the Gemini verdict (tier A, which gives an analytic proof of the residue identity), and the coordinator's checks. New section 9 in `code/particle_reflection.py` uses the book's own `K1nm`, `S`, `Prod`; all 47 checks PASS.

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 1, 2 (residue signs; sign change at $B=\frac1{2n-1}$) | **Verified** | Section 9: $-i\,\mathrm{Res}>0$ below and $<0$ inside the window for $n=2,3,4$; $-1.18749$ at $B=\frac9{20}$ (Gemini's contour integral agrees). The sign convention $K\simeq\frac i2g^2/(\theta-iv)$ was confirmed in ch. 28. | ch. 35 paragraph rewritten |
| 3 (threshold, virtual state) | **Accepted** (worker and Gemini, independent code) | not in the book script | ch. 35 bullet |
| 4 (D1 realizable in the window; order $2\to1$) | **Accepted** (geometry re-derived by Gemini) | | ch. 35 bullet |
| 5 (exact identity $\mathrm{Res}=\frac12g^2\Gamma^2K_2'(\theta_R)$) | **Verified**, and it is a consequence of the bootstraps (Gemini's proof; the coordinator checked the zero-derivative step $2i/g^2$) | Section 9 at generic $B$, $\le3\times10^{-21}$. At $B=\frac12$, a level crossing, the formula degenerates; the first sampling hit it. | ch. 35 |
| 6 (Lee–Yang calibration) | **Not integrated** | | — |
| 7 ("naive rule incomplete; collapsed configurations") | **Not integrated as a mechanism** (Gemini: unverified hypothesis; the sign of the worker's loop evaluation is not cross-checked). The book says only that the naive rule has not been derived. | | open problem 8 |
| 8 (coincidence at $B=\frac14$) | **Not integrated** | | — |
| 9 (DG spectrum stands) | **Integrated** | | ch. 35; Appendix D DG §6.2 entry updated (still *Incomplete*: their attribution of the pole) |
| 10 (text) | Integrated as modified (Gemini's refined open-problem text, shortened) | | |
