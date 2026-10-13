# Brief: deep referee pass on chapter 42 (defects meet boundaries)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/integrable-defects-framework-a54a84`, branch `claude/integrable-defects-framework-a54a84`).

**Rules.** Other agents are editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Write scratch scripts only under `referee_scratch/defects_42/`.
- Write the report to `briefs/referee_defects_42_report.md`.

## Context

Chapter 42 (`part8-defects/42-defects-and-boundaries.qmd`, `sec-defect-boundary`) was drafted on 2026-10-12 by a Claude agent (commit 581bed0). Its verification script is `textbook_code/defect_doubling.py` (17 checks, about 30 s). The drafting report is `briefs/draft_ch42_report.md`.

The chapter treats a defect as a boundary of the doubled theory (the doubling trick of ch. 36, `eq-doubled-k`, `prp-doubling`) and asks which quantum-group symmetry that boundary has. It also fuses defects onto boundaries, and relates the results to tensor K-matrices and non-invertible boundary scattering.

Read these first:
- `briefs/part8_defects_conventions.md`: the conventions, verified by script, and the source dictionary. The coproduct is Δe = e⊗1 + k⊗e, Δf = f⊗k⁻¹ + 1⊗f. 𝕋: V⊗𝒲 → 𝒲⊗V.
- The chapter-42 section of `briefs/part8_defects_status.md`, which includes a scope flag for item 6.
- Ch. 36 `sec-doubling-trick`, ch. 37 `eq-sg-transmission-matrix` and `eq-sgd-scalar`, ch. 40 `sec-borel`, `prp-two-borels` and `tbl-boundary-defect-algebra`, ch. 34 (`part7-half-line/34-quantum-symmetric-pairs.qmd`), and ch. 29–30 (`eq-bsg-k-full`, `eq-bsg-oscillator`, `sec-dynamical-boundary`).
- Source texts in `~/referee-copy/papers/*.txt`, notably:
  - 1202.6016 (cz2012);
  - bajnok_simon2008;
  - bajnok_george2006;
  - 1610.09271 (bao_wang2018b; find it by grep if the file name differs);
  - 2402.16676 (Appel–Vlaar 2024);
  - 2504.08375 (Shimamori–Yamaguchi 2025).

To find an object: `grep -rn '{#<id>}' part*/ index.qmd`.

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_36_37_deep_report.md` shows the expected depth.

## Items (most important first)

1. **The refutation of `cnj-twisted-diagonal`, and its logical scope.** The chapter claims that no twisted diagonal (1⊗ψ)Δ(U), with ψ a torus twist or τ∘torus, is a symmetry of the doubled boundary of a sine-Gordon type-I₋ or I₊ defect, "for any θ-independent action on 𝒲" (residuals 0.93–0.97, stable on windows N = 5, 7, 9).
   - Does the test really range over all θ-independent actions on 𝒲, i.e. are Y = χ(x) and κ solved for as unknowns? Or does it fix them?
   - Is restricting to torus and τ∘torus twists justified, i.e. are these all the relevant ψ?
   - Could a truncation artefact (window edges, interior masks) produce a residual near 1?
   - Build your own doubled K from `eq-sg-transmission-matrix` and test independently.
   - Is the label "[Conjecture; refuted as stated by the test in …]" accurate?
2. **`prp-diagonal-borel` [Theorem (1), (4); Numerical (2), (3), (5)].**
   - Check (1): Δ(U_q(𝔟^∓)) is a symmetry, with 𝒲 the Weyl module.
   - Check (2)–(3): the second module ρ̃ (`eq-dbd-second-module`) satisfies the U_q(𝔟^±) relations, and (1⊗τ)Δ(U_q(𝔟^±)) holds on it.
   - Check (4): the obstruction. Any coideal containing both halves contains e_j⊗1, which fails.
   - Check (5): the θ_d → ±∞ limits restore the untwisted and the τ-twisted diagonal respectively, with decay like e^{−λ|θ_d|}.
   - Are the [Theorem] parts proved in the text? Does the claim refine ch. 40's statement (in `sec-perturbed-defects` / `prp-two-borels`) consistently, or contradict it?
3. **`prp-defect-fusion-k`.**
   - Check that 𝕋̄(1⊗K_D)𝕋 for a type-I₋ defect fused onto Dirichlet equals the Ghoshal–Zamolodchikov K-matrix with (η, ϑ) = (η_D, λθ_d), sector by sector. Check also that ϑ = −λθ_d and ϑ = θ_d fail.
   - Check the scalar-factor identity: the chapter says no CDD factor appears; the script agrees to 3e-9.
   - Check the infinite-dimensional K_f against cz2012 (2.18)–(2.19). Work in their own conventions first, then translate.
   - Check the sinh-Gordon identity of Bajnok–Simon.
   - Check the consequence λθ_d = (λ+1)θ_d^cl, and its consistency with ch. 37's `cnj-defect-parameter`.
4. **`prp-reflected-copy` and `prp-diagonal-coideal`.**
   - Check that the reflected copy is acted on through Δ^op, so that the doubled bulk algebra is U⊗U^cop.
   - Check the claim V(1/y)∘ω_q = V̄(y).
   - Check that the largest left coideal in the commutant of the doubled K of an invertible line is exactly the twisted diagonal.
   - Check the map ω_q⊗id from the book's left coideal to Bao–Wang's right coideal (ω⊗1)Δ(U).
   - Check the source claims in the Notes: Bao–Wang Rem. 4.10(2) is a finite-type remark; Kolb 2014 and Appel–Vlaar 2022 assume an indecomposable Cartan matrix; Kolb–Yakimov's "diagonal type" is unrelated.
5. **The RTT algebra inside the doubled reflection-equation algebra.** The chapter says the RTT algebra of ch. 40 is the off-diagonal part of the doubled reflection-equation algebra, with 𝕃⁻¹ adjoined. Is this precise, and correct?
6. **`prp-fused-bdg` and the classical dressing (`eq-dbd-dressed-k`).**
   - Check by your own finite-difference computation that K_φ(λ) = [L_d(1/λ)⁻¹]ᵀ K_ψ(λ) L_d(λ)⁻¹ satisfies the classical K-matrix equation on the fused system, for a₁ and a₂.
   - Check that the sine-Gordon defect potential equals BDG's H_osc at the stated complex (p̂, q̂).
   - Check that the fused system is not the BDG boundary.
   - **Scope:** the chapter derives new classical dynamical a_n boundaries from the dressing. Is the claim stated modestly, with its status label matching what is proved (conservation only; involution and quantization not examined)?
7. **`sec-tensor-k` and `sec-noninvertible-boundaries`.**
   - Check the Appel–Vlaar citations (Thm 5.6.1, 6.4.2, 6.5.2, Rem. 6.5.3) against the source.
   - Check the claim that Shimamori–Yamaguchi's modified boundary crossing (their (3.29), `eq-dbd-sy-crossing`) becomes the usual Ghoshal–Zamolodchikov crossing after the kink renormalization (`exr-dbd-boundary-crossing`). Work it out yourself.
8. **Status labels, exercises and bibliography.**
   - Is every [Theorem] backed by a proof in the text?
   - Spot-check at least three exercises.
   - Check that every citation says what the source says.

## Report format (`briefs/referee_defects_42_report.md`)

A header with the date, the model and the commit. Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, corrected text ready to paste, with file and line.

End with other observations, labelled as unverified.
