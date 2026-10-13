# Brief: referee pass on the Part VIII statements upgraded after their referee pass

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/integrable-defects-framework-a54a84`, branch `claude/integrable-defects-framework-a54a84`).

**Rules.** Other agents are editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Write scratch scripts only under `referee_scratch/defects_upgrades/`.
- Write the report to `briefs/referee_defects_upgrades_report.md`.

## Context

Chapters 38, 39 and 41 were refereed by you (Codex) on 2026-10-12/13; the reports are `briefs/referee_defects_38_40_report.md`, `briefs/referee_defects_39_report.md` and `briefs/referee_defects_41_report.md`, and each ends with an "Integration" section that says what was changed. During integration, Claude agents **upgraded some labels to [Theorem] with new proofs that no referee has seen**. This pass checks only those new proofs and the few new claims listed below. Do not re-referee what your earlier reports already covered.

Commits: ch. 39 integration ad9342c; ch. 41 integration df973fb; ch. 41 triangle relation 9441742; ch. 38 a₃ spin-two script e748b4f; ch. 42 integration 90061c8; ch. 43 integration fd38d74. Use `git show <commit> -- <file>` to see exactly what changed.

Items 8–11 were added on 2026-10-13, after the chapter-42 and chapter-43 integrations. Chapter 43 needs extra care. You ran out of quota during its referee pass. Its report, `briefs/referee_defects_43_report.md`, has items 1–4 by you, and items 5–8 written by Claude from your handoff file. So the claims of items 5–8 have had no independent referee; item 11 lists them.

**Conventions** (as in the ch. 41 brief): ω = 4π/β² − 1, q = −e^{−iπω}, X = e^{ω(θ−θ_d−iπ/2)}, coproduct Δe = e⊗1 + k⊗e, Δf = f⊗k⁻¹ + 1⊗f; 𝕋: V⊗𝒲 → 𝒲⊗V. Read `briefs/part8_defects_conventions.md` §§10–11 first.

To find an object: `grep -rn '{#<id>}' part*/ index.qmd`.

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts and the integration scripts (`briefs/referee_defects_{39,41}_scripts/integration_*.py`, `textbook_code/defect_spin_two_a3.py`) for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear".

## Items (most important first)

1. **All-rank uniqueness, `prp-an-transmission` item 3 (ch. 41), now [Theorem].**
   - Check the proof line by line. The key step is that Ω_j = Δ(f_j)Δ(f_{j+1}) − q⁻¹Δ(f_{j+1})Δ(f_j) has no term acting on 𝒲 alone, and that on every charge sector of V⊗𝒲 its image is the line of v_{j+2} (for V̄⊗𝒲: its kernel is the hyperplane i ≠ j+2). From this, any charge-conserving endomorphism that commutes with the Borel subalgebra is diagonal, and then a scalar.
   - Check the type-I₊ version, with Δ(e_j)Δ(e_{j+1}) − qΔ(e_{j+1})Δ(e_j).
   - Check the hypotheses: the sector determinants are nonzero and independent of the sector and of s; the charge sectors are h-dimensional for every offset s.
   - Check the claim that q is a root of unity exactly when ω is rational, and the q⁵ = 1 counterexample (ω = 3/5) to uniqueness without charge conservation.
   - Is the step "E = (T⁰)⁻¹T commutes with the Borel action" valid with 𝒲 infinite-dimensional, given that T⁰ is invertible block by block?
2. **Type-I₊ intertwiner for V, `prp-an-transmission` item 2 (ch. 41), now [Theorem].**
   - Check `eq-dt-reversed`, the reversed-orientation two-entry matrix, as a U_q(𝔟⁻)-intertwiner for every r. Is it proved for every r, or only checked for h = 3–6?
   - Check the Chevalley-involution argument, (ϖ⊗ϖ)Δ = Δ^op ϖ, and the relabelling N + s → −(N + s), which together turn the inverse of `eq-dt-reversed`, with X → −X₊, into the type-I₊ intertwiner for V.
   - Check the claimed normalization: with the same μ the reversed module gives CZ2007 (5.10) with x̂ = qX; keeping x̂ = X needs μ/q.
3. **The selection rule, `prp-an-transmission` item 5 (ch. 41), now [Theorem].** It is read off from the closed forms. Check the statement as written, including the bound min(|X|^j, |X|^{j−h+1}) and the type-I₊ exchange of V and V̄.
4. **`prp-an-linear-transmission` (ch. 39), upgraded from [Numerical] to [Theorem].** The proof is an exact sympy elimination, `textbook_code/defect_solitons_an.py`. Check that the elimination proves the statement for every r and both types, and is not a check at sample values.
5. **New corrections to the literature from the ch. 41 integration** (`appendices/D-corrections.qmd`, entry "Corrigan and Zambon 2011", and the extended "Corrigan and Zambon 2009a" entries):
   - CZ2011 (3.12) as printed pairs Gamma functions at different k and drifts like K^{−2z}; grouped at equal k it solves their (3.8) and (3.10) with z = i(θ−ν)/ξ. Hence the slip is in the printed pole positions, not in (3.12). Source: `~/referee-copy/papers/1012.4186.txt` and the PDF.
   - The CZ2009a (6.10) entry now states the normalization at z = 0.

6. **The triangle relation for every rank, `prp-an-transmission` item 4 (ch. 41), upgraded from [Numerical] to [Theorem] in commit 9441742.** The proof evaluates the universal R-matrix of U_q(a_r^(1)) on V(y)⊗𝒲 as a formal series, uses item 3 over the field of formal Laurent series to identify it with a scalar series times the displayed matrix, and gets the triangle relation from the Yang–Baxter equation. Hypothesis (ii) of `thm-t-triangle` is no longer used. Check in particular:
   - which universal R goes with which type: ℛ ∈ U_q(𝔟⁻)⊗̂U_q(𝔟⁺) for type I₊, and ℛ₂₁⁻¹ for type I₋. Does this agree with the Part VI coproduct, with ch. 40 (`prp-two-borels`, `sec-universal-r-affine`) and with the ch. 41 matrices? Do the terms of order zero and one of the evaluated R-matrix reproduce the matrix parts exactly, including X = −(q−q⁻¹)μ y/y_d and X₊ = (q−q⁻¹)μ y_d/y?
   - that the intertwining identity and the Yang–Baxter equation are valid when 𝒲 is only a Borel module (all third tensor factors lie in that Borel subalgebra);
   - the formal-series bookkeeping: the variables u, uv, v, the finiteness of every coefficient, and invertibility of the scalar ψ_ab(u)φ_a(uv)φ_b(v) in ℂ((u))((v));
   - the transfer of uniqueness (item 3, and the Part VI uniqueness of Ř) to formal Laurent series;
   - the dropped factor q^{c⊗d+d⊗c} at level zero, the treatment of the non-integral offset s in the Cartan factor q^t, and the extension to roots of unity by rationality;
   - whether the cited sources ([@khoroshkin_tolstoy1991], [@tolstoy_khoroshkin1992], [@frenkel_hernandez2015]) state what the proof needs. The proof gives no pinpoint references.

7. **`prp-spin-two` item (3) for a₃^(1) (ch. 38), upgraded from [Numerical] to [Theorem] in commit e748b4f.** The new script `textbook_code/defect_spin_two_a3.py` (9 checks, about 3 min) builds the a_r root data exactly, takes the bulk spin-2 tensor A_abc = 2 Σ_k μ_ka μ_kb μ_kc with B = 3(1 − D), and constructs 𝓑_{±2} by the homotopy formula of check 11 of `textbook_code/defect_classical_an.py`. It then finds defect flux − d𝓑_{±2}/dt = 0 exactly in sympy, for both types, both chiralities and symbolic σ. Check in particular:
   - that the flux expression (`boundary_G`) is the correct defect flux of the spin-±2 charges for a₃: the use of the defect conditions to replace x-derivatives, the second derivative d_ε²φ via the field equation, and the sign conventions for spin (the book's spin s goes with e^{sθ}; see `sec-qg-gradation`);
   - that B = 3(1 − D), together with the A_abc formula, gives the bulk spin-2 current of a₃ in the normalization the flux uses. The script's check 4 tests conservation on shell, but is the scale of B relative to A, which the docstring ties to a factor of 3 between check 5 and check 11, consistent?
   - that the exact zero is a real identity and not an artefact, e.g. of simplification or of the homotopy construction producing G by design. The control with a non-integrable D fails; is that control sensitive enough?
   - whether the label change in `prp-spin-two` is now justified as stated.

8. **`prp-diagonal-borel` (ch. 42), items (2) for type I₋, (3) and (5) upgraded to [Theorem] in 90061c8; new proof of (4) for type I₊.** All four rest on comparing powers of y in one generator equation, with 𝕃 = P𝕋 = 𝒟 + y𝒪 (type I₋) or 𝒟₊ + y⁻¹𝒪₊ (type I₊), where the blocks are θ-independent. Check in particular:
   - (3): that the y² coefficient (untwisted e₁) and the y⁻¹ coefficient (τ-twisted f₁) are contradictory for **arbitrary** θ-independent operators Y, κ on 𝒲, with a torus factor t kept explicit; and the type-I₊ analogues. Does [@kolb2014, Thm. 3.2] say that every Hopf automorphism of the algebra in question is a torus rescaling, possibly composed with τ? Is that algebra the one the book's U is (derived algebra, generic q, level-zero modules)?
   - (5): the forced operators χ(Δk₁) = t⁻¹AD⁻¹ and χ(Δk₁⁻¹) = t⁻¹D⁻¹A in the diagonal limit, the charge-shifting pair in the off-diagonal limit, and the exchange of the two directions for type I₊. Is "no twisted diagonal acts on all of 𝒲" proved for every Hopf twist?
   - (2), type I₋: is the proof block complete, i.e. does it check both e-generators, both k-equations and all relations of ρ̃ (`eq-dbd-second-module`)?
   - (4), type I₊: conjugation by Δ(k₁) separates the two terms of f₁⊗k₀⁻¹ + 1⊗f₀, which isolates 1⊗f₀, whose y⁻² coefficient is uncancellable.
   - Comparison scripts: `briefs/referee_defects_42_scripts/integration_obstructions.py`, `textbook_code/defect_doubling.py` (20 checks).
9. **`prp-an-particle-transmission` (ch. 43), upgraded from [Numerical, h = 2,…,8] to [Theorem] for every h in fd38d74.** Items 1–4 are now proved by hand in the proof block (rational identities in x = e^{−(θ−θ_d)}, u = e^{iγ_a}, v = e^{iγ_b}, for both fusings and both types); item 5 is `prp-an-linear-transmission`, which is item 4 above. Check the hand proof, in particular the shifts for the fusing c = a + b − h and the type-I₊ version, which the proof gives only in outline.
10. **`prp-dpt-doubled` (ch. 43): the new explicit second module `eq-dpt-second-module` (fd38d74).** The label stays [Numerical], but the proof now claims by hand:
    - on V, that the y² terms force κ to be the plain shift by α_j − α_{j−1};
    - that the y terms give Y = ρ̃(e_j) + c_j q^{(N+s)_j} times the shift;
    - that the Cartan relation forces c_j = 0, after which ρ̃ is a U_q(𝔟⁺)-module, including the q-Serre relations for every r;
    - that the h = 2 case is ch. 42's `eq-dbd-second-module` after the basis change d_n = τ^{−n}q^{n(n−1)/2+ns/2}, τ = −ie^{iπω/2};
    - that the failing twists fail by an exact y-power argument for every r, and that the two Borel halves do not lie in one coideal for q³ ≠ 1.
    
    The V̄ half is numerical only. Is the V half proved, and is there a short proof of the V̄ half?
11. **Claims in ch. 43 from the Claude-written items 5–8 of its referee report (integrated in fd38d74).** For each, an independent computation:
    - the residue of 𝕋_V̄ at the species-r pole has rank 1. At that pole X^h = q, and on a sector of fixed charge M_V̄ = q u wᵀ. Check also that the pole of ḡ is not cancelled. `tbl-dpt-spectrum` and problem 10 now say this;
    - the soliton CDD condition `eq-dpt-soliton-cdd`, ∏_{k=0}^{h−1} f(θ − iπ + 2πik/h) = 1, derived from the crossing relations of ch. 41 and `eq-dt-bootstrap`. Is it necessary and sufficient? Is the example f = 𝕋^cl_1, with its coupling-dependent breather multiplier, right?
    - the exceptions at θ_d = iv with v ∈ πℤ (Θ unbroken, 𝕋_a a phase, σ = −1 inside the window for a < h/2; for sine-Gordon only for ½ < λ < 1), in `tbl-dpt-spectrum`, "Reading the table", `exr-dpt-imaginary` and ch. 37 `sec-sg-defect-spectrum`;
    - the Lee–Yang poles are in 𝕋̄_b (table caption);
    - the TBA paragraph of `sec-defect-finite-volume` (fugacity form, F_ā = conj(F_a), |F_a| < 1, F_{h/2} = −tanh((θ−θ_d)/2), reality only for a conjugation-symmetric solution), and the rewritten Jiang bullet against 1703.03562 (not in `~/referee-copy/papers`; your download is in `referee_scratch/defects_43/sources/`);
    - the new statements in `prp-dpt-vacuum`'s proof: the explicit unbounded direction φ = ψ = tv for σ < 0, and the linear-layer estimate for the total static energy;
    - the citation and novelty changes in the Notes, Corrigan–Zambon (7.11) and (6.12).
    
    Comparison scripts: `briefs/referee_defects_43_scripts/` (Codex's, `second_module_checks.py`, and `integration_*`), `textbook_code/defect_particle_transmission.py` (11 checks).

## Report format (`briefs/referee_defects_upgrades_report.md`)

A header with the date, the model and the commit. Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, corrected text ready to paste, with file and line.

End with other observations, labelled as unverified.
