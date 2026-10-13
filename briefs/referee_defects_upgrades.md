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

Commits: ch. 39 integration ad9342c; ch. 41 integration df973fb; ch. 41 triangle relation 9441742; ch. 38 a₃ spin-two script e748b4f. Use `git show <commit> -- <file>` to see exactly what changed.

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

## Report format (`briefs/referee_defects_upgrades_report.md`)

A header with the date, the model and the commit. Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, corrected text ready to paste, with file and line.

End with other observations, labelled as unverified.
