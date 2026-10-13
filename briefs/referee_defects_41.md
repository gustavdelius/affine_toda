# Brief: deep referee pass on chapter 41 (soliton transmission matrices)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/integrable-defects-framework-a54a84`, branch `claude/integrable-defects-framework-a54a84`).

**Rules.** Other agents are editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Write scratch scripts only under `referee_scratch/defects_41/`.
- Write the report to `briefs/referee_defects_41_report.md`.

## Context

Chapter 41 (`part8-defects/41-soliton-transmission.qmd`, `sec-soliton-transmission`) was drafted on 2026-10-12 by Claude agents (commit 4258b88). Its verification script is `textbook_code/defect_an_transmission.py` (21 checks, 100–180 s). It uses helpers appended to `textbook_code/defect_common.py` after the "ch. 41" marker.

Documents:
- the drafting report, `briefs/draft_ch41_report.md`;
- the exploration checkpoint, `briefs/draft_ch41_checkpoint.md`;
- the chapter-41 section of `briefs/part8_defects_status.md`.

**Convention warning.** In this chapter ω = 4π/β² − 1 is the coupling (Part VI's ω, Corrigan–Zambon's γ), and X = e^{ω(θ−θ_d−iπ/2)} is CZ's x̂. In chapter 39, ω is a root of unity. The physical q is q = −e^{−iπω}, and y = e^{ωθ}. The coproduct is Δe = e⊗1 + k⊗e, Δf = f⊗k⁻¹ + 1⊗f; 𝕋: V⊗𝒲 → 𝒲⊗V.

Read these first:
- `briefs/part8_defects_conventions.md`;
- ch. 36 (`eq-triangle-braid`, `eq-defect-unitarity`, `eq-defect-crossing`, `eq-defect-bootstrap`);
- ch. 37 (the sine-Gordon template and its q-sign paragraph);
- ch. 38 (`eq-dc-ar-delay`);
- ch. 39 (`eq-ds-class-shift`, `prp-ds-absorption`);
- ch. 40 (`sec-borel`, `thm-t-triangle`);
- Part VI ch. 21 (`sec-an-example`, `sec-unitarity-crossing`, `sec-scalar-factor`);
- the source texts in `~/referee-copy/papers/*.txt`: 0705.1066 (cz2007), 0902.1307 (cz2009a), 1012.4186 (cz2011), 1902.10690 (Penati–Polvara).

To find an object: `grep -rn '{#<id>}' part*/ index.qmd`.

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

## Items (most important first)

1. **The a_r Weyl module and the intertwiners (`eq-dt-weyl`, `eq-dt-tv`, `eq-dt-tvbar`, `eq-dt-tplus`, `prp-an-transmission`).**
   - Build the U_q(𝔟⁻) module yourself and check the q-Serre relations, including the claim that they hold iff the commutation factor is q^{±1}.
   - Solve the intertwining equation V⊗𝒲 → 𝒲⊗V and V̄⊗𝒲 → 𝒲⊗V̄ for a₂ and a₃. Check uniqueness and the closed forms.
   - Check the gauge match with cz2009a (5.8) and (5.11), working in their conventions.
   - Is the triangle relation with the Part VI Ř satisfied?
   - Does "orientation −1 = the charge-conjugate defect (cz2007's Dᵀ)" hold?
2. **Crossing and the q-sign (`eq-dt-crossing`, `prp-transmission-crossing`).**
   - Check C = 1, c_I = 1 + (−q)^r X^h, c_II = 1 − qX^h, and their independence of the offset s.
   - Check that cz2009a (6.8) follows from them.
   - Check the q-sign claims at −q:
     - for h odd, no conjugation exists at the physical crossing point, so crossing fails with both C and C′ = diag((−1)^i);
     - for h even, it holds with C′.
   - Is the statement properly convention-qualified? Could a Klein-type redefinition relate ±q for h odd, as it does for sine-Gordon?
   - Check the claim that the matrices are not unitary for r ≥ 2, and what replaces unitarity.
3. **The scalar factor (`eq-dt-g`) and the correction to cz2009a (6.10).**
   - Derive the minimal solution of the functional equations yourself.
   - Check that cz2009a (6.10) as printed pairs Gamma functions across different k, so that its partial products drift like K^{−z}.
   - Check that the book's equal-k pairing converges and solves the equations.
   - Check the bootstrap ḡ = ∏ g (`eq-dt-bootstrap`).
4. **Poles, windows and the kernel bootstrap (`eq-dt-pole`, `eq-dt-window`, `prp-a2-unstable`).**
   - Check that the block determinant is q^r(1 + (−q^{−1})^r X^h), so that the residue has rank h−1 per sector.
   - Check the window max(0, a − h/2) ≤ 1/ω < a against cz2009a (6.13)–(6.15). Check the correction of their text ("1/γ < h/2 − 1 or (h−1)/2").
   - Check the species relabelling against cz2007: their species 1 is the book's species 2, and the book's species 1 is unstable for ω > 1.
   - Check the kernel bootstrap: on the kernel of M_V(θ₁) the transmission becomes the class-1 module with θ_d → θ_d + 2πi/h.
   - Check the claim that classical real absorptions become resonances at θ_d − iπ/(hω).
5. **The selection rule and the semiclassical test.**
   - Check the tunnelling suppression min(|X|^j, |X|^{j−h+1}) in M_V̄, and its classical limit against ch. 39's `prp-stored-charge`.
   - Check the Wigner test against `eq-dc-ar-delay` (`prp-transmission-semiclassical`), including the prefactor 2π/(h(ω+1)).
6. **a₂^(2) (`prp-a22-transmission`).**
   - Check the conversion of cz2011's algebra to the Part VI form.
   - Check the module, the uniqueness, the triangle relation, and crossing at x_c = −q^{−6} with the antidiagonal C.
   - Check the claim that the q-sign is seen only relative to C.
7. **cz2007 (6.14).** Does the corrected (−q)^{+1/2} reproduce the sine-Gordon matrix quoted before it?
8. **The Penati–Polvara caveat and the status labels.**
   - Is every a_r claim correctly labelled: established as a solution of the axioms, conjectural as the Lagrangian defect's transmission matrix?
   - Is "[Established; numerical re-derivation]" justified for `prp-an-transmission`?
   - Spot-check at least three exercises.
   - Check that every citation says what the source says.

## Report format (`briefs/referee_defects_41_report.md`)

A header with the date, the model and the commit. Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, corrected text ready to paste, with file and line.

End with other observations, labelled as unverified.
