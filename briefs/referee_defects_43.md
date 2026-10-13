# Brief: deep referee pass on chapter 43 (breather and particle transmission amplitudes)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/integrable-defects-framework-a54a84`, branch `claude/integrable-defects-framework-a54a84`).

**Rules.** Other agents are editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Write scratch scripts only under `referee_scratch/defects_43/`.
- Write the report to `briefs/referee_defects_43_report.md`.

## Context

Chapter 43 (`part8-defects/43-particle-transmission.qmd`, `sec-particle-transmission`) is the last chapter of Part VIII. It was drafted on 2026-10-12 by a Claude agent (commit f677c84). Its verification script is `textbook_code/defect_particle_transmission.py` (8 checks, about 60 s). The drafting report is `briefs/draft_ch43_report.md`, and the chapter-43 section of `briefs/part8_defects_status.md` lists its verification record and the drafter's own flags for the referee.

**Notation.** ω = 4π/β² − 1 is the coupling, as in chapter 41. Chapter 39's root of unity is written e^{2πi/h} here. The physical q is q = −e^{−iπω}. The coproduct is Δe = e⊗1 + k⊗e, Δf = f⊗k⁻¹ + 1⊗f; 𝕋: V⊗𝒲 → 𝒲⊗V.

Read these first:
- `briefs/part8_defects_conventions.md`;
- ch. 37 (`sec-sg-breather-transmission`, `eq-shg-linear-transmission`);
- ch. 39 (`eq-an-linear-transmission`);
- ch. 41 (`eq-dt-tv`, `eq-dt-tvbar`, `eq-dt-g`, `eq-dt-crossing`, `prp-a2-unstable`; helpers at the end of `textbook_code/defect_common.py`);
- ch. 42 (`prp-diagonal-borel`, `sec-jump-coideal`);
- the book's a_r S-matrix (`eq-an-smatrix`; find it by grep);
- source texts in `~/referee-copy/papers/*.txt`: 0705.1066 (cz2007 §7), 0712.4292 (Bajnok–Simon), 1902.10690 (Penati–Polvara), and 1703.03562 (Jiang) if present.

To find an object: `grep -rn '{#<id>}' part*/ index.qmd`.

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

## Items (most important first)

1. **`prp-lightest-breather-universal` [Numerical, 2 ≤ r ≤ 4].**
   - Build the lowest breather B^(a)_1 yourself from the chapter-41 soliton matrices and scalar factors: the U_q singlet of species h−a at θ−iū and species a at θ+iū, with ū = π/2 − π/(hω).
   - Check that it crosses the type-I₋ defect as 𝕋^cl_a(θ−θ_d) times the identity on 𝒲, independent of ω, the defect label and the offset s. Do a₂ and a₃ at least.
   - Check the claim that with the shifts reversed there is no singlet.
   - Check the identification with cz2007 (7.5) and the claimed label swap (their T^{b1}_1 is the book's particle 2).
   - Is the label justified, and is "for every r" stated only as an expectation?
2. **`prp-an-particle-transmission`.** For the factors `eq-dpt-particle` (both types), check:
   - unitarity, crossing (a ↔ ā) and Hermitian analyticity;
   - the bootstrap at both kinds of a_r fusing (c = a+b and c = a+b−h);
   - the absence of physical-strip poles for real θ_d;
   - agreement with ch. 39's linearized factors.

   Is "there is no S-factor in the transmission bootstrap" correct?
3. **Bajnok–Simon and the continuation conjecture.**
   - Check against 0712.4292 (their (16), (23) and Sec. 4.2) that their one-loop sinh-Gordon factor is exactly the O(b²) term of the exact factor with θ_d = (1 − B/2)θ_d^cl.
   - Check that `cnj-transmission-continuation` lists its assumptions (minimality, vacuum, parameter renormalization, quantum integrability), and that the a_r one-loop prediction is a falsifiable consequence of the conjecture rather than a claim.
   - Check the conditional duality remark.
4. **`prp-dpt-vacuum` [Theorem].** Prove or refute: at real coupling with σ > 0 the defect potential ℬ_d (book normalization, `eq-defect-potential`) is ≥ 0, with equality only at φ = ψ = 0, for both types and all r; for σ < 0 it is unbounded below.
5. **`prp-dpt-doubled` [Numerical].**
   - Check that the doubled boundary of the a_r type-I₋ defect (a₂, a₃) keeps Δ(U_q(𝔟⁻)) and (1⊗π)Δ(U_q(𝔟⁺)) with π(e_j) = e_{j−1}, and that the opposite rotation fails.
   - Check the logical scope: are θ-independent Y and κ solved for, or fixed?
   - Check consistency with ch. 42's sine-Gordon result (`prp-diagonal-borel`).
6. **`tbl-dpt-spectrum` and `sec-defect-spectrum`.** Check every row against chapters 36, 37, 41 and 42 and their sources. Example: ch. 41 says the a₂ species-1 unstable state exists for all ω > 1, not only in a window.
7. **`sec-defect-finite-volume` and `sec-defect-open-problems`.**
   - Does each problem state accurately what is known and what is missing?
   - Are there claims of novelty for results already in the literature? Cf. items 10 and 17 of the ch. 38–40 report on Rossi–Weston, BLZ and Weston.
   - Is the Jiang paraphrase faithful? It was written from the abstract only.
   - Check the claim that a soliton CDD factor generically destroys the coupling independence. It is argued, not computed.
8. **Status labels, exercises and bibliography.**
   - Is every [Theorem] proved in the text?
   - Spot-check at least three exercises.
   - Check that every citation says what the source says.

## Report format (`briefs/referee_defects_43_report.md`)

A header with the date, the model and the commit. Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, corrected text ready to paste, with file and line.

End with other observations, labelled as unverified.
