You are integrating independent referee reports into a mathematical-physics textbook (a Quarto book on affine Toda field theory). Read `AGENTS.md` at the repository root first; its rules apply throughout. Work in the git worktree

  `/home/gustav/Git/affine_toda/.claude/worktrees/part3-int-ch11`   (branch `part3-int-ch11`)

Run every command from that directory, and never `cd` to the main repository. Never use `git stash`. Never use `git add -A` or `git add .`; stage files by name. Do not push. Do not spawn subagents.

Scratch directory for your own check scripts: `/tmp/claude-1000/-home-gustav-Git-affine-toda/7fa1324b-155d-4b76-8703-4c3aaf6d6c0d/scratchpad/part3_int_B/`. Before you commit, copy the scripts the record relies on into `briefs/referee_classical_09_11_integration_scripts/`, with the prefix `integration_B_`.

## The task

Part III of the book (chapters 9–11) had a two-phase referee pass by two referees, Gemini and Codex. **You integrate everything that concerns chapter 11**, `part3-classical/11-reality.qmd` (`sec-classical`), together with what it entails elsewhere: the Part III item of Appendix E, new bibliography entries, gap H6 in Part IX, and the chapter-11 figure. A second agent is integrating chapters 9 and 10 at the same time, in another worktree, and owns `part3-classical/09-lagrangian.qmd`, `part3-classical/10-solitons.qmd`, `textbook_code/toda_classical.py` and the chapter 9–10 figure scripts. Do not edit those.

The reports, all committed unchanged in 45955bb, are in `briefs/`:
- `referee_classical_09_11_phase1_report.md` (Gemini, phase 1) and `referee_classical_09_11_report.md` (Gemini, phase 2). Reviewed commit a03edd4; line numbers match the current files.
- `referee_classical_09_11_codex_phase1_report.md` and `referee_classical_09_11_codex_report.md` (Codex). Reviewed commit bcd82b9; chapter 11 line numbers match. Match passages by content.
- `referee_classical_09_11_phase1_combined.md`: the coordinator's merge of the two phase-1 reports, with spot checks. Read it first; items M1–M4 are yours, and it records which referee claims the coordinator already found false.
- Scripts: `briefs/referee_classical_09_11_scripts/` (Gemini) and `briefs/referee_classical_09_11_codex_scripts/` (Codex, with outputs). Codex's copies of the four Fring-group papers (PDF and text) are in `/home/gustav/Git/affine_toda/.claude/worktrees/referee-briefs-part-iii/referee_scratch/classical_09_11/codex/codex_source_*.{pdf,txt}`; you may read them as sources, but check quotations against the PDF or arXiv.

Calibration: Codex's report has been reliable so far. Gemini located real problems, but several of its claims were false. Known false or overstated Gemini claims relevant to you: that a single soliton at Im θ = π has "real negative energy" (Codex's counterexample, θ=iπ, ξ=iπ/2 in sine-Gordon, is a regular real antikink with E=+8, verified by the coordinator); the label **[Theorem, given additivity…]** for the selection rule; "real and positive" as part of the general rule; "closes gap H6 for Hirota soliton solutions". Treat every finding from either referee as a claim to check.

## Rule 1: verify before you edit

Do not trust any finding of the reports. For every item, re-derive it yourself before changing the book, by hand or with your own short sympy/numpy/mpmath script written from the chapter's definitions. Do NOT import the referees' scripts. You may import the book's helpers (e.g. `textbook_code/toda_classical.py`, read-only for you). Rerunning a referee's script is a second check, not verification. For claims about papers, read the source: texts are in `~/referee-copy/papers/*.txt` (named by arXiv id with `/` replaced by `_`); fetch anything else from arXiv. Mark each item verified, modified or rejected, with the reason.

## Read first

1. `01-introduction.qmd`, section `sec-conventions`; `README.md`, section on status labels.
2. Chapter 11 in full, chapter 10 for context, and `part9-foundations/44-twelve-gaps.qmd` (H6, lines ~33, 60, 81).
3. The combined report, then the four reports' sections on chapter 11.
4. As the model for procedure and record: commit df973fb (`git show df973fb`) and the `## Integration` section at the end of `briefs/referee_defects_41_report.md`.

## Checkpoint (mandatory)

Create `briefs/part3_integration_B_checkpoint.md`. After every completed item, append one line: the item, the verdict, the files changed, what is half done, and the next step. A fresh agent must be able to resume from this file and `git diff`.

## Items (coordinator's notes on the judgment calls)

1. **`prp-hirota-energy`** (11:32–41; combined report M2). E = Σ M cosh θ_i needs Re cosh θ_i > 0 (equivalently Re e^{θ_i}>0, |Im θ_i|<π/2), a finite nonzero leading interaction product, and no zero of τ_j on the integration line. Counterexample: sine-Gordon at θ=iπ, ξ=iπ/2 gives a regular antikink with E=+8, formula −8. Olive–Turok–Underwood (hep-th/9305160 §8.3) state the condition ("providing that its real part is positive") and that the contribution is independent of the constants Q_i; Takács–Watts eq. 4.6 is the a₂ formula. Restrict the statement, cite precisely, and keep the label **[Established; reformulated]**. Check the uses of `prp-hirota-energy` in other chapters (`17-thimbles.qmd`, `30-boundary-conditions.qmd:319`, `39-defect-solitons.qmd`, `44-twelve-gaps.qmd`, `49-unitarity.qmd`) and make minimal edits only where they now overclaim.
2. **The Θ selection rule** (11:43–49; M1). Make it a labelled statement with a proof. Hypotheses, combined from both referees: (i) a local conserved charge Q, normalized to have definite Θ-parity ε_Q (an antilinear involution allows a basis of definite parity, but not every current has one: (1+i)·J has none, Codex phase 2 §1); (ii) on the regular solutions considered, Q's value depends only on {(a_i,θ_i)}, not on the ξ_i (additivity alone does not give this; Codex phase 2 §1 sketches why phase independence holds within a regular domain); (iii) Θ-invariance of the reduced data. Topological charges are excluded: the sine-Gordon kink has Θ-invariant data and topological charge 1 with parity −1 (Codex P1 §3). Distinguish a Θ-invariant *field* from Θ-invariant *data* (11:26, 11:54). Gemini's claim that Freeman (hep-th/9408092, eq. 3.49) gives q_ā = (−1)^{s−1} q_a, hence ε_Q = (−1)^{s−1} after real normalization, is plausible (for a_n, sin(sπ(h−a)/h) = (−1)^{s−1} sin(sπa/h)); verify it from the source before using it.
   **Label: the coordinator's default.** The general mechanism (antilinear symmetry plus asymptotic additivity) is prior art: Cen–Fring arXiv:1602.05465, Cen–Correa–Fring arXiv:1608.01691 (higher charges, KdV), Fring–Taira arXiv:2007.15425, Correa–Fring–Taira arXiv:2110.06825 (Bullough–Dodd = a₂⁽²⁾ solitons and breathers); energy/momentum reality of conjugate pairs is in Takács–Watts (a₂, N=2, p. 11 before eq. 4.8) and OTU. So label the new statement **[Established; reformulated]**, not [Theorem], and cite these. Read the four papers yourself and say exactly what each proves (Codex phase 2 §1 and Gemini phase 2 Item 1 give competing summaries; Codex also found Cen–Fring "Multicomplex solitons", J. Nonlinear Math. Phys.; check it).
3. **Appendix E** (`appendices/E-original-results.qmd`, the first item, "Classical solitons (chapter 11)"). Default: remove the [Theorem, …] item, since the rule is now an established reformulation. If, after reading the prior art, something genuinely new and checked remains, for example that C-odd higher charges (a₂⁽¹⁾ spin 2 and spin 4) are imaginary on breathers and additive there, which both referees confirmed numerically and which you must verify with your own script, then replace the item by an honest **[Numerical]** item citing a new statement that carries that label. Do not keep any [Theorem] label for this rule. `check-refs.py` must pass. Record your reasoning; the editor will review the choice.
4. **Additivity test** (11:46, coordinator's phase-2 item 3). Both referees integrated the a₂⁽¹⁾ spin-2 charge (Codex also spin 4) on explicit breathers and found Q₂/q₁ = 2i e^{2θ₀} sin 2u. Construct the current yourself (or via the Lax pair), integrate it on a breather, and if it holds, state it as checked, with a script, instead of "assumed".
5. **Khastgir–Sasaki** (11:47; hep-th/9507001). Their c-solitons (single solitons at complex rapidity) have complex energy and develop singularities in finite time; their finite negative-energy solutions start from φ=0, ∂_tφ = i f(x) (eq. 5.13) and are Θ-invariant as fields. So Θ gives reality, not positivity, and "not in the physical sector" is an interpretation, not a classical theorem. Rewrite. The chapter-10 agent rewrites the matching sentence at 10:180 so that it does not call Θ-invariant solutions "the physical solutions".
6. **H6** (11:49; `44-twelve-gaps.qmd:33, 60, 81`). H6 asks which solutions count as physical. A sufficient reality condition on regular solutions does not close it. Change "closes gap H6" and the table cell to "partly addressed" (or similar), with what remains open (necessity, singular solutions, physical-state selection).
7. **Singular solutions** (11:58–62 and `exr-singular-points` 11:90–92; M3). At a simple zero of τ_j on a time slice the energy density has a double pole with zero residue, but the real-line integral diverges (static sine-Gordon τ₀,₁ = 1 ± e^{2x}: density −8 csch² 2x; with (−ε, ε) removed the integral is 8(1−coth 2ε) ≈ −4/ε). Only a contour or finite-part prescription gives a finite value. When a zero crosses the real line, the topological charge jumps by a coroot. "The Hirota data continue analytically": the τ polynomials do, the field and the real-line integrals do not. "Generically isolated points": true when the real and imaginary parts have independent gradients; a real-rapidity single soliton has lines. Rewrite the paragraph and the exercise accordingly. At 11:62 the identification of Takács–Watts's singular final states with a Θ-broken quantum sector is an interpretation; say so.
8. **Breathers and the figure** (11:23, caption 11:28, `fig-breather-theta`, `textbook_code/figures/ch11_f16_breather_theta.py`; M4). The coordinator verified: a neutral (type B) breather of species (a, ā=h−a) at θ₀±iu with ξ₂ = ξ̄₁+iπ has A_{a,ā}(2iu) = 1 − sin²(πa/h)/cos²u; it is regular iff A<0, i.e. |π/2 − πa/h| < u < π/2, so its mass 2M_a cos u fills (0, 2M_a sin(πa/h)); the full (0, 2M_a) only for a = h/2. For 0<A<1 no choice of phases avoids a zero (Codex P1 §4). Harder–Iskandar–McGhee hep-th/9409035 §3.2. The chapter-10 agent states exactly this in `sec-breathers`; use the same statement. Redo the figure's right panel to show the regular band for each species of one algebra (e.g. a₃⁽¹⁾ or a₅⁽¹⁾), with the excluded part marked, following `briefs/figures_for_the_book.md`; update the caption and the script's checks (band, sign of A, a τ zero for some 0<A<1 case). Do not use Gemini's type-A formula A_aa(2iu) = −tan²(πa/h)/tan²u; it is wrong.
9. **Missing charges** (11:72): "no static saddle of the right energy" → "no regular static solution approaching vacua at both ends" (Codex P1 item 19 proved this scope). Check `prp-static-saddles` (ch. 17) for consistency; edit there only if it overclaims.
10. **Notes** (11:96): replace "The formulation as a selection rule from Θ is new here" by an accurate account of the prior art (item 2) and of what, if anything, the chapter adds. Learning goals at 11:5–6 (Codex phase 2 §1). Also check `index.qmd:43` ("the antilinear selection rule that guarantees real energy and momentum") and `44-twelve-gaps.qmd:60`.
11. Other chapter-11 points in the reports: `i φ` versus `−i φ` (both referees: +i is right), `eq-theta-data`, `exr-theta-interaction`, `exr-breather-real` (correct per both; check), `exr-theta-solutions`. Verify the substantive ones; record the rest.

## Bibliography

Add entries to `references.bib` for the papers you cite, with metadata checked on arXiv/INSPIRE (journal, volume, year, pages, arXiv id). Follow the existing key style (lower case, surname plus year, underscores for several authors, e.g. `fring_koberle1994a`); for example `cen_fring2016`, `cen_correa_fring2017`, `fring_taira2020`, `correa_fring_taira2022`. Check that no key collides. In prose, spell out author names.

## Scripts

Write a new check script for chapter 11, `textbook_code/classical_reality.py` (one PASS/FAIL line per check, a few minutes at most), covering what the chapter newly asserts: the domain of `prp-hirota-energy` (including the θ=iπ antikink), the selection rule's examples, the breather additivity check, the divergence near a τ zero, the breather band. Cite it at each claim, as a full repository path in inline code. Update the figure script `textbook_code/figures/ch11_f16_breather_theta.py`; it writes its SVG into the repository; `git restore` any other regenerated SVG whose diff is only ids and timestamps. Do not edit `textbook_code/README.md`; put the rows it needs (new script, changed figure row) in your final message.

## Files you may edit

- `part3-classical/11-reality.qmd`.
- `appendices/E-original-results.qmd` (the Part III item only).
- `references.bib` (new entries only).
- `part9-foundations/44-twelve-gaps.qmd` (H6 only), `index.qmd` (line 43 only), and minimal consistency edits in other chapters that repeat a claim you correct (find them by grepping ids and formulas). Do not edit `part3-classical/09-lagrangian.qmd`, `part3-classical/10-solitons.qmd`, `textbook_code/toda_classical.py` or the chapter 9–10 figure scripts; record needed edits there as pending.
- `appendices/D-corrections.qmd`, only if a verified correction to a published paper arises (unlikely; follow AGENTS.md exactly).
- `textbook_code/classical_reality.py` (new), `textbook_code/figures/ch11_f16_breather_theta.py` and its SVG.
- Your checkpoint, your scripts, and `briefs/part3_integration_B.md`: write the integration record there (a table of item, source in the reports, verdict, how you checked it, files changed; then what was not integrated and why). The coordinator will merge it into the reports' `## Integration` sections. Do not edit the reports themselves.

## Style

Adapt the referees' text to the book's style; do not paste it. Short declarative sentences; define every symbol; LaTeX math, not backticks. Statements are `::: {#prp-x}` divs with a bold status label first; keep `@` references out of the label's bracket. Never cite the reports, `briefs/` or `referee_scratch/` in the book, and never mention a referee.

## Final checks, then commit

- `textbook_code/classical_reality.py` and the figure script print only PASS.
- `python3 check-refs.py` passes. Do not run a Quarto render.
- Commit the files by name, in one commit, with a message that lists the verdicts and ends with the line `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- In your final message: the verdict per item, every status label you changed or created (and whether its proof has been seen by a referee), what you did with Appendix E and why, pending edits for files you may not touch, the `textbook_code/README.md` rows, and anything the coordinator must decide.
