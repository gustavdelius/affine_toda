You are integrating independent referee reports into a mathematical-physics textbook (a Quarto book on affine Toda field theory). Read `AGENTS.md` at the repository root first; its rules apply throughout. Work in the git worktree

  `/home/gustav/Git/affine_toda/.claude/worktrees/part3-int-ch9-10`   (branch `part3-int-ch9-10`)

Run every command from that directory, and never `cd` to the main repository. Never use `git stash`. Never use `git add -A` or `git add .`; stage files by name. Do not push. Do not spawn subagents.

Scratch directory for your own check scripts: `/tmp/claude-1000/-home-gustav-Git-affine-toda/7fa1324b-155d-4b76-8703-4c3aaf6d6c0d/scratchpad/part3_int_A/`. Before you commit, copy the scripts the record relies on into `briefs/referee_classical_09_11_integration_scripts/`, with the prefix `integration_A_`.

## The task

Part III of the book (chapters 9–11) had a two-phase referee pass by two referees, Gemini and Codex. **You integrate everything that concerns chapters 9 and 10**: `part3-classical/09-lagrangian.qmd` (`sec-theories`) and `part3-classical/10-solitons.qmd` (`sec-classical-solitons`). A second agent is integrating chapter 11 at the same time, in another worktree, and owns `part3-classical/11-reality.qmd`, Appendix E, `references.bib`, Part IX and the chapter-11 figure. Do not edit any of those.

The reports, all committed unchanged in 45955bb, are in `briefs/`:
- `referee_classical_09_11_phase1_report.md` (Gemini, phase 1) and `referee_classical_09_11_report.md` (Gemini, phase 2). Reviewed commit a03edd4; line numbers match the current files.
- `referee_classical_09_11_codex_phase1_report.md` and `referee_classical_09_11_codex_report.md` (Codex). Reviewed commit bcd82b9; chapter 9 line numbers after line 117 are 2 lower than now. Match passages by content.
- `referee_classical_09_11_phase1_combined.md`: the coordinator's merge of the two phase-1 reports, with spot checks. Read it first; it ranks the findings (M1–M11) and records which referee claims the coordinator already found false.
- Scripts: `briefs/referee_classical_09_11_scripts/` (Gemini) and `briefs/referee_classical_09_11_codex_scripts/` (Codex, with outputs).

Calibration: Codex's report has been reliable so far. Gemini located real problems, but three of its claims were false (a type-A breather formula, a misreading of Olive–Turok–Underwood eq. 6.3, a mislabelled d4 coupling table). Treat every finding from either referee as a claim to check.

## Rule 1: verify before you edit

Do not trust any finding of the reports. For every item, re-derive it yourself before changing the book, by hand or with your own short sympy/numpy/mpmath script written from the chapter's definitions. Do NOT import the referees' scripts. You may import the book's helpers (`textbook_code/toda_classical.py`, `textbook_code/lie_algebra_checks.py`). Rerunning a referee's script is a second check, not verification. For claims about papers, read the source: texts are in `~/referee-copy/papers/*.txt` (named by arXiv id with `/` replaced by `_`); fetch anything else from arXiv. Mark each item verified, modified or rejected, with the reason.

## Read first

1. `01-introduction.qmd`, section `sec-conventions`.
2. Chapters 9 and 10 in full, and chapter 11 for context.
3. The combined report, then the four reports' sections on chapters 9–10.
4. As the model for procedure and record: commit df973fb (`git show df973fb`) and the `## Integration` section at the end of `briefs/referee_defects_41_report.md`.

## Checkpoint (mandatory)

Create `briefs/part3_integration_A_checkpoint.md`. After every completed item, append one line: the item, the verdict, the files changed, what is half done, and the next step. A fresh agent must be able to resume from this file and `git diff`.

## Items (coordinator's notes on the judgment calls)

Chapter 9:
1. **Spin of `eq-spin3-current`** (09:96, `exr-spin3`). Both referees: built from ∂₊ = ∂_t+∂_x, the charge scales as e^{−3θ} (spin −3); the spin-+3 charge comes from exchanging ∂₊ and ∂₋. Both computed Q = −16 m₀³ e^{−3θ}/(3β_SG²) on the kink. Check against the book's rapidity convention (`eq-charge-eigenvalue`, ch. 2) and the boost convention at 09:77. Also: "replace β² by −β²" for sine-Gordon. The coordinator judged Gemini's "false" overstated (cosh βφ is even in β), but make the substitution explicit (β→iβ, cosh→cos). The example's sinh-Gordon m, β are not the root-normalized ones of `exr-lax`; say which dictionary it uses (`sec-sg-dictionary`).
2. **Conservation-law sign** (09:41, Codex phase 2 §5): the generic form "∂₋T_{s+1}=∂₊X_{s−1}" versus the displayed current's sign. Make them consistent.
3. **Quantum conservation citations** (09:32, Codex phase 2 §5 and §9): feigin1993, dgz1992, niedermaier1994 have different scopes. Check the quoted restrictions in the sources yourself and soften the sentence only as far as the sources require.
4. **"Minimum"** of the complex potential (09:30, Codex P1 item 17): say what is minimal (the real part, on real fields) or rephrase.
5. **Lax-pair generators for twisted algebras** (09:57, Codex P1 item 5, phase 2 §5): the description of e_j, f_j is untwisted only. Verify with an explicit a₂⁽²⁾ check of your own, and fix the wording. Codex also checked the "another gradation needs a constant gauge transformation" remark; verify or leave.
6. **Area rule** (09:130–138, Codex P1 item 7; Gemini P1 rank 11): the displayed equality must be restricted to nonzero couplings (a₂: C₁₁₂=0 though the triangle exists; the book itself says at 09:142 that many triangles carry no coupling). "(with m=1)": make the normalization of Δ consistent with C's m². Gemini's separate point: for d_{2m}^{(1)} the centre characters are ±1, so the charge basis is real; "the complex eigenvectors of the centre" and "in a real orthonormal basis … neither holds" need rewording (an *arbitrary* real basis fails). Check e.g. a₄ and d₄ yourself.
7. Spin lists for twisted algebras (09:88) and the simply-laced qualifier at 09:98: Codex's phase-2 centralizers support the lists. Spot-check one (e.g. a₄⁽²⁾ or (g₂⁽¹⁾, d₄⁽³⁾)) yourself; change nothing unless wrong.

Chapter 10:
8. **Vertex-operator formula** (10:72): the displayed τ_j has no x, t dependence; chapter 6 (`06-affine-algebras.qmd:168`) has the factor W_k(x,t), as does Olive–Turok–Underwood (hep-th/9305160, eqs. 2.28–2.29). Gemini's defence via OTU eq. (6.3) is a misreading (6.3 is the level-x tensor power). Fix 10:72–75 consistently with chapter 6.
9. **Degrees of folded τ-functions** (10:77, Gemini P1 rank 6): `tbl-foldings` (ch. 6) includes a₂ₙ⁽²⁾ and defines n_j^∨ = k n_j α_j²/2, which gives (2,…,2,1); a coincident pair from a₂ₙ⁽¹⁾ has degree 2 at every node, so "|O| n_j^∨" fails there. Ch. 6 already calls a₂ₙ⁽²⁾ "an exception to both statements". Codex checked g₂ orbit degrees (3,3,6)=3n^∨. Fix the scope; also say that "τ₁=τ₃" uses the parent a₃⁽¹⁾ labels. Do not upgrade anything that depends on @hyp-D.
10. **Position of the soliton** (10:59): −Re ξ/m_a is the position at t=0 only for θ=0; in general −Re ξ/(m_a cosh θ). Verify and fix.
11. **Time delays** (10:122 and the caption of `fig-two-soliton-shift`): "Δt = −Δx/v₁ < 0" and "shifted forward" need v₁>0; the caption's "each soliton is displaced forward along its own direction of motion" holds in a frame with v₁>0>v₂ (true in the figure). Use "displacement magnitude" or state the frame. Rename the energy E₁ in this paragraph (clash with the exponential E₁; the notation note at 10:61 forbids it).
12. **Wigner relation versus the quantum S-matrix** (10:122). The two referees disagree. Codex read Hollowood (hep-th/9203076, eqs. 5.2–5.8) and found that the leading semiclassical derivative is −(2h/β²) ln A, the opposite sign, *before* matching in/out and amplitude conventions, and refused to call it an error. Gemini's "confirmation" used an integral representation it wrote itself. Resolve this in the book's own conventions: take the book's a_n^{(1)} soliton S-matrix in `sec-an-example` (`part6-soliton-smatrices/21-construction.qmd`) and the book's sine-Gordon S₀ (`eq-sg-s0`, ch. 3), and determine the sign of the leading semiclassical phase derivative for a₂ (and a₃) with your own script. Then either (a) confirm the sentence, (b) restrict the claim to sine-Gordon, as Codex proposes, or (c) if the book's chapter-21 S-matrix and chapter-10 classical delay genuinely disagree, do not edit chapter 21: record it as pending, with your evidence, in the integration record, and tell the coordinator in your final message. Do not propose a correction to a published paper without the full AGENTS.md procedure.
13. **Fusion** (10:126 and `exr-fusion`, both referees): needs a+b<h; a+b=h is not fusion to a species (vacuum, double pole); for a+b>h the channel is species a+b−h with rapidity difference iπ(2h−a−b)/h. The limit must be joint in θ₁−θ₂ and ξ₁,₂. Verify and fix both places.
14. **Breathers** (`sec-breathers`, 10:128–136, both referees; HIM hep-th/9409035 §3.2). The coordinator verified: a neutral (type B) breather of species (a, ā=h−a) at θ₀±iu with ξ₂ = ξ̄₁+iπ has A_{a,ā}(2iu) = 1 − sin²(πa/h)/cos²u; it is regular iff A<0, i.e. |π/2 − πa/h| < u < π/2, so its mass 2M_a cos u fills (0, 2M_a sin(πa/h)), the full (0, 2M_a) only for a = h/2. For 0<A<1 no choice of phases avoids a zero of τ (Codex P1 §4 gives a proof; check it). HIM: "it is not possible to have type B breathers with positive interaction coefficient". Mention that same-species (type A, charged) breathers also exist (HIM §3.1). **Do not use Gemini's type-A formula A_aa(2iu) = −tan²(πa/h)/tan²u; it is wrong.** From `eq-interaction-coefficient`, A_aa(2iu) = −sin²u/(sin²(πa/h) − sin²u), which changes sign at u = πa/h. The chapter-11 agent will state the same band in chapter 11 and redo its figure; use exactly this statement so the two chapters agree. Check that the quantum sentence after the mass formula ("the lowest breather state … mass ratio sin(πa/h)") is still consistent.
15. **Notes** (`sec-classical-literature`): (a) 10:174, "Reality.": the total-derivative stress tensor holds in the standard normalized τ gauge (Codex phase 2 §2: C is not a gauge-invariant local function of φ); replace "and the energies are real" by a statement restricted to regular solutions with convergent integrals. (b) 10:180, Khastgir–Sasaki: they found single solitons at complex rapidity with complex energy that blow up in finite time, and finite negative-energy solutions from initial data φ=0, ∂_tφ = i f(x) (hep-th/9507001 eq. 5.13; these data are Θ-invariant). Rewrite the sentence so that it does not say that Θ-invariant solutions are "the physical solutions"; say that chapter 11 gives a sufficient condition for real charges of regular Hirota solutions. The chapter-11 agent writes the matching text in chapter 11.
16. Anything else in the reports about chapters 9–10 (e.g. Lie-duality citation `kneipp1996`, which Codex found supported; McGhee's count; exercises in Codex phase 2 §9). Verify the substantive ones; record the rest.

## Scripts

`textbook_code/toda_classical.py` must reflect the corrected claims and check what the chapters newly assert (for example the spin −3 scaling, a vanishing coupling with an existing triangle, the fusion range, the breather band and the zero for 0<A<1, the position formula). Keep its runtime reasonable (now about 3 minutes). Note that `smatrix_simply_laced.py` and `algebra_data.py` import it: rerun both afterwards, in the background, waiting on PIDs (never with `pgrep -f` loops), with `python3 -u`. The figure scripts `textbook_code/figures/ch9_f13_dorey_d4.py`, `ch10_f14_topological_charges.py` and `ch10_f15_two_soliton_shift.py` may need changes only if a caption changes what they check; they write SVGs into the repository; `git restore` any SVG whose diff is only ids and timestamps. Cite a script at each new claim, as a full repository path in inline code. Do not edit `textbook_code/README.md`; put the row changes it needs in your final message.

## Files you may edit

- `part3-classical/09-lagrangian.qmd` and `part3-classical/10-solitons.qmd`.
- Minimal consistency edits in other chapters that repeat a claim you correct (find them by grepping ids and formulas), except: `part3-classical/11-reality.qmd`, `appendices/E-original-results.qmd`, `references.bib`, `part9-foundations/44-twelve-gaps.qmd`, `part6-soliton-smatrices/21-construction.qmd` (record needed edits as pending instead). If you need a new bibliography entry, put it in your final message instead.
- `appendices/D-corrections.qmd`, only if a verified correction to a published paper arises (unlikely; follow AGENTS.md exactly).
- `textbook_code/toda_classical.py` and the three figure scripts above, with their SVGs.
- Your checkpoint, your scripts, and `briefs/part3_integration_A.md`: write the integration record there (a table of item, source in the reports, verdict, how you checked it, files changed; then what was not integrated and why). The coordinator will merge it into the reports' `## Integration` sections. Do not edit the reports themselves.

## Style

Adapt the referees' text to the book's style; do not paste it. Short declarative sentences; define every symbol; LaTeX math, not backticks. Never cite the reports, `briefs/` or `referee_scratch/` in the book, and never mention a referee. No status label changes are expected in chapters 9–10; do not upgrade any.

## Final checks, then commit

- `textbook_code/toda_classical.py`, `smatrix_simply_laced.py`, `algebra_data.py` and any figure script you touched print only PASS.
- `python3 check-refs.py` passes. Do not run a Quarto render.
- Commit the files by name, in one commit (or one per chapter), with a message that lists the verdicts and ends with the line `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- In your final message: the verdict per item, any label change (none expected), pending edits for files you may not touch, the `textbook_code/README.md` row change, and anything the coordinator must decide (especially item 12).
