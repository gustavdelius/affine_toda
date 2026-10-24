# Agent briefs

Self-contained task briefs for follow-up work on the book. On 2026-10-10 a new Part VII (the half-line, chapters 28–35, `part7-half-line/`) was inserted, and the Foundations part became Part VIII with chapters 36–41 (`part8-foundations/`, formerly `part7-foundations/28–33`); dated reports keep the numbering they were written with. On 2026-10-12 a new Part VIII (integrable defects, chapters 36–43, `part8-defects/`) was inserted: Foundations moved from Part VIII (ch. 36–41) to Part IX (ch. 44–49, `part9-foundations/`) and the outlook from ch. 42 to ch. 50 (`outlook/50-outlook.qmd`); dated reports from before that date use the old numbers. The first two were written when the book had two parts; their "Part I" is now Parts V and IX and their "Part II" is Part VI, and the files they name under `foundations/` and `smatrices/` have moved (use the id table below). On 2026-10-14 the code directories were merged: `foundations_code/X` is now `code/X`, and `soliton_rmatrix_code/` is now `code/rmatrix/`; dated briefs and their `*_scripts/` folders keep the old paths. Also on 2026-10-14 an unnumbered preface became the home page `index.qmd`, and the introduction (chapter 1, including `sec-new-results` and `sec-conventions`) moved from `index.qmd` to `01-introduction.qmd`; chapter numbers did not change. The list of original results, formerly in `sec-new-results`, is now Appendix E (`appendices/E-original-results.qmd`). Each brief can be given verbatim as the prompt of a background agent, or pasted into a fresh Claude Code session started in this repository.

| Brief | Task | Paper location |
|---|---|---|
| `exchange_statistics.md` | Decide which exchange of soliton labels the Bethe–Yang equations require, and hence which multi-soliton sectors are $\Theta$-broken. **Answered 2026-10-10:** a graded exchange is a Bloch twist, and the global form selects the twists; integrated as Propositions `prp-exchange-twist` and `prp-bloch-sectors` (scripts in `foundations_code/sectors/bloch/`). Refereed on 2026-10-10 (below); the fixes are in PR #2, merged. | Section `sec-sector-status` (paragraph "Exchange statistics and Bloch sectors"), open problem 6 (see below) |
| `referee_exchange_statistics.md` | Adversarial referee pass on the answer to `exchange_statistics.md`: the proofs of `prp-exchange-twist` and `prp-bloch-sectors`, the scope claim "in every global form", the numerics in `foundations_code/sectors/bloch/`, and prior art. **Run 2026-10-10** (items 1–5, partly 6); its findings are integrated in PR #2 | Section `sec-sector-status`, open problem 6, `sec-crossing-sign` |
| `referee_exchange_statistics_report.md` | Report of the referee pass on PR #2 (2026-10-10), with the integration notes and the open items; scripts in `referee_exchange_statistics_scripts/`. Not a brief | `part9-foundations/48-scattering.qmd`, `foundations_code/sectors/bloch/` |
| `referee_textbook.md` | Referee the textbook chapters 1–16 and 20, and the textbook framing of 17–27: what was verified, and a prioritized list of unverified formulas, conventions, claims and references | `part1-background/`, `part2-algebra/`, `part3-classical/`, `part4-real-coupling/`, `part5-semiclassics/`, `part6-soliton-smatrices/`, `code/` |
| `referee_textbook_report.md` | Report of the referee pass of 2026-10-10 on `referee_textbook.md`: a verdict per item, with replacement text for each error; scripts in `referee_textbook_scripts/`. Not a brief: it records findings to be integrated. Item 2 (the sign of $q$, crossing sign) was integrated on 2026-10-10, with the dependent results recomputed at the physical $q$; the other items on branch `textbook-referee-fixes`, except those listed at the top of the report | all parts, `references.bib` |
| `referee_half_line_28_29.md` | Referee chapters 28–29 of Part VII: scattering axioms with a boundary, boundary Lee–Yang, boundary sine-Gordon and sinh-Gordon; the corrections to Ghoshal–Zamolodchikov, the status of the UV–IR and $E(\varepsilon,\beta)$ relations, the unreproduced BPTT factor | `part7-half-line/28-boundary-scattering.qmd`, `29-boundary-sine-gordon.qmd`, `code/boundary_axioms.py`, `boundary_sine_gordon.py`, `shg_boundary_breathers.py` |
| `referee_half_line_30_31.md` | Referee chapters 30–31 of Part VII: the boundary convention, Sklyanin's formalism, the classification table (corrected $g_2$, $c_2$ rows), the soliton-preserving condition, the method of images, boundary breathers and real-coupling vacua; corrections to Bowcock et al. and Delius 1998 | `part7-half-line/30-boundary-conditions.qmd`, `31-half-line-solitons.qmd`, `code/boundary_classical.py`, `boundary_images.py`, `boundary_vacua.py` |
| `referee_half_line_32_33.md` | Referee chapters 32–33 of Part VII: the boundary coideal, $K$ as intertwiner, the Neumann/uniform identification that overturns Delius–MacKay, the $a_n$ and $d_n$ reflection matrices, XXZ of any spin; corrections to Delius–MacKay, Delius–Nepomechie, Delius–Gandenberger and Delius–George | `part7-half-line/32-boundary-quantum-groups.qmd`, `33-soliton-reflection.qmd`, `code/boundary_coideal.py`, `reflection_matrices.py` |
| `referee_half_line_34.md` | Referee chapter 34 of Part VII: quantum symmetric pairs, the dictionary with Kolb's coideals, the Serre obstruction of the boundary coideal and the answer that the $s=\pm1$ (Neumann) coideals are not quantum symmetric pairs in Kolb's sense, the q-Onsager algebra, universal and trigonometric K-matrices, generic irreducibility, generalized Satake diagrams | `part7-half-line/34-quantum-symmetric-pairs.qmd`, `code/qsp_checks.py` |
| `referee_half_line_35.md` | Referee chapter 35 of Part VII: breather and particle reflection amplitudes of $a_n$ (uniform, Neumann, solitonic boundaries), the $a_2$ boundary spectrum, the corrected Delius–Gandenberger (6.11) and (6.26), the $n=1$ discrepancy, and the open problems of the part | `part7-half-line/35-particle-reflection.qmd`, `code/particle_reflection.py` |
| `referee_half_line_28_29_report.md`, `referee_half_line_30_31_report.md`, `referee_half_line_32_33_report.md`, `referee_half_line_34_report.md`, `referee_half_line_35_report.md` | Reports of the referee passes of 2026-10-11 on the Part VII briefs (chapters 34 and 35: same date, later), with a verdict per item, replacement text, and a closing section recording what was integrated and what was not. Not briefs | `part7-half-line/28–35`, `code/boundary_*`, `reflection_matrices.py`, `shg_boundary_breathers.py`, `particle_reflection.py`, `qsp_checks.py`, `references.bib` |
| `referee_part8_report.md` | Adversarial referee report of 2026-10-11 on the foundations part (then Part VIII, chapters 36–42; now Part IX, chapters 44–50), with a closing Integration section of 2026-10-13: no universal $\Theta$-breaking for four $a_2^{(1)}$ solitons (window $\operatorname{Re}q<-\frac12$ and the parity rule near $q=-1$), the open non-degeneracy of the vacuum, the scope of "no sign problem", the coupling profile, $\beta^2_{\rm UV}=4\pi$ for all non-simply-laced algebras, a rejected $q^{-4}$ correction. Not a brief | `part9-foundations/44–49`, `outlook/50-outlook.qmd`, `foundations_code/scan_n4.py`, `uv_thresholds.py`, `zero_mode.py`, `referee_part8_scripts/` |
| `rigorous_heat_trace.md` | Turn Proposition `prp-net-count` part 2 and Proposition `prp-one-loop-mass` from sketches into theorems. Done 2026-10-10: proofs integrated (PR #3) and refereed (below) | Sections `sec-counting-rule-transmission-factors` and `sec-one-loop-masses-folded`, open problem 3 (see below) |
| `referee_heat_trace.md` | Adversarial referee of those proofs (`eq-one-loop-mellin`, `prp-net-count` parts 2–3, `lem-born`, `prp-one-loop-mass`) before their [Theorem] labels are kept. Done 2026-10-10: no fatal errors; all [Theorem] labels kept after fixes made in place (proof details, applicability and status statements, the Mellin-form gap stated in Chapter 17, prior art credited), plus the independent $\phi^4$ check `phi4_check.py` | Chapters 17–19, `foundations_code/oneloop/heat_trace_checks.py` |
| `halfline_open_problems.md` and `halfline_open_NN_*.md` | One brief per open problem of Part VII (all except 10, defects), each with a literature survey and plan of attack, written 2026-10-11; the index file ranks them by feasibility (3, 12, 6, 1, 11, 7, 8, 2, 9, 5, 4) and lists the new references found | `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`) and the sections each brief names |
| `halfline_consolidation_common.md`, `halfline_consolidation_A<k>_prompt.md`, `_report.md`, `_scripts/` | Consolidation round of 2026-10-13 (section A of `halfline_open_continuation.md`): prompts, agent reports and scripts for items A1–A6, and `halfline_consolidation_referee_prompt.md` (GPT6-Astra referee brief for the new theorems, literature corrections and numerical results), each verified by the coordinator before integration; outcomes in the status block of `halfline_open_problems.md` and the Integration tables of the round-1 reports | `part7-half-line/30, 31, 33, 35`, `appendices/D-corrections.qmd`, `code/` (`solitonic_linear`, `sp_scalar`, `boundary_dirichlet`, `dn_nondiagonal_scalar` and extended scripts) |
| `figures_for_the_book.md` | How to add figures to chapters that have none: the decisions made for the introduction and Part I (committed SVGs from matplotlib scripts that check what they plot, mermaid only for simple graphs, 7 inch canvas), the recipe and verification steps, the pitfalls met, and figure ideas for Parts II–VII. Written 2026-10-11 after nine figures were added to the introduction and chapters 2–4 | `figures/`, `part1-background/figures/`, `code/figures/` |

Recommended protocol, which worked in the session that wrote these briefs:

1. Run one worker per brief; at most two or three agents at a time, because long parallel runs hit usage limits.
2. The worker reports in its final message; it does not edit the repository.
3. Spot-check its key numbers by re-running its scripts, then integrate the results into the paper.
4. Send the integrated text to a fresh adversarial referee agent before keeping any [Theorem] label.

## Locating sections and results

The papers are chapters of a Quarto book, and Quarto renumbers sections and theorems automatically. The briefs therefore cite **ids**, which do not change. Find an id with, for example, `grep -rn '{#prp-net-count}' part*/`. Id prefixes: `sec-` section, `thm-`/`prp-`/`lem-`/`cnj-` theorem, proposition, lemma, conjecture, `eq-` equation.

| Id | File (current) | Number when the briefs were written | Former number |
|---|---|---|---|
| `sec-lagrangians` | `part3-classical/09-lagrangian.qmd` | 2.2 | 2.2 |
| `sec-counting-rule` | `part5-semiclassics/17-thimbles.qmd` | 8.2 | 8.2 |
| `sec-counting-rule-transmission-factors` | `part5-semiclassics/18-jordan-chains.qmd` | 8.5 | 8.5 |
| `sec-hirota-form-dn1-en1` | `part5-semiclassics/18-jordan-chains.qmd` | 8.8 | 8.8 |
| `sec-one-loop-masses-folded` | `part5-semiclassics/19-one-loop-masses.qmd` | 8.9 | 8.9 |
| `sec-krein-unitarity` | `part9-foundations/48-scattering.qmd` | 9.2 | 9.2 |
| `sec-sector-status` | `part9-foundations/48-scattering.qmd` | 9.3 | 9.3 |
| `prp-net-count` | `part5-semiclassics/18-jordan-chains.qmd` | Proposition 8.3 | Proposition 8.5 |
| `lem-wronskian` | `part5-semiclassics/18-jordan-chains.qmd` | Lemma 8.1 | Lemma 8.4 |
| `thm-hirota-factorization` | `part5-semiclassics/18-jordan-chains.qmd` | Theorem 8.2 | Theorem G |
| `prp-cn-threshold-blocks` | `part5-semiclassics/18-jordan-chains.qmd` | Proposition 8.5 | Proposition 8.7 |
| `prp-hirota-form` | `part5-semiclassics/18-jordan-chains.qmd` | Proposition 8.6 | Proposition 8.8 |
| `prp-one-loop-mass` | `part5-semiclassics/19-one-loop-masses.qmd` | Proposition 8.7 | Proposition 8.9 |
| `lem-born` | `part5-semiclassics/19-one-loop-masses.qmd` | new (2026-10-10) | — |
| `eq-one-loop-mass` | `part5-semiclassics/17-thimbles.qmd` | (8.1) | (8.1) |
| `eq-net-factorization` | `part5-semiclassics/18-jordan-chains.qmd` | (8.2) | (F) |
| `eq-one-loop-mellin` | `part5-semiclassics/17-thimbles.qmd` | new (2026-10-10) | — |
| `prp-krein-unitarity` | `part9-foundations/48-scattering.qmd` | Proposition 9.1 | Proposition 9.1 |
| `prp-krein-fundamental` | `part9-foundations/48-scattering.qmd` | Proposition 9.3 | Proposition 9.3 |
| `prp-scalar-breathers` | `part9-foundations/48-scattering.qmd` | Proposition 9.4 | Proposition 9.4 |
| `prp-exchange-twist` | `part9-foundations/48-scattering.qmd` | (new, 2026-10-10) | — |
| `prp-bloch-sectors` | `part9-foundations/48-scattering.qmd` | (new, 2026-10-10) | — |
| `sec-crossing-sign` | `part6-soliton-smatrices/21-construction.qmd` | 14.6 | Section 2.6 of the S-matrix note |

The numbers in the third column may drift as chapters change; the ids will not.

Open problems are items of numbered lists, not labelled objects. Use the first words of the item:

- Former Part I open problems, now in `outlook/50-outlook.qmd`: item 3 is the one that begins with a reference to `cnj-counting-rule` ("… in general"); item 6 is "The sector classification for self-conjugate theories" (its former sub-item "Exchange statistics" is now "Bloch sectors").
- Former Part II open problems, now in `part6-soliton-smatrices/26-f4-solitons.qmd` (section `sec-sm-discussion`): item 4 is "Crossing convention".

## Part VIII: integrable defects (2026-10-11 to 2026-10-12)

Part VIII (`part8-defects/`, chapters 36–43) was planned, drafted, refereed and integrated with the files below. Drafting is done by Claude only; other models referee (Gemini through `run_referee_agy.sh`, Codex for the deep passes).

| File | What it is |
|---|---|
| `part8_defects_conventions.md` | The conventions, names and id registry of Part VIII, with the source dictionary and the conventions added by each chapter (§§1–11). Read it first |
| `part8_defects_status.md` | Working record: the chapter status table, the verification record of every chapter (by script, by hand, unchecked), the corrections to the sources, the pending edits elsewhere (each marked done or not needed) and the obligations between chapters |
| `part8_defects_sources.md` | Digests of the defect literature (2026-10-11), with equation numbers of the sources; working notes, not checked results |
| `draft_ch39.md` | Drafting brief for chapter 39 (the other chapters were briefed in their sessions) |
| `draft_ch39_report.md`, `draft_ch41_report.md`, `draft_ch42_report.md`, `draft_ch43_report.md` | Drafting reports: summary, script output, verification record, corrections, requests for edits elsewhere and proposed conventions text. Not briefs |
| `draft_ch41_checkpoint.md`, `draft_ch42_checkpoint.md`, `draft_ch43_checkpoint.md`, `integration_38_40_checkpoint.md`, `global_integration_part8_checkpoint.md` | Checkpoints that let a fresh agent resume a drafting or integration task after an interruption |
| `referee_defects_36_37.md`, `referee_defects_36_37_report.md` | Light referee brief for chapters 36–37 (run on Gemini) and its report |
| `referee_defects_36_37_deep.md`, `referee_defects_36_37_deep_report.md` | Deep referee brief for the high-risk claims of chapters 36–37 (Codex) and its report; scripts in `referee_defects_36_37_deep_scripts/` |
| `referee_defects_38_40.md`, `referee_defects_38_40_report.md` | Deep referee brief for chapters 38 and 40 and its report, with the integration notes of 2026-10-12 at the end; scripts and outputs in `referee_defects_38_40_scripts/` |
| `referee_defects_39.md`, `referee_defects_41.md`, `referee_defects_42.md`, `referee_defects_43.md` | Deep referee briefs for chapters 39, 41, 42 and 43; their reports are `referee_defects_<n>_report.md` when present |
| `referee_defects_38_l6.md` | Referee brief (Codex, 2026-10-23) for the three [Theorem] labels that PR #40 added to chapter 38 from Lecture 6 of the defect lecture notes: the $a_1$ case of `prp-defect-involution` and the uniqueness of its r-matrix, `prp-sg-spin-three` (spin-$\mp3$ defect charges of sine-Gordon), and `prp-dc-gauge` (the $a_r$ defect matrix as a gauge transformation in both light-cone directions). Report: `referee_defects_38_l6_report.md` |
| `run_referee_agy.sh` | Runs a referee brief with a non-Claude model through the Antigravity CLI, in a throwaway copy of the repository; writes `briefs/<brief>_report.md` |

The book cites only the scripts in `code/` (`defect_*.py`, listed in `code/README.md`); the referee scripts are kept for reproducibility.

## Lecture notes on defects (from 2026-10-15)

The lecture notes "Windows in spacetime" (`lecture-notes/defects/`, a separate Quarto book) tell the story of Part VIII in 9 lectures. Lectures 1–5 were merged with PR #30, and the book improvements found while writing Lecture 5 with PR #33.

| File | What it is |
|---|---|
| `defects_lectures_handoff.md` | Handoff of 2026-10-17 for the agent that writes Lectures 6–9: the editor's decisions, the infrastructure of the notes, the voice and structure of a lecture, the notation registry, a ready plan for Lecture 6 with every verified formula and the items of its check script, and the promises Lectures 1–5 make to Lectures 7–9. Read it first |

## Part II: model-diverse referee pass (from 2026-10-15)

Part II (`part2-algebra/`, chapters 5–8) had one referee pass, on 2026-10-10, by Claude agents only (items 2–4 and 9–12 of `referee_textbook.md`). Much of its current text is that pass's replacement wording or was added by the reading pass (`part2_reading_fixes.md`), and has had no referee since. The two briefs below ask for a non-Claude pass, weighted toward that text. The chapters were reviewed as on main at 9d85ab7.

| File | What it is |
|---|---|
| `referee_algebra_05_06.md` | Referee brief for chapters 5–6: duality of the dual pairs and its normalization against ch. 14, the two foldings to $a_{2n}^{(2)}$, the folding and label tables, exponents and vertex operators, affine-diagram symmetries, Perron–Frobenius and conserved charges, citations, exercises. Report: `referee_algebra_05_06_report.md` |
| `referee_algebra_07_08.md` | Referee brief for chapters 7–8: twisted normalizations, the central element and the $S^2$ shift (with $a_{2n}^{(2)}$), the sign of $q$ and the Part VI dictionary, the twisted tensor-product graph for the $U_q(d_{n+1}^{(2)})$ spinor, where the graph fails, the history of the R-matrices (DGZ, GMW), roots of unity, fusion and IRF, citations, exercises. Report: `referee_algebra_07_08_report.md` |

**Files under review: do not edit until the reports are in.** `part2-algebra/05-root-systems.qmd`, `06-affine-algebras.qmd`, `07-quantum-groups.qmd`, `08-r-matrices.qmd`, and the captions of the figures in `part2-algebra/figures/`.
Both reports arrived on 2026-10-15: chapters 5–6 by GPT-Astra (high), chapters 7–8 by Gemini 3.8 Flash (high). They were integrated the same day on branch `part2-referee-integration`. Each report's `## Integration` section gives the verdict per item, and its scripts are in `briefs/<report>_scripts/`. The new check script is `code/lie_algebra_checks.py`. No labels changed, and the files are no longer frozen. Main corrections:
- The $S^2$ shift for $a_{2n}^{(2)}$ is $q^{h^\vee}$ in Kac's labelling (ch. 7), with a matching exception in ch. 20.
- MacKay–Watts use non-minimal twisted labels.
- The twisted $2T$ formula in ch. 6 needs $\beta_{\rm lit}$.
- The Coxeter-plane radii need not differ for long and short roots.
- Diagram symmetries permute species.
- The two parents of $a_{2n}^{(2)}$ need $m_F^2=m_A^2=2m_D^2$.
- The $q$ in the twisted tensor-product-graph formula is normalized on the long roots.
- The fusion formula needs no projector.
- The vertex–IRF relation is stated precisely.

Added afterwards, at the editor's request: the duality sentence for the self-dual $a_{2n}^{(2)}$ in ch. 6 ($\beta^\vee=8\pi/\beta$, $B=\frac{\beta^2/4\pi}{1+\beta^2/8\pi}\to2-B$). Sources: Dorey's lectures (`dorey1998`, §self-dual theories) and, for $a_2^{(2)}$, the Bullough–Dodd action and S-matrix of Alekseev, arXiv:1106.4758, eqs. (3.1) and (3.5). The book states no $a_{2n}^{(2)}$ particle S-matrix elsewhere, so ch. 13 and 14 needed no change.

## Part III: two-phase referee pass (from 2026-10-15)

An experiment in briefing. A brief that lists what to check carries the drafting model's view of where the risks are, so a referee from another model tends to find the errors Claude would have found anyway. This pass therefore has two phases. In phase 1 the referee gets only context (rules, conventions, history, sources, standard of evidence), writes its own ranked list of risky claims, and checks it. In phase 2 it gets the coordinator's items, checks what its own list missed, and compares the two lists. The chapters are reviewed as on main at 2d36c9c.

| File | What it is |
|---|---|
| `referee_classical_09_11.md` | Phase-1 brief for chapters 9–11. Phase-1 report: `referee_classical_09_11_phase1_report.md` |
| `referee_classical_09_11_items.md` | Phase-2 items, written before phase 1 and held outside the repository until phase 1 is complete. SHA-256 `7e46e33bf78567a291bb5350625a8a67b9139c75665e02d0a9de3252de1f41b6`. Report: `referee_classical_09_11_report.md`, with the overlap table |
| `referee_classical_09_11_phase1_report.md`, `referee_classical_09_11_report.md` | Gemini 4 Argon (high), phase 1 and phase 2 |
| `referee_classical_09_11_codex_phase1_report.md`, `referee_classical_09_11_codex_report.md` | Codex GPT-6 Astra (high), phase 1 and phase 2, with checkpoints and `referee_classical_09_11_codex_scripts/` |
| `referee_classical_09_11_phase1_combined.md` | The coordinator's merge of the two phase-1 reports (M1–M11), and at its end the `## Integration` section: verdict per item, the coordinator's checks, decisions for the editor, and the overlap table of the experiment |
| `part3_integration_coordinator.md` | Coordinator's state, restart procedure and decisions for the integration; the agents' prompts are `part3_integration_A_prompt.md` (ch. 9–10) and `part3_integration_B_prompt.md` (ch. 11), their checkpoints `part3_integration_A_checkpoint.md` and `part3_integration_B_checkpoint.md`, their scripts `referee_classical_09_11_integration_scripts/` |

Protocol for the editor:
1. Start the referee in this worktree with the phase-1 brief.
2. When it reports "Phase 1 complete", commit the phase-1 report as it stands. Then copy the held items file into `briefs/`, check its SHA-256 against the value above, and commit it.
3. Send the referee, in the same session, "Phase 2: follow `briefs/referee_classical_09_11_items.md`."
4. When the report is in, integrate it with the `integrate-referee-report` skill, and record here which errors each phase found.

Both referees, Gemini 4 Argon (high) and Codex GPT-6 Astra (high), reported both phases on 2026-10-15. The reports were integrated on 2026-10-24 on branch `part3-referee-integration`. The verdicts are in the `## Integration` section of `referee_classical_09_11_phase1_combined.md`, and each report ends with a pointer there. The new check script is `code/classical_reality.py`; `code/toda_classical.py` was extended. **The Part III files are no longer frozen.** Main changes:
- The Θ selection rule is now `prp-theta-reality`, **[Established; reformulated]**, with explicit hypotheses and the Fring-group prior art. The Appendix E item [Theorem, …] is replaced by a [Numerical] item for `prp-breather-charges`.
- `prp-hirota-energy` needs $\operatorname{Re}\cosh\theta_i>0$, a finite nonzero leading coefficient and no zero of $\tau_j$.
- Regular neutral breathers have masses in $(0,2M_a\sin\frac{\pi a}h)$ (ch. 10, 11 and the figure).
- Θ gives reality, not positivity; singular integrals diverge; H6 is only partly addressed.
- Chapter 9: spin $-3$ (already on main), the area rule for nonzero couplings, twisted Lax generators, citation scopes. Chapter 10: $W_k(x,t)$ in the vertex formula, the $a_{2n}^{(2)}$ degrees, the frame of the time advance, the fusion range.

**Outcome of the experiment.** Of 29 corrected errors, 28 were found by a referee, 23 of them already in the blind phase 1, including every major one. Codex found the prior art in phase 1 without being pointed to it. Phase 2 added five minor errors, all of scope or attribution, and settled the open questions with evidence (breather charges of spin 4, Freeman's parity, the Hollowood sign). Codex made no false claim that the integration found; Gemini made about a dozen. Recommendation: keep the blind phase first, and aim the directed phase at prior art, citation scope and interpretive sentences.
