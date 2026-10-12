# Global integration of Part VIII: checkpoint

Brief: the task message of 2026-10-12 (global integration of Part VIII, ch. 36–43). Sources of truth: `briefs/part8_defects_status.md` (pending-edit lists), the drafting reports `briefs/draft_ch39/41/42/43_report.md`, the Integration section of `briefs/referee_defects_38_40_report.md`, and `briefs/part8_defects_conventions.md`.

Never touch `referee_scratch/` or `briefs/referee_defects_*`. Do not push. Do not run `quarto render`.

## Done
- Task 1 (pending edits in chapters) and task 2 (theorem titles in ch. 36–39; ω clash decided: Parts III and VI already use both meanings, so one clarifying sentence in ch. 39 and its script docstring, no rename; "type II" naming in scripts: not needed, the scripts already say I_-/I_+ and every remaining "type II" is the Corrigan–Zambon defect or Weston's module). Commit 1834767 "Part VIII integration: pointers and pending edits".
  Edits: ch. 35 problems 9 (→ prp-fused-bdg) and 10 (rewritten: purely transmitting, prp-rt-nogo, sec-defect-open-problems); ch. 48 twist/seam paragraph after prp-bloch-sectors; ch. 28 LY → sec-defect-lee-yang; ch. 27 notes → sec-topological-lines (cck2024); ch. 3 → sec-backlund; ch. 10 → prp-bianchi and prp-stored-charge; ch. 16 → eq-an-linear-transmission; ch. 29 Dirichlet bullet → sec-sg-breather-transmission, sec-bsg-dynamical → prp-fused-bdg, cnj-shg-e bullet → cnj-defect-parameter; ch. 34 → sec-diagonal-qsp; ch. 36 unstable a₂ window + π-line limits in sec-topological-lines; ch. 37 a₂ window, → prp-lightest-breather-universal, → prp-defect-fusion-k, → prp-defect-frozen-soliton; ch. 38 absorption restriction (→ prp-ds-absorption, prp-stored-charge); ch. 39 → sec-defect-open-problems, ω sentence; ch. 40 sec-perturbed-defects → prp-diagonal-borel, tbl-boundary-defect-algebra third column + status + "settled"; ch. 41 breather result; ch. 42 → prp-dpt-doubled.

- Tasks 3–4: index.qmd (status, Part VIII sentence, "Defects" original-results entry, Part VIII bullet, mermaid node VIII with edges from II, III, VI, VII and no edge to IX, Part II needed for VI–VIII, VIII/IX meet at ch. 48, reading path "Integrable defects: chapters 2–3, 7–11 and 20–21, then 28–32, then Part VIII", conventions sentence, alt text); intro_f1_reading_paths.py (six paths, VIII colour #8f9a1f, regex) 7 PASS, svg regenerated. Commit 418d192 "Part VIII integration: index and reading paths".

- Tasks 5–7 (uncommitted at time of writing, commit "Part VIII integration: appendices, READMEs and conventions brief"): outlook defects bullets; Appendix D Part VIII entries (+ kind "clarification"); Appendix C section sec-app-code-defects; textbook_code/README rows; briefs/README Part VIII section; root README; conventions brief §5/§6/§8/§9 edits and new §10–§13; status file marks + "Global integration" section; defect_common.sg_tmatrix docstring trap (verified numerically); boundary_classical.tau_bracket sign note (43 PASS).

## In progress
- Final checks: all defect_*.py rerun in parallel (defect_common docstring changed); ref checker 0 problems.

## Remaining
- Optional: refresh the stale comment/label about BCZ (3.44) in defect_sg_transmission.py (the r of (3.44) diverges on its own; the ratio r/rbar with common cutoff converges), then rerun it.
- Final report.

## Original task list
1. Pending edits in chapters (pointers, ch. 35 problems 9–10, a₂ window, ch. 38 absorption restriction, ch. 42 → prp-dpt-doubled, ch. 37 → prp-lightest-breather-universal, ch. 40 sec-perturbed-defects, tbl-boundary-defect-algebra third column, optional pointers).
2. Global: theorem titles in ch. 36–39; ω clash; "type II" naming in defect_*.py.
3. index.qmd: part bullet, original results, mermaid, reading path, conventions sentence, alt text.
4. intro_f1_reading_paths.py + svg.
5. Outlook defects paragraph.
6. Appendix D, Appendix C, textbook_code/README.md, briefs/README.md, root README.md.
7. Conventions brief §10, §11, ch. 41, ch. 38–40, ch. 39 items.
8. Final checks (refs, duplicates, citations; scripts; git status).
