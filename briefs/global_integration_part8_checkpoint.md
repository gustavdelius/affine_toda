# Global integration of Part VIII: checkpoint

Brief: the task message of 2026-10-12 (global integration of Part VIII, ch. 36–43). Sources of truth: `briefs/part8_defects_status.md` (pending-edit lists), the drafting reports `briefs/draft_ch39/41/42/43_report.md`, the Integration section of `briefs/referee_defects_38_40_report.md`, and `briefs/part8_defects_conventions.md`.

Never touch `referee_scratch/` or `briefs/referee_defects_*`. Do not push. Do not run `quarto render`.

## Done
- Task 1 (pending edits in chapters) and task 2 (theorem titles in ch. 36–39; ω clash decided: Parts III and VI already use both meanings, so one clarifying sentence in ch. 39 and its script docstring, no rename; "type II" naming in scripts: not needed, the scripts already say I_-/I_+ and every remaining "type II" is the Corrigan–Zambon defect or Weston's module). Commit: see git log "Part VIII integration: pointers and pending edits".
  Edits: ch. 35 problems 9 (→ prp-fused-bdg) and 10 (rewritten: purely transmitting, prp-rt-nogo, sec-defect-open-problems); ch. 48 twist/seam paragraph after prp-bloch-sectors; ch. 28 LY → sec-defect-lee-yang; ch. 27 notes → sec-topological-lines (cck2024); ch. 3 → sec-backlund; ch. 10 → prp-bianchi and prp-stored-charge; ch. 16 → eq-an-linear-transmission; ch. 29 Dirichlet bullet → sec-sg-breather-transmission, sec-bsg-dynamical → prp-fused-bdg, cnj-shg-e bullet → cnj-defect-parameter; ch. 34 → sec-diagonal-qsp; ch. 36 unstable a₂ window + π-line limits in sec-topological-lines; ch. 37 a₂ window, → prp-lightest-breather-universal, → prp-defect-fusion-k, → prp-defect-frozen-soliton; ch. 38 absorption restriction (→ prp-ds-absorption, prp-stored-charge); ch. 39 → sec-defect-open-problems, ω sentence; ch. 40 sec-perturbed-defects → prp-diagonal-borel, tbl-boundary-defect-algebra third column + status + "settled"; ch. 41 breather result; ch. 42 → prp-dpt-doubled.

## In progress
- (next) index.qmd, figure script.

## Remaining
1. Pending edits in chapters (pointers, ch. 35 problems 9–10, a₂ window, ch. 38 absorption restriction, ch. 42 → prp-dpt-doubled, ch. 37 → prp-lightest-breather-universal, ch. 40 sec-perturbed-defects, tbl-boundary-defect-algebra third column, optional pointers).
2. Global: theorem titles in ch. 36–39; ω clash; "type II" naming in defect_*.py.
3. index.qmd: part bullet, original results, mermaid, reading path, conventions sentence, alt text.
4. intro_f1_reading_paths.py + svg.
5. Outlook defects paragraph.
6. Appendix D, Appendix C, textbook_code/README.md, briefs/README.md, root README.md.
7. Conventions brief §10, §11, ch. 41, ch. 38–40, ch. 39 items.
8. Final checks (refs, duplicates, citations; scripts; git status).
