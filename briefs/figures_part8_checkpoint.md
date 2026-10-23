# Part VIII figures: checkpoint (2026-10-13)

Task: figures F37–F4x for Part VIII (`part8-defects/`), following `briefs/figures_for_the_book.md`.
Only edit ch. 36, 37, 41, 42 (figure + one referring sentence), new `code/figures/chNN_fMM_*.py`, SVGs in
`part8-defects/figures/`, README figure rows, this file. Do not touch ch. 38, 39, 40, 43 or `referee_scratch/`.
Commit message: "Part VIII figures F37–F4x" + Co-Authored-By line; stage named paths only.

| # | Figure | Script | Id | Status |
|---|---|---|---|---|
| F37 | ch. 37 soliton delay through jump defect | `ch37_f37_sg_defect_delay.py` | `fig-sg-defect-delay` | DONE: script (7 PASS), SVG, inserted after proof of prp-sg-delay, README row |
| F38 | ch. 36 doubling trick | `ch36_f38_doubling_trick.py` | `fig-doubling-trick` | DONE: script (6 PASS), SVG, inserted before "The doubled theory is not parity invariant", README row |
| F39 | ch. 36 Lee–Yang defect poles | `ch36_f39_ly_defect_poles.py` | `fig-ly-defect-poles` | DONE: script (5 PASS), SVG, inserted before "**Non-unitarity.**" in sec-defect-lee-yang, README row |
| F40 | ch. 41 a2 unstable-state windows | `ch41_f40_a2_unstable_windows.py` | `fig-a2-unstable-windows` | DONE: script (4 PASS), SVG, inserted before "**Corrected windows.**", README row |
| F41 | ch. 42 doubled-boundary symmetries (optional) | `ch42_f41_doubled_symmetries.py` | `fig-dbd-symmetries` | DONE: script (5 PASS), SVG, inserted before "**Beyond sine-Gordon.**", README row |

## Final checks (2026-10-13)

- All five scripts exit 0 with 27 PASS and no FAIL (F37 7, F38 6, F39 5, F40 4, F41 5).
- Reference check over all 64 book `.qmd` files: no duplicate ids, no unresolved references; each new `fig-` id
  is defined once and referenced once.
- Viewed every SVG as PNG (light at 1000 px, dark contact sheet at 520 px).
- Next step if interrupted: commit the named paths with message "Part VIII figures F37–F41".
