# Scripts of the lecture notes

These scripts check the calculations and draw the figures of the lecture notes in `lecture-notes/`. Each is written from the definitions of the lectures, independently of the book's scripts, and prints one `PASS` or `FAIL` line per check; the exit status is 1 on any `FAIL`. The figure scripts import `../figures/style.py` and write their SVGs into `lecture-notes/<notes>/figures/`. Run them from any directory.

## Windows in spacetime: lectures on integrable defects (`lecture-notes/defects/`)

| Script | Lecture (ids) | What it checks | Runtime |
|---|---|---|---|
| `defects_l1_mirrors_windows.py` | 1: `sec-l1-delta`, `sec-l1-algebra`, `thm-l1-nogo`, `tbl-l1-histories`, `eq-l1-triangle-index`, `eq-l1-triangle`, `sec-l1-particle-window`, `sec-l1-unitarity-crossing`, `eq-l1-ly-window`, `sec-l1-operators`, `exr-l1-ising`, `exr-l1-circle` | the $\delta$-defect amplitudes from the matching conditions, their unitarity relations, bound state and pole regimes; the unitarity relations of the exchange algebra; the two histories of the no-go argument by symbolic rewriting, coefficient by coefficient; index and map forms of the triangle relation; a particle as a window (six-vertex $\check R$), with a failing control; windows from shifted bulk amplitudes; the Lee–Yang window (axioms, excited window, energies, unitary lines); the Verlinde lines of $\mathcal M(2,5)$; the Ising defect; the twist window; the $\delta$-defect on a circle | 27 checks, a few seconds |
| `defects_l1_f1_mirror_window.py` | 1: `fig-mirror-window` | $\lvert R\rvert^2+\lvert T\rvert^2=1$ at the plotted point; widths and slopes of the world lines | seconds |
| `defects_l1_f2_delta_poles.py` | 1: `fig-delta-poles` | the plotted poles solve $2im\sinh\theta=g$; a pole is in the physical strip exactly for $-2m<g<0$; the bound-state energy | seconds |
| `defects_l1_f3_histories.py` | 1: `fig-histories` | the order of events, the rapidity difference of each collision and its side of the defect, computed from the world lines | seconds |
