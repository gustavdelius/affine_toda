# Prompt for the agents that integrate the half-line lecture referee reports

Three agents (A, B, C) run this prompt at the same time, each with its own assignment, in one worktree. They edit disjoint files and do not commit; the coordinator reviews and commits.

---

You are integrating an independent referee report into the lecture notes of a mathematical-physics project (a Quarto book on affine Toda field theory with companion courses). Read `AGENTS.md` at the repository root first; its rules apply throughout. Work in the git worktree

  `/home/gustav/Git/affine_toda/.claude/worktrees/referee-halfline-lectures`   (branch `claude/referee-halfline-lectures`)

Run every command from that directory, and never `cd` to the main repository. Do not spawn subagents. Never use `git stash`, `git add`, `git commit`, `git checkout` or `git reset`: two other agents are editing other files in this worktree at the same time, and the coordinator commits. You may run `git diff` and `git restore` on **your own** SVG files only.

## The task

The course `lecture-notes/half-line/` (six lectures and `index.qmd`, drafted by Gemini) was refereed by Codex in two phases. Read, in this order:

1. `briefs/referee_lectures_halfline_combined.md`: the coordinator's merged checklist. Your items are in the sections for your agent letter, plus the general items G1–G4 in your files.
2. The two reports it points to: `briefs/referee_lectures_halfline_phase1_report.md` (P1) and `briefs/referee_lectures_halfline_report.md` (P2). They hold the referee's evidence and proposed replacement text.
3. Your lecture files in full, and the book chapters they accompany (`part7-half-line/29-…` to `36-…`), as far as you need them.

Line numbers in the reports refer to the lectures at commit 035351d; the files have not changed since, but match passages by content.

## Rule 1: verify before you edit

Do not trust any finding of the report, and do not paste its replacement text unchecked. For every item, re-derive it yourself before changing the lecture, by hand or with your own short sympy/numpy/mpmath script written from the lecture's and the book's definitions. Do NOT import the referee's scripts (`briefs/referee_lectures_halfline_scripts/`); rerunning them is a second check, not verification. You may import the book's helpers in `code/` (for example `code/boundary_common.py`, `code/particle_reflection.py`, `code/an_nondiagonal_scalar.py`).

For claims about published papers, read the source. Texts are in `~/referee-copy/papers/*.txt` (named by arXiv id, `/` → `_`), and the referee's downloaded arXiv TeX is in `referee_scratch/lectures_halfline/sources/<id>/` in this worktree. If you cannot read a source, keep the lecture's claim modest rather than reconstructing the source.

Use generic couplings (for example $\omega=\sqrt2$ or $e/2$), never short decimals, which make $q$ a root of unity.

Mark each item **verified** (edited as the report proposes, up to wording), **modified** (edited differently, say how and why) or **rejected** (not edited, say why).

## How to edit

- Keep the lecture's teaching voice (second person, blackboard calculations), but write short declarative sentences, define every symbol, and add no new exclamation marks or superlatives. Remove drafting debris ("Wait—", "= ⋯", "No, …").
- Fix the mathematics and the status of claims; do not rewrite passages the report does not touch.
- When a statement rests on a book result labelled [Conjecture] or [Numerical], say so in the lecture's words ("conjectured", "proposed", "checked numerically"), with a `[text](book:<id>)` link. Never type a book chapter or statement number (item G3).
- Never cite `briefs/`, `referee_scratch/` or the reports, and never mention a referee.
- Cite literature with `[@key]` keys from `references.bib` (for example `corrigan_delius1999`, `gandenberger1999b`).
- Do not edit book files (`part*/`, `appendices/`, the book's `code/` scripts). If a finding needs a book change, record it as pending in your checkpoint; the coordinator collects these for a separate book PR.
- Do not edit `code/lectures/README.md`; report the new state of your scripts (check names, count, runtime) to the coordinator instead.

## Scripts and figures

Your lectures' scripts are `code/lectures/halfline_lN_f1_*.py` (each also draws the lecture's figure). Update them so that they test the corrected claims, and add a check for each new or corrected formula that the lecture now asserts (for example a corrected identity or residue). Each check prints one PASS or FAIL line. Keep each script's runtime at seconds to about a minute. A check must test the claim as stated in the lecture; it must not, for example, choose whichever sign pattern fits best.

The scripts hard-code their SVG output paths under the repository root, so running one overwrites its SVG in `lecture-notes/half-line/figures/` in this worktree. That is intended when you fix a figure. If a rerun changes only ids and timestamps, `git restore` that SVG. Run scripts with `python3 -u`; never use `pgrep -f` or `pkill -f`.

Your scratch directory is `/tmp/claude-1000/-home-gustav-Git-affine-toda/af9e27ec-122e-43f3-a547-2dbab370a65e/scratchpad/integ_<letter>/`. Before you finish, copy the scripts your verdicts rely on into `briefs/referee_lectures_halfline_scripts/` with the prefix `integration_<letter>_`.

## Checkpoint (mandatory)

Create `briefs/referee_lectures_halfline_integration_<letter>_checkpoint.md`. After every completed item, append one line: the item id, the verdict, the files changed, the evidence (script name or "by hand"), anything half done, and the next step. A fresh agent must be able to resume from this file and `git diff`.

## Final checks

- Your scripts print only PASS.
- `python3 check-refs.py lecture-notes/half-line` and `python3 check-refs.py` pass. Do not run `quarto render`.
- At the end of your checkpoint file, add a section `## Table` with one row per item: id, verdict, how checked, files changed. Then a section `## Pending` listing book findings, items you could not finish, and the new state of your scripts (check names, count, runtime) for the README.

In your final message, summarize the verdicts (counts of verified / modified / rejected), list rejected and modified items with one-line reasons, and list anything pending.
