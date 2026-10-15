# Template: prompt for a Claude agent that integrates a Codex referee report

This is the Part VIII version, with its worktree and files hard-coded. The general procedure is the Claude Code skill `.claude/skills/integrate-referee-report/`, whose `agent-prompt.md` generalizes this template.

Used for chapter 41 on 2026-10-13: the referee report was committed first as 3f7883a, then integrated as df973fb. Before each integration:

1. Commit the report unchanged.
2. Copy the referee's check scripts into `briefs/referee_defects_<N>_scripts/`, rerun them from there, and compare their outputs with the saved ones. Leave out report-format validators and rendered pages of papers.
3. Fill in the placeholders below: `<N>`, the chapter file and its id, the report's main points, the chapters that repeat corrected claims, and the files that must not be touched because a referee is reading them.

Launch it as a general-purpose agent on Opus, in the background. Afterwards, the coordinator re-derives every label upgrade, reruns the chapter script, pushes, and lists label upgrades without a referee pass in the PR.

---

You are integrating an independent referee report into a mathematical-physics textbook (Quarto book on affine Toda field theory). Work in the git worktree

  /home/gustav/Git/affine_toda/.claude/worktrees/integrable-defects-framework-a54a84   (branch claude/integrable-defects-framework-a54a84)

Run every command from that directory. Never `cd` to the main repository. Never use `git stash`. Never `git add -A` or `git add .`: the untracked `referee_scratch/` belongs to the Codex referee and must not be committed or edited. Do not push.

Scratch directory for your own check scripts: `<scratchpad>/int<N>/`. Copy the scripts the record relies on into `briefs/referee_defects_<N>_scripts/`, with the prefix `integration_`, before you commit.

## The task

Chapter <N>, `part8-defects/<file>.qmd` (id `<sec-id>`), was refereed by Codex. The report is `briefs/referee_defects_<N>_report.md`, committed unchanged in <commit>. The referee's scripts are in `briefs/referee_defects_<N>_scripts/` with their outputs. Line numbers in the report refer to the review copy under `referee_scratch/defects_<N>/`, which is read-only for you. Match passages by content, because the live file may have moved on.

## Rule 1: verify before you edit

Do not trust any finding of the report. For every item, re-derive it yourself before changing the book, by hand or with your own short sympy/numpy/mpmath script. Write the script from the chapter's definitions. Do NOT import the referee's scripts; you may import `textbook_code/defect_common.py` and the book's own helpers. Rerunning the referee's scripts is a second check, not verification.

For claims about published papers, read the source in `~/referee-copy/papers/*.txt` (and the PDFs there). When a paper contradicts itself, decide which side is wrong, for example by checking which formula solves the paper's own equations.

Mark each item verified, modified or rejected, with the reason.

## Read first

1. `briefs/part8_defects_conventions.md`.
2. The chapter-<N> section and "Global pending edits" of `briefs/part8_defects_status.md`.
3. The whole chapter and the report.
4. As the model for procedure and record: commit df973fb (`git show df973fb`), the "Integration" section at the end of `briefs/referee_defects_41_report.md`, and `briefs/integration_41_checkpoint.md`.

## Checkpoint (mandatory)

Create `briefs/integration_<N>_checkpoint.md` and append one line after every completed item: the item, the verdict, the files changed, what is half done, and the next step. A fresh agent must be able to resume from this file and `git diff`.

## Items

<one entry per report section, with the coordinator's notes on the judgment calls: which corrections touch other chapters, which labels may change, which new Appendix D entries to expect>

## Script

`textbook_code/<script>.py` must reflect the corrected claims and check whatever the chapter newly asserts. Keep the runtime reasonable.

If you add helpers to `textbook_code/defect_common.py`, only append. Then rerun every `textbook_code/defect_*.py` that imports it. Run them in the background and wait on their PIDs, never with `pgrep -f` loops.

Update the script's row in `textbook_code/README.md` if its check count or runtime changes.

## Files you may edit

- The chapter (main target).
- Minimal consistency edits in other `part8-defects/` chapters that repeat a corrected claim. Not <files under review by a referee right now>; record those edits as pending in the status file instead.
- `appendices/D-corrections.qmd` and `textbook_code/README.md` (only the script's row).
- The chapter's script and its figure scripts and SVGs. Figure scripts write into the repo regardless of OUT_DIR, so run `git restore` on any other regenerated SVG whose diff is only ids and timestamps.
- `textbook_code/README.md`.
- `briefs/part8_defects_status.md` and `briefs/part8_defects_conventions.md`.
- The report: append an "## Integration" section only.
- Your checkpoint and scripts.

## Style of the book text

- Adapt the referee's text to the book's style rather than pasting it.
- Use short declarative sentences and define every symbol.
- Write cross-references with `@id` and `(-@eq-…)`.
- Never put a status label and a cross-reference in one bracket: write `[Theorem] (@sec-x)`.
- Never cite `referee_scratch/` or the report in the book, and never mention "the referee". Corrections to published papers go in the chapter Notes and in Appendix D.

## Final checks, then commit

- The chapter script prints only PASS, as do any other scripts you touched.
- Run a reference check over all book `.qmd` files: every `@sec-/@eq-/@prp-/@fig-/@exr-/@thm-/@cnj-/@tbl-…` resolves, no id is duplicated, and every citation key is in `references.bib`. Do not run a full Quarto render.
- Append the Integration section to the report.
- Update the status file and the conventions brief.
- Commit with explicit paths. Title the commit "Integrate the Codex deep referee report on chapter <N>" and end it with the line:
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>

## Your final message

Report:
- the verdict per item and what changed;
- what was rejected, and why;
- the script results;
- the commit hash;
- anything pending;
- every label upgraded without an independent referee pass.
