# Prompt template: an agent that integrates a referee report

Fill in every `<…>` and delete the lines that do not apply. Launch it as a background general-purpose agent on Opus. Afterwards the coordinator re-derives every label upgrade, reruns the chapter scripts, and records the upgrades without a referee pass for the PR. This template generalizes `briefs/integration_prompt_template.md`, which was used for chapters 41–43.

---

You are integrating an independent referee report into a mathematical-physics monograph (a Quarto book on affine Toda field theory). Read `AGENTS.md` at the repository root first; its rules apply throughout. Work in the git worktree

  `<absolute worktree path>`   (branch `<branch>`)

Run every command from that directory, and never `cd` to the main repository. Never use `git stash`. Never use `git add -A` or `git add .`; stage files by name. Do not push. Do not spawn subagents.

Scratch directory for your own check scripts: `<scratchpad>/<tag>/`. Before you commit, copy the scripts the record relies on into `<briefs/…_scripts/>`, with the prefix `integration_`.

## The task

`<chapter file(s) and their section ids>` were refereed by `<Codex / Gemini / Claude>`. The report(s) `<briefs/…_report.md (and briefs/…_phase1_report.md / briefs/…_phase1_combined.md for a two-phase pass)>` were committed unchanged in `<commit>`. The referee's scripts are in `<briefs/…_scripts/>`, with their outputs. Line numbers in the report refer to `<the reviewed commit or review copy>`. Match passages by content, because the live files may have moved on.

## Rule 1: verify before you edit

Do not trust any finding of the report. For every item, re-derive it yourself before changing the book, by hand or with your own short sympy/numpy/mpmath script. Write the script from the chapter's definitions. Do NOT import the referee's scripts. You may import the book's helpers (`<e.g. code/defect_common.py>`). Rerunning the referee's scripts is a second check, not verification.

For claims about published papers, read the source (`<path to paper texts, e.g. ~/referee-copy/papers/*.txt>`). When a paper contradicts itself, decide which side is wrong, for example by checking which formula solves the paper's own equations.

Mark each item verified, modified or rejected, with the reason.

## Read first

1. `<the part's conventions brief, e.g. briefs/part8_defects_conventions.md>`.
2. `<the relevant section of the part's status record>`.
3. The whole chapter and the report(s).
4. As the model for procedure and record: commit df973fb (`git show df973fb`), the `## Integration` section at the end of `briefs/referee_defects_41_report.md`, and (for a two-phase pass) `briefs/referee_classical_09_11_phase1_combined.md`.

## Checkpoint (mandatory)

Create `briefs/<tag>_checkpoint.md`. After every completed item, append one line: the item, the verdict, the files changed, what is half done, and the next step. A fresh agent must be able to resume from this file and `git diff`.

## Items

<One entry per report section, with the coordinator's notes on the judgment calls: which corrections touch other chapters, which labels may change, and which new Appendix D entries to expect.>

## Scripts

`<chapter script>` must reflect the corrected claims and check whatever the chapter newly asserts. Keep its runtime reasonable. If you append helpers to a shared module, rerun every script that imports it, in the background, waiting on the PIDs (never with `pgrep -f` loops). Update the script's row in its directory README if its check count or runtime changes.

## Files you may edit

- The chapter(s) above (the main target).
- Minimal consistency edits in other chapters that repeat a corrected claim. Do not edit `<files currently under review by a referee>`; record the edits they need as pending in the status record instead.
- `appendices/D-corrections.qmd`, and the script's row in its directory README.
- `appendices/E-original-results.qmd`: this part's items in the list of original results, when a verdict changes what they say.
- The chapter's check script, its figure scripts and their SVGs. Figure scripts hard-code output paths under the repository root (ignoring any `OUT_DIR` environment variable), so run `git restore` on every regenerated SVG whose diff is only ids and timestamps.
- The script directory's README.
- `<the part's status record and conventions brief>`.
- The report (or combined report): append an `## Integration` section only.
- Your checkpoint and your scripts.

## Final checks, then commit

- The chapter script, and any other script you touched, prints only PASS.
- `python3 check-refs.py` passes. Do not run a full Quarto render.
- The `## Integration` section is appended to the report, and the status record is updated.
- Commit the files by name, in one commit (or one per chapter), with a message that lists the verdicts.
- In your final message, list every status label you changed. Mark which upgrades rest on a proof that no referee has seen.
