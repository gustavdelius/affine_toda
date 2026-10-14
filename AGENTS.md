# Working rules for AI agents

This repository is a Quarto book on affine Toda field theory, with the scripts that check it. `README.md` describes the layout, the build and the scripts. This file holds the rules that every agent working here must follow, whether Claude, Codex or Gemini. Claude Code reads it through `CLAUDE.md`. A task brief can override a rule here only by saying so explicitly.

## Roles

- **Claude drafts.** All book text and its check scripts are written by Claude, so that the style stays consistent.
- **Codex and Gemini only referee.** If you are a referee, follow your brief and the rules under "Refereeing" below. Do not edit book files.
- **The editor decides.** The editor (Gustav Delius) decides what is merged, and agents commit or push only when asked. `briefs/` holds the record of what was drafted, refereed and checked; `briefs/README.md` is its index.

## Rule 1: verify before you edit

Do not trust any finding in a brief, a referee report or a review, including findings from earlier Claude passes. Before you change the book, check each one yourself: re-derive it, recompute it with your own short script written from the book's definitions, or read the source. Rerunning the author's or the referee's script is a second check, not verification. Record in the report's `## Integration` section which items were verified, which were modified and which were rejected, and why. Referees, including Claude ones, regularly make convention slips and overclaim.

## Writing book text

- **Statements.** Write `::: {#prp-x}` (likewise `thm-`, `lem-`, `cnj-`), with a bold status label first, e.g. `**[Theorem]**`. Proofs go in `::: {.proof}`. A `##` heading inside a theorem div is the theorem's name; do not demote it when moving text.
- **Ids and references.** Ids are `sec-`, `eq-`, `prp-`, `thm-`, `lem-`, `cnj-`, `fig-`, `tbl-`, `exr-`, and `hyp-` for the named hypotheses (R, H, D, X1–X3) of chapter 18, which are floats `::: {#hyp-D}` around the hypothesis. Quarto takes a float's last block as its caption, so end each `hyp-` float with an empty `[]{}` line; refer with `@hyp-D`, which `hyp-names.lua` prints as "Hypothesis D". Never type a number. Refer with `@id`; write `[-@id]` for the bare number and `(-@eq-x)` for an equation. On unnumbered pages (part pages, `citation.qmd`) Quarto prints the target's title instead of its number.
- **Status labels.** `README.md` gives their meanings. Keep `@` references out of the label's bracket: write `**[Numerical]** (@sec-x)`. Pandoc reads `[Numerical, @sec-x]`, and any bracket whose `;`-separated parts all contain an `@key`, as a citation. The brackets are then lost and the label is not coloured. A new status word must be added to `status-tags.lua`, which matches on the first word.
- **Changing a label.** Upgrade to [Theorem] only with a complete proof in the text. Keep it only after an independent referee pass. Every upgrade made without one must be listed in the PR description.
- **Style.** Use short declarative sentences and define every symbol. Never cite `briefs/`, `referee_scratch/` or a report in the book, and never mention "the referee".
- **Prior art.** Before calling a result new, search the literature, including the editor's own papers. Several "new" R-matrices turned out to be in Delius–Gould–Zhang.
- **Corrections to published papers.** Read the source line yourself, quote it exactly, and recompute with the book's code. A misquoted prefactor once produced a false "correction". Each correction goes in the chapter's Notes and in Appendix D (`appendices/D-corrections.qmd`), in the same change. Update or drop the Appendix D entry whenever a chapter gains or loses a correction.
- **The list of original results.** Appendix E (`appendices/E-original-results.qmd`) lists the book's original results by part, in a `::: {.original-results}` div. Each item opens with a status label and ends with the ids of the statements it summarises: `- **[Theorem]** … (@prp-x).` The label may shorten a long one, but must keep its first word and its colour. Keep the list in step with the chapters, in the same change. Add an item for a new original result. Update or drop an item when its result is proved, refuted or relabelled, when it turns out to be in earlier literature, and when a correction to a paper is added or dropped. `check-refs.py` fails when an item's label no longer matches its statement's; it cannot notice a result that is missing from the list.
- **Moving or renumbering chapters.**
  - File names carry global chapter numbers, figure scripts are named `chNN_fM_*.py`, and `textbook_code/figures/style.py` has one `OUT_PARTn` per part.
  - Regenerate the chapter list in `_quarto.yml`.
  - Make a pure-rename commit first, then a text commit.
  - Dated briefs keep the numbering of their date; `briefs/README.md` maps the old numbers.

## Conventions and known traps

The book's conventions are in `01-introduction.qmd`, section `sec-conventions`. Part-specific conventions and clash tables are in briefs such as `briefs/part8_defects_conventions.md`. Each of these traps has caused a wrong claim before:

- **The sign of q.** The physical q is fixed in `eq-q-physical` (ch. 20): q = e^{−4π²i/β²} = −e^{−iπω}, with ω = 4π/β² − 1. Part VI's coproduct is opposite to chapter 7's, so R_VI(x;q) = R_7(x;1/q). The map q → −q is a Klein-factor twist, not a gauge. "Unchanged under q → −q" usually holds only for moduli or for even soliton number.
- **Twisted algebras.** Some literature gives the longest root of X_n^(k) squared length 2k, so β²_lit = β²/k.
- **The Bethe–Yang operator.** It is R in some sources and P·Ř in others; check which one a formula uses.
- **Bibliography keys.** `dgz1992` is Delius–Grisaru–Zanon and `dgz1996` is Delius–Gould–Zhang. In prose, spell out the names.
- **Half-line notation.** Ghoshal–Zamolodchikov's boundary ξ clashes with the book's ξ = π/λ; use the (η, ϑ) of `eq-bsg-eta-theta`. In sine-Gordon the coideal and boundary labels cross: ε̂_1 goes with C_0 and ε̂_0 with C_1 (`sec-bsg-parameters`).
- **Roots of unity.** q is a root of unity exactly when ω is rational. Short decimals such as 0.3, 0.7 or 1.37 are therefore roots of unity, where representations degenerate and R-matrix code silently returns garbage. Use generic values such as √2 or e/2.

## Check scripts and numerics

- **Every new claim gets a check.** A claimed formula or [Numerical] statement needs a check in `textbook_code/` or one of its subdirectories (see `README.md`). Each new script prints one `PASS` or `FAIL` line per check. Add or update the script's row in `textbook_code/README.md` (or its subdirectory README) and in `appendices/C-code.qmd`.
- **Cite the script at the claim.** Name it in the text where the claim is made, as a full repository path in inline code, e.g. `` `textbook_code/hirota/fold.py` ``. `script-links.lua` turns such paths into GitHub links, and `check-refs.py` fails on a path that does not exist. Appendix C indexes every script, one section per part.
- **Shared helpers.** Only append to shared helpers such as `textbook_code/defect_common.py`, then rerun every script that imports them.
- **Figures.** Follow `briefs/figures_for_the_book.md`. Figure scripts write SVGs into the repository whatever `OUT_DIR` says. After a rerun, `git restore` every SVG whose diff is only ids and timestamps.
- **Waiting on processes.** Never use `pgrep -f` or `pkill -f` with a pattern that appears in your own command line: it matches the shell running it, and wait loops then never end. Wait on a PID (`while kill -0 $pid; do sleep 5; done`) and kill by PID.
- **Output buffering.** Python output redirected to a file is buffered, so an empty log does not mean no progress; use `python3 -u`.
- **Caches** must be keyed on every parameter that changes the result. A cache keyed on ω alone once silently reused the other sign of q.
- **Poles.** Small residues hide poles. Count pole orders from the structure, then confirm with 1/ε scaling, not with a print threshold.
- **Large tensor products.** For dimensions of 64 and up, a dense intertwiner solve is unusable; restrict the unknowns to weight-preserving entries.
- **ODE scans.** Start from vectors normalized to 1 with atol of about 1e-12. Exponentially scaled starts with atol 1e-30 once ran for hours and were killed for memory.

## Building and checking the book

- **`python3 check-refs.py`** checks that every id is defined once, every cross-reference resolves, every citation key is in `references.bib`, no status label is parsed as a citation, and every item in the list of original results carries the label of the statement it cites. It takes under a second. Run it before every commit that touches a `.qmd` file.
- **`./render-scratch.sh`** (HTML) or **`./render-scratch.sh --pdf`** (HTML and PDF, about 7 minutes) renders in a scratch copy and lists unresolved references. Never run `quarto render` in a checkout where `quarto preview` may be running; the editor usually has one open. `--output-dir` does not prevent the race.
- **LaTeX fixes.** Test them by running `lualatex` on the `index.tex` left in the scratch copy. The math uses `\rm` and `\bf`, so the KOMA class option `enabledeprecatedfontcommands` must stay. Mermaid figures need `%%| fig-width:` of 6 inches or less.
- **Website styling.**
  - `theme.scss` serves both themes.
  - Theme variables such as `$body-bg` exist only in `scss:rules`, not in `scss:defaults`.
  - `scroll-padding-top` keeps link targets clear of the title bar. Do not add `scroll-behavior: smooth`, which broke fragment links.
- **Bibliography backlinks** come from the post-render script `add-backlinks.py`. A Lua filter or citetools cannot do it, because citeproc runs last and refills the references page.
- **Visual checks.** Serve the scratch `_site` with `python3 -m http.server`. Pages opened with `file://` cannot be screenshotted, and the first screenshot after navigating is often blank.
- **CI and the PDF.** CI renders HTML only, on pushes to `main`. The PDF is built locally by the editor with `publish-pdf.sh`.

## Git and coordination

- **Where to work.** Start new work in a git worktree under `.claude/worktrees/` (git-excluded). The main checkout often holds another session's uncommitted files and a running preview.
- **Shared worktrees.** Several agents, and a referee, may share one worktree:
  - never use bare `git stash`, because the stash stack is shared by all worktrees;
  - never use `git add -A` or `git add .`; stage files by name;
  - do not edit a file that a referee is currently reviewing.
- **Merging main** into a branch: do it only when no agent is editing that worktree, then rerun `check-refs.py` and render.

## Running agents (for coordinating sessions)

- **Concurrency.** Run at most two or three agents at a time; more hit the usage limits.
- **Every agent prompt says "Do not spawn subagents."** Referees that spawned helpers once exhausted the quota overnight.
- **Checkpoints.** Long agents append one line per finished item to `briefs/<task>_checkpoint.md`, so that a fresh agent can resume from it and `git diff`.
- **Scratch space.** Give each worker its own scratch subdirectory. Copy its scripts to `briefs/<task>_scripts/` as soon as it reports, because scratchpads are temporary.
- **Skills.** Claude Code has the repository skills `integrate-referee-report` and `write-referee-brief` in `.claude/skills/`.

## Refereeing (Codex, Gemini, or a Claude referee)

- **Write access.** Edit no existing file. Write scratch work only under `referee_scratch/<tag>/` and the report only where the brief says (`briefs/<name>_report.md`).
- **Commands.** Run no git command that changes the index or the working tree (commit, checkout, stash, reset, rebase, restore), and no `quarto render`.
- **Evidence.** For every item, do one independent computation: a derivation written out, or a script of your own. The book's scripts are not evidence. A verdict without its own evidence counts as "unclear".
- **Literature claims.** Quote the exact source line.
- **Verdicts.** Give a verdict per item (correct / wrong / unclear), your computation and its output, and, for errors, replacement text with file and line.
- **Labels.** Recommend [Theorem] only for a complete proof. Identities between conjectured formulas, and convention derivations, are not theorems.
- **Report header.** Give the date, your model and the commit you reviewed.
