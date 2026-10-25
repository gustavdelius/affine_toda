# Working rules for AI agents

This repository holds a Quarto book on affine Toda field theory, three companion courses of lecture notes in `lecture-notes/`, and the scripts that check them. `README.md` describes the layout, the build and the scripts. This file holds the rules that every agent working here must follow, whether Claude, Codex or Gemini. Claude Code reads it through `CLAUDE.md`. A task brief can override a rule here only by saying so explicitly.

## Roles

- **Claude drafts the book.** All book text and its check scripts are written by Claude, so that the style stays consistent. The companion courses in `lecture-notes/` and their scripts in `code/lectures/` follow their course briefs (`lecture-notes/defects/` is drafted by Claude, `lecture-notes/half-line/` and `lecture-notes/defects-gemini/` by Gemini).
- **Codex and Gemini referee the book.** If you are a referee, follow your brief and the rules under "Refereeing" below. Do not edit book files.
- **The editor decides.** The editor (Gustav Delius) decides what is merged, and agents commit or push only when asked. `briefs/` holds the record of what was drafted, refereed and checked; `briefs/README.md` is its index.

## Rule 1: verify before you edit

Do not trust any finding in a brief, a referee report or a review, including findings from earlier Claude passes. Before you change the book, check each one yourself: re-derive it, recompute it with your own short script written from the book's definitions, or read the source. Rerunning the author's or the referee's script is a second check, not verification. Record in the report's `## Integration` section which items were verified, which were modified and which were rejected, and why. Referees, including Claude ones, regularly make convention slips and overclaim.

## Writing book text

- **Statements.** Write `::: {#prp-x}` (likewise `thm-`, `lem-`, `cnj-`), with a bold status label first, e.g. `**[Theorem]**`. Proofs go in `::: {.proof}`. A `##` heading inside a theorem div is the theorem's name; do not demote it when moving text.
- **Ids and references.** Ids are `sec-`, `eq-`, `prp-`, `thm-`, `lem-`, `cnj-`, `fig-`, `tbl-`, `exr-`, and `hyp-` for the named hypotheses (R, H, D, X1–X3) of chapter 19 (`part5-semiclassics/19-jordan-chains.qmd`), which are floats `::: {#hyp-D}` around the hypothesis. Quarto takes a float's last block as its caption, so end each `hyp-` float with an empty `[]{}` line; refer with `@hyp-D`, which `hyp-names.lua` prints as "Hypothesis D". Never type a number. Refer with `@id`; write `[-@id]` for the bare number and `(-@eq-x)` for an equation. On unnumbered pages (part pages, `index.qmd`, `preface.qmd`, `citation.qmd`) Quarto prints the target's title instead of its number.
- **Status labels.** `README.md` gives their meanings. Keep `@` references out of the label's bracket: write `**[Numerical]** (@sec-x)`. Pandoc reads `[Numerical, @sec-x]`, and any bracket whose `;`-separated parts all contain an `@key`, as a citation. The brackets are then lost and the label is not coloured. A new status word must be added to both `status-tags.lua` (which matches on the first word) and `STATUS_WORDS` / `status_kind()` in `check-refs.py`.
- **Changing a label.** Upgrade to [Theorem] only with a complete proof in the text. Keep it only after an independent referee pass. Every upgrade made without one must be listed in the PR description.
- **Style.** Use short declarative sentences and define every symbol. Never cite `briefs/`, `referee_scratch/` or a report in the book, and never mention "the referee".
- **What the book is.** The book is a review, with new results, not a textbook or a course; do not describe it as one. Step-by-step teaching belongs in the lecture notes in `lecture-notes/`.
- **Prior art.** Before calling a result new, search the literature, including the editor's own papers. Several "new" R-matrices turned out to be in Delius–Gould–Zhang.
- **Corrections to published papers.** Read the source line yourself, quote it exactly, and recompute with the book's code. A misquoted prefactor once produced a false "correction". Each correction goes in the chapter's Notes and in Appendix D (`appendices/D-corrections.qmd`), in the same change. Update or drop the Appendix D entry whenever a chapter gains or loses a correction.
- **The list of original results.** Appendix E (`appendices/E-original-results.qmd`) lists the book's original results by part, in a `::: {.original-results}` div. Each item opens with a status label and ends with the ids of the statements it summarises: `- **[Theorem]** … (@prp-x).` The label may shorten a long one, but must keep its first word and its colour. Keep the list in step with the chapters, in the same change. Add an item for a new original result. Update or drop an item when its result is proved, refuted or relabelled, or when it turns out to be in earlier literature (corrections to papers go in Appendix D, not here, unless the item itself states a new result that overturns a paper's claim). `check-refs.py` fails when an item's label no longer matches its statement's; it cannot notice a result that is missing from the list.
- **Moving or renumbering chapters.**
  - File names carry global chapter numbers, figure SVGs and figure scripts are named `chNN_fM_*`, and `code/figures/style.py` has one `OUT_PARTn` per part (see `briefs/superalgebra_toda_renumber.py`).
  - Regenerate the chapter list in `_quarto.yml`.
  - Make a pure-rename commit first, then a text commit.
  - Dated briefs keep the numbering of their date; `briefs/README.md` maps the old numbers.
- **Lecture notes (`lecture-notes/`).**
  - Each course (`lecture-notes/defects/`, `lecture-notes/defects-gemini/`, `lecture-notes/half-line/`) is a standalone Quarto project with its own `_quarto.yml`.
  - Refer to the book with `[text](book:<id>)`, which `lecture-notes/book-links.lua` turns into a link to the published chapter and `python3 check-refs.py lecture-notes/<name>` checks; never type a book chapter or statement number.
  - Put lecture check and figure scripts in `code/lectures/`, written from the lectures' own definitions, and index them in `code/lectures/README.md`.
  - Bibliography entries used only by a course go in its own `.bib` file (e.g. `references-notes.bib`), because the book's `references.bib` uses `nocite: @*`. Keep book edits in a separate PR from lecture notes.

## Conventions and known traps

The book's conventions are in `01-introduction.qmd`, section `sec-conventions`. Part-specific conventions and clash tables are in briefs such as `briefs/part8_defects_conventions.md`. Each of these traps has caused a wrong claim before:

- **The sign of q.** The physical q is fixed in `eq-q-physical` (ch. 21): q = e^{−4π²i/β²} = −e^{−iπω}, with ω = 4π/β² − 1. Part VI's coproduct is opposite to chapter 7's, so R_VI(x;q) = R_7(x;1/q). The map q → −q is a Klein-factor twist, not a gauge. "Unchanged under q → −q" usually holds only for moduli or for even soliton number.
- **Twisted algebras.** Some literature gives the longest root of X_n^(k) squared length 2k, so β²_lit = β²/k (and in the twisted crossing period 2T = 4πkh/β² − h^∨). MacKay–Watts (`mackay_watts1995`) use non-coprime twisted Kac labels with ∑ n_a = kh, whereas the book divides out the common factor c so ∑ n_a = h (rescaling the folded mass parameter by m_F² = c m_P², with c = 2 for d_{2n+2}^{(1)} → a_{2n}^{(2)}). The algebra a_{2n}^{(2)} is an exception in several places: in the book's Toda tables α_0 has squared length 2 (n_0 = 1, n_0^∨ = 2) and α_1, …, α_n have squared lengths 1 and 1/2, whereas in Kac's labelling (used in ch. 7) α_0 has squared length 1 (n_0 = 2, n_0^∨ = 1) and the S² shift is q^{h^∨}, not q^{2h^∨}; and because k = k' = 2, self-duality of a_{2n}^{(2)} maps β → β^∨ = 8π/β, not 4π/β.
- **The Bethe–Yang operator.** It is R in some sources and P·Ř in others; check which one a formula uses.
- **Bibliography keys.** `dgz1992` is Delius–Grisaru–Zanon and `dgz1996` is Delius–Gould–Zhang. In prose, spell out the names.
- **Half-line notation.** Ghoshal–Zamolodchikov's boundary ξ clashes with the book's ξ = π/λ; use the (η, ϑ) of `eq-bsg-eta-theta`. In sine-Gordon the coideal and boundary labels cross: ε̂_1 goes with C_0 and ε̂_0 with C_1 (`sec-bsg-parameters`).
- **Roots of unity.** q is a root of unity exactly when ω is rational. Short decimals such as 0.3, 0.7 or 1.37 are therefore roots of unity, where representations degenerate and R-matrix code silently returns garbage. Use generic values such as √2 or e/2.

## Check scripts and numerics

- **Every new claim gets a check.** A claimed formula or [Numerical] statement needs a check in `code/` or one of its subdirectories (see `README.md`). Each new script prints one `PASS` or `FAIL` line per check. Add or update the script's row in `code/README.md` (or its subdirectory README, or `code/lectures/README.md` for lecture scripts), which is the only index of the scripts.
- **Cite the script at the claim.** Name it in the text where the claim is made, as a full repository path in inline code, e.g. `` `code/hirota/fold.py` ``. `script-links.lua` turns such paths into GitHub links, and `check-refs.py` fails on a path that does not exist. Appendix C (`appendices/C-code.qmd`) describes the code for the reader and points to the README; it lists no scripts, so a new script needs no edit there.
- **Shared helpers.** Only append to shared helpers such as `code/defect_common.py`, then rerun every script that imports them.
- **Figures.** Follow `briefs/figures_for_the_book.md`. Figure scripts hard-code output paths under the repository root (ignoring any `OUT_DIR` environment variable), so running one overwrites tracked SVGs in place; referees must not run them in a shared worktree. After a rerun, `git restore` every SVG whose diff is only ids and timestamps.
- **Waiting on processes.** Never use `pgrep -f` or `pkill -f` with a pattern that appears in your own command line: it matches the shell running it, and wait loops then never end. Wait on a PID (`while kill -0 $pid; do sleep 5; done`) and kill by PID.
- **Output buffering.** Python output redirected to a file is buffered, so an empty log does not mean no progress; use `python3 -u`.
- **Caches** must be keyed on every parameter that changes the result. A cache keyed on ω alone once silently reused the other sign of q.
- **Poles.** Small residues hide poles. Count pole orders from the structure, then confirm with 1/ε scaling, not with a print threshold.
- **Large tensor products.** For dimensions of 64 and up, a dense intertwiner solve is unusable; restrict the unknowns to weight-preserving entries.
- **ODE scans.** Start from vectors normalized to 1 with atol of about 1e-12. Exponentially scaled starts with atol 1e-30 once ran for hours and were killed for memory.

## Building and checking the book

- **`python3 check-refs.py`** checks that every id is defined once, every cross-reference resolves, every citation key is in `references.bib`, no status label is parsed as a citation, every item in the list of original results carries the label of the statement it cites, and every script path cited in inline code exists. Pass a lecture-notes directory (`python3 check-refs.py lecture-notes/<name>`) to check a course and its `[text](book:<id>)` links into the book. It takes under a second. Run it before every commit that touches a `.qmd` file.
- **`./render-scratch.sh`** (HTML), **`./render-scratch.sh --pdf`** (HTML and PDF, about 7 minutes) or **`./render-scratch.sh --project lecture-notes/<name>`** renders in a scratch copy and lists unresolved references. Never run `quarto render` in a checkout where `quarto preview` may be running; the editor usually has one open. `--output-dir` does not prevent the race.
- **LaTeX fixes.** Test them by running `lualatex` on the `affine-toda.tex` left in the scratch copy. The math uses `\rm` and `\bf`, so the KOMA class option `enabledeprecatedfontcommands` must stay. Mermaid figures need `%%| fig-width:` of 6 inches or less.
- **Website styling.**
  - `theme.scss` serves both themes; `theme-dark.scss` overrides `$link-color` and `$mermaid-label-fg-color` in `scss:defaults` for the dark theme.
  - Theme variables such as `$body-bg` exist only in `scss:rules`, not in `scss:defaults`.
  - `scroll-padding-top` keeps link targets clear of the title bar. Do not add `scroll-behavior: smooth`, which broke fragment links.
- **Bibliography backlinks** come from the post-render script `add-backlinks.py`. A Lua filter or citetools cannot do it, because citeproc runs last and refills the references page.
- **Visual checks.** Serve the scratch `_site` with `python3 -m http.server`. Pages opened with `file://` cannot be screenshotted, and the first screenshot after navigating is often blank.
- **CI and the PDF.** CI renders HTML only (the book and each `lecture-notes/<name>/`), on pushes to `main`. The PDFs are built locally by the editor with `publish-pdf.sh`.

## Git and coordination

- **Where to work.** Start new work in a git worktree (under `.claude/worktrees/` for Claude Code, or Antigravity's worktrees directory). The main checkout often holds another session's uncommitted files and a running preview.
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
- **Skills.** Repository skills `integrate-referee-report` and `write-referee-brief` live in `.claude/skills/` and are symlinked in `.agents/skills/` for Codex and Gemini.

## Refereeing (Codex, Gemini, or a Claude referee)

- **Write access.** Edit no existing file, and do not run figure scripts in a shared worktree. Write scratch work only under `referee_scratch/<tag>/` and the report only where the brief says (`briefs/<name>_report.md`).
- **Commands.** Run no git command that changes the index or the working tree (commit, checkout, stash, reset, rebase, restore), and no `quarto render`.
- **Two-phase passes.** A referee pass may be split into a blind phase 1 and a directed phase 2 (see `briefs/referee_classical_09_11.md`). In phase 1 the brief gives only context, conventions, sources and the evidence rule; the referee writes its own ranked list of risky claims, checks them, and writes `briefs/<name>_phase1_report.md`. Only after phase 1 is committed is the held items list given for phase 2, which aims at prior art, citation scope, interpretive sentences, and claims that a script checks something.
- **Evidence.** For every item, do one independent computation: a derivation written out, or a script of your own. The book's scripts are not evidence. A verdict without its own evidence counts as "unclear".
- **Literature claims.** Quote the exact source line.
- **Verdicts.** Give a verdict per item (correct / wrong / unclear), your computation and its output, and, for errors, replacement text with file and line.
- **Labels.** Recommend [Theorem] only for a complete proof. Identities between conjectured formulas, and convention derivations, are not theorems.
- **Report header.** Give the date, your model and the commit you reviewed.
