---
name: write-referee-brief
description: Write a self-contained referee brief on chapters or results of the affine Toda book, which the editor can hand to Codex (GPT6-Astra) or Gemini as it stands. Use when asked for a referee brief or prompt, or when new chapters, label upgrades or literature corrections need an independent check before their labels are kept.
---

# Writing a referee brief

Claude drafts and other models referee (`AGENTS.md`, "Roles"). The referee can be started interactively or through `briefs/run_referee_agy.sh <brief-name> [codex|gemini|both|<model>]`, which runs `codex exec` or `agy -p` in a throwaway copy under `~/referee-copy/`. In either case the brief must work with no context but the repository.

## Choosing the pass structure and the referee model

- **Two-phase pass (default for chapter or part reviews).** A brief that lists what to check carries the drafting model's view of where the risks are, so a referee tends to find only what Claude already suspected. Split chapter/part passes into:
  - **Phase 1 (blind):** `briefs/referee_<topic>.md` gives only context, conventions, sources and the evidence rule. The referee reads the chapters cold, writes its own ranked risk list (15–30 entries) before checking, checks at least the top 10, writes `briefs/referee_<topic>_phase1_report.md`, and stops with `"Phase 1 complete"`. Model: `briefs/referee_classical_09_11.md`. Run **Codex** (`codex`), or **Codex and Gemini in parallel** (`both`) for maximum coverage (Gemini's blind pass catches complementary geometric and analytic issues).
  - **Phase 2 (directed):** `briefs/referee_<topic>_items.md` gives the coordinator's pre-written items, aimed at what blind passes miss: prior art, Appendix E novelty claims, citation scope/fidelity, interpretive sentences, exercises, and claims that a script checks something. Hold this file outside the repository until the phase-1 report is committed (recording its SHA-256 in `briefs/README.md` if desired). Output goes to `briefs/referee_<topic>_report.md`, with an overlap table comparing the two lists. Model: `briefs/referee_classical_09_11_items.md`. Run **Codex** (`briefs/run_referee_agy.sh referee_<topic>_items codex`), which has the highest precision on literature, prior art and label discipline.
- **Single-phase pass (for narrow targeted checks).** When only a few specific `[Theorem]` upgrades, literature corrections or consolidation items need checking, put the items directly in one brief `briefs/referee_<topic>.md`. Good models are `briefs/referee_defects_38_l6.md` (three theorem upgrades) and `briefs/referee_defects_43.md` (one chapter). Run **Codex** by default, **Gemini** (`gemini-4-argon-high` or `gemini-3.8-flash-high`) for fast triage when Codex quota is scarce, or **Claude (Opus)** only when both are unavailable.

## Where things go

- **Brief:** `briefs/referee_<topic>.md` (and `briefs/referee_<topic>_items.md` for phase 2 of a two-phase pass), or `briefs/<round>_referee_prompt.md` within a round of briefs.
- **Reports:** `briefs/referee_<topic>_phase1_report.md` (phase 1) and `briefs/referee_<topic>_report.md` (phase 2 or single-phase). The brief names them exactly.
- **Referee scratch work:** `referee_scratch/<tag>/` (untracked, never committed).

## Structure of the brief

1. **Role and place.** "You are an adversarial mathematical-physics referee", plus the worktree path and branch the referee works in. For a two-phase pass, explain the two phases up front.
2. **Rules block.** Copy this verbatim, filling in the tag and name:
   - Other agents may be editing files in this worktree at the same time.
   - Do not edit any existing file, and do not run figure scripts in a shared worktree (they overwrite tracked SVGs in place).
   - Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
   - Do not run `quarto render`.
   - Do not spawn subagents.
   - Write scratch scripts only under `referee_scratch/<tag>/`.
   - Write the report to `briefs/<name>_phase1_report.md` (or `briefs/<name>_report.md`).
   - *(Two-phase phase 1 only)* Until phase 1 is complete, do not read earlier referee briefs, reports or reading-fix briefs on these chapters (name the files explicitly), and do not look for the phase-2 list.
3. **Context.**
   - What is under review: files, section ids, the drafting commit (and `git diff <last-refereed-commit> HEAD -- <files>` if applicable).
   - Statements carrying status labels, and any items in Appendix E (`appendices/E-original-results.qmd`) or Appendix D (`appendices/D-corrections.qmd`) for this part.
   - The check scripts, with their check counts and runtimes (noting that they are the drafter's checks, not evidence).
   - *(Single-phase or phase 2)* The drafting report, and the drafter's own flags for the referee from the status record.
4. **Notation block.** The conventions in force: the coupling, the physical q, the coproduct, the tensor order, and any symbols that clash with other parts of the book. The referee will not know the book's conventions, and convention slips are the most common false findings.
5. **Read first.**
   - `01-introduction.qmd` (`sec-conventions`), `README.md` ("Status labels"), the conventions brief, and the relevant background sections by id.
   - The source papers by arXiv number. On the editor's machine their texts are in `~/referee-copy/papers/*.txt`, and `run_referee_agy.sh` copies them to `papers/`. If a source cannot be read, the referee must mark that citation "unclear" rather than reconstructing what it probably says.
   - How to find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.
6. **Evidence rule.**
   - For each item, one independent computation: a derivation written out, or a script of the referee's own.
   - The book's scripts may be run for comparison but are not evidence.
   - A verdict without its own evidence counts as "unclear".
   - Name a past report that shows the expected depth, e.g. `briefs/referee_defects_38_40_report.md`.
7. **What to check.**
   - **In a two-phase pass (Phase 1):** Instruct the referee to read the chapters in full first, write a ranked risk list (15–30 entries with location, quoted claim, why it is at risk, and how to check it) as section 1 before checking anything, check at least the top 10 in rank order, list unchecked entries and why, and stop with `"Phase 1 complete"`.
   - **In a single-phase brief or Phase 2 (`*_items.md`):** List the items, most important first. Each item states:
     - the exact claim and its current label;
     - the computation that would confirm or refute it;
     - "is the label justified?".

     Always include an item on prior art and novelty claims, covering the part's items in the list of original results (Appendix E, `appendices/E-original-results.qmd`), and one on whether each citation says what the source says. Weight the items toward what changed or what blind passes miss: label upgrades without a referee pass, new [Theorem]s, new Appendix D entries, citation scope, interpretive sentences, and claims that a book script checks something. In Phase 2, ask the referee to note for each item whether its Phase 1 risk list already covered it, check what Phase 1 missed, and include an overlap table.
8. **Report format.**
   - A header with the date, the model and the commit reviewed.
   - For Phase 1: section 1 (the risk list as written before checking), section 2 (one subsection per checked entry: verdict correct / wrong / unclear, computation and output, and replacement text with file and line for errors), section 3 (entries not checked and why), section 4 (other observations, labelled as unverified).
   - For single-phase or Phase 2: one section per item (with the Phase 1 overlap note in Phase 2), the overlap table and summary of what each list found (in Phase 2), claims not checked and why, and other observations labelled as unverified.

## Before handing it over

- **Commit the brief** (for a two-phase pass, commit the Phase 1 brief while holding `*_items.md` outside the repository until Phase 1 is committed), and put the reviewed commit in the brief.
- **Record the files under review** in the part's status record or `briefs/README.md`, and tell the editor. Nobody edits them until the reports are in.
- **Gemini briefs.** Gemini has been over-positive: blanket "CONFIRMED" verdicts, and label upgrades. Ask it explicitly to list the claims it could not check. The coordinator, not Gemini, chooses which items get a follow-up Opus pass.
- **If Codex is out of quota** and the editor asks Claude to run the brief, do so. Say in the report header that the pass is not model-diverse.

When the report arrives, use the `integrate-referee-report` skill.
