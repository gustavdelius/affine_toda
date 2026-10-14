---
name: write-referee-brief
description: Write a self-contained referee brief on chapters or results of the affine Toda book, which the editor can hand to Codex (GPT6-Astra) or Gemini as it stands. Use when asked for a referee brief or prompt, or when new chapters, label upgrades or literature corrections need an independent check before their labels are kept.
---

# Writing a referee brief

Claude drafts and other models referee (`AGENTS.md`, "Roles"). The editor starts the referee by hand: Codex by default, Gemini through `briefs/run_referee_agy.sh <brief-name>` on request. So the brief must work with no context but the repository. Good models are `briefs/referee_defects_43.md` (short, deep) and `briefs/referee_defects_38_40.md` (two chapters).

## Where things go

- **Brief:** `briefs/referee_<topic>.md`, or `briefs/<round>_referee_prompt.md` within a round of briefs.
- **Report:** `briefs/referee_<topic>_report.md`. The brief names it exactly.
- **Referee scratch work:** `referee_scratch/<tag>/` (untracked, never committed).

## Structure of the brief

1. **Role and place.** "You are an adversarial mathematical-physics referee", plus the worktree path and branch the referee works in.
2. **Rules block.** Copy this verbatim, filling in the tag and name:
   - Other agents may be editing files in this worktree at the same time.
   - Do not edit any existing file.
   - Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
   - Do not run `quarto render`.
   - Do not spawn subagents.
   - Write scratch scripts only under `referee_scratch/<tag>/`.
   - Write the report to `briefs/<name>_report.md`.
3. **Context.**
   - What is under review: files, section ids, the drafting commit.
   - The check scripts, with their check counts and runtimes.
   - The drafting report, and the drafter's own flags for the referee from the status record.
4. **Notation block.** The conventions in force: the coupling, the physical q, the coproduct, the tensor order. The referee will not know the book's conventions, and convention slips are the most common false findings.
5. **Read first.**
   - The conventions brief, and the relevant sections by id.
   - The source papers by arXiv number. On the editor's machine their texts are in `~/referee-copy/papers/*.txt`, and `run_referee_agy.sh` copies them to `papers/`.
   - How to find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.
6. **Evidence rule.**
   - For each item, one independent computation: a derivation written out, or a script of the referee's own.
   - The book's scripts may be run for comparison but are not evidence.
   - A verdict without its own evidence counts as "unclear".
   - Name a past report that shows the expected depth, e.g. `briefs/referee_defects_38_40_report.md`.
7. **Items, most important first.** Each item states:
   - the exact claim and its current label;
   - the computation that would confirm or refute it;
   - "is the label justified?".

   Always include an item on prior art and novelty claims, covering the part's items in the list of original results (Appendix E, `appendices/E-original-results.qmd`), and one on whether each citation says what the source says. Weight the items toward what changed: label upgrades without a referee pass, new [Theorem]s, and new Appendix D entries.
8. **Report format.**
   - A header with the date, the model and the commit reviewed.
   - One section per item: the verdict (correct / wrong / unclear); the referee's own computation (script inline or path, and its output); and, for errors, replacement text ready to paste, with file and line.
   - End with other observations, labelled as unverified.

## Before handing it over

- **Commit the brief,** so that the referee reviews a fixed state, and put the reviewed commit in the brief.
- **Record the files under review** in the part's status record, and tell the editor. Nobody edits them until the report is in.
- **Gemini briefs.** Gemini has been over-positive: blanket "CONFIRMED" verdicts, and label upgrades. Ask it explicitly to list the claims it could not check. The coordinator, not Gemini, chooses which items get a follow-up Opus pass.
- **If Codex is out of quota** and the editor asks Claude to run the brief, do so. Say in the report header that the pass is not model-diverse.

When the report arrives, use the `integrate-referee-report` skill.
