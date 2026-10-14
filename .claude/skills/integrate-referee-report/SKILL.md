---
name: integrate-referee-report
description: Verify and integrate a referee report on the affine Toda book (from Codex, Gemini or a Claude referee) into the chapters. Use when a briefs/*_report.md has arrived and its findings should go into the book, or when asked to integrate a referee report, review or reader's report.
---

# Integrating a referee report

A referee report is a list of claims to check, not a list of edits to make. Rule 1 of `AGENTS.md` applies to every item. This procedure was used for Parts VII–IX in October 2026. Good models of the record it leaves are commit df973fb and the `## Integration` section of `briefs/referee_defects_41_report.md`.

## 1. Freeze the evidence

1. Make sure the referee has finished. A Codex report ends with its last item. A Gemini coordinator ends with "Referee pass complete."
2. Note the commit the referee reviewed, from the report header. Line numbers in the report refer to that commit or to its review copy, so match passages by content.
3. Commit the report unchanged ("… as received"), staging it by name.
4. Copy the referee's check scripts into `briefs/<report-name>_scripts/`. Leave out report-format validators and rendered pages of papers. Rerun the scripts from there and compare their outputs with those in the report. Some scripts write relative to their cwd (for example into `referee_scratch/`); run those from a scratch directory with symlinks.

## 2. Triage

- **Small reports** (a few items, no proofs) can be integrated in the main session.
- **Large reports** go to one background agent per report, on Opus. Fill in `agent-prompt.md` (next to this file). An integration agent typically runs 30–50 minutes. It must keep a checkpoint file, so that a fresh agent can resume after a usage cutoff; resuming a large agent with SendMessage also works.
- **Files under review.** Name any files another referee is reviewing right now. Nobody edits them; edits they need are recorded as pending in the status file.

## 3. Verify each item

Re-derive by hand, or with a short sympy/mpmath script written from the book's definitions. Import the book's helpers, never the referee's scripts. Patterns seen in past referee output:

- **Reruns are not independent.** A rerun of the worker's own script was presented as an independent check. Once, that "confirmed" an error in a published paper which came from the worker's mislabelled S-matrix.
- **Blanket praise is weak evidence.** Blanket "CONFIRMED" or "no errors" verdicts, especially on qualitative claims, need your own check.
- **Literature corrections.** One was built on a misquote of the source. Before any Appendix D entry, read the source line yourself and recompute with book code.
- **Label inflation.** Referees upgrade convention derivations, and identities between conjectured formulas, to [Theorem]. Don't follow.
- **Superseded items.** Some "missing" results are already proved elsewhere in the book. Search the ids before adding anything.
- **Referee arithmetic.** A referee's own numbers can contradict its conclusion. The rejected q⁻⁴ finding of the Part IX report was one such case.

Mark each item **verified**, **modified** or **rejected**, with the reason.

## 4. Edit

- **Style.** Adapt the referee's text to the book's style rather than pasting it. Never cite the report or `referee_scratch/`, and never mention the referee.
- **Corrections to papers** go in the chapter's Notes and in `appendices/D-corrections.qmd`.
- **Repeated claims.** Fix every chapter that repeats a corrected claim. Find them by grepping the id and the key formula.
- **Original results.** If a verdict changes a result in the list of original results (Appendix E, `appendices/E-original-results.qmd`), update its item there too (see `AGENTS.md`). `check-refs.py` catches a changed label, but not a new result that belongs in the list. Integration agents need `appendices/E-original-results.qmd` in their list of editable files: one agent that lacked it left the entry for Part VIII weaker than a new theorem, and the coordinator had to fix it.
- **Check scripts.** Update the chapter's check script to test the corrected claims, and its row in the directory README and `appendices/C-code.qmd` if the check count or runtime changes.
- **Label upgrades.** Re-derive every upgrade yourself, even when an agent made it and says it verified it.

## 5. Record and finish

1. Append `## Integration` to the report. Give a table of item, verdict, how it was checked, and files changed, and list what was not integrated and why.
2. Update the part's status record if it has one, for example `briefs/part8_defects_status.md`, or the status block of `briefs/halfline_open_problems.md`.
3. Run the chapter scripts (all PASS) and `python3 check-refs.py`. Before a PR, also run `./render-scratch.sh`.
4. Commit by named files. Push or open a PR only when the editor asks. The PR description lists every label upgraded without a second referee pass.
