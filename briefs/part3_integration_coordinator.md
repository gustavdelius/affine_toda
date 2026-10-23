# Part III referee integration: coordinator's state and restart procedure

Written 2026-10-15 by the coordinating Claude session, so that the integration survives a usage-limit stop of the coordinator and of its agents. If you are a coordinator restarting after a stop, read this file first.

## State

- **Branch `part3-referee-integration`** (worktree `.claude/worktrees/part3-referee-integration`), from `origin/main` a03edd4.
  - 45955bb: the four referee reports as received, the phase-2 items, both referees' scripts and the combined phase-1 report.
  - 3b38470: reruns of the referee scripts. All of Codex's reproduce exactly, and Gemini's run cleanly.
- **Agent A** (chapters 9–10, `toda_classical.py`, figures ch9/ch10): worktree `.claude/worktrees/part3-int-ch9-10`, branch `part3-int-ch9-10`. Checkpoint: `briefs/part3_integration_A_checkpoint.md`. Scratch: `<scratchpad>/part3_int_A/`.
- **Agent B** (chapter 11, Appendix E, `references.bib`, H6 in ch. 44, `classical_reality.py`, figure ch11_f16): worktree `.claude/worktrees/part3-int-ch11`, branch `part3-int-ch11`. Checkpoint: `briefs/part3_integration_B_checkpoint.md`. Scratch: `<scratchpad>/part3_int_B/`.
- `<scratchpad>` = `/tmp/claude-1000/-home-gustav-Git-affine-toda/7fa1324b-155d-4b76-8703-4c3aaf6d6c0d/scratchpad`. It is temporary. If it has gone, the agents' scratch scripts are lost, and a restarted agent must recreate the ones it needs.
- Both agents were launched as background general-purpose agents on Opus. They commit only at the end, on their own branches.

## Restart procedure (after a usage-limit stop)

1. **Find the agents.** Use `ListAgents`; their descriptions are "Integrate Part III referee: ch 9-10" and "Integrate Part III referee: ch 11". Note which ones finished. A finished agent has a commit on its branch beyond 45955bb and a final message.
2. **Check for leftover processes** from the agents, such as `toda_classical.py`, `smatrix_simply_laced.py` or `algebra_data.py` reruns. Use `ps -ef | grep python3`, and kill by PID only, never with `pkill -f`.
3. **For each unfinished agent**, read its checkpoint and `git -C <worktree> status --short; git -C <worktree> diff --stat`. Then choose:
   - **Resume** with SendMessage if the agent stopped in the middle of an item that needs much context (e.g. A on item 12, the Hollowood sign, or B on the prior-art reading). Message: "The usage limit has reset. Continue where you stopped. First re-read your checkpoint and `git diff`, and check that no process of yours is still running."
   - **Start a fresh agent** if it stopped between items or after a long run (resuming replays the whole context, which is expensive). Use the same prompt (below), with the resume preamble prepended.
4. Do not run more than these two agents at once.
5. When both have finished, do the coordinator steps at the end of this file.

### Resume preamble for a fresh agent

> A previous agent started this task and was stopped by a usage limit. Its work is in your worktree: read `briefs/part3_integration_<A|B>_checkpoint.md`, run `git status` and `git diff`, and list your scratch directory. Do not redo items marked complete in the checkpoint. Re-verify the item that was half done before you build on it. Rule 1 still applies: a previous agent's conclusions are claims to check, unless the checkpoint records your own verification script and its output. Continue appending to the same checkpoint.

## Prompt for agent A (chapters 9–10)

The prompt as launched is reproduced in `part3_integration_A_prompt.md` next to this file.

## Prompt for agent B (chapter 11)

The prompt as launched is reproduced in `part3_integration_B_prompt.md` next to this file.

## Coordinator decisions the editor may reverse

- The Θ selection rule becomes a labelled statement **[Established; reformulated]**, not [Theorem]. The general mechanism is prior art (Cen–Fring; Cen–Correa–Fring; Fring–Taira; Correa–Fring–Taira), and the energy case is in Olive–Turok–Underwood and Takács–Watts.
- The Appendix E item for chapter 11 is dropped. The exception: if agent B verifies something genuinely new, it becomes a [Numerical] item. The candidate is that the a₂ spin-2 and spin-4 charges are imaginary and additive on breathers.
- Open question: the sign of the semiclassical Hollowood S-matrix phase against $(2h/\beta^2)\ln A$ (agent A, item 12). If chapters 10 and 21 genuinely disagree, ask the editor before editing chapter 21.

## Coordinator steps after both agents finish

1. Read each agent's final message and `briefs/part3_integration_<A|B>.md`.
2. **Rule 1 on the agents.**
   - Re-derive every label change yourself (expected: the new Θ statement, and possibly a [Numerical] item).
   - Spot-check the key edits: the breather band in ch. 10 and ch. 11 must agree; so must the spin −3 sign, the domain of `prp-hirota-energy` and the Appendix E entry.
   - Run `python3 check-refs.py` in each worktree.
3. **Merge** `part3-int-ch9-10` and `part3-int-ch11` into `part3-referee-integration`. Their file sets are disjoint; resolve any conflict by hand. Apply the pending edits each agent recorded for files it could not touch.
4. **`textbook_code/README.md`:** apply the row changes from both final messages (`toda_classical.py`, the new `classical_reality.py`, the figure rows).
5. **Integration record.**
   - Merge `part3_integration_A.md` and `part3_integration_B.md` into one `## Integration` section, appended to `briefs/referee_classical_09_11_phase1_combined.md`.
   - Append a short `## Integration` section to each of the four referee reports that points there.
   - Record which errors each phase and each referee found (overlap table) for the two-phase experiment.
6. **`briefs/README.md`:**
   - add rows for the four reports, the combined report, the integration record and this file;
   - say that the Part III files are no longer frozen;
   - record the outcome of the experiment.
7. **Check and commit.** Run the chapter scripts (all PASS), `python3 check-refs.py` and `./render-scratch.sh`. Commit by named files. Do not push or open a PR until the editor asks. The PR description must list every label change made without a second referee pass.
8. **Memory.** Update `two-phase-referee-experiment.md` with the outcome. If the experiment worked, consider updating the `write-referee-brief` skill (blind phase first, then the coordinator's items).
9. **Clean up.** The main checkout holds untracked copies of Gemini's reports and the items file in `briefs/`. Delete them after the integration merges into main; they would block a pull.
