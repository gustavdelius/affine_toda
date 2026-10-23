# Brief: coordinate the referee passes for the half-line open problems (Gemini, Antigravity)

Written 2026-10-11. You are the referee coordinator. You run in Antigravity on Gemini. A separate coordinator, Claude Code, runs the worker agents that attack the open problems of Part VII (`briefs/halfline_open_problems.md`). Whenever a worker report arrives, Claude saves it and writes a referee prompt file. You watch for those prompt files, run an independent adversarial referee pass on each one, and save the verdicts. The two coordinators communicate only through files in `briefs/`.

The repository is `/home/gustav/Git/affine_toda`, a Quarto textbook on affine Toda field theory with Python scripts in `code/`. The half-line material is Part VII, `part7-half-line/` (chapters 28–35). The open problems are §35.7, `sec-boundary-open-problems` in `part7-half-line/35-particle-reflection.qmd`.

## Why you, and what you are for

The workers are Claude Opus. Earlier referee passes on this book were Claude agents too, and every one of them found real errors that the author pass had missed: a doubled logarithm in a lemma, an unproved orbit case, overclaimed status labels, and prior art that went uncredited. You are a different model family, so your errors should not correlate with the author's. That independence is the point. Do not defer to the report's confidence, its tone, or its claims to have checked something. Recheck every step yourself, especially the ones that look plausible.

## The file protocol

All paths are relative to the repository root. `<tag>` is the problem number, zero-padded, with an optional suffix, for example `03`, `12a`, `12b` or `06_01`.

| File | Written by | Meaning |
|---|---|---|
| `briefs/halfline_open_NN_<name>.md` | (exists) | the original brief for problem NN |
| `briefs/halfline_open_<tag>_report.md` | Claude | the worker's report, verbatim |
| `briefs/halfline_open_<tag>_scripts/` | Claude | the worker's scripts, if any |
| `briefs/halfline_open_<tag>_referee_prompt.md` | Claude | **your trigger**: what to referee, and the tier |
| `briefs/halfline_open_<tag>_referee.md` | **you** | your verdicts; its existence means "done" |
| `briefs/halfline_open_<tag>_referee_opus.md` | Claude | a later Opus pass on tier-B items |
| `briefs/halfline_referees_done` | Claude | sentinel: no more prompts will come; stop |

A prompt file has this shape:

```
# Referee prompt: half-line open problem <tag>
Tier: A | B
Brief: briefs/halfline_open_NN_<name>.md
Report: briefs/halfline_open_<tag>_report.md
Scripts: briefs/halfline_open_<tag>_scripts/   (or "none")
## Claims to referee
1. [Theorem] ...   2. [Numerical] ...   ...
## Coordinator notes
(anything specific: suspected weak points, scope, things already known)
```

**Tier A**: your verdict is the referee pass. **Tier B**: the report claims a [Theorem] or a correction to published literature. You referee everything first. Claude then runs an Opus pass only on the tier-B claims you leave standing, so your list of "items for an Opus pass" decides what Claude spends its quota on.

## Your loop

1. **Find pending prompts.** A prompt is pending when its `_referee.md` does not exist and it is not claimed. Keep claim markers outside the repository, in `/tmp/halfline-referee/claimed/`:

   ```bash
   cd /home/gustav/Git/affine_toda
   mkdir -p /tmp/halfline-referee/claimed
   pending() {
     for p in briefs/halfline_open_*_referee_prompt.md; do
       [ -e "$p" ] || continue
       t=${p#briefs/halfline_open_}; t=${t%_referee_prompt.md}
       [ -e "briefs/halfline_open_${t}_referee.md" ] || [ -e "/tmp/halfline-referee/claimed/$t" ] || echo "$t"
     done
   }
   # wait up to 10 minutes for work, then return (rerun if empty)
   for i in $(seq 60); do out=$(pending); [ -n "$out" ] && { echo "$out"; break; }; [ -e briefs/halfline_referees_done ] && { echo DONE; break; }; sleep 10; done
   ```

   If your shell tool has a shorter timeout, shorten the loop and call it repeatedly. Do not use `pgrep -f` or `pkill -f` in wait loops: they match their own command line.

2. **Claim it.** `touch /tmp/halfline-referee/claimed/<tag>`.

3. **Referee it.** If Antigravity lets you start separate agents, start one per prompt with a fresh context. Give it the prompt file, the "Standard referee instructions" below verbatim, and the paths of the brief, the report and the scripts. Run at most three referees at once. Otherwise referee the prompts yourself, one at a time, and start each in a fresh context if you can. Use the strongest Gemini setting you have (High) for every referee.

4. **Check and save the result.** Make sure the verdict file has every section of the template below and a verdict for every numbered claim. Write it first to `/tmp/halfline-referee/<tag>/referee.md`, then copy it to `briefs/halfline_open_<tag>_referee.md` in one step. The file appears only when it is complete, and Claude treats its appearance as completion.

5. **Repeat** until `briefs/halfline_referees_done` exists or the user tells you to stop. Before you stop, finish or explicitly abandon every claimed prompt. For an abandoned one, write a verdict file that says so and lists what was done, and remove its claim marker.

## What you must not do

- Do not edit anything in the repository except the `briefs/halfline_open_<tag>_referee.md` files you write. In particular, leave the book chapters (`*.qmd`), `references.bib`, `code/`, the reports, the prompts and the briefs untouched.
- No git operations at all: no commit, branch, stash or checkout.
- Do not contact, stop or modify the Claude agents or their files.
- Work in `/tmp/halfline-referee/<tag>/`. Copy scripts there before you run them.

## Standard referee instructions

*(Give this section verbatim to every referee agent.)*

You are an independent adversarial referee for a research report on affine Toda field theory on the half-line. The author is another AI model. Treat the report as a manuscript submitted by a capable but fallible colleague who has an incentive to declare success. Your job is to find what is wrong, unproved, overclaimed, or already published. A pass that finds nothing is acceptable only if you actually redid the work.

**Setup.** Work in `/tmp/halfline-referee/<tag>/`. Copy the worker's scripts from `briefs/halfline_open_<tag>_scripts/` and any `code/` modules they import into that directory, and run them from there with `python3`, which has numpy, scipy, mpmath and sympy. If a package is missing, make a venv in your scratch directory; do not install system-wide. Do not edit the repository.

**Read first.** The prompt file, the original brief, the report, and the book sections the brief names. Find them by id with `grep -rn "#sec-..." part*/`. The book's conventions are the reference: where the report and the book disagree, decide which is right and say so.

**For every numbered claim:**

1. **Re-derive it.** Redo analytic steps by hand or with sympy or mpmath. Do not reread the report's derivation and nod along. For a closed form, test it at a few random complex parameter values to at least 30 digits with mpmath, and at the special points the brief names, such as $n=1$, $n=2$ and $B\to0$.
2. **Recompute it.** Rerun the worker's scripts. Then compute at least one key number per [Numerical] claim with code you write yourself, without importing the worker's functions. A rerun alone only shows that the script is deterministic.
3. **Check the inputs.** When the report quotes an equation from a paper ("DG (3.11)", "AKR eq. (10)", "Gandenberger (5.16)"), open the paper on arXiv and check that the equation says what the report claims, in the stated conventions. Several earlier errors in this project came from misquoted equations and prefactors.
4. **Check prior art.** Search arXiv and INSPIRE for the result. Include Delius, Gandenberger, Delius–Gould–Zhang, Bowcock–Perkins, Corrigan and collaborators, Bajnok–Palla–Takács, Fateev, and Regelskis–Vlaar. If it is already published, say where.
5. **Check the status label.** The book uses [Theorem] (complete proof), [Numerical] (computed, with stated precision and parameter range), [Conjecture], and similar labels. Is the proposed label earned? A [Theorem] proved only for $n\le5$ is a [Numerical] statement for general $n$.

**Pitfalls that recur in this project:**

- The sign of $q$, and the coproduct convention ($q$ versus $q^{-1}$, $\pm q$). The book's physical value is $q=-e^{-i\pi\omega}$; at $-q$ some K-matrices fail the reflection equation.
- $R$ versus $P\check R$ as the operator in the Bethe–Yang or crossing equations.
- Scalar-factor prefactors and signs, such as $f(\theta)=-(1-i\sinh\theta)/(1+i\sinh\theta)$, and CDD ambiguities. One earlier "correction" of a published formula came from a misquoted prefactor.
- Boundary-parameter dictionaries: Bowcock's $A_i=-C_i$; Delius 1998's $C\equiv-\epsilon$; in sine-Gordon the node crossing $\hat\epsilon_1\leftrightarrow C_0$, $\hat\epsilon_0\leftrightarrow C_1$; GZ's boundary $\xi$ versus the book's $\xi=\pi/\lambda$; the boundary-energy normalization $\mathcal B(0)=0$.
- Coupling normalizations: imaginary versus real coupling ($\omega\to-2/B$); $b^2$ versus $\beta^2$; twisted algebras with long roots of length$^2=2k$.
- Caches keyed on too few parameters, which silently reuse stale results. Check this whenever a parameter changes.
- Python output redirected to a file is buffered, so an empty log does not mean no progress.

**Verdict for each claim**, one of:
- **CONFIRMED**: you re-derived or recomputed it independently.
- **CONFIRMED (numerics only)**: your independent numbers agree, but you could not check the analytic argument.
- **CORRECTED**: true after a fix; give the fix.
- **REFUTED**: give the counterexample or the broken step.
- **UNCHECKED**: say why, for example a paper you could not access.

**Write your report** in exactly this form:

```
# Referee report: half-line open problem <tag>
Date: YYYY-MM-DD. Referee: Gemini <version> (<setting>), via Antigravity. Tier: A|B.

## Summary
(3–6 sentences: does the report's main conclusion stand?)

## Verdicts
| # | Claim (short) | Label claimed | Verdict | Label you recommend |
|---|---|---|---|---|

## Details
### Claim 1
What you did, what you found, numbers with precision, scripts used.
(one subsection per claim)

## Errors and overclaims not tied to a single claim

## Prior art found

## Proposed book text
(Comment on the report's proposed replacement text: wrong statements, wrong labels, missing citations.)

## Items for an Opus pass
(Tier B only: every [Theorem] or literature-correction claim that you left standing, and every one you could not check, each with the specific step a second referee should attack. Write "none" for tier A.)

## Scripts
(Each script you wrote, with a one-line description and where it is in /tmp/halfline-referee/<tag>/.)

Referee pass complete.
```

## Final message to the user

When you stop, give one table with the tag, the tier, the summary verdict (stands / stands with corrections / fails / incomplete), the number of claims in each verdict category, and the path of the verdict file. Then list the tier-B items you flagged for an Opus pass.
