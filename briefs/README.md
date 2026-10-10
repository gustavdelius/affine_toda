# Agent briefs

Self-contained task briefs for follow-up work on the book. The first two were written when the book had two parts; their "Part I" is now Parts V and VII and their "Part II" is Part VI, and the files they name under `foundations/` and `smatrices/` have moved (use the id table below). Each brief can be given verbatim as the prompt of a background agent, or pasted into a fresh Claude Code session started in this repository.

| Brief | Task | Paper location |
|---|---|---|
| `exchange_statistics.md` | Decide which exchange of soliton labels the Bethe–Yang equations require, and hence which multi-soliton sectors are $\Theta$-broken. **Answered 2026-10-10:** a graded exchange is a Bloch twist, and the global form selects the twists; integrated as Propositions `prp-exchange-twist` and `prp-bloch-sectors` (scripts in `foundations_code/sectors/bloch/`). Still needs a referee pass. | Section `sec-sector-status` (paragraph "Exchange statistics and Bloch sectors"), open problem 6 (see below) |
| `referee_exchange_statistics.md` | Adversarial referee pass on the answer to `exchange_statistics.md`: the proofs of `prp-exchange-twist` and `prp-bloch-sectors`, the scope claim "in every global form", the numerics in `foundations_code/sectors/bloch/`, and prior art. **Run 2026-10-10** (items 1–5, partly 6); its findings are integrated in PR #2 | Section `sec-sector-status`, open problem 6, `sec-crossing-sign` |
| `referee_exchange_statistics_report.md` | Report of the referee pass on PR #2 (2026-10-10), with the integration notes and the open items; scripts in `referee_exchange_statistics_scripts/`. Not a brief | `part7-foundations/32-scattering.qmd`, `foundations_code/sectors/bloch/` |
| `referee_textbook.md` | Referee the textbook chapters 1–16 and 20, and the textbook framing of 17–27: what was verified, and a prioritized list of unverified formulas, conventions, claims and references | `part1-background/`, `part2-algebra/`, `part3-classical/`, `part4-real-coupling/`, `part5-semiclassics/`, `part6-soliton-smatrices/`, `textbook_code/` |
| `referee_textbook_report.md` | Report of the referee pass of 2026-10-10 on `referee_textbook.md`: a verdict per item, with replacement text for each error; scripts in `referee_textbook_scripts/`. Not a brief: it records findings to be integrated. Item 2 (the sign of $q$, crossing sign) was integrated on 2026-10-10, with the dependent results recomputed at the physical $q$; the other items on branch `textbook-referee-fixes`, except those listed at the top of the report | all parts, `references.bib` |
| `rigorous_heat_trace.md` | Turn Proposition `prp-net-count` part 2 and Proposition `prp-one-loop-mass` from sketches into theorems. Done 2026-10-10: proofs integrated (PR #3) and refereed (below) | Sections `sec-counting-rule-transmission-factors` and `sec-one-loop-masses-folded`, open problem 3 (see below) |
| `referee_heat_trace.md` | Adversarial referee of those proofs (`eq-one-loop-mellin`, `prp-net-count` parts 2–3, `lem-born`, `prp-one-loop-mass`) before their [Theorem] labels are kept. Done 2026-10-10: no fatal errors; all [Theorem] labels kept after fixes made in place (proof details, applicability and status statements, the Mellin-form gap stated in Chapter 17, prior art credited), plus the independent $\phi^4$ check `phi4_check.py` | Chapters 17–19, `foundations_code/oneloop/heat_trace_checks.py` |

Recommended protocol, which worked in the session that wrote these briefs:

1. Run one worker per brief; at most two or three agents at a time, because long parallel runs hit usage limits.
2. The worker reports in its final message; it does not edit the repository.
3. Spot-check its key numbers by re-running its scripts, then integrate the results into the paper.
4. Send the integrated text to a fresh adversarial referee agent before keeping any [Theorem] label.

## Locating sections and results

The papers are chapters of a Quarto book, and Quarto renumbers sections and theorems automatically. The briefs therefore cite **ids**, which do not change. Find an id with, for example, `grep -rn '{#prp-net-count}' part*/`. Id prefixes: `sec-` section, `thm-`/`prp-`/`lem-`/`cnj-` theorem, proposition, lemma, conjecture, `eq-` equation.

| Id | File (current) | Number when the briefs were written | Former number |
|---|---|---|---|
| `sec-lagrangians` | `part3-classical/09-lagrangian.qmd` | 2.2 | 2.2 |
| `sec-counting-rule` | `part5-semiclassics/17-thimbles.qmd` | 8.2 | 8.2 |
| `sec-counting-rule-transmission-factors` | `part5-semiclassics/18-jordan-chains.qmd` | 8.5 | 8.5 |
| `sec-hirota-form-dn1-en1` | `part5-semiclassics/18-jordan-chains.qmd` | 8.8 | 8.8 |
| `sec-one-loop-masses-folded` | `part5-semiclassics/19-one-loop-masses.qmd` | 8.9 | 8.9 |
| `sec-krein-unitarity` | `part7-foundations/32-scattering.qmd` | 9.2 | 9.2 |
| `sec-sector-status` | `part7-foundations/32-scattering.qmd` | 9.3 | 9.3 |
| `prp-net-count` | `part5-semiclassics/18-jordan-chains.qmd` | Proposition 8.3 | Proposition 8.5 |
| `lem-wronskian` | `part5-semiclassics/18-jordan-chains.qmd` | Lemma 8.1 | Lemma 8.4 |
| `thm-hirota-factorization` | `part5-semiclassics/18-jordan-chains.qmd` | Theorem 8.2 | Theorem G |
| `prp-cn-threshold-blocks` | `part5-semiclassics/18-jordan-chains.qmd` | Proposition 8.5 | Proposition 8.7 |
| `prp-hirota-form` | `part5-semiclassics/18-jordan-chains.qmd` | Proposition 8.6 | Proposition 8.8 |
| `prp-one-loop-mass` | `part5-semiclassics/19-one-loop-masses.qmd` | Proposition 8.7 | Proposition 8.9 |
| `lem-born` | `part5-semiclassics/19-one-loop-masses.qmd` | new (2026-10-10) | — |
| `eq-one-loop-mass` | `part5-semiclassics/17-thimbles.qmd` | (8.1) | (8.1) |
| `eq-net-factorization` | `part5-semiclassics/18-jordan-chains.qmd` | (8.2) | (F) |
| `eq-one-loop-mellin` | `part5-semiclassics/17-thimbles.qmd` | new (2026-10-10) | — |
| `prp-krein-unitarity` | `part7-foundations/32-scattering.qmd` | Proposition 9.1 | Proposition 9.1 |
| `prp-krein-fundamental` | `part7-foundations/32-scattering.qmd` | Proposition 9.3 | Proposition 9.3 |
| `prp-scalar-breathers` | `part7-foundations/32-scattering.qmd` | Proposition 9.4 | Proposition 9.4 |
| `prp-exchange-twist` | `part7-foundations/32-scattering.qmd` | (new, 2026-10-10) | — |
| `prp-bloch-sectors` | `part7-foundations/32-scattering.qmd` | (new, 2026-10-10) | — |
| `sec-crossing-sign` | `part6-soliton-smatrices/21-construction.qmd` | 14.6 | Section 2.6 of the S-matrix note |

The numbers in the third column may drift as chapters change; the ids will not.

Open problems are items of numbered lists, not labelled objects. Use the first words of the item:

- Former Part I open problems, now in `part7-foundations/33-outlook.qmd`: item 3 is the one that begins with a reference to `cnj-counting-rule` ("… in general"); item 6 is "The sector classification for self-conjugate theories" (its former sub-item "Exchange statistics" is now "Bloch sectors").
- Former Part II open problems, now in `part6-soliton-smatrices/26-f4-solitons.qmd` (section `sec-sm-discussion`): item 4 is "Crossing convention".
