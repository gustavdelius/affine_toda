# Agent briefs

Self-contained task briefs for follow-up work on the book. The first two were written when the book had two parts; their "Part I" is now Parts V and VII and their "Part II" is Part VI, and the files they name under `foundations/` and `smatrices/` have moved (use the id table below). Each brief can be given verbatim as the prompt of a background agent, or pasted into a fresh Claude Code session started in this repository.

| Brief | Task | Paper location |
|---|---|---|
| `exchange_statistics.md` | Decide which exchange of soliton labels the Bethe–Yang equations require, and hence which multi-soliton sectors are $\Theta$-broken | Section `sec-sector-status` (caveat on exchange statistics), open problem 6 (see below) |
| `referee_textbook.md` | Referee the textbook chapters 1–11: what was verified, and a prioritized list of unverified formulas, conventions, claims and references | `part1-background/`, `part2-algebra/`, `part3-classical/`, `textbook_code/` |
| `rigorous_heat_trace.md` | Turn Proposition `prp-net-count` part 2 and Proposition `prp-one-loop-mass` from sketches into theorems | Sections `sec-counting-rule-transmission-factors` and `sec-one-loop-masses-folded`, open problem 3 (see below) |

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
| `eq-one-loop-mass` | `part5-semiclassics/17-thimbles.qmd` | (8.1) | (8.1) |
| `eq-net-factorization` | `part5-semiclassics/18-jordan-chains.qmd` | (8.2) | (F) |
| `prp-krein-unitarity` | `part7-foundations/32-scattering.qmd` | Proposition 9.1 | Proposition 9.1 |
| `prp-krein-fundamental` | `part7-foundations/32-scattering.qmd` | Proposition 9.3 | Proposition 9.3 |
| `prp-scalar-breathers` | `part7-foundations/32-scattering.qmd` | Proposition 9.4 | Proposition 9.4 |
| `sec-crossing-sign` | `part6-soliton-smatrices/21-construction.qmd` | 14.6 | Section 2.6 of the S-matrix note |

The numbers in the third column may drift as chapters change; the ids will not.

Open problems are items of numbered lists, not labelled objects. Use the first words of the item:

- Former Part I open problems, now in `part7-foundations/33-outlook.qmd`: item 3 is the one that begins with a reference to `cnj-counting-rule` ("… in general"); item 6 is "The sector classification for self-conjugate theories".
- Former Part II open problems, now in `part6-soliton-smatrices/26-f4-solitons.qmd` (section `sec-sm-discussion`): item 4 is "Crossing convention".
