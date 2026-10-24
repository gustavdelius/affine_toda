# Codex phase-2 evidence

Report: `briefs/referee_classical_09_11_codex_report.md`.
Reviewed worktree commit: `bcd82b9ee9667975eed02109f08a257d50ba38b1`.

These are independent referee scripts, not book scripts. Run from the reviewed repository root with Python 3, SymPy, NumPy and mpmath installed. All output names below should remain Codex-specific. The scripts print their results to stdout, except the current constructor, which also writes the JSON file described below.

- `codex_phase2_currents_v3.py`: exact polynomial conservation-law solve for A2 spins 2 and 4. Writes `referee_scratch/classical_09_11/codex/codex_phase2_currents_v3.json`; create that directory first if scratch has been cleaned. The JSON is also archived here.
- `codex_phase2_integrals_v2.py`: direct canonical stress/current integrations at 45 digits on [-90,90]. Reads the above JSON from the scratch path. If not rerunning the constructor, restore the archived JSON to that exact path before running. All reported comparisons pass; spin-2 breather residuals are below 8e-34.
- `codex_phase2_integrals.py`: original 32-digit, [-40,40] run, retained with its output to document two spin-2 tail-truncation FAILs. Use v2 for the converged results.
- `codex_phase2_centralizers.py`: graded matrix-loop centralizers, including explicit triality construction. Nullspaces are numerical, with singular-value gaps reported.
- `codex_phase2_couplings.py`: A4/E6 potential cubic tensors, charge-basis area/orbit tests, and real-basis counterexamples.
- `codex_phase2_extra.py`: symbolic twisted A2 Lax pair, folded G2 tau polynomials, embedded SG root calculations, and a complex-soliton singularity.
- `codex_phase2_phase.py`: leading Hollowood scalar-phase integral and fused A2/A3 comparisons. PASS means the independent source-side integral and product identities agree. It does **not** mean the book/source phase convention has been reconciled; the report marks that comparison unclear.

The `*_output.txt` files record the executions used in the report. None of these checks imports an author's helper. Source PDFs, chapter snapshots and preliminary unsuccessful scripts remain under the Codex scratch directory. The phase-1 evidence already in this directory is unchanged.
