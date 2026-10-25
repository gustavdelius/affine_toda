# Referee scripts for the half-line lecture notes

The Codex referee's own checks for `briefs/referee_lectures_halfline_phase1_report.md` and `briefs/referee_lectures_halfline_report.md`, copied unchanged from `referee_scratch/lectures_halfline/`. They were written in a throwaway copy of the repository and expect to be run from the repository root.

- `check*.py` and `verification.log`: phase 1, numbered by risk-list entry.
- `phase2_checks.py`, `phase2_additional.py` and their logs: phase 2. `phase2_fetch.py` downloaded the arXiv sources.
- `source_checksums.txt`: SHA-256 of the arXiv TeX sources the referee read (`sources/<arXiv id>/`) and of its PNG renders of the six lecture figures (`figures/`). Those files are not archived here.
