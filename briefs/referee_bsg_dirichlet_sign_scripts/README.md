# Evidence for the boundary sine-Gordon Dirichlet-sign report

These are the independent referee computations supporting [the report](../referee_bsg_dirichlet_sign_report.md), originally run against commit `ed0585ae78572d7ada8986337f89ec685c0b6462`. They are archived referee evidence, not maintained book check scripts. No book module is imported.

## Reproduce the independent computations

Requires Python 3, NumPy and SciPy. From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 briefs/referee_bsg_dirichlet_sign_scripts/independent_checks.py
PYTHONDONTWRITEBYTECODE=1 python3 briefs/referee_bsg_dirichlet_sign_scripts/algebra_checks.py
PYTHONDONTWRITEBYTECODE=1 python3 briefs/referee_bsg_dirichlet_sign_scripts/analytic_flux_check.py
```

| Script | Evidence | Recorded output |
|---|---|---|
| `independent_checks.py` | Exact three-soliton solution, PDE and boundary residuals, nine flux-integral cases, Gamma-product phase derivative and cutoff dependence | `independent_checks.out` |
| `algebra_checks.py` | Static energies, pole scaling, scalar-factor evenness at finite cutoff, node Taylor coefficients, fusion/continuation identities, charge conjugation and strip winding | `algebra_checks.out` |
| `analytic_flux_check.py` | Closed boundary derivative and large-lambda convergence; imports the sibling `independent_checks.py` | `analytic_flux_check.out` |

These scripts print diagnostic values; they are not a comprehensive automated assertion suite. In particular, finite Gamma products have a cutoff error and the common phase agrees with the classical answer only in the large-lambda limit. Interpret the results using the derivations and limits in the report.

## Secondary checks of the book's scripts

The six `book_*.out` files record executions of the book scripts at **ed0585a**, not the later merged head. They contain 132 PASS lines and no FAIL lines altogether. `refs_complete.out` records the reference check against the original pinned snapshot, after extracting the configuration and code README files needed by the checker. These are secondary comparisons, not independent evidence.

## Sources and provenance

The report links to the source papers and identifies TeX source locations. `source-archives.sha256` records the SHA-256 hash and download URL of every arXiv source archive actually read. The downloaded articles, extracted TeX, full repository snapshot and broad search outputs are not included in this PR. Their original local location is:

```text
.claude/worktrees/eager-bun-c1cfa6/referee_scratch/bsg_dirichlet_sign/
```

The computations above do not need that directory or the downloaded sources. For source verification, download the listed archives into a fresh referee scratch directory and compare the hashes; the unversioned arXiv endpoints may eventually serve a revised archive. The report's historical source line numbers refer to the downloaded versions.

## Integration

Pending Claude's independent verification under Rule 1. Follow the report's post-PR-#63 checklist and record decisions in an Integration section of the report. Do not infer a proof of the finite-coupling UV–IR relation or a priority claim from these numerical checks.
