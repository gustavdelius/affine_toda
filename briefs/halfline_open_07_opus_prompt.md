# Opus referee prompt: half-line open problem 07 (tier-B follow-up)

Written 2026-10-12 by the Claude coordinator, after reading the Gemini verdict `briefs/halfline_open_07_referee.md`.

**What the Gemini pass established.** It independently established the following, so do not redo them:
- $\check R^\dagger(C\otimes C)\check R=C\otimes C$ with its own code (48 cases);
- the coideal $\tau$-stability computation;
- the reality of the full problem-3 scalar factors;
- the UN VV threshold bracket (unbroken at $\omega=0.60$, broken at $0.80$).

**Gemini's own Opus items** are research problems and go to the book's open-problem list: an analytic explanation of the thresholds, TCSA, and part (c).

Read:
- the brief `briefs/halfline_open_07_nonunitary_half_line.md`;
- the report `briefs/halfline_open_07_report.md`;
- the Gemini verdict;
- the scripts in `briefs/halfline_open_07_scripts/` (copy `textbook_code/boundary_common.py` alongside).

## Items

1. **Claim 3, the left-wall convention against the book.**
   - The worker says the left wall must be $K^L_{\bar a}(-\theta)=K^R_a(\theta)$, the image under $e\leftrightarrow f$, $k\to k^{-1}$, and that exchanging $K$ and $K^{\rm T}$ on the left breaks commutativity.
   - Find what the book's double-row transfer matrix and boundary Bethe–Yang equations use (`sec-boundary-bethe-yang`, `sec-boundary-finite-size` in ch. 28, and any use in ch. 35).
   - Say whether the book is consistent with the worker's convention, silent, or wrong, and draft the correction if needed.
2. **Are the $\zeta\ne1$ sectors physical?** This decides how strongly "a single soliton between a uniform and a Neumann wall has complex energies" may be stated.
   - Analyse topological charge on the interval: which vacua each wall allows (uniform: $\vec\phi=0$; Neumann: the degenerate family of `sec-neumann-duality`, at imaginary coupling the lattice of vacua), which soliton charges can interpolate between them, and how $\mathcal Z$ acts on those states.
   - Decide whether $\zeta\ne1$ sectors are excluded, are Fourier modes of degenerate boundary vacua, or are physical.
   - Give book text that states the conclusion at the right strength. Gemini's proposal says "formal Bethe–Yang diagonalization".
3. **Consistency of $C$ with `prp-krein-fundamental`.** Check that the book's `prp-krein-fundamental` (`part8-foundations/40-scattering.qmd`) uses the same $C$ (basis reversal) in the same gradation and coproduct, so that claim 6 ("untwisted $C$ works for $V\bar V$") is a statement about the book's $C$. If the book's $C$ differs, translate.

For each item give the verdict (CONFIRMED / CORRECTED / REFUTED / UNCHECKED), the evidence with file and line references, and book text. Give script paths.
