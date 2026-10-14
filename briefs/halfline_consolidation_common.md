# Common instructions for the consolidation agents (round 2, 2026-10-12)

Prepended verbatim to each `halfline_consolidation_A*_prompt.md` when the agent is launched.

You are verifying, independently, results that a worker agent obtained in round 1 on the open problems of Part VII (the half-line) of a Quarto textbook on affine Toda field theory. The repository is `/home/gustav/Git/affine_toda` (branch `half-line-consolidation`). The book's Part VII is in `part7-half-line/` (chapters 28–35); its check scripts are in `textbook_code/` (index: `textbook_code/README.md`).

**Rules.**
- Do not edit, create or delete any file in the repository. Do not run git commands. Do not spawn subagents.
- Work only in your own scratch directory, given below. Run book code from there with `PYTHONPATH=/home/gustav/Git/affine_toda/textbook_code`, or copy the book scripts you need into your directory.
- **Independence.** The point of this pass is a check that does not reuse the worker's code. You may read the worker's scripts to understand what was done, but your checks must be built on the book's own code (`boundary_common.py`, `particle_reflection.py`, `smatrix_simply_laced.py`, `an_nondiagonal_scalar.py`, `boundary_classical.py`, `boundary_vacua.py`, …) or written from scratch. Say explicitly, for every claim, which code your check uses.
- Distrust steps that merely look plausible. Earlier referee passes on this project found real errors every time (swapped S-matrix entries, misquoted papers, fitted constants presented as derived). If a claim is wrong, say so plainly; a refutation is as useful as a confirmation.
- When you quote a paper, read it at the source (arXiv TeX or PDF; download into your scratch directory) and give equation numbers. Never build a correction to a published paper on a paraphrase.
- Usage limits may stop you; you will then be resumed. Run long computations as background jobs writing to log files in your directory, and on resumption check those logs before recomputing. Python output redirected to a file is buffered: use `python3 -u`. Never use `pkill -f` or `pgrep -f` (they match your own shell); kill by PID.

**Technical pitfalls known from round 1.**
- `smatrix_simply_laced.py`, `boundary_coideal.py` and `reflection_matrices.py` run all their checks on import, which is slow. To get only the definitions, `exec` the source up to the first check (e.g. split at `'\nB = 0.37'`).
- `boundary_common.rcheck` is the generic intertwiner solver and is far too slow for $d_n$ tensor products.
- `boundary_coideal.k_closed` takes the wrong branch for odd $N$ when $\cos\frac{\pi\omega}2<0$. Use `kmat` from `an_nondiagonal_scalar.py` (analytic branch).
- Avoid the special couplings $B=\frac1{2k}$ in residue checks (level crossings). Normalize eigenvectors and intertwiners by a fixed entry, not by `max|·|`.
- Infinite Gamma products: sum them with Stirling/Hurwitz tails as in `an_nondiagonal_scalar.logprod`, at 30 digits, rather than by truncation.

**Book notation.** $h=n+1$, $T=\frac{h\omega}2$, $\mu=-\frac{iT\theta}\pi$, blocks $(x)_I=\sin\frac{\pi(\mu+x)}{2T}/\sin\frac{\pi(\mu-x)}{2T}$, real-coupling blocks $(x)$ as in ch. 35, $\omega\to-\frac2B$ for the continuation. The book's $\Phi$ in ch. 33 is taken; the crossing function is $\Psi_s$. Write the fundamental range as $\frac\pi2(\lambda+1)$, not $\eta_N$. Status tags (**[Theorem]**, **[Numerical]**, **[Conjecture]**) are bold labels at the start of a statement; [Theorem] needs a proof, not numerics. No status tag inside prose that contains a cross-reference.

**Deliverable.** Your final message is your full report (it is saved verbatim; files you write are not read unless listed). It contains:
1. For every claim in your task: verdict **VERIFIED / CORRECTED / REFUTED / UNCHECKED**, the evidence (precision, parameters, which code), and what remains open.
2. A **check script for `textbook_code/`**, written in your scratch directory, in the style of `textbook_code/an_nondiagonal_scalar.py`: a docstring naming the book ids it supports, the conventions and each numbered check; one `PASS`/`FAIL` line per check via `boundary_common.report`; imports only from `textbook_code/` and standard packages; runtime a few minutes at most. Run it and include its output. (If your task says to extend an existing book script, write the extended copy in your scratch directory instead.)
3. **Proposed book text** (Quarto markdown, book notation, cross-reference ids from the chapters), with the place where it goes. Keep it short and in the register of the surrounding chapter.
4. Any statement in the book or in Appendix D (`appendices/D-corrections.qmd`) that your results contradict.
5. The full paths of all your scripts and logs.
