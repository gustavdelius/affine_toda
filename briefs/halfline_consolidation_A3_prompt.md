# A3. Problem 11: the two classical claims (coweights outside $2\Lambda_W^\vee$; real-coupling Dirichlet)

Scratch directory: `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a3/`. Model: Opus, effort max.

**Background.** While solving open problem 11, a worker made two classical claims by hand (report `briefs/halfline_open_11_report.md`, §1 "Classical side"). Neither has been checked. Claim 4 contradicts the book and an entry in Appendix D, so it needs special care.

Read:
- the report §1 and its Integration table (claims 4, 5);
- the coordinator's follow-up prompt `briefs/halfline_open_11_opus_prompt.md`, items 1 and 2 (written in round 1 but never run). **Those two items are your task.** Item 3 has been settled; skip it.
- the Gemini verdict `briefs/halfline_open_11_referee.md`, whose treatment of claim 5 drops the constants $e^{\beta\alpha_j\cdot\phi_0/2}$;
- ch. 30 `part7-half-line/30-boundary-conditions.qmd`, in particular the compatibility condition (`eq-boundary-compatibility`), `sec-durham-classification` and `sec-soliton-preserving` (with `prp-soliton-preserving`);
- `textbook_code/boundary_classical.py` (it checks the soliton-preserving condition);
- Appendix D, entry "Delius 1998b" in `appendices/D-corrections.qmd`, item (1.4): it says $\vec\lambda$ must lie in $2\Lambda_W^\vee$ "because for any other coweight $\operatorname{Re}\vec\phi(0)=-\frac\pi\beta\vec\lambda$ is not a vacuum". Read Delius's paper (hep-th/9804030 or the arXiv number in `references.bib`, key `delius1998b`) at (1.4) and around it.

**Claim 4** (report §1): with $\sigma=\dagger$ and constant diagonal $K$, compatibility gives $\operatorname{Re}\vec\phi(0)\in\frac\pi\beta\Lambda_W^\vee$ and $K_+=e^{i\pi\mu\cdot H}$, so coweights outside $2\Lambda_W^\vee$ are allowed, contrary to the book. Decide:
- Is the compatibility condition satisfied for $\vec\mu\notin2\Lambda_W^\vee$ with $K_+=e^{i\pi\mu\cdot H}$? Check by hand from the book's Lax pair and numerically with `boundary_classical.py`'s machinery.
- What exactly does the book claim (integrability, or that the boundary value is a vacuum, or that the ground state is the trivial one)? Separate "classically integrable" from "the boundary value is a vacuum" and from "the vacuum next to the boundary is $\vec\phi=$ const". The sine-Gordon Dirichlet boundary at $\phi_0=\frac\pi\beta$ is integrable although $\phi_0$ is a maximum of the potential; is the $a_n$ case the analogue?
- Is the Appendix D entry right, wrong, or right for a different reason? Draft the corrected book sentence and, if needed, the corrected Appendix D item.

**Claim 5**: real-coupling Dirichlet $\vec\phi(0)=\vec\phi_0$ is not classically integrable for $n\ge2$ (field-independent $K$, linear $\sigma$, fundamental representation). Re-derive carefully, keeping the constants $e^{\beta\alpha_j\cdot\phi_0/2}$ and allowing $\sigma$ to permute nodes (diagram automorphisms) and to be antilinear. Say whether the restrictions leave a loophole (field-dependent $K$, other representations, a nonlocal or dynamical boundary). Check prior art: Bowcock–Corrigan–Dorey–Rietdijk hep-th/9501098 (they classify integrable boundaries of real-coupling Toda; do they find Dirichlet conditions, and for which algebras?), Corrigan's reviews (e.g. hep-th/9612138), Corrigan–Dorey–Rietdijk–Sasaki hep-th/9404108, and anything on Dirichlet conditions in $a_n$ Toda you find. Is the statement known, or contradicted? Note the connection to open problem 7(c) ("Dirichlet walls") and say what your result implies for it.

The check script: a new `textbook_code/boundary_dirichlet.py` (or extended copy of `boundary_classical.py`, if that is cleaner) that tests both claims numerically: the compatibility residual for $\vec\mu\in\Lambda_W^\vee\setminus2\Lambda_W^\vee$ with and without $K_+=e^{i\pi\mu\cdot H}$, and the non-existence of a constant $K$ for real Dirichlet at generic $\vec\phi_0$ for $n=2,3$ (with the sinh-Gordon case $n=1$ as the control, where it must exist).
