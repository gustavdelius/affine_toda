# A2. Problem 11: scalar factors, CDD factor and breather amplitudes of the soliton-preserving matrices

Scratch directory: `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a2/`. Model: Opus, effort xhigh.

**Background.** Open problem 11 (§35.7, `sec-boundary-open-problems` in `part7-half-line/35-particle-reflection.qmd`) asks for the scalar factors and particle amplitudes of the soliton-preserving (SP) reflection matrices of $a_n^{(1)}$. A worker solved part of it; the book already contains the verified part (the matrices in book conventions, the antisoliton partner with $\bar\eta=\eta+i\pi\omega(l-\frac N2)$, scalar unitarity, the $\mathbb Z_N$ obstruction, $B_1\to\bar B_1$; see ch. 33 `sec-preserving-conjugating` and `textbook_code/soliton_preserving.py`). The rest was never checked independently.

Read:
- the brief `briefs/halfline_open_11_soliton_preserving_scalar.md`;
- the report `briefs/halfline_open_11_report.md`, especially §2–§4 and its Integration table;
- the Gemini verdict `briefs/halfline_open_11_referee.md` (it reran the worker's scripts; that is not independent);
- `textbook_code/soliton_preserving.py`, `textbook_code/an_nondiagonal_scalar.py` (the template: same method for the conjugating $s=+1$ matrix), and ch. 33 `sec-preserving-conjugating`, `sec-an-nondiagonal-k`;
- for the bootstrap, the method of `briefs/halfline_open_03_scripts/sol_boot.py` (read it; write your own version on book code).

**Tasks, in order of priority.**
1. **Crossing function.** Recompute $\Phi_{1,2}$ of report §2 with the book's $\check R_{V\bar V}$, $\check R_{\bar VV}$ (`boundary_common`) for $n=1,2,3$, all $l$, two values of $\omega$. Check the claimed closed form $-(-q)^{\pm(1-l)}\sin\pi(2\mu-T)/\sin\pi(2\mu+T)$, its independence of $\eta$, and the claim that the soliton–antisoliton factor is $S^T(T-2\mu)$ in book conventions (ch. 26 / `smatrix` conventions; state the formula you use).
2. **Scalar factors.** Derive the functional equations for $a$ (on $K_V$) and $b$ (on $\lambda K_{\bar V}$) yourself, then test the worker's Gamma-product solution (report §3) with your own Stirling/Hurwitz implementation at 30 digits, for $(n,l)=(1,1),(2,0),(2,1),(3,0),(3,1),(3,2)$. The worker fixed a constant $k$ "numerically on the minimal branch": find what it is in closed form, or show the solution is ambiguous. Verify the $n=1$ statement: GZ Dirichlet $=$ minimal solution $\times(-\frac T2)_I^2$, for all $\eta$.
3. **Breather amplitude.** Compute the $B_1\to\bar B_1$ singlet projection with book code and obtain the amplitude in closed form. The worker's $K_0$ is a *fit* (six blocks, residual $5\times10^{-6}$): derive it (from the Gamma products, as `an_nondiagonal_scalar.py` check 3 does for $s=+1$) or confirm it at 30 digits for $n=1,2,3$ and, if feasible, $n=4$. Check $E_\eta$ and the $l=0,N$ statement that $\eta$ drops out.
4. **CDD factor.** Settle it for $n=2$ and $n=3$ with the SP soliton bootstrap $V^{\otimes n}\to\bar V$ (fusion of $n$ solitons into an antisoliton, then reflection, compared with the reflection of the antisoliton by $\lambda bK_{\bar V}$), measured relative to a reference so that bulk factors cancel, as problem 3 did. Decide between $\sigma_A=(-\frac T2)_I^2$, $\sigma_B$ of report §3, or neither (then find what the bootstrap requires). Also test the $n=1$ case as a control.
5. **Continued amplitudes** (only if 1–4 are done). With the CDD factor fixed, continue $\omega\to-\frac2B$ and check unitarity, crossing-unitarity and the bulk–boundary bootstrap for the particle amplitudes $K_{a\to h-a}$ with `particle_reflection.py`'s conventions and the book's S-matrix; the $B\to0$ limit.

The check script should be a new `textbook_code/sp_scalar.py` (name it so), or an extension of `soliton_preserving.py` if that is cleaner. The proposed text replaces the parts of open problem 11 that your results settle, and adds a paragraph to ch. 33 (`sec-preserving-conjugating`) or ch. 35, as you judge.
