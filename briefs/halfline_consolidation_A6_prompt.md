# A6. Problem 2 leftovers: the classical boundary of the $d_n$ $\hat\epsilon=0$ amplitude, and the $d_n$ $s=\pm1$ scalar factors

Scratch directory: `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a6/`. Model: Opus, effort xhigh.

**Background.** Open problem 2 of Part VII (§35.7, `sec-boundary-open-problems` in `part7-half-line/35-particle-reflection.qmd`) asks for reflection amplitudes beyond $a_n^{(1)}$. Round 1 put the $d_n$ $\hat\epsilon=0$ scalar factor, breather amplitude and particle amplitude $K_1$ (= Castro-Alvaredo–Fring's $\hat R_1$) into the book: ch. 33 `prp-dn-scalar`, `eq-dn-scalar`, `eq-dn-breather`; ch. 35 `sec-dn-uniform`, `eq-dn-uniform-amplitude`; check script `textbook_code/dn_uniform.py`. The book says the identification of the classical boundary is open and expected by analogy.

Read:
- the brief `briefs/halfline_open_02_beyond_an.md`;
- the report `briefs/halfline_open_02_report.md`, especially §3 (classical identification: $d_4$ Tzitzéica and $d_5$ sinh-Gordon reductions, the invariant static kink as the vacuum of $C\equiv+1$, the linearized reflection) and §4–§5, and its Integration table (claims 8, 9);
- the Gemini verdict `briefs/halfline_open_02_referee.md` (it checked the reductions and the boundary condition, but the linearization numbers come from the worker's script only);
- `textbook_code/dn_uniform.py`, `textbook_code/reflection_matrices.py` (`prp-dn-r`, `prp-dn-k`: the $d_4$, $d_5$ vector $\check R$ and $K$, including the $s=\pm1$ matrices), `textbook_code/boundary_vacua.py` and `textbook_code/an_nondiagonal_scalar.py` (the template for a non-diagonal scalar factor: unitarity with the antisoliton matrix, crossing function, Stirling/Hurwitz products, breather projection, continuation, Kim's one loop);
- ch. 30 `sec-durham-classification` (the boundary potential and $C_j$ for $d_n$), ch. 31 (`sec-classical-reflection-factors`, `eq-linear-reflection`), ch. 33 (`sec-dn-vector-k`).

**Tasks, in order of priority.**
1. **Classical identification (verification).** For $d_4$ and $d_5$ with $C\equiv+1$ (all boundary coefficients at the value the book uses for the "uniform" boundary; check what that means for $d_n$, whose $\vec\phi=0$ is not a vacuum when $|C_j|=1$), recompute independently, with your own code:
   - the static boundary vacuum (the worker's invariant static kink, via the Tzitzéica / sinh-Gordon reductions, or by direct minimization of the energy with the book's boundary potential), its regularity and energy, and whether it is the minimum among the candidates;
   - the linearized reflection matrix about it (exact if possible, e.g. by the Hirota/dressing method; otherwise a carefully conditioned ODE), in the particle basis, and compare with the $B\to0$ limit of the book's $d_n$ amplitude $K_1$ (`eq-dn-uniform-amplitude`) and, if the book has them, the other $K_a$.
   - Say whether "$\hat\epsilon=0$ is $C\equiv+1$ for $d_4,d_5$" holds, at what precision, and what remains (other $a$, $n\ge6$, other sign patterns).
2. **The $d_n$ $s=\pm1$ scalar factors (new work).** For the $s=\pm1$ vector reflection matrices of `prp-dn-k`, derive the unitarity and crossing-unitarity functional equations in book conventions (antisoliton matrix, crossing function, as in ch. 33 for $a_n$), find their minimal solutions as Gamma products, compute the lowest breather amplitude and its continuation $\omega\to-\frac2B$. Test against Kim's $d_4$ one-loop Neumann amplitude (Kim, hep-th/9506031, read at the source) and say which $s$ (if any) gives the Neumann amplitude; compare with the self-dual Neumann conjecture mentioned in report §2. Do $n=4$ first, then $n=5$ if time allows.
3. **The $c_2$ step** (only if 1–2 are done). Report §4 says folding $a_3$ quantum amplitudes cannot give $c_2$ amplitudes; check that argument, and if a meaningful classical check exists (the folded $B\to0$ targets $-(1)(3)$, $-(2)^2$ at $h=4$), do it.

The check script: an extension of `dn_uniform.py` (write the extended copy in your directory) for task 1, and a new `textbook_code/dn_nondiagonal_scalar.py` for task 2 if it succeeds. Proposed text for `sec-dn-uniform`, `sec-dn-vector-k` and open problem 2.
