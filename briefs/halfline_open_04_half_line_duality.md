# Brief: weak–strong duality on the half-line — solitonic boundaries and coideals (half-line open problem 4)

You are a mathematical physicist with expertise in boundary conformal field theory, Toda theories, the boundary TBA and quantum groups. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 4 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Duality on the half-line.** Duality exchanges $C\equiv+1$ and Neumann, and maps the $s=-1$ amplitude to Corrigan et al.'s $C\equiv-1$ conjecture. What does it do to the solitonic boundaries, whose amplitudes (`eq-solitonic-factor`) have duals that satisfy every axiom? Is there a statement at the level of coideals?

## What the book already has

`sec-neumann-duality` (ch. 35): bulk S-matrices are self-dual block by block; reflection amplitudes are not. Duality maps $\omega\to\omega'=-\omega/(\omega+1)$, and how the coideal parameters transform is not known. Since $\hat\epsilon=0$ and $\hat\epsilon=\epsilon_*$ are exchanged, the map cannot be a symmetry of the coideals at fixed $\omega$. Also `sec-solitonic-boundaries` (the amplitude $K^{(b)}_a=K_aM^{(b)}_a$ and its quantum equivalence classes), `exr-kim-expansion`, `sec-strong-weak`, and the sinh-Gordon case in `sec-shg-boundary-breathers` ("Self-duality"). Script: `textbook_code/particle_reflection.py`.

## Literature

- Corrigan (`corrigan1998`) and Corrigan–Delius (`corrigan_delius1999`): sinh-Gordon reflection factors are not self-dual; duality relates different boundaries.
- **Ahn–Kim–Rim (`akr2002`, hep-th/0110218).** The boundary action is $\mu_B\sum\eta_ie^{be_i\cdot\phi/2}$ with $\mu_B^2=\frac\mu2\cot\frac{\pi b^2}2$ (their (10)). The Neumann amplitude is conjectured to be $R^{(f)}_j(b)=R^{(+)}_j(1/b)$ (their (14)). In the FZZ parametrization of the boundary cosmological constant, the two boundaries sit at $s^{(+)}=ib/2$ and $s^{(f)}=i/(2b)$, related by $b\to1/b$ (machine summary; check). Their TBA confirms both for $a_n$. The mixed (+f) boundaries are self-dual.
- **Fateev–Onofri, hep-th/0207152 (not in `references.bib`):** for non-simply-laced ATFTs, quite different classical integrable boundary theories $G$ and $G^\vee$ are one quantum theory related by $b\to1/b$. Boundary parameters $w_{0,r}=s_{0,r}b/(1+b^2)$ transform so that the dual classical limits coincide, and a "dual to Neumann" boundary appears (their (24)). This is the most explicit statement of how boundary parameters transform under duality.
- Fateev–Zamolodchikov–Zamolodchikov, boundary Liouville I, hep-th/0001012 (not in `references.bib`): the boundary cosmological constant $\mu_B$ and its dual $\tilde\mu_B$, the relation $\cosh^2\pi bs\propto\mu_B^2\sin\pi b^2/\mu$, and its dual under $b\to1/b$. This is the CFT origin of boundary duality.
- Fateev, hep-th/0103014: boundary one-point functions and duality in deformed Toda theories.
- Faddeev's modular double and Bytsko–Teschner: the bulk quantum group with dual parameter $\tilde q=e^{i\pi/b^2}$. These are the natural algebraic setting for a coideal-level statement (search for coideals of the modular double; we found none in a quick search).

## Suggested plan of attack

1. **The CFT mechanism.** Write the $a_n$ boundary action in the UV as a boundary Toda CFT with both the boundary exponentials $e^{be_i\cdot\phi/2}$ and their duals $e^{e_i\cdot\phi/(2b)}$. Following FZZ and AKR, identify $(+)$ and Neumann as points of a single boundary-parameter space exchanged by $b\to1/b$. Check that AKR's (10) and the FZZ relation reproduce `sec-uniform-boundary`'s statement for $n=1$ (see also the brief for problem 12).
2. **Solitonic duals.** Compute the dual $\widetilde{K^{(b)}_a}=K^{(b)}_a|_{B\to2-B}$ and its classical limit $B\to0$. That limit is the $B\to2$ limit of $K^{(b)}_a$, the strong-coupling limit of the solitonic boundary. Find a real-coupling boundary with a vacuum whose linearized reflection factors reproduce it. Candidates: Neumann with a stationary soliton (a "dual solitonic" boundary), and boundaries of the dual exponentials with mixed signs.
3. **Boundary energies.** Run the boundary TBA (diagonal, ground state) with $K^{(b)}_a$ and with its dual for $a_2$ and $a_3$. Compare the UV boundary entropy $\log g$ and the boundary energy with candidates. The $g$-function (`sec-boundary-finite-size`) is a duality-covariant check: dual boundaries should have equal $g$ if they are the same quantum boundary.
4. **Coideals.** For $n=1$ the map is known through the UV–IR relation. Write $\hat\epsilon_{0,1}$ as functions of $(\tilde\mu,\phi_0)$ (`prp-dm-first-order`) and of the dual parameters, and find the transformation $\hat\epsilon\to\hat\epsilon'$ at $\omega\to\omega'$. Then test whether the $a_n$ map $\{0\}\leftrightarrow\{s=+1\}$ and $\{s=-1\}\leftrightarrow\{C\equiv-1\text{ minimal}\}$ is the restriction of a natural map, for example an exchange of the coideal of $U_q$ with a coideal of $U_{\tilde q}$ in the modular double.

## Known traps

- Duality acts on $B$ by $B\to2-B$, on $\omega$ by $\omega\to-\omega/(\omega+1)$, and on $\beta$ by $\beta\to4\pi/\beta$. Check every step on the bulk blocks first.
- Kim's self-dual Neumann conjectures agree with DG only at one loop (`exr-kim-expansion`); one loop cannot test duality.
- The solitonic amplitudes depend on the class $\{b,h-b\}$ only (`prp-solitonic-amplitude` (ii)), so their duals do too.

## Deliverable

- The CFT-level statement of half-line duality for $a_n$, with the parameter map.
- A proposal (with evidence) for the boundary dual to each solitonic class, or a proof that no real-coupling boundary with the given classical limit exists among the known types.
- The $n=1$ coideal transformation and the status of a general coideal-level statement.
- Proposed text for open problem 4 and for the paragraph "Duality on the half-line" of `sec-neumann-duality`.

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium-low.** Steps 1 and 3 are concrete (literature translation and TBA numerics). Step 2 may produce a clean answer or nothing. Step 4 is open-ended research with no existing algebraic framework for coideals of the modular double.
