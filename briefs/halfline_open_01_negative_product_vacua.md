# Brief: vacua and particle amplitudes for the $a_n^{(1)}$ boundaries with $\prod_jC_j=-1$ (half-line open problem 1)

You are a mathematical physicist with expertise in classical and quantum affine Toda theory on the half-line. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 1 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Boundaries with $\prod_jC_j=-1$.** For $a_n^{(1)}$ with these sign patterns no regular real vacuum is known, except $\vec\phi=0$ for $C\equiv-1$ with $h$ odd (`sec-real-coupling-vacua`), and no particle reflection amplitudes have been derived [DG sec. 7]. Are there vacua at all, and if so, which soliton configuration do they continue?

Related loose ends: the quantum point $s=-1$ has no known boundary condition (`sec-epsilon-constraint`, `sec-neumann-duality` "The third matrix"). The CDRS $C\equiv-1$ amplitudes built on $\vec\phi=0$ differ from DG's (`sec-solitonic-boundaries`, "What is assumed"). For $h$ even, $C\equiv-1$ is solitonic and has the self-conjugate vacuum family (`prp-self-conjugate-vacuum`).

## What the book already has

Read `sec-images-boundary`, `sec-real-coupling-vacua`, `sec-classical-reflection-factors`, `sec-boundary-types` (`tbl-boundary-types`, last two rows) in ch. 31; `sec-epsilon-constraint` in ch. 32; `sec-neumann-duality` and `sec-solitonic-boundaries` in ch. 35. Scripts: `code/boundary_vacua.py`, `boundary_images.py`, `boundary_classical.py`, `particle_reflection.py`. Facts to use:

- At imaginary coupling, for $h$ odd, every pattern with $\prod C_j=-1$ is a coweight shift of $C\equiv-1$; for $h$ even no image solution is known (`tbl-boundary-types`).
- At real coupling, $\vec\phi=0$ satisfies the boundary condition only when all $C_j$ are equal (or Neumann). For $C\equiv-1$ the linearized factor (`eq-linear-reflection`) has classical boundary bound states for every $a\ne h/2$.
- The $s=-1$ breather amplitude for $a_2$, $K^{(-)}_1=-(1)(-1+\frac B2)(3-\frac B2)$, tends to 1 as $B\to0$. Its dual is CDRS's minimal $C\equiv-1$ proposal.

## Literature

- **Bowcock–Perkins, hep-th/9909174 (`bowcock_perkins2003`).** By a machine summary: a Bogomolny bound on the half-line and stability iff $A_i\le1$ (Bowcock's sign; $C_i\ge-1$ in the book's), so all $C=\pm1$ boundaries are bounded below; tables of vacua for $A_r^{(1)}$, $r\le5$, by sign pattern, including flat, single-pair, multi-soliton and "exceptional" backgrounds; and an observation that numerically the non-flat solutions occur "when a maximal number of tau functions vanish at the boundary". **Start by extracting their entries for $\prod C_j=-1$ patterns.**
- Bowcock, hep-th/9609233 (`bowcock1998`), and Perkins's thesis (Durham etheses 4560).
- Corrigan–Dorey–Rietdijk–Sasaki (`cdrs1994`) and Corrigan–Dorey–Rietdijk (`cdr1995`): the $C\equiv-1$ conjectures on $\vec\phi=0$ and the classical bound states.
- Fateev–Onofri, hep-th/0207152 (not in `references.bib`): classical vacua from the classical limit of boundary one-point functions. Only positive signs are treated; negative signs are "more subtle and will be considered in a separate publication". If such a publication exists, find it.
- Fateev, hep-th/0103014 (not in `references.bib`): boundary one-point functions and conjectured boundary ground-state energies in affine Toda theory.
- Ahn–Kim–Rim (`akr2002`): their conclusion lists boundaries with $\eta_i=-1$ as open.

## Suggested plan of attack

1. **Literature extraction.** Translate the Bowcock–Perkins tables for $\prod C_j=-1$ into book conventions (`tbl-boundary-dictionary`: Bowcock's $A_i=-C_i$). Recompute every claimed vacuum with `boundary_vacua.py` conventions: regularity on $x\le0$, the boundary condition, the energy in the book's normalization $\mathcal B(0)=0$.
2. **Numerical minimization.** For every $\prod C_j=-1$ pattern of $a_2,a_3,a_4$ up to the symmetries that survive at real coupling (the diagram automorphisms; there is no coweight shift), minimize the static energy on $[-L,0]$ from many random starts, as in the brief for problem 6. Classify the outcome: a regular minimizer; a flat vacuum $\vec\phi=0$ when all $C_j=-1$; or no minimizer (runaway to a boundary singularity). Shared code with problem 6 is natural.
3. **Identify the continuation.** For each regular minimizer, find the imaginary-coupling static configuration it continues. Candidates: a coweight-shifted image of the $C\equiv-1$ boundary with a stationary soliton (for $h$ odd), multi-soliton backgrounds, and solutions with $\tau_j(0)=0$ for a maximal set of $j$.
4. **Is $\vec\phi=0$ the vacuum of $C\equiv-1$ for $h$ odd?** Compare its energy with any other regular static solution found in step 2. The linearized spectrum about $\vec\phi=0$ has bound states but no negative modes; check this too.
5. **Particle amplitudes (stretch).** If a vacuum is found, compute its linearized reflection factors. Then test candidates for the quantum amplitude: DG's soliton route (`cnj-reflection-continuation`) using the $C\equiv-1$ soliton reflection matrix with a stationary soliton, the $s=-1$ breather amplitude and its dual (needs the scalar factor of problem 3 for $n\ge3$), and the CDRS proposals. Check their $B\to0$ limits against the linearized factors.
6. **The $s=-1$ question.** Report whether any $\prod C_j=-1$ boundary has a classical limit compatible with $K^{(-)}\to1$ as $B\to0$ (a trivial classical reflection, like Neumann). This would suggest assigning $s=-1$ to it.

## Known traps

- Real-coupling theories with patterns related by coweight shifts are different; the shift is an imaginary-coupling symmetry only.
- $C\equiv-1$ with $h$ even is solitonic ($\prod=+1$) and belongs to problem 6 and `prp-self-conjugate-vacuum`, not here.
- Energies depend on the normalization of $\mathcal B$; compare in one normalization.
- A configuration with a zero of some $\tau_j$ exactly at $x=0$ satisfies several patterns at once.

## Deliverable

- A table, for $a_2$–$a_4$ (and $a_5$ from Bowcock–Perkins), of $\prod C_j=-1$ patterns, vacuum (formula or numerics), energy, continuation, and linearized reflection factors.
- A verdict on "are there vacua at all", with the evidence.
- Any amplitude candidates and their checks.
- Proposed replacement text for the last row of `tbl-boundary-types` and for the closing paragraph of `sec-real-coupling-vacua`.
- Scripts.

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium-high for the classical part.** It is the same numerics as problem 6, and Bowcock–Perkins may have done much of it. **Low for the quantum amplitudes:** they need either a new soliton reflection matrix with a background or the $s=-1$ scalar factor, and there are no independent checks beyond the classical limit.
