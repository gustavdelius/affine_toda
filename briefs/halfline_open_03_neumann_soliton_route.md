# Brief: derive the Neumann particle amplitudes of $a_n^{(1)}$ from the $s=+1$ soliton reflection matrix (half-line open problem 3)

You are a mathematical physicist with expertise in exact S-matrices, quantum groups and the boundary bootstrap. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 3 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`) reads:

> **Neumann from the soliton side.** Derive (`eq-neumann-amplitude`) for $n\ge3$ from the $s=+1$ reflection matrix (`eq-an-k-matrix`). This needs its scalar factor, which requires the eigenvalue of $\check R_{\bar V\bar V}$ in crossing-unitarity, and the singlet projection (`sec-an-nondiagonal-k`). Agreement would confirm both `cnj-neumann-coideal` and `cnj-neumann-amplitude`; disagreement would break the duality picture.

The target is

$$K^N_a(\theta)=\prod_{c=1}^a\frac{(c-1)\,(h-c+\frac B2)}{(h-c)\,(c-1+\frac B2)},$$

the dual ($B\to2-B$) of the uniform amplitude. It is supported by duality, the classical limit, Kim's one loop and the boundary TBA, but it has been derived from solitons only for $n=1$ (Ghoshal–Zamolodchikov at the Neumann point) and $n=2$ (Gandenberger's $K^{(+)}_1=(-2)(-\frac B2)(2+\frac B2)$).

## What the book already has

Read `sec-k-procedure`, `sec-an-diagonal-k`, `sec-an-nondiagonal-k` (ch. 33), `sec-epsilon-constraint` (ch. 32), and `sec-breather-reflection`, `sec-continuation-assumptions`, `sec-neumann-duality` (ch. 35).

- The matrix (`eq-an-k-matrix`) in book conventions, with $X=(-q)^{1/2}y=e^{\omega(\theta-i\pi/2)}$; `prp-an-nondiagonal` shows it is the unique intertwiner, solves the twisted reflection equation, and that **its matrix part satisfies boundary crossing-unitarity** (checked for $n=2,3$). So crossing-unitarity reduces to one scalar function.
- The diagonal case as a template: the functional equations (`eq-an-diagonal-scalar-eqs`), the Gamma-product solution (`eq-an-diagonal-scalar`), the prefactor ambiguity $f(\theta)=-(1-i\sinh\theta)/(1+i\sinh\theta)$, and the breather projection $K_{B_1}=A\,F_{11}(2\mu)\,A$ that gives DG (3.11).
- Scripts: `code/boundary_coideal.py` (the $s=\pm1$ intertwiners), `reflection_matrices.py` (book $\check R$, $F_{11}$, the diagonal scalar factor, the breather projection for $n=2,3$), `particle_reflection.py` (`K_neumann`, the axioms checker, block functions). Copy them to your scratch directory before running.

## Literature

- Gandenberger, hep-th/9806003 (`gandenberger1999a`): $a_2^{(1)}$ non-diagonal $K^\pm$, their scalar factors from boundary unitarity and crossing, breather amplitudes (5.16)–(5.18), and (6.4). This is the $n=2$ template; read its appendix for how the scalar factor was obtained.
- Gandenberger, hep-th/9911178 (`gandenberger1999b`): the $a_n^{(1)}$ non-diagonal solutions. It states explicitly that the scalar factor is **not** determined there and is left open (abstract and conclusion checked).
- Regelskis–Vlaar (`regelskis_vlaar2016`, App. D, (D.2)): the same dense matrix from the q-Onsager-type constraint. No scalar factor.
- Appel–Vlaar (`appel_vlaar2025`): trigonometric K-matrices $K_V(z)=g_V(z)\mathcal K_V(z)$ with a scalar series $g_V$ from the universal K-matrix. This covers quantum symmetric pairs only; the $s=\pm1$ coideals are **not** QSPs (`tbl-physics-qsp`). Their method of fixing $g_V$ through crossing may still be a useful model.
- Ahn–Kim–Rim, hep-th/0110218 (`akr2002`): boundary TBA checks of the uniform and Neumann amplitudes for $a_n$. This is non-perturbative support, independent of the soliton route.
- Ghoshal–Zamolodchikov (`gz1994`) and Ghoshal (`ghoshal1994`): the $n=1$ case, which every general-$n$ formula must reduce to.

## Suggested plan of attack

1. **The scalar crossing function.** For $n=2,3,4,5$ compute, with the book's $\check R$ and $K^{(s=+1)}$, the matrix that appears in boundary crossing-unitarity (`sec-boundary-unitarity-crossing`). By `prp-an-nondiagonal`(iii) it is a scalar $\Phi_n(\theta)$ times the identity. Fit $\Phi_n$ to a ratio of products of $\sin$ or $\sinh$ of linear functions of $\theta$ (its zeros and poles are at rapidities fixed by $\omega$ and $h$). Guess the closed form for general $n$, and prove it if possible, for example from the eigenvalue of $\check R_{\bar V\bar V}$ on the relevant tensor component.
2. **Boundary unitarity.** Compute $K(\theta)K(-\theta)$ for (`eq-an-k-matrix`). It should be a scalar $U_n(\theta)$, since the diagonal entries are $\mp i\sinh\Phi$ or $i\cosh\Phi$ over $\sin\frac{\pi\omega}2$.
3. **Minimal solution.** Solve $A(\theta)A(-\theta)=U_n^{-1}$ and $A(\frac{i\pi}2-\theta)/A(\frac{i\pi}2+\theta)=\Phi_n(2\theta)\times$(the bulk factor), by the standard infinite Gamma-product method used for (`eq-an-diagonal-scalar`). Fix the ambiguous prefactor (the analogue of $f(\theta)$) by demanding the $n=1$ limit (GZ at the Neumann point) and the $n=2$ case (Gandenberger). Record the CDD freedom.
4. **The breather projection.** The lowest breather $B_1$ is the singlet in $V\otimes\bar V$ at the fusing angle (`sec-breather-reflection`). With a non-diagonal $K$ the reflected pair is no longer a sum of identical-soliton states, so the projection is a genuine matrix computation. Do it numerically for $n=2,\dots,5$ with the full $\check R$, at random $\omega$ and $\theta$, and check that the singlet is mapped to a multiple of itself.
5. **Continuation.** Continue the breather amplitude by $\omega\to-2/B$ and compare with $K^N_1$ of (`eq-neumann-amplitude`). Then bootstrap to $K^N_a$ (`prp-uniform-amplitude` (i) machinery) or, better, compute the higher breathers directly.
6. **Semiclassics (optional).** With the scalar factor in hand, test $s=+1$ against the Neumann time advance $\frac1{m_a\sinh\theta}\ln(1-s^2)$ (`sec-boundary-semiclassics`, last paragraph). This is the analogue of `prp-semiclassical-uniform`.
7. **$s=-1$ (bonus).** Repeat steps 1–5 for $s=-1$. For $n=2$ the result must be $K^{(-)}_1=-(1)(-1+\frac B2)(3-\frac B2)$; its dual is the CDRS $C\equiv-1$ "minimal" proposal. A general-$n$ formula would feed problems 1 and 5.

## Known traps

- The sign of $q$: the physical $q=-e^{-i\pi\omega}$. At $-q$ the same matrix fails the reflection equation (`prp-an-nondiagonal`(iv)). Use the book's conventions throughout and check them on the diagonal case first.
- The prefactor of the diagonal scalar factor was once misquoted in the book as $-(\frac h4\omega)_I$. That misquote led to a false "correction" of DG (3.11). An analogous sign ambiguity will appear here, and only the breather limits ($n=1$ GZ, $n=2$ Gandenberger) fix it.
- Gandenberger's free parameters are a torus gauge (`sec-an-nondiagonal-k`). Use the parity-symmetric representative.
- Disagreement is a possible, publishable outcome: it would break the duality picture of `sec-neumann-duality`. Check CDD factors before concluding it.

## Required checks

- $n=1$: Ghoshal's $K_{B_1}$ at $\eta=\frac\pi2(\lambda+1)$, $\vartheta=0$ (book's `eq-bsg-b1-reflection`).
- $n=2$: Gandenberger's (5.16)–(5.18) and (6.4) in book blocks.
- All axioms of `prp-uniform-amplitude` (i)–(iii) for the result, and $K\to1$ as $B\to0$.
- Kim's one-loop $a_3^{(1)}$ Neumann amplitudes (`kim1995`, eqs. (33), (39)) at $O(\beta^2)$.

## Deliverable

In your final message:

- closed forms for $\Phi_n$, $U_n$ and the scalar factor $A_n^{(+)}$ (and $A_n^{(-)}$ if reached), with the precision of every check;
- the breather amplitude and its continuation, and the verdict on `cnj-neumann-amplitude` for $n=3,4,5$ (general $n$ if proved);
- proposed replacement text for the paragraph "Which soliton reflection matrix?" of `sec-neumann-duality` and for `sec-an-nondiagonal-k`, in the book's style (short declarative sentences, LaTeX in `$...$`, status labels);
- scripts, with a one-line description each.

## Feasibility (assessment when the brief was written, 2026-10-11)

**High.** Every input exists in the repository, the matrix part of crossing-unitarity is already known to factor, and there are two exact checks ($n=1,2$). The main work is the closed form of $\Phi_n$ and the Gamma-product solution. Expected effort: one focused agent session.
