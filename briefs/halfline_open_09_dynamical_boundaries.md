# Brief: integrable dynamical boundaries for $a_n^{(1)}$, $n\ge2$ (half-line open problem 9)

You are a mathematical physicist with expertise in classical and quantum integrable boundaries, q-Onsager algebras and integrable defects. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 9 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Higher-rank dynamical boundaries.** The sine-Gordon boundary with an oscillator (`sec-dynamical-boundary`) has no $a_n^{(1)}$ counterpart yet: which dynamical boundary systems are integrable, and what are their coideals and reflection matrices? Baseilhac and Belliard's generalized q-Onsager algebras give candidates.

This overlaps with the planned Part on integrable defects (worktree `integrable-defects-framework-a54a84`, see the memory notes). Coordinate rather than duplicate: this brief is about **boundaries** with degrees of freedom.

## What the book already has

- `sec-dynamical-boundary` (ch. 30): the Baseilhac–Delius–George oscillator boundary for sine-Gordon, $C_0=\cos\hat p$, $C_1=\cos\hat q$, with $\{\hat p,\hat q\}=\beta^2/2$, and its $\Theta$-invariance.
- `sec-bsg-dynamical` (ch. 29).
- `sec-q-onsager` (ch. 34): the q-Dolan–Grady relations, Baseilhac–Koizumi's dynamical reflection matrices, and the statement that a dynamical boundary is a $B$-module $W$ at the right end, with $\hat E_j=\pi_W(b_j)$, which must represent the q-Onsager relations.
- `prp-an-k-existence`: for fixed boundaries of $a_n$, $n\ge2$, $\hat\epsilon\in\{0\}\cup\{\pm\epsilon_*\}^{n+1}$. For a dynamical boundary the $\hat\epsilon_j$ become operators. The discreteness suggests the operators must satisfy $\hat E_j^2=\epsilon_*^2$-type relations, or that the relevant algebra is a different coideal.

## Literature

- Baseilhac–Belliard, arXiv:0906.1215 (`baseilhac_belliard2010`): generalized q-Onsager algebras $\mathcal O_q(\hat g)$ and boundary (soliton non-preserving) affine Toda theories. All integrable boundary conditions are classified by algebraic arguments, and the dynamical case is described through representations of $\mathcal O_q(\hat g)$.
- **Belliard–Fomin, arXiv:1106.1317 (J. Phys. A 45 (2012) 025201; not in `references.bib`):** a procedure that builds **dynamical** (operator-valued) K-matrices from the generalized q-Onsager algebra, with $\mathcal O_q(a_2^{(1)})$ and $\mathcal O_q(a_2^{(2)})$ treated in detail. This is very likely the quantum answer for $a_2^{(1)}$. Check its conventions against the book's coideal (`eq-coideal-generators`).
- Baseilhac–Belliard(–Crampé?), "Higher order relations for ADE-type generalized q-Onsager algebras", arXiv:1312.5897 (not in `references.bib`; check its authors).
- Baseilhac–Koizumi (`baseilhac_koizumi2003`), Baseilhac–Delius–George (`bdg2002`).
- Avan–Doikou, "Boundary Lax pairs for the $A_n^{(1)}$ Toda field theories", arXiv:0809.2734 (NPB 821 (2009) 481; not in `references.bib`): a systematic construction of boundary Lax pairs, with two distinct boundary types for $A_2^{(1)}$. Also their general framework with dynamical boundary matrices in an ultralocal Poisson algebra (arXiv:1110.4235, "Integrable boundary conditions and modified Lax equations"; check).
- **Hills, PhD thesis, University of York (2016), "Generating boundary conditions for integrable field theories using defects"** (etheses.whiterose.ac.uk/16379): generalized, infinite-dimensional reflection matrices for $a_1^{(1)}$, $a_2^{(1)}$ and $a_2^{(2)}$ ATFT, built by dressing an integrable boundary with a defect.
- Corrigan–Zambon, arXiv:1202.6016 (JHEP 06 (2012) 050): infinite-dimensional sine-Gordon reflection matrices from boundary topological charge, related to defects. Bajnok–Simon, arXiv:0712.4292: new integrable boundaries by fusing defects to boundaries.
- Bowcock–Corrigan–Zambon, `bcz2004b`: $a_r^{(1)}$ type-I defects (classical).

## Suggested plan of attack

1. **Classical construction by dressing.** Place an $a_2^{(1)}$ type-I defect (`bcz2004b`) at $x=-d$ in front of a fixed integrable boundary (uniform, Neumann or solitonic), and let $d\to0$. The defect field becomes a boundary degree of freedom with its own Poisson bracket. Derive the resulting boundary Lagrangian and equations of motion, and check classical integrability: the energy and spin-$\pm2$ charges, or a Sklyanin dynamical $K_+$ in the style of `sec-sklyanin`. Compare with Hills's construction and with Avan–Doikou's boundary Lax pairs.
2. **Quantum algebra.** Take Belliard–Fomin's dynamical K-matrix for $\mathcal O_q(a_2^{(1)})$ and translate it to book conventions. Check that it intertwines the book's coideal with operator-valued $\hat E_j=\pi_W(b_j)$, and solves the twisted reflection equation (numerically on a truncated $W$). Identify the module $W$ (q-oscillator? finite-dimensional?).
3. **Match.** Compare the classical limit of step 2 (the Poisson algebra of the $\hat E_j$) with the boundary degrees of freedom of step 1. Check that freezing the boundary degrees of freedom reduces both to the fixed boundaries of `prp-an-k-existence`. Two candidates for this reduction: the one-dimensional modules ($\hat E_j=\pm\epsilon_*$ or $0$), and the defect sent to infinity.
4. **Physical content.** For the simplest finite-dimensional $W$, compute the reflection of the vector soliton including boundary-state transitions, and the analogue of boundary bound states. State whether the boundary is $\Theta$-invariant (`prp-boundary-theta`).

## Known traps

- $H_{\rm osc}$ is complex in sine-Gordon (BDG); expect the same for $a_n$. Use $\Theta$-invariance rather than Hermiticity.
- The sign convention $\dot f=\{H,f\}$ of BDG versus $\dot f=\{f,H\}$.
- Type-I $a_r$ defects keep $Q_{+2}$ or $Q_{-2}$, not both. Higher-spin charges have unremovable quantum anomalies (Penati–Polvara, arXiv:1902.10690). A dressed boundary may inherit these, so quantum statements built on Lagrangian dressing are conjectural.
- Baseilhac–Koizumi's (19) is stronger than the q-Onsager relations (`sec-q-onsager`).

## Deliverable

- A classical $a_2^{(1)}$ dynamical boundary with a proof (or a check) of integrability.
- The Belliard–Fomin dynamical K in book conventions with checks, and the match to the classical boundary.
- A short classification statement of what is known for $a_n$, $n\ge3$.
- Proposed text replacing open problem 9 and a paragraph for `sec-dynamical-boundary` or `sec-q-onsager`.
- Scripts.

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium.** The quantum algebra for $a_2$ seems to exist already (Belliard–Fomin), and the dressing method has been developed in Hills's thesis. The work is translation, verification and matching. **Low** for a general classification. Best run after, or together with, the defects Part.
