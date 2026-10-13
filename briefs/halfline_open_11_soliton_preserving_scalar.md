# Brief: scalar factors and breather amplitudes for soliton-preserving reflection in $a_n^{(1)}$ (half-line open problem 11)

You are a mathematical physicist with expertise in exact S-matrices, the reflection equation and quantum groups. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 11 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Scalar factors for soliton-preserving reflection.** The soliton-preserving matrices (`sec-preserving-conjugating`) have no scalar factors and no breather amplitudes yet, so the soliton-preserving boundary has no known particle amplitudes.

## What the book already has

Read `sec-soliton-preserving`, `sec-dynamical-boundary` and `sec-boundary-theta` (ch. 30); `sec-images-preserving` (ch. 31); `sec-reflection-equation` and `sec-boundary-unitarity-crossing` (ch. 28); `sec-preserving-conjugating` with `tbl-k-types` (ch. 33); and the last row of `tbl-physics-qsp` (ch. 34). Facts to use:

- Classically (`prp-soliton-preserving`, [Theorem]): with $\sigma$ = Hermitian conjugation and $K_+=1$, the boundary condition at imaginary coupling is Dirichlet on $\operatorname{Re}\vec\phi$ (at a vacuum, $\frac{2\pi}\beta\Lambda_W^\vee$) and Neumann on $\operatorname{Im}\vec\phi$. There is no free parameter. Images are $\Theta P$ images (`eq-theta-p-reflection`).
- Quantum: the reflection matrices map $V_\mu\to V_\mu$ and obey the braid-form reflection equation. The candidates are the diagonal de Vega–González-Ruiz matrices $K_a(\theta)\propto\sinh(\xi\pm\theta)e^{\pm\theta}$ with a discrete parameter $l$, which are **unchecked in book conventions**. Their coideal is expected to be a QSP with $\tau\ne{\rm id}$ (Regelskis–Vlaar family A.3 with $M_2=0$, matched only by shape).
- For $a_1^{(1)}$ the soliton-preserving condition is Dirichlet $\phi(0)\in\frac{2\pi}\beta\Lambda_W^\vee$: the GZ Dirichlet matrix, whose breather amplitude continues to the sinh-Gordon Dirichlet factor $K_D=(2-\frac B2)(1+\frac B2)/(1)$, with $K_D\to-1$ as $B\to0$.

## Literature

- de Vega–González-Ruiz (`devega1993`): diagonal solutions for $a_{n-1}$ vertex models.
- Regelskis–Vlaar, arXiv:1803.06491 (J. Phys. A 51 (2018); not in `references.bib`): **the complete set of invertible solutions of the untwisted and twisted reflection equations for the $A^{(1)}_{N-1}$ Bazhanov–Jimbo R-matrix**, all obtained by affinization of constant solutions. This is the classification against which "which soliton-preserving matrix" should be decided.
- Regelskis–Vlaar (`regelskis_vlaar2016`) and Appel–Vlaar (`appel_vlaar2025`, secs. 7.7–7.8): the QSP origin of the untwisted matrices, and trigonometric K-matrices with a scalar series $g_V(z)$ from the universal K-matrix, with crossing.
- Arnaudon, Avan, Crampé, Doikou, Frappat, Ragoucy, "General boundary conditions for the $sl(N)$ and $sl(M|N)$ open spin chains", math-ph/0406021 (not in `references.bib`): reflection equations for both SP and SNP, classification, transfer matrices and analytical Bethe ansatz. Their eigenvalue formulas contain the crossing properties of SP matrices.
- Doikou (`doikou2000`) and earlier Doikou–Nepomechie work on SP/SNP: dressing and crossing in spin-chain language.
- Delius, hep-th/9809140 (`delius1998b`): the classical soliton-preserving condition.
- Lima-Santos and Malara, nlin/0412058 (not in `references.bib`): regular solutions of the (untwisted) reflection equation for non-exceptional vertex models; useful for $d_n$ and $c_n$ extensions.

## Suggested plan of attack

1. **Which matrix?** Solve the soliton-preserving intertwining problem in book conventions (step 2 of `sec-k-procedure`, $V(\theta)\to V(-\theta)$) for coideals of the A.3 type, for $a_2$ and $a_3$. Match the solutions with Regelskis–Vlaar 2018's list and with de Vega–González-Ruiz. Use classical data to decide which one is the boundary of `prp-soliton-preserving`: it has no free parameter, and its classical limit should be the $\Theta P$-image time delays of `sec-images-preserving`.
2. **Crossing-unitarity for SP.** Derive the soliton-preserving form of boundary crossing-unitarity. The crossed channel involves $\check R_{V\bar V}$, because the crossed soliton is the antisoliton (GZ's derivation, adapted). Reduce it to a scalar function, as for the diagonal soliton-conjugating case (`eq-an-diagonal-scalar-eqs`).
3. **Minimal scalar factor.** Solve the unitarity and crossing equations by Gamma products. Fix the prefactor and CDD freedom by the $n=1$ limit (GZ Dirichlet scalar factor at $\phi_0=0$) and by the soliton bootstrap condition $\sigma(\mu)=\prod_k\sigma(\mu+\frac h2\omega-k\omega)$ (`sec-an-diagonal-k`).
4. **Breather amplitudes.** Bootstrap the lowest breather (the singlet in $V\otimes\bar V$) and continue by $\omega\to-2/B$. Check every axiom of `prp-uniform-amplitude` (i)–(iii).
5. **Which real-coupling boundary?** Decide which real-coupling boundary the continued amplitudes describe. There is a tension to resolve. Classically, at real coupling with real fields, Hermitian conjugation equals transposition, so $K_+=1$ with $\sigma=\dagger$ is the Neumann condition. On the other hand, for $n=1$ the continuation of the soliton-preserving (Dirichlet) breather amplitude is the sinh-Gordon **Dirichlet** factor $K_D$, not Neumann. Determine which statement generalizes, using the $B\to0$ limit of the continued amplitude. Candidates: $-1$ for every particle (Dirichlet $\vec\phi(0)=0$ at real coupling; check whether this is classically integrable for $a_n$), or $1$ (Neumann). If the limit is 1, compare the full amplitude with (`eq-neumann-amplitude`); this would give a second soliton route to the Neumann amplitudes and a cross-check of problem 3.

## Known traps

- The braid-form reflection equation and the twisted one involve different R-matrices (`tbl-k-types`); the crossing relation differs accordingly.
- Physical $q=-e^{-i\pi\omega}$ (the sign matters in the scalar factor, not always in the reflection equation).
- The de Vega–González-Ruiz matrices are written in a homogeneous-gradation, additive-rapidity convention; convert to the principal gradation before comparing.
- The Dirichlet value must be a vacuum: $\vec\lambda\in2\Lambda_W^\vee$ in Delius's form (`sec-soliton-preserving`, image form).

## Deliverable

- The soliton-preserving reflection matrix in book conventions, with its identification in Regelskis–Vlaar 2018.
- Closed-form scalar factor and breather amplitudes, all axiom checks, and the $n=1$ reduction.
- A resolution of the real-coupling interpretation (step 5).
- Proposed text for `sec-preserving-conjugating` and a new row or paragraph for ch. 35.
- Scripts.

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium-high** for steps 1–4: standard machinery with an exact $n=1$ check, and Regelskis–Vlaar 2018 removes the guesswork about which matrices exist. **Medium** for step 5, which needs a classical computation (integrability of real-coupling Dirichlet for $a_n$) and may change the book's reading of the soliton-preserving boundary.
