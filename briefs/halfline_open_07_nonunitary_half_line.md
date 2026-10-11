# Brief: the non-unitary half-line theory — Krein structure, reality of boundary Bethe–Yang spectra, Dirichlet walls (half-line open problem 7)

You are a mathematical physicist with expertise in non-unitary integrable QFT, Krein spaces and constructive field theory. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 7 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`) asks three questions:

> (a) With a boundary: is there a Krein structure (`sec-krein-unitarity`) for which boundary unitarity (`eq-boundary-unitarity`) is pseudo-unitarity?
> (b) Are the boundary Bethe–Yang spectra (`sec-boundary-bethe-yang`) real for $\Theta$-invariant boundaries (`prp-boundary-theta`)?
> (c) Are the Dirichlet walls of the constructive continuum limit (`sec-collapse-thresholds`, `thm-window`) integrable boundaries? On real field configurations a Dirichlet wall $\vec\phi=0$ satisfies the soliton-preserving condition (`eq-soliton-preserving`), but that condition was derived for complex classical solutions, and whether the Euclidean theory with such a wall is the soliton-preserving integrable boundary is open.

## What the book already has

Read `sec-krein-unitarity`, `prp-krein-fundamental`, `sec-sector-status` (`part8-foundations/40-scattering.qmd`); `sec-boundary-bethe-yang` and `sec-boundary-finite-size` (ch. 28); `sec-boundary-theta` and `prp-boundary-theta` (ch. 30); `sec-an-diagonal-k`, `sec-an-nondiagonal-k`, `sec-dn-vector-k` (ch. 33); `sec-collapse-thresholds` (ch. 38) and `thm-window` (ch. 39). Facts to use:

- In the bulk, $T_1^\dagger=\mathcal C\,T_1^{-1}\mathcal C$ for the fundamental-soliton transfer matrices (`prp-krein-fundamental`, [Theorem] for the pairing). So eigenvalues are phases or pairs $(s,1/\bar s)$, and the sector table records where reality is broken.
- All rigid boundaries $C_j\in\{0,\pm1\}$ are $\Theta$-invariant (`prp-boundary-theta`).
- The lattice of Part VIII has Dirichlet rows with site-by-site Wick ordering; `thm-window` covers the boundary layer.
- Scripts: `foundations_code/krein_bethe.py` (bulk Krein and Bethe–Yang checks), `textbook_code/boundary_coideal.py`, `reflection_matrices.py`.

## Literature

- Takács–Watts, "Non-unitarity in quantum affine Toda theory and perturbed CFT", hep-th/9810006 (`tw1999`?), and "RSOS revisited", hep-th/0203073 (check the bib keys): real finite-volume spectra for all $a_1^{(1)}$-based theories, generally complex for $a_2^{(2)}$, with exceptions. This is the bulk precedent for (b).
- Exact solutions of non-Hermitian $\mathcal{PT}$-symmetric spin chains with integrable boundary fields: arXiv:2301.06004 and arXiv:2406.10334 (a dissipative Kondo effect at the edge; not in `references.bib`). There, $\mathcal{PT}$ symmetry with complex-conjugate boundary parameters gives real ground-state rapidities, with $\mathcal{PT}$-broken phases where energies are complex. This is the lattice analogue of (b) and a model for the phase diagram to expect.
- Bender et al. on $\mathcal{PT}$-symmetric QFT, for the language of (a); Mostafazadeh's pseudo-Hermiticity.
- Bajnok–Böhm–Takács, hep-th/0309119: the boundary reduction formula and the boundary analogue of unitarity (Cutkosky rules), useful for stating pseudo-unitarity of $K$ on the boundary Fock space.
- For (c): Ang–Remy–Sun(–Zhu) on probabilistic boundary Liouville (cited in the outlook, ch. 42) as the model of a constructive boundary theory that is matched to integrable data.

## Suggested plan of attack

1. **(a) Algebra.** For the soliton reflection matrices of the book (the $a_n$ diagonal $K$, the $s=\pm1$ matrices, the $d_n$ matrices, the sine-Gordon GZ matrix at real $\hat\epsilon$), test whether $K(\theta)^\dagger=\mathcal C_{\bar V}\,K(\theta)^{-1}\,\mathcal C_V$ (or a variant with the conjugation of the bulk Krein form) at real $\theta$, up to the scalar factor. Then prove it from the coideal: $\Theta$-invariance of the boundary should make $B_{\hat\epsilon}$ stable under the antilinear map that defines the Krein form, when $\hat\epsilon$ is real (or satisfies a reality condition). Identify that condition, and compare it with `prp-boundary-theta` ($C_j$ real).
2. **(a) Double-row transfer matrix.** Show that Sklyanin's double-row transfer matrix built from Krein-unitary $\check S$ and $K$ is Krein-unitary on the boundary Fock space. This gives the conclusion of `prp-krein-unitarity` on the interval.
3. **(b) Numerics.** For $a_2^{(1)}$ with two boundaries (uniform/uniform, Neumann/Neumann, uniform/$s=+1$), solve the boundary Bethe–Yang equations (nested) for two and three fundamental solitons on an interval at several couplings. Record where eigenvalues leave the unit circle. Compare with the bulk sector table (`sec-sector-status`): does the boundary restore or break reality relative to the periodic case?
4. **(c) Classical.** Determine whether the Dirichlet wall $\vec\phi=0$ is classically integrable as a boundary of the **Euclidean** imaginary-coupling theory with real fields. At imaginary coupling, $K_+=1$ with $\sigma=\dagger$ gives Dirichlet on $\operatorname{Re}\vec\phi$; real Euclidean fields have $\operatorname{Im}\vec\phi=0$. Check whether the charges built from $\tau(\lambda)$ are conserved on real configurations, and whether the identification survives quantization (counterterms at the wall: compare the site-by-site Wick ordering of ch. 37–39 with boundary Wick ordering in boundary CFT).
5. **(c) Quantum, a first test.** If a lattice computation is cheap with existing code, compare the boundary free energy of the Part VIII lattice with Dirichlet rows (extrapolated) against the boundary ground-state energy predicted by the boundary TBA with the soliton-preserving amplitudes of problem 11, once those exist. Otherwise, state the test precisely.

## Known traps

- The bulk theory is not unitary in the QFT sense; only Krein-unitarity holds (`sec-unitarity`). Do not impose ordinary boundary unitarity $K^\dagger K=1$.
- `prp-krein-unitarity` is [Proposition, sketch]; `prp-krein-fundamental` is the proved statement. Distinguish them.
- Part (c) mixes Euclidean real fields with Minkowski complex classical solutions; state which one each claim refers to.

## Deliverable

- (a) A precise statement with proof (or counterexample) of boundary Krein-unitarity, and the reality condition on $\hat\epsilon$.
- (b) Tables of boundary Bethe–Yang spectra and their reality, with code.
- (c) A verdict on the classical integrability of Dirichlet walls in the Euclidean theory, and a precise formulation of the quantum question.
- Proposed text for open problem 7 (split into what is solved and what is open).

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium-high for (a)**: the algebra is a small extension of `prp-krein-fundamental`. **Medium for (b)**: nested boundary Bethe–Yang numerics for $a_2$ are standard but fiddly. **Low for (c)**: it touches the constructive theory and needs the amplitudes of problem 11.
