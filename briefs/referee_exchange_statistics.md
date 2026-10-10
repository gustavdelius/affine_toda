# Brief: referee the exchange-statistics results

You are an adversarial mathematical-physics referee. Your job is to find errors, gaps, overclaimed status labels and uncredited prior art. Do not confirm what you have not checked.

The material is uncommitted work in the git worktree `/home/gustav/Git/affine_toda/.claude/worktrees/exchange-statistics` (branch `worktree-exchange-statistics`). The main checkout at `/home/gustav/Git/affine_toda` does **not** contain it, so read every file in the worktree. Do not edit any file there. Work in your own scratch directory, copy scripts there before running them, and deliver your full report in your final message: subagents may be unable to write report files.

## Context

The brief `briefs/exchange_statistics.md` asked whether the Bethe–Yang equations of imaginary-coupling affine Toda solitons use the ordinary flip of soliton labels or a $\mathbb Z_2$-graded flip. The answer now in the book is that the question dissolves:

1. a graded flip is exactly the ordinary flip with a Bloch-angle twist (`prp-exchange-twist`);
2. on a circle, the vacuum shifts $\vec\phi\mapsto\vec\phi+\frac{2\pi}\beta\vec\lambda$ split the states into Bloch sectors $\alpha$, the Bethe–Yang equations in sector $\alpha$ carry the twist $e^{-i\alpha\cdot\mu_j}$, and the global form of the theory decides which $\alpha$ occur (`prp-bloch-sectors`);
3. the sector $Q=0$, $\alpha=0$ exists in every global form, and the bosonic breaking found earlier sits there, so the bosonic statements hold in every form.

See all changes with `git -C /home/gustav/Git/affine_toda/.claude/worktrees/exchange-statistics diff` plus the new folder `foundations_code/sectors/bloch/`. Locations (line numbers as of 2026-10-10):

| Location | What |
|---|---|
| `part7-foundations/32-scattering.qmd:185–240` | paragraph "Exchange statistics and Bloch sectors": `prp-exchange-twist` [Theorem] with proof, `prp-bloch-sectors` [Proposition, sketch] with sketch, three paragraphs of consequences, results at general $\alpha$, Krein pairing |
| `32-scattering.qmd:143–145, 153–155` | revised and new rows of the sector table |
| `32-scattering.qmd:118, 179, 181, 183` | one-sentence additions ($a_n^{(1)}$ colour sector, excited solitons, $\psi=\pi$ exception, spinor sectors) |
| `part7-foundations/33-outlook.qmd:55, 69, 115, 168–169` | "with bosonic exchange" replaced by "in the Bloch sector $\alpha=0$"; open problem 6 sub-item "Bloch sectors" |
| `part6-soliton-smatrices/21-construction.qmd:72`, `26-f4-solitons.qmd:161` | the crossing sign $(-1)^n$ is neither a grading nor the symmetry type of $C$ |
| `references.bib` | new key `bajnok2000` |
| `foundations_code/sectors/bloch/` | scripts `t1`–`t13`, outputs in `out/`, `README.md` with arguments and conventions |

Conventions: $q^2=-e^{i\psi}$, $\psi=2\pi^2/3\xi$ reduced to $(-\pi,\pi]$, in the conventions of Takács–Watts (hep-th/9810006); the code uses $q_{\rm code}=\bar q_{\rm TW}$. The Bethe–Yang operator is the particle-labelled $R=P\check R$, $T_1=R_{1N}\cdots R_{12}$. Spinor amplitudes use $q=e^{-i\pi\omega}$. $Q$ is the total topological charge, $\mu_j$ that of soliton $j$, and $\Omega_\alpha=e^{i\alpha\cdot\mu_1}$. "Unbroken" means $\max\lvert\lvert s\rvert-1\rvert<10^{-6}$ over the rapidities sampled.

## Already checked by the author (redo only if you find a reason)

- `t1`, `t5`: $T_1^{\rm gr}=(-1)^{g_1(G-1)}T_1$ as an operator identity, residual $\le3\times10^{-15}$, for $a_2^{(2)}$ kinks ($N\le4$) and the $c_2^{(1)}$, $c_3^{(1)}$, $a_5^{(2)}$ spinors ($N\le3$).
- `t3`: two kinks, $Q=0$: the bisected critical twist agrees with $\max(0,\min(2\lvert\psi\rvert-\pi,\pi-\lvert\psi\rvert))$ to $10^{-5}$ at 13 couplings.
- `t4`: for $N=3,4,5$ kinks at $\alpha=0$, just above $\lvert\psi\rvert=\pi/N$ only $Q=0$ breaks.
- `t11`: at $\psi=\pi$ the particle-labelled kink R-matrix is diagonal to $10^{-15}$.
- `t6`: Bloch-torus counts for the spinor sector $Q=0$, reproduced identically on a re-run.
- `t7`, `t8`: $R_{12}(x)^\dagger=R_{12}(1/x)$, $\mathcal CR_{12}\mathcal C=R_{21}$ for Izergin–Korepin and spinors; $R_{12}^{\mathsf T}=R_{21}$ for spinors only.
- `t12`, `t13`: excited-soliton breaking at $\alpha=0$ lies in $Q=0$; two-kink $Q=\pm1$ unbroken for all $\alpha\in[-\pi,\pi]$.
- The book renders with no unresolved cross-references.

## To be checked, in priority order

1. **The proof of `prp-exchange-twist`.** Re-derive it independently, without reading the proof first.
   - Check the graded embedding of $R^{\rm gr}_{12}$ into the factors $(1,j)$ and the sign $(-1)^{g_1^{(j)}g'_j+(g_j+g'_j)s_{j-1}}$.
   - Check the telescoping, whether "outgoing labels" and left multiplication are stated correctly, and whether "preserves the total grading" is the right hypothesis.
   - Check the second sentence, the twist $\alpha=\pi(G-1)\gamma$ "up to a phase common to all states". Is that phase really common within each sector?
   - Verify symbolically or numerically for $N=2,3$ with a random even matrix, not only the R-matrices of the book.
   - Search for prior art: the equivalence of graded and ungraded Bethe ansatz up to twists may be folklore (Jordan–Wigner; supersymmetric $t$-$J$ and Hubbard chains; graded quantum inverse scattering). If it is known, the label can stay but the result must be credited.
2. **`prp-bloch-sectors` and its sketch.**
   - (a) Can one-soliton states always be chosen covariant under the $U_\lambda$? Is the R-matrix used in the code the covariant one, with no vacuum-dependent phases?
   - (b) Check the form and sign of the twist, and that $\Theta=PK$ preserves $\alpha$: $\Theta U_\lambda\Theta^{-1}=U_\lambda^{-1}$ with $\Theta$ antilinear.
   - (c) Read Bajnok–Palla–Takács–Wágner (hep-th/0004181), especially eqs. 2.9–2.10, 3.8–3.9 and section 5.3. Record exactly what they establish and in which conventions, and whether the book's statement attributes to them only what they show.
   - (d) Check the global-form statement: $Q\in k\Lambda_W^\vee$, the allowed characters, the $k=1$ case, and "not identified, every $\alpha$, only $Q=0$ with periodic $\vec\phi$". Does the last clause ignore twisted-boundary sectors that a non-compact theory might include?
   - (e) In rank $\ge2$: is the Bloch torus correctly parametrised for the spinor charges, given that they generate a lattice larger than $\mathbb Z^n$?
   - (f) Check the Ramond/Neveu–Schwarz remark.
3. **The convention argument in `32-scattering.qmd:230`.** A change of coproduct by a bicharacter $F$ is said to conjugate $\check R$ by $F$, with the symmetric part acting by a similarity on $T_1$ and the antisymmetric part shifting $\alpha$ linearly in $Q$.
   - Derive this, for example with $G=\prod_{i<j}F(\mu_i,\mu_j)$ and the cyclic shift.
   - Are there other convention changes that move the $Q=0$, $\alpha=0$ sector, such as a different coproduct ordering, $q\to q^{-1}$, or parity?
4. **The scope of "in every global form" and "in every case computed".**
   - Is there a natural global form without the sector $Q=0$, $\alpha=0$, for example a fermionic form with only the Ramond sector, or boundary conditions twisted by charge conjugation or the diagram automorphism?
   - Is the list of cases in line 230 complete and accurate?
   - Does the book's lattice, a Dirichlet strip (`thm-krein-transfer`), fix a global form at all?
5. **Numerical claims: re-run and stress them.**
   - The two-kink threshold formula: more couplings, both signs of $\psi$, finer rapidity grids including $\lvert\log x\rvert<10^{-2}$ and large $\log x$.
   - The $\pi/N$ thresholds by sector for $N\le5$.
   - The $\psi=\pi$ exception for $N\ge4$.
   - The graded rule "unbroken iff $\lvert\psi\rvert\le2\pi/N$". It was only spot-checked at $N=4,5$, and at $N=5$ no coupling above $2\pi/5$ was tested.
   - The spinor torus counts: grid resolution, and rapidity sampling near $\theta=0$ given the window $\lvert\theta\rvert<\theta_c$.
   - The $(s,1/\bar s)$ pairing within each $(Q,\alpha)$ sector to $10^{-9}$.
   - The table row "[Numerical], $N\le3$, spot checks at $N=4,5$".
6. **The spinor Krein-pairing argument at the end of the new block.**
   - Check that the three identities plus Yang–Baxter give pairing within $(Q,\alpha)$. In particular, check the reordering similarity $R_{12}\cdots R_{1N}\sim R_{1N}\cdots R_{12}$: it must commute with $\Omega_1$ and preserve $Q$.
   - Can the identities be proved rather than checked numerically?
   - For Izergin–Korepin, find a basis or gradation in which $R^{\mathsf T}=R_{21}$, or another proof (`t9` tried diagonal $x^{am}$ gauges only).
7. **The crossing-sign remarks.**
   - Check that $C^{-1}C^{\mathsf T}$ is diagonal and that $\epsilon_n=-1,-1,+1,+1$.
   - Check that the finite part of the $c_n^{(1)}$ spinor's quantum group is $b_n$, and that these are its Frobenius–Schur indicators.
   - Check that "the factor is a scalar, so it is not a grading" is not misleading.
8. **The quantum-trace coincidence (`32-scattering.qmd:238`).** Is $q^{\pm2m}$, in these conventions, the twist implementing the quantum trace of the $U_q(a_1)$ subalgebra of the $(1,2)$ restriction? Check $q$ versus $q_{\rm TW}$ and the power. If not, say what the boundary $\lvert\alpha\rvert=\lvert\arg q^2\rvert$ is. The other branch of the two-kink threshold, $\lvert\alpha\rvert=2\lvert\psi\rvert-\pi$ for $\pi/2\le\lvert\psi\rvert\le2\pi/3$, equals $\lvert\arg(-q^4)\rvert$. Is it the quantum-trace twist of the $U_{q^4}(a_1)$ subalgebra of the $(1,5)$ restriction, up to a sign $(-1)^m$? The book does not state this; propose text only if you can confirm an interpretation.
9. **Labels and wording of every edit.**
   - Check the status labels against the evidence.
   - Check consistency between the table, the paragraphs and `33-outlook.qmd`.
   - Check for overclaiming: for example "every global form", "the S-matrix does not select the statistics", and the Notes line "Krein-unitarity and the sector classification are new here".
   - Check the Pandoc pitfall: a status label and a citation in one bracket get swallowed.
10. **Wider prior art.** Look for:
    - Feverati–Ravanini–Takács on sine-Gordon versus massive Thirring and twisted sectors;
    - non-linear integral equations or TBA for $a_2^{(2)}$ or Izergin–Korepin with a twist;
    - Pasquier–Saleur quantum-group-invariant twists;
    - any earlier treatment of soliton statistics as a boundary condition in affine Toda theory.

## Known traps

- The Bethe–Yang operator is $R=P\check R$, not $\check R$, whose eigenvalues are trivially phases.
- Gradations act by diagonal similarities, which commute with every twist $\Omega_\alpha$. They never change moduli.
- $(Q,\alpha)$ and $(-Q,-\alpha)$ have related spectra (charge conjugation plus the Krein relation). The scans used $Q\ge0$ and mostly $\alpha\ge0$.
- Couplings are offset by $0.0013\pi$ in $\xi$, because roots of unity (for example $\xi=2\pi$) make the projectors singular.
- Twisted-algebra literature normalises the longest root to squared length $2k$, so its coupling is $\beta^2/k$ in the book's normalisation.
- Usage is limited. Do items 1, 2, 4 and 5 first, run computations incrementally, and save outputs.

## Deliverable

In your final message:

- a verdict for each of `prp-exchange-twist`, `prp-bloch-sectors`, the convention argument, the scope claim, each new or revised table row, the Krein-pairing argument and the crossing-sign remark: keep, downgrade (to which label) or fix;
- findings ranked most severe first, each with its location (file and line, or id), a concrete failure scenario, and proposed replacement text in the book's style (short declarative sentences, LaTeX in `$...$`, no hype, labels as **[Theorem]**, [Proposition, sketch], [Established], [Numerical], [Conjecture]);
- the numbers you re-ran, with precision and the scripts used;
- prior art found, with exact citations and what each establishes;
- every doubt, stated plainly.
