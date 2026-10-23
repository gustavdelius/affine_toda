<!-- Worker report, 2026-10-12. Model: Claude Opus 5.5, effort xhigh. Scope: problem 7 parts (a) and (b) only (brief halfline_open_07_nonunitary_half_line.md); part (c) deferred until problem 11 is done. Interrupted twice by the usage limit and resumed. Saved verbatim by the coordinator. Scripts and outputs copied to briefs/halfline_open_07_scripts/; they import boundary_common.py, an unchanged copy of code/boundary_common.py, which must be copied alongside. -->

## Report: half-line open problem 7, parts (a) and (b)

**Scope.** I did parts (a) and (b) for the $a_n^{(1)}$ vector solitons, with side checks for sine-Gordon and $d_4$. Part (c) was not started, as instructed. Nothing in the repository was edited. All scripts and outputs are in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/p07/` (section 6).

### Verdict

- **(a) Yes.**
  - There is a Krein structure for which boundary unitarity is pseudo-unitarity. It uses the same colour reversal $C$ as the bulk (`prp-krein-fundamental`).
  - The reality condition is $\hat\epsilon_{\nu(j)}=\overline{\hat\epsilon_j}$, where $\nu$ is the diagram automorphism $j\mapsto n+1-j$ that fixes node 0.
  - All three $a_n$ solutions satisfy it ($\hat\epsilon=0$, and $s=\pm1$ with real $\epsilon_*=1/(2\cos\frac{\pi\omega}2)$).
  - The double-row transfer matrices are therefore exactly Krein-unitary. Their eigenvalues are phases or pairs $(s,1/\bar s)$.
- **(b) No, not in general.**
  - Every $\mathcal T_j$ commutes with the $\mathbb Z_{n+1}$ cyclic shift $\mathcal Z=Z^{\otimes N}$ (twisted, $Z^{n+1}=s$, for $s=-1$). The Krein form maps $\mathcal Z$ to $\mathcal Z^{-1}$, so it pairs sector $\zeta$ with sector $\bar\zeta$. Only the sectors with $\zeta=\pm1$ can be unbroken.
  - For $a_2$, a single soliton between a uniform wall and a Neumann ($s=+1$) wall, or between two Neumann walls, already has complex energies in the $\zeta=e^{\pm2\pi i/3}$ sectors. This is proved in closed form.
  - In the $\Theta$-invariant sector $\zeta=1$:
    - two solitons are unbroken below coupling thresholds;
    - between uniform walls they are unbroken at every coupling sampled, including the soliton–antisoliton pair, which is broken in the bulk;
    - three solitons are broken for every pair of walls tested.

### 1. Part (a): statement and proof

Conventions are those of the code: Part VI coproduct, $q=-e^{-i\pi\omega}$, principal gradation with $y=e^{\omega\theta}$, and $V$, $\bar V$ as in `an_vector`. $C$ is the basis reversal $v_i\mapsto v_{n+2-i}$ on both $V$ and $\bar V$.

**Theorem (boundary Krein-unitarity, $a_n^{(1)}$).** Let $|q|=1$. Let $\tau$ be the antilinear anti-automorphism with $\tau(e_j)=e_{\nu j}$, $\tau(f_j)=f_{\nu j}$, $\tau(k_j)=k_{\nu j}$.

1. $\tau$ respects the relations and satisfies $(\tau\otimes\tau)\Delta=\Delta\tau$. For real $y$, $\pi_{V(y)}(a)^\dagger=C\,\pi_{V(y)}(\tau a)\,C$, and the same holds for $\bar V(y)$.
2. For real rapidities and all $a,b\in\{V,\bar V\}$, the R-checks normalized to 1 on highest weights satisfy $\check R_{ab}^\dagger(C\otimes C)\check R_{ab}=C\otimes C$.
   - The physical scalar factors are phases: $F_{11}$, and for $V\bar V$ the crossed transmission $S^T(T-\mu)$ of `prp-an-diagonal-scalar`.
   - So every $\check S_{ab}$ is Krein-unitary, including soliton–antisoliton.
3. $\tau(b_j)=e_{\nu j}+q^{-1}f_{\nu j}k_{\nu j}+\bar\epsilon_jk_{\nu j}$. In general $\tau$ also sends $c_j$ to $c_{\nu j}=q^{-2}\bar c_j$, which the parity-symmetric choice $c_j=q^{-1}$ satisfies. Hence:
   - $B_{\hat\epsilon}$ is $\tau$-stable if and only if $\hat\epsilon_{\nu(j)}=\overline{\hat\epsilon_j}$.
   - In that case, by uniqueness of $K$ (Schur), $K^\dagger CK=c_KC$ with $c_K$ real.
4. Explicitly:
   - for $\hat\epsilon=0$, $c_K=1$;
   - for (`eq-an-k-matrix`), $c_K=\dfrac{\cosh2T\theta-s\cos\pi T}{2\sin^2\frac{\pi\omega}2}\ge0$.
5. The physical amplitudes are $\mathbb K_V=AK$ and $\mathbb K_{\bar V}=AK^{\rm T}=C\,\mathbb K_V\,C$.
   - Boundary unitarity $\mathbb K_{\bar V}(-\theta)\mathbb K_V(\theta)=1$ becomes pseudo-unitarity $\mathbb K_V^\dagger C\,\mathbb K_V=C$ if and only if $\overline{A(\theta)}=s\,A(-\bar\theta)$ (with $s=1$ for $\hat\epsilon=0$).
   - Equivalently, $C\,\mathbb K_V(\theta)^\dagger C=\mathbb K_{\bar V}(-\theta^*)$, a boundary Hermitian analyticity twisted by $C$.
   - Ordinary unitarity $\mathbb K^\dagger\mathbb K=1$ fails, because $\mathbb K_V^{-1}=C\,\mathbb K_V(-\theta)\,C$, not $\mathbb K_V(-\theta)$.
6. Every double-row transfer matrix $\mathcal T_j$ on $\bigotimes_kV_{a_k}$ built from these amplitudes satisfies $\mathcal T_j^\dagger\mathcal C\,\mathcal T_j=\mathcal C$, with $\mathcal C=C^{\otimes N}$. This is the conclusion of `prp-krein-unitarity` on the interval, without asymptotic completeness.

**Proof.**

- **Part 1.** Entrywise: $CE_jC=F_{\nu j}$ and $Cq^{H_j}C=q^{-H_{\nu j}}$, and both are compatible with $\bar q=q^{-1}$.
- **Part 2.** Take the adjoint of $\check R\,\Delta_{12}(a)=\Delta_{21}(a)\,\check R$ and use part 1. Then $(C\otimes C)\check R^\dagger(C\otimes C)$ intertwines like $\check R^{-1}$, so it equals $c\,\check R^{-1}$ with $c$ real.
  - Evaluating on the highest-weight vector gives $c$ as the conjugate of the transmission element of the colour-reversed pair.
  - For $VV$ that element is the identical-soliton element, equal to 1 by the $\mathbb Z_N$ rotation.
  - For $V\bar V$ all transmissions equal the crossed $VV$ transmission. I checked this crossing: the spectrum of the partial transpose of $S_{VV}(i\pi-\theta)$ equals $t_{VV}(i\pi-\theta)$ times the spectrum of $S_{V\bar V}(\theta)$ to $5\times10^{-15}$, and fails with $-y$.
  - Hence $c=1$; numerically $|c-1|<10^{-14}$ for all four pairs, $n=2,3$ and $\omega\in\{0.37,1.43,2.6\}$.
- **Parts 3 and 4.**
  - The Schur step is the same argument as in part 2.
  - For $s=\pm1$ there is a direct proof. Three identities hold:
    - $K(\theta)^\dagger=s\,K(-\theta)$ (entrywise);
    - $CKC=K^{\rm T}$;
    - $K(-\theta)^{\rm T}K(\theta)=U=s\,c_K$, which is the unitarity function proved in the problem-3 report and rechecked here.
  - Together they give $K^\dagger CK=sK(-\theta)K(\theta)^{\rm T}C=c_KC$.
- **Part 5.** The scalar factors have the required reality.
  - DG's $A_d$ consists of phases.
  - In the problem-3 $A^{(\pm)}$, the Gamma products and the $(a)_I$ blocks with real $a$ are phases on the real axis. The prefactor $\sin\frac{\pi\omega}2/\sin\pi(\mu-\frac T2)$ (with $\cos$ in the denominator for $s=-1$) has modulus squared exactly $1/c_K$. The extra $-i$ for $s=-1$ supplies $\overline{A}=-A(-\theta)$.
  - I checked this prefactor ($A_{\min}$ in the code) to $10^{-16}$. I did not re-evaluate the full infinite products of problem 3.
- **Part 6.** A product of Krein-unitary maps between Krein spaces is Krein-unitary.

**Left wall.** The left wall is the image of the right wall under the parity map $\omega$ ($e\leftrightarrow f$, $k\to k^{-1}$). This map is an algebra automorphism and a coalgebra anti-automorphism, and $\pi_{V(y)}\circ\omega=\pi_{\bar V(1/y)}$ exactly in the code bases. So $K^L_{\bar a}(-\theta)=K^R_a(\theta)$ as matrices.

- With this convention $[\mathcal T_i,\mathcal T_j]=0$ to $10^{-14}$ for $N=2,3$ and all wall pairs.
- The alternative, with $K$ and $K^{\rm T}$ exchanged on the left, gives $\|[\mathcal T_1,\mathcal T_2]\|\sim0.2$–$7$ whenever the left wall is non-diagonal. So the choice is forced, not a convention, and it matters: for $N=1$ between two Neumann walls the alternative would give $K^{\rm T}K$, which is unimodular.

**Numerical checks** (all to $10^{-14}$ or better):

| Check | Range |
|---|---|
| part 1 | $n=2,3,4$ |
| $K$ (`eq-an-k-matrix`) equals the coideal intertwiner | $n=2,3,4$; $\omega=0.37,1.43,2.6$; $s=\pm1$ |
| $K^\dagger CK=c_KC$ with the formula of part 4 | same |
| $K_{\bar V}=K^{\rm T}$ | same |
| $K^\dagger=sK(-\theta)$ | same |
| part 3, for random complex $\epsilon$ | $n=2$ |

**Comparison with `prp-boundary-theta`.**

- The classical condition is "all $C_j$ real", with no $\nu$. The quantum form $C$ corresponds to field parity composed with $\nu$, so the condition is $\Theta\nu$-invariance: $\Theta$-invariance together with charge-conjugation symmetry.
- For $a_n$ vectors no Krein form with $\nu=\mathrm{id}$ exists: $V$ composed with the untwisted $\tau$ is isomorphic to $\bar V$, so the untwisted form pairs a sector with its charge conjugate.
- The gauge-fixed $a_n$ solutions satisfy both conditions. A Cartan gauge applied to both walls is a similarity. A relative gauge between the walls acts like a twist, which I did not study.

**Sine-Gordon** ($\nu=\mathrm{id}$, $\eta=\sigma^x$):

- $K^\dagger\sigma^xK\propto\sigma^x$ holds if and only if $\epsilon_0,\epsilon_1$ are real. Through the book's $\epsilon(\mu_\pm)$ this is $\phi_0\in\frac{2\pi}\beta\mathbb Z$, which matches `prp-boundary-theta`.
- The second structure ($\nu$ exchanging nodes 0 and 1, $\eta=1$) gives ordinary unitarity if and only if $\epsilon_0=\bar\epsilon_1$.
- With $K_{12}=1$, $c=1+4\sin^2(\pi\omega)\dfrac{\epsilon_0\epsilon_1(y^2+y^{-2})+\epsilon_0^2+\epsilon_1^2}{(y^2-y^{-2})^2}$. This can be negative when $\epsilon_0\epsilon_1<0$ and $|\epsilon|>1/|\sin\pi\omega|$, and pseudo-unitarity then holds only up to sign. For $a_n$, $c_K\ge0$ always.

**$d_4$ vector.**

- With $\nu=\mathrm{id}$, $\eta$ is unique: a signed anti-diagonal ($v_i\leftrightarrow v_{-i}$) with $\eta^2\propto1$.
- $R^\dagger(\eta\otimes\eta)R=\eta\otimes\eta$ with $c=1.000000$.
- All the $K$ of `prp-dn-k` ($\hat\epsilon=0$ and all $\pm\epsilon_*$ patterns) satisfy the relation, with $c>0$ in the samples.
- With $\nu=(n{-}1\ n)$, patterns with $\epsilon_{n-1}\ne\epsilon_n$ fail (residual 1.0), as the reality condition predicts.

**Bulk by-product.** This settles the remark after `prp-krein-fundamental`: soliton–antisoliton amplitudes are Krein-unitary with the *untwisted* $C$ in the principal gradation. Bulk $V\bar V$ eigenvalues come in pairs off the unit circle, as Takács and Watts found. This also covers the crossed case of outlook problem 5 for $a_n$.

### 2. Part (b): $a_2^{(1)}$ results

Method and notation:

- I diagonalized $\mathcal T_j$ jointly with $\mathcal Z$ at real rapidities, using the deviation $\max\big||\Lambda|-1\big|$.
- Bulk scalar phases are omitted (they are unimodular). The boundary uses $A_{\min}$ ($|A_{\min}|^2c_K=1$).
- Walls: U is uniform ($K=1$), N is $s=+1$, M is $s=-1$; written left/right.

**One soliton, closed form.**

- $K$ is a polynomial in $Z$, with eigenvalue $\kappa(\zeta)\propto(X\zeta+q)/(X\zeta-1)$.
- On the $\zeta$-mode, $|\Lambda_{\rm UN}|=\rho^{\pm1}$ and $|\Lambda_{\rm NN}|=\rho^{\pm2}$, where $\rho=\dfrac{|e^{\omega\theta}\zeta-w|}{|e^{\omega\theta}\zeta-\bar w|}$ and $w=e^{-i\pi\omega/2}$.
- $\rho=1$ only for real $\zeta$. The ratio of moduli between the $\zeta$ and $\bar\zeta$ sectors is independent of the scalar factor, so the breaking is robust.
- Bethe–Yang solutions at $mL=10$, mode $k$; energies $E/m$:

| Walls | $\omega$ | $\zeta=1$, $k=1,2,4,8$ | $\zeta=e^{\mp2\pi i/3}$, $k=1,2,4,8$ |
|---|---|---|---|
| UN | 0.37 | 1.01391, 1.11808, 1.52688, 2.63090 (real) | $1.04737\mp0.00513i$, $1.17857\mp0.00901i$, $1.59991\mp0.01281i$, $2.69407\mp0.01410i$ |
| NN | 1.43 | real | $1.09058\mp0.04041i$, …, $2.79452\mp0.01190i$ |

**Table 1.** Maximum deviation by $\mathbb Z_3$ sector, random rapidities in $[0.05,2.5]$:

| Walls, content | $\omega=0.37$ | $0.83$ | $1.43$ | $2.6$ |
|---|---|---|---|---|
| UU, $N=1$ | 0 | 0 | 0 | 0 |
| UU, VV ($\zeta=1$ / $\zeta\ne1$) | 0 / 1.1 | 0 / 1.6 | 0 / 2.2 | 0 / 0.42 |
| UU, VW ($\zeta=1$ / $\zeta\ne1$) | 0 / 1.2 | 0 / 1.7 | 0 / 4.4 | 0 / 4.8 |
| UU, VVV (all sectors) | 0.79 | 0.25 | 0.46 | 0.46 |
| UN, $N=1$ ($\zeta=1$ / $\zeta\ne1$) | 0 / 0.42 | 0 / 1.5 | 0 / 5.3 | 0 / 6.6 |
| NN, $N=1$ ($\zeta=1$ / $\zeta\ne1$) | 0 / 1 | 0 / 5.3 | 0 / 20 | 0 / 14 |
| MM, $N=1$ ($\zeta=-1$ / other) | 0 / 8.9 | 0 / 49 | 0 / 2 | 0 / 1.3 |
| NM, $N=1$ (no $\mathbb Z_3$) | 1.3 | 3 | 3.2 | 1 |

Here 0 means below $10^{-12}$.

**Table 2.** Sector $\zeta=1$, fine logarithmic rapidity grids from 0.001 to 10 (the breaking first appears at $\theta_2\to0$ with $\theta_1$ large):

| Walls, content | Unbroken | Broken |
|---|---|---|
| UU, VV and VW | all $\omega$ sampled, 0.1–3 | none sampled |
| UN, VV | $\omega\le0.68$ | $\omega\ge0.70$ (threshold near $2/3$) |
| UN, VW | $\omega\le0.65$ | $\omega\ge0.68$ |
| NN, VV | $\omega\le0.53$ | $\omega\ge0.54$ |
| NN, VW | $\omega\le0.50$ | $\omega\ge0.51$ |
| UU, VVV | — | every $\omega$ except the special points $\omega=1,2,3$ ($q=\pm1$) |
| UN and NN, VVV | — | every $\omega\ge0.2$ tested |
| MM, VV | $1.4\le\omega\le2.5$ (coarse sampling only) | elsewhere |

- At $\omega=2$ every case is unbroken. NN VV is also unbroken at the isolated point $\omega=0.67$.
- A first coarse scan suggested NN VVV was unbroken at small $\omega$; that was a sampling artefact, which the fine grid removed.

**Comparison with the periodic bulk** (maximum deviation of $T_1$, same four $\omega$):

| Content | Deviation |
|---|---|
| VV | 0 |
| VW | 1.2–5.1 |
| VVV | 0.14–0.34 |
| VVW | 0.3–2.2 |

- The boundary *restores* reality for the soliton–antisoliton pair in the $\zeta=1$ sector of UU, and below threshold for N walls.
- It *breaks* reality for one soliton ($\zeta\ne1$, whenever a wall is non-diagonal), and for two solitons with N walls above threshold.
- Three solitons stay broken, as in the bulk three-colour sector.

**Caveats.**

- I did not settle whether every $\zeta$-sector is physical. For Neumann walls they should be the Fourier modes of the degenerate boundary vacua; this needs TCSA or the lattice.
- The thresholds are numerical brackets with no analytic explanation.
- The $s=\pm1$ scalar factors come from the unrefereed problem-3 report. Only their moduli matter here, and I checked those analytically.

### 3. Proposed text for open problem 7

"7. **The non-unitary half-line theory.** For the $a_n^{(1)}$ vector solitons, (a) and (b) are now settled; see @prp-boundary-krein [new: the theorem of section 1, with its proof].

*Answer to (a).* The antilinear anti-automorphism $\tau(e_j)=e_{\nu j}$ ($\nu$: $\alpha_j\to\alpha_{n+1-j}$) gives the colour-reversal Krein form of @prp-krein-fundamental on $V$ and $\bar V$. It maps $B_{\hat\epsilon}$ to itself exactly when $\hat\epsilon_{\nu(j)}=\overline{\hat\epsilon_j}$, which all three solutions of @prp-an-k-existence satisfy. Boundary unitarity then becomes pseudo-unitarity $\mathbb K^\dagger C\mathbb K=C$, and the double-row transfer matrices are Krein-unitary. The quantum condition is $\Theta$-invariance together with charge conjugation; @prp-boundary-theta has no $\nu$.

*Answer to (b): no.* The cyclic $\mathbb Z_{n+1}$ commutes with the transfer matrices, and the Krein form pairs its sector $\zeta$ with $\bar\zeta$. For $a_2$ one soliton between a uniform and a Neumann wall already has complex levels for $\zeta\ne1$. In $\zeta=1$, two solitons are unbroken between uniform walls at every coupling sampled, and with Neumann walls below $\omega\approx0.5$–$0.7$. Three solitons are broken [Numerical].

*Still open:*
- (b′) An explanation of the thresholds.
- Whether all $\zeta$-sectors are physical (boundary vacua; TCSA).
- The $d_n$, $s=-1$ and mixed-wall cases beyond the checks above.
- Relative gauges between the two walls.
- (c) as before: are the Dirichlet walls of the constructive continuum limit (@sec-collapse-thresholds, @thm-window) integrable boundaries? On real field configurations a Dirichlet wall satisfies the soliton-preserving condition (-@eq-soliton-preserving), but that condition was derived for complex classical solutions, and whether the Euclidean theory with such a wall is the soliton-preserving integrable boundary is open."

Related edits:
- `sec-krein-unitarity`: replace "We have not checked whether part 2 extends…" with the untwisted-$C$ result for $V\bar V$.
- `eq-boundary-unitarity`: needs the $K^{\rm T}$ form, as problem 3 already flags.

### 4. Not done

- Part (c), as instructed.
- Full Newton solutions for $N\ge2$; for $N=1$ they are in the table above.
- An analytic proof of the thresholds.

### 5. Brief premises that turned out wrong

- No twisted $C$ is needed for soliton–antisoliton amplitudes.
- "Real for $\Theta$-invariant boundaries" is false as stated. The relevant condition is $\Theta\nu$-invariance, and even then reality holds only in self-conjugate sectors.

### 6. Scripts

All in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/p07/`. Run them from that directory.

- `boundary_common.py` — unchanged copy of `code/boundary_common.py`.
- `p07_common.py` — $C$, $\nu$, cached R-checks, `k_book` (analytic branch, $\epsilon_*=1/(2\cos\frac{\pi\omega}2)$), $c_K$, $A_{\min}$.
- `a_krein_algebra.py` — part (a) checks 1–5 (reps, R-checks, K-matrices, coideal $\tau$-stability, sine-Gordon).
- `a_crossing_vw.py` — the $V\bar V$ amplitude is the crossed $VV$ amplitude.
- `a_dn.py` — the $d_4$ Krein form, R-check and K-matrices.
- `b_transfer.py` — double-row transfer matrices, with the $\omega$ and swap left-wall conventions.
- `b_checks.py` — commutativity, Krein-unitarity and $[\mathcal Z,\mathcal T]$.
- `b_scan.py` → `b_scan.out` — Table 1 and the bulk comparison.
- `b_threshold.py` → `b_threshold.out` — coarse $\omega$ grid.
- `b_refine.py` → `b_refine.out`.
- `b_grid2.py` → `b_grid2.out` — two-soliton rapidity grids.
- `b_bisect.py` → `b_bisect.out`, and `b_bisect2.py` — Table 2.
- `b_by_n1.py` — the one-soliton closed form and the Bethe–Yang energies.

The scratch directory is temporary, so copy anything worth keeping.

## Integration (branch half-line-open-problems-integration)

The Opus follow-up prompt `halfline_open_07_opus_prompt.md` was not launched. Its three items were handled by the coordinator during integration (rows "claim 3", "$\zeta\ne1$ sectors" and "claim 6" below). New script: `code/boundary_krein.py` (6 checks, all PASS, seconds), using only book functions.

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 1, parts 1–4 (boundary Krein-unitarity) | **Verified, integrated** as [Theorem] `prp-boundary-krein` | Checks 1–4: $\pi(x)^\dagger=C\pi(\tau x)C$ (exact); Krein relation for all four $\check R_{ab}$ ($9\times10^{-15}$); $\tau(b_j)=b_{\nu j}(\bar\epsilon_j)$ (exact), with failure for non-real patterns; $K^\dagger CK=c_KC$ and $CKC=K^{\rm T}$. The proof was redone by hand. | New `sec-boundary-krein` in ch. 35, before the open problems |
| 2 (reality of the scalar factors) | **Accepted** for $\hat\epsilon=0$ (DG's $A_d$ is a phase) and $s=+1$ (Gemini checked the full products; $|A^{(+)}|^2c_K=1$ follows from `prp-an-nondiagonal-scalar`) | | Part 3 of the theorem is stated as an iff |
| 3 (left-wall convention forced) | **Consistent with the book**: `sec-boundary-bethe-yang` already says the left wall is the mirror image of a right wall; no correction needed. The commutativity numbers (worker, Gemini) were not rechecked. | | One sentence in `sec-boundary-krein` |
| 4, 5 (sine-Gordon, $d_4$) | **Not integrated** | Not rechecked | — |
| 6 (bulk $V\bar V$ Krein-unitary with untwisted $C$) | **Verified** (check 2), and `prp-krein-fundamental` uses the same colour reversal | | ch. 40 remark after `prp-krein-fundamental` rewritten |
| 7 ($\mathbb Z_{n+1}$ sectors; one soliton UN) | **Verified** | Check 6 ($[Z,K]=0$, $[Z\otimes Z,\check R]=0$, $CZC=Z^{-1}$); check 5 (moduli $1$, $\rho^{\pm1}$ with the worker's $\rho$; the coordinator's first transcription used $X$ instead of $y$ and failed, then was corrected) | `sec-boundary-krein`, "Reality of the spectrum", with the caveat that the $\zeta\ne1$ sectors may not be physical (the coordinator's analysis of this question was not done; it is stated as open) |
| 8 (two- and three-soliton tables, thresholds) | **Verified, integrated** (consolidation, branch `half-line-consolidation`, 2026-10-13), with one item **refuted** | Coordinator's own double-row transfer matrices from `rcheck` and the analytic-branch $K$ (scratch `a4/dr.py`; book check 7 of `boundary_krein.py`): $[\mathcal T_1,\mathcal T_2]=0$ to $10^{-13}$, and 0.92 with the unmirrored left wall (so the convention is forced); UU unbroken at $\omega=0.1$–$1.43$; UN, NN unbroken at $\omega=0.1,0.2,0.37$ and below threshold; thresholds bisected at $\theta_1=12$: UN VV $(0.6981,0.6984)$, UN VW $(0.6659,0.6669)$, NN VV $(0.5361,0.5363)$, NN VW $(0.5006,0.5009)$, stable for $\theta_1=10$–$25$ (below threshold the deviation is rounding, $\le2\times10^{-9}$ at $\theta_1=25$); VVV broken for UU, UN, NN at $\omega=0.3,0.5,0.83,1.43$. **Refuted:** "NN VV is unbroken at the isolated point $\omega=0.67$" (deviation 2.1 at $0.6667$ and at $\frac23$). MM not rechecked. | `sec-boundary-krein`: new paragraph "Two and three solitons" and a sentence on the mirror convention; open problem 7 updated |
| 9 (premise of (b) false as stated) | **Integrated** in substance | | Open problem 7 rewritten |
