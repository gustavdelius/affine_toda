<!-- Worker report, 2026-10-12. Model: Claude Opus 5.5, effort xhigh. Scope: brief halfline_open_02_beyond_an.md with coordinator priorities (step 1, then step 4, then step 3, then step 2); given problem 3's report and the 06_01 vacuum code as templates. Interrupted once by the usage limit and resumed. Saved verbatim by the coordinator. Scripts copied to briefs/halfline_open_02_scripts/ (dn_common.py imports code/boundary_common.py read-only from the repository). -->

## Report: half-line open problem 2, the d_n ε̂=0 boundary and a map of the rest

**Status.**
- **Step 1 ($d_n$, $\hat\epsilon=0$): done.** This covers the scalar factor, the breather, the continued particle amplitudes, the axiom checks, the comparison with CAF and FO, and the classical identification for $d_4$ and $d_5$. Two items stay open: the spinor amplitudes for odd $n$, and $d_n$ vacua for $n\ge6$.
- **Step 4 (the map): done.** It rests on the literature, with items marked where I did not verify them.
- **Step 3 ($c_2$): not done.** I only analysed its premises; two of them are wrong (see §5).
- **Step 2 ($d_n$, $s=\pm1$): not started.**
- Nothing in the repository was edited. Scripts are listed in §7.

### 1. Soliton level, book conventions
Notation: $h=2n-2$, $T=(n-1)\omega=h\omega/2$, $\mu=-iT\theta/\pi$, blocks $(a)_I$ of DG, $g(\mu,c)=\Gamma(\mu+c)/\Gamma(-\mu+c)$.

**The $\hat\epsilon=0$ matrix and its unitarity.**
- $K(\theta)=A(\theta)\,J$. Here $J$ is the intertwiner $V(y)\to V(1/y)$ at $\hat\epsilon=0$. It is antidiagonal with entries $\pm1$, independent of $y$, and $J^2=1$ ($n=4$, three values of $y$, residual 2e-15).
- Since $\bar V(y)\cong V(y)$ and $K$ is unique, every member of the multiplet, the antiparticle included, is reflected by the same $K$. So unitarity is simply $A(\theta)A(-\theta)=1$. There is no $K^{\rm T}$ analogue.

**Crossing.**
- The crossing matrix $C$ ($V(y)^*\cong V(-y/q)$) is antidiagonal, symmetric, with signs $\pm1$, and $C^2=1$.
- The crossing factor of the R-matrix, from $(R^{-1})^{t_1}=c_u\,(C^{-1}\otimes1)R(x x_c)(C\otimes1)$, is
  $c_u=\dfrac{\sin\pi\mu\,\sin\pi(\mu+T-\omega)}{\sin\pi(\mu+\omega)\sin\pi(\mu+T)}$
  ($n=4,5$, $\omega=0.37,1.31$, 1e-14).
- The boundary-state vector $(C\otimes1)\mathrm{vec}(J^{\rm T})=\sum\pm v_{-i}\otimes v_{-i}$ has weights $-2\epsilon_i$, so it lies in $L(2\omega_1)$. There $\check R=1$ (checked, 1e-13), so $\Phi\equiv1$.
- Hence the equations are DG's (3.8) with the $d_n$ factor: $A(\mu)A(-\mu)=1$ and $A(\tfrac T2-\mu)/A(\tfrac T2+\mu)=F(2\mu)$.

**The vector–vector scalar factor** (new; normalized to 1 on $L(2\omega_1)$, $F(0)=-1$):

$$F(\mu)=\prod_{k\ge0}\frac{g(\mu,1{+}T{+}2Tk)\,g(\mu,1{+}\omega{+}2Tk)\,g(\mu,T{+}2Tk)\,g(\mu,2T{-}\omega{+}2Tk)}{g(\mu,1{+}2T{+}2Tk)\,g(\mu,1{+}\omega{+}T{+}2Tk)\,g(\mu,2Tk)\,g(\mu,T{-}\omega{+}2Tk)}$$

- It is unitary, and $F(T-\mu)=F(\mu)\,c_u(-\mu)$ ($n=4,5,6$, $\omega=0.37,1.31$, ≤1e-7).
- Physical-strip poles: $\mu=\omega$ (the bound state $V_2$, from $\check R$), $T-\omega$ (crossed), $T-m$ (breathers, from $\check R$) and $m$ (crossed breathers). Zeros at $1+\omega+m$ cancel the unphysical poles of $\check R$.
- **Not done:** a comparison with Johnson 1997, and the bulk breather test (that the lowest breathers reproduce $S_{11}$ of $d_n$).

**The boundary scalar factor:**

$$A(\mu)=\prod_{k\ge0}\frac{g(2\mu,1{+}3T{+}4Tk)\,g(2\mu,1{+}\omega{+}2T{+}4Tk)\,g(2\mu,T{+}4Tk)\,g(2\mu,2T{-}\omega{+}4Tk)}{g(2\mu,1{+}2T{+}4Tk)\,g(2\mu,1{+}\omega{+}T{+}4Tk)\,g(2\mu,2T{+}4Tk)\,g(2\mu,3T{-}\omega{+}4Tk)}$$

- $A(0)=1$. It has a simple zero at $\theta=i\pi/2$ and no other poles or zeros in the strip.
- There is **no** $(-T/2)_I$ prefactor. Adding one, as in DG's $a_n$ formula, makes crossing fail by the sign $-1$; here the sign comes from $g(\mu,0)$ in $F$.
- Both equations hold for $n=4,5,6$ and $\omega=0.37,1.31,2.2$ (≤1e-6).

**Breather.**
- The singlet projection $(1\otimes J)\check R(1\otimes J)$ was computed as a full matrix product. It maps the $B_1$ singlet to exactly 1 times itself ($n=4,5$, two values of $\omega$, residual 1e-14).
- The breather amplitude is

$$K_{B_1}=A(\mu-\tfrac T2+\tfrac12)\,F(2\mu)\,A(\mu+\tfrac T2-\tfrac12)=\Big(\tfrac\omega2-T\Big)_I\Big(\tfrac12-T\Big)_I\Big(-\tfrac\omega2-\tfrac12\Big)_I\,\frac{(\frac T2-\frac\omega2)_I\,(\frac T2+\frac\omega2+\frac12)_I}{(\frac T2)_I\,(\frac T2+\frac12)_I}.$$

- This is DG (3.11) with $T=(n-1)\omega$, times four extra blocks.
- I derived it by pairing the Gamma progressions by hand. It agrees numerically to ≤1e-6 ($n=4,5,6$; three values of $\omega$).

### 2. Particle level ($\omega\to-2/B$)

$$K_1(\theta)=\frac{(h-\frac B2)\,(\frac h2-1)\,(\frac h2+1-\frac B2)}{(h-1)\,(1-\frac B2)\,(\frac h2)\,(\frac h2-\frac B2)},\qquad h=2n-2.$$

**Checks.** The continuation was checked numerically. The following hold for $n=4..7$, in exact block arithmetic to 1e-30:
- unitarity, and crossing-unitarity $K_a(\theta)K_a(\theta+i\pi)=S_{aa}(2\theta)$;
- the bulk boundary bootstrap at every fusing $a+b\to a+b$ with $a,b\le n-2$;
- the closure $1+2\to1$ at $u=\pi-\pi/h$.

The amplitude is not self-dual.

**Relation to Castro-Alvaredo–Fring (hep-th/0308004).**
- $K_a=\hat R_a=\tilde R_a(\theta,B)\prod_jS_{aj}(\theta,1-\frac B2)$. This is exactly their (Rhat): their self-dual "universal" solution dressed by their ambiguity (am3).
- It holds for all vector-like $a$ and $n=4..7$, and $a_2$'s uniform amplitude is the same dressing.
- So the soliton route **derives** CAF's $\hat R$ for $d_n$.
- Their $\tilde R(\theta+i\pi)$ is Kim's self-dual Neumann conjecture, which is relevant to step 2, not here.

**Spinors.**
- $\hat R_s=\hat R_{s'}$ satisfies all axioms and the bootstraps $1+s\to s'$ and $s+s'\to1$ for even $n$ (4, 6).
- **For odd $n$ (5, 7)** it satisfies crossing-unitarity only with $S_{ss'}$ in place of $S_{ss}$ (error 1e-30 versus 0.1). The same holds for the $s+s\to1$ bootstrap, which applies for odd $n$ (I checked that this is the correct bulk fusion).
- This failure is in CAF's (dn2) $\tilde R_s$ itself. For odd $\ell$ it violates their own (cuR). This is a **candidate entry for Appendix D**.
- The d_5 linearization (§3) shows that the spinors reflect diagonally, so the odd-$n$ spinor amplitudes are open. A correction factor $X$ with $X(\theta)X(\theta+i\pi)=S_{ss}/S_{ss'}$ and $X\to1$ classically is needed.

**Classical limit.**
$$K_1^{\rm cl}=-\frac{1}{(1)(h-1)}\cdot\frac{(\frac h2-1)(\frac h2+1)}{(\frac h2)^2}=\prod_jS_{1j}(\theta,B=1).$$

### 3. Which classical boundary: $\hat\epsilon=0$ is $C\equiv+1$ (d_4, d_5) [Numerical]

**Vacua.** Units $\beta=m=1$, book $V$ and $\mathcal B$. Here $\vec\phi=0$ violates the boundary condition. In the sector invariant under the diagram automorphisms the problem reduces to one field:

| | $d_4$ | $d_5$ |
|---|---|---|
| field | $\vec\phi=f(e_2-e_3)$ | $\vec\phi=f(e_2-e_4)$ |
| reduction | Tzitzéica, $f''=2e^{2f}-2e^{-f}$ | sinh-Gordon, $f''=4\sinh f$ |
| solution | $u=e^f$, $u'=(u-1)\sqrt{2u(u+2)}$ | $\tanh(f/4)=e^{2(x-x_0)}$ |
| boundary value | $e^{f(0)}=9/8$ | $e^{f(0)}=4-2\sqrt2$, so $e^{f(0)/2}=m_1$ |
| energy ($m/\beta^2$) | $-0.069521$ | $-0.092709$ |

- Both solutions are regular on $x\le0$ and satisfy the boundary condition exactly.
- **Not proved to be the vacuum.** A lattice energy minimization for $d_4$ (8 random starts, poorly converged) found nothing lower; its best value is $-0.0692$.

**Linearized reflection about these vacua.**
- $d_5$, particle 1: $K=K_1^{\rm cl}$ to 1e-8 at $\theta=0.3$ and $0.6$.
- $d_4$: the three degenerate light particles (1, 3, 4) reflect diagonally, each with $K_1^{\rm cl}$.
- $d_5$, spinors: diagonal, equal to $\hat R_s(B{=}0)$.
- The $a_n$-type factor $-1/((1)(h-1))$ is clearly different.

So $\hat\epsilon=0$ is the uniform boundary, as for $a_n$. The difference is that the vacuum is a static kink rather than $\vec\phi=0$, and this adds the factor $(\frac h2-1)(\frac h2+1)/(\frac h2)^2$.

**Remaining ambiguity.** Possible CDD factors that tend to 1 as $B\to0$, such as DG's $k_c$, are not excluded. No one-loop data exist for this boundary. The minimal $A$, without $(T/2)_I^2$, is fixed by the classical limit.

**Fateev–Onofri observation.**
- $K_1$ equals FO (r1c) exactly as a block identity ($n=4..7$, with $h\to2n-2$, $x\to B/2$, $w_0=w_r=x/2$).
- This is their $B_r$ specialization at $w_r=x/2$.
- $b_r^{(1)}$ is the folding of $d_{r+1}^{(1)}$, and $C\equiv+1$ folds to $C_{\rm short}=1$, which is $w_r\to0$ classically. So the match is suggestive.
- The FO $x$ ↔ book $B$ map for $b_r$ has not been matched (the brief's trap).

### 4. Corrections to the brief
- hep-th/9907072 and hep-th/0002213 are **bulk** (UV, Liouville-type) reflection amplitudes, by Ahn–Fateev–Kim–Rim–Yang and Ahn–Baseilhac–Fateev–Kim–Rim. They are not boundary amplitudes. Only akr2002 is boundary work, and its "D, E blocked by CDD" claim is unverified.
- FO's $C_r$ family is the branch with $C_{\rm short}=0$ and $C_0,C_n=\cos\pi w_{0,n}$. My translation is book $C_j=A_j\sqrt{\alpha_j^2/2}$. Only positive signs are treated. It is **not** the all-$\pm1$ branch that contains the $c_2$ $C\equiv-1$ vacuum. So step 3's comparison does not apply as planned.
- Folding $a_3$ *quantum* amplitudes cannot give $c_2$ quantum amplitudes, because the $c_2$ S-matrix (DGZ) is not the restricted $a_3$ one. Only the $B\to0$ limits fold. The step-3 targets would be the folded $-(1)(3)$ and $-(2)^2$ at $h=4$, which I did not recompute.

### 5. Step 4: proposed table to replace the prose of problem 2

| algebra | soliton K | scalar factor | particle amplitudes | vacuum | identification | missing |
|---|---|---|---|---|---|---|
| $a_n^{(1)}$ | $\hat\epsilon=0$, $s=\pm1$ | all (DG; problem 3, $n\le5$) | uniform, Neumann, $s=-1$, solitonic | $\vec\phi=0$; continued pairs; BP | $C\equiv+1$, Neumann | boundary for $s=-1$; proofs for all $n$ |
| $d_n^{(1)}$ | $\hat\epsilon=0$, $s=\pm1$ (prp-dn-k) | $\hat\epsilon=0$: this report; $s=\pm1$: none | $\hat\epsilon=0$: CAF $\hat R$ (vector-like; spinors, even $n$); Kim self-dual Neumann conjecture $=\tilde R(\theta+i\pi)$ | $C\equiv+1$: invariant static kink ($d_4$, $d_5$) | $\hat\epsilon=0$ is $C\equiv+1$ ($d_4$, $d_5$) | $s=\pm1$ scalar factors and breathers (Kim $d_4$ one-loop test); odd-$n$ spinors; vacua for $n\ge6$ and other patterns |
| $e_{6,7,8}^{(1)}$ | none | – | CAF $\tilde R$, $\hat R$; Riva (bound states, $E_n$); Kim | $\vec\phi=0$ not a solution; unknown | none | everything on the soliton side; vacua |
| $b_n^{(1)}$ | none (GMW bulk) | – | CAF; FO $B_r$ ($C_{\rm short}=\cos\frac{\pi w_r}2$) | FO τ-functions (unchecked) | FO classical limits; $d_{n+1}$ match (§3) | soliton K; FO↔book coupling map |
| $c_n^{(1)}$ | none (Part VI spinor bulk) | – | CAF; FO $C_r$ (short $C=0$ branch only) | images by folding; $c_2$, $C\equiv-1$ family | none | spinor K; FO translation; anything for the all-$\pm1$ branch |
| $g_2^{(1)}$, $f_4^{(1)}$ | none ($f_4$ bulk open) | – | CAF | unknown | none | all |
| $a_{2n}^{(2)}$ ($BC_r$) | none | – | FO conjecture (BC1)→(BC2) | FO | flows between two classical boundaries | K; checks |
| $a_{2n-1}^{(2)}$, $d_{n+1}^{(2)}$ | none | – | FO as $B_r^\vee$, $C_r^\vee$ (dual of $b_n$, $c_n$); CAF as dual pairs | FO | FO: dual pairs are one quantum theory | K |
| $e_6^{(2)}$, $d_4^{(3)}$ | none | – | CAF (as dual pairs with $f_4$, $g_2$) | unknown | none | all |

The CAF coverage is from their section list. I verified it only for $d_n$, and for $a_2$ via their text.

### 6. Proposed text

**Open problem 2 (replacement):**
> **Beyond $a_n^{(1)}$.** For $d_n^{(1)}$ the $\hat\epsilon=0$ vector reflection matrix now has its scalar factor and breather amplitude (@sec-dn-vector-k). The continued particle amplitudes are Castro-Alvaredo and Fring's $\hat R$, and their classical limit is the linearized reflection about the uniform vacuum of $d_4$ and $d_5$, a static kink (@sec-dn-uniform). Open: the scalar factors of the $s=\pm1$ matrices, with Kim's $d_4$ one-loop Neumann amplitude as the test; the spinor amplitudes for odd $n$; and the vacua for $n\ge6$ and for the other sign patterns. For the other algebras see @tbl-beyond-an. In no case beyond $a_n$, $d_n$ and $a_1$ is a soliton reflection matrix known. Fateev and Onofri's $B_r$, $C_r$ and $BC_r$ families have continuous parameters and belong to branches other than the all-$\pm1$ ones of @sec-images-cn.

**New ch. 33 subsection.** Append to sec-dn-vector-k a paragraph "Scalar factor of the $\hat\epsilon=0$ matrix" ([Numerical], $n=4,5,6$). It should contain:
- $K=AJ$ with $J^2=1$, and $\Phi\equiv1$ (with the $L(2\omega_1)$ argument);
- the formulas for $F$ and $A$ of §1;
- the note that $A$ has no $(-T/2)_I$ prefactor.

**New ch. 35 subsection `sec-dn-uniform`.** It should contain:
- the breather formula, and $K_1$ of §2 with its checks ([Numerical]);
- $K_a=\hat R_a$ (CAF);
- the $d_4$/$d_5$ vacua and the linearization match ([Numerical]), stating explicitly that lowest energy is assumed, not proved;
- the odd-$n$ spinor caveat.

### 7. Scripts
All are in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/p02/`. They import `code/boundary_common.py` read-only from the repository. Copy them soon, since this is a temp directory.

**Libraries and checks:**
- `dn_common.py`: weight-restricted Ř, K intertwiner.
- `dn_scalar.py`: $F$, $A$, $K_B$ numerics and blocks.
- `proj_dn.py`: singlet projection.
- `dn_particle.py`: $d_n$ S-matrix, blocks; bulk checks.
- `dn_boot.py`: vector-like bootstrap.
- `caf.py`, `caf_all.py`: CAF $\tilde R$ and $\hat R$; spinors.
- `fo_cmp.py`: the FO identity.

**Classical:**
- `dn_classical.py`: masses, static shooting, `linear_K`.
- `sym_vac2.py`: vacua and the linearization match.
- `spin_cl.py`: spinor classical check.
- `dn_min.py`: lattice minimization.

**Small tests:**
- `t1.py`: $J$.
- `t2.py`: failed invariant-vector search.
- `t3.py`: $c_u$ and $C$.
- `t4.py`, `t5.py`: $F$, $A$, $K_B$.
- `t6.py`: $K_1$ and continuation.
- `t7.py`: axioms.
- `t10.py`: spinor crossing.
- `t11.py`: $\Phi\equiv1$.

**Superseded:** `sym_vac.py`, `t8.py`, `t9.py`, `fix_prof.py`, `fix2.py`.

**Downloaded sources:** `dl_caf/` (CAF TeX) and `dl_fo/` (FO TeX, arXiv abstracts).

## Integration (branch half-line-open-problems-integration)

The Opus follow-up prompt `halfline_open_02_opus_prompt.md` was not launched; the coordinator verified instead. New script `code/dn_uniform.py` (5 checks, all PASS, about 30 s).

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 1 ($K=AJ$, $J$ antidiagonal $\pm1$, unitarity $A(\theta)A(-\theta)=1$) | **Verified** | Check 1 (coideal intertwiner, $n=4,5$) | ch. 33 new paragraph and `prp-dn-scalar` |
| 2 (crossing factor 1 on $L(2\omega_1)$) | **Verified by argument**: the boundary vector has only extremal weights $-2\epsilon_i$, which lie in $L(2\omega_1)$ where the book's `eq-dn-vector-r` is 1. The generic `rcheck` was too slow for $d_n$ to test numerically. | | ch. 33 |
| 3, 4 ($F$, $A$, no $(-T/2)_I$ prefactor) | **Verified** (check 2, 3; own Stirling/Hurwitz products, $\le10^{-26}$). The crossing factor $c(\mu)$ is the worker's formula; the chain is validated end to end by check 4. The "with prefactor crossing fails" statement was not tested; the text says only that there is no prefactor. | | `prp-dn-scalar`, `eq-dn-scalar` |
| 5 (breather amplitude closed form) | **Verified** (check 3). The singlet coefficient 1 was not recomputed; it is stated as an assumption ("with the matrix part acting trivially on the singlet") and supported by check 4. | | `eq-dn-breather` |
| 6 ($K_1$; axioms; equals CAF $\hat R_1$) | **Verified**: check 4 uses the book's $d_n$ S-matrix (`smatrix_simply_laced.setup`), independent of the worker's hand-coded `S_dn`; check 5 recomputes CAF's $\hat R_1$ from their block definition (blockw). Only $a=1$ checked, not all vector-like $a$. | | ch. 35 new `sec-dn-uniform` |
| 7 (CAF spinor (dn2) violates (cuR) for odd $n$; Appendix D candidate) | **REFUTED** | CAF's (dn2) with their own block satisfies (cuR) with $S_{ss}$ for $n=4,\dots,7$ to $10^{-15}$ (check 5). The worker's `dn_particle.S_dn` swaps $S_{ss}=\prod\{4k-3\}$ and $S_{ss'}=\prod\{4k-1\}$ for odd $n$. The book's S-matrix has $S_{ss}=\prod\{4k-3\}$ for all $n$, and for odd $n$ particle $n-1$ is the antiparticle of $n-2$ ($S_{33}(i\pi-\theta)=S_{34}(\theta)$). Gemini's "confirmation" ran the worker's script, so it was not independent. | Nothing in Appendix D. The odd-$n$ "spinor caveat" is dropped; open problem 2 just lists spinor amplitudes as open. |
| 8 (classical identification with $C\equiv+1$ via $d_4$, $d_5$ kinks) | **Not integrated**: no independent check (Gemini checked the reductions and the boundary condition, but the linearization numbers come from the worker's script only) | | `sec-dn-uniform` says the identification is open and expected by analogy |
| 9 (FO match) | **Not integrated** | | — |
| 10 (corrections to the brief) | Accepted for the brief; the folding remark is in open problem 2 | | open problem 2 |
| 11 (step-4 table) | **Not integrated** (literature map, CAF coverage unverified) | | Open problem 2 rewritten more briefly |

## Consolidation A6 (branch half-line-consolidation, 2026-10-13)

Agent report `briefs/halfline_consolidation_A6_report.md`, scripts in `briefs/halfline_consolidation_A6_scripts/`. The coordinator re-derived the key claims with its own code (scratch `a6v/`: full-space boundary-value solution of the vacuum without the symmetry ansatz, own multichannel linearization, own one-loop expansion against the Kim–Cho TeX, the proposed scalar-factor formula typed from the text, the sign rule from the corner entry of the book's intertwiner).

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 2 (crossing factor 1), 5 (singlet coefficient 1) | **Now computed** | `dn_uniform.py` check 6 (agent's; read and rerun) | hedge in `prp-dn-scalar` removed |
| 8 (classical identification with $C\equiv+1$) | **Verified, extended to $d_6$ and to every particle** | Coordinator: BVP vacuum in the full field space gives $e^{f(0)}=9/8$, $4-2\sqrt2$ and the closed-form energies to $10^{-13}$; $d_6$ energy $-0.1032515$; linearized $K_a=\prod_bS_{ab}(\theta,1)$ for all particles to $10^{-11}$, off-diagonal $<10^{-11}$; $K_1$ = $B\to0$ limit of `eq-dn-uniform-amplitude`. Convexity proof checked by hand, including the $d_4$, $d_5$ reductions and $f(0)$. Pole positions $p=\frac i2 m_2$ checked from the S-matrix blocks. No-bound-state claim: agent's check 8 only. | `prp-dn-uniform-vacuum` [Theorem], `prp-dn-uniform-classical` [Numerical] |
| $d_n$ $s=\pm1$: $\varsigma=(-1)^ns$, $U_\varsigma$, $\Psi_\varsigma$, $c_\varsigma$ | **Verified** | Coordinator: corner entry of the book intertwiner vanishes at $y^h=\varsigma$; $K(-\theta)K(\theta)=U_\varsigma$ to $10^{-13}$ ($n=4,5,6$). $\Psi$, $c$: agent's checks 2–3 (matrix level, read and rerun); $c_\pm$ match the book's $a_n$ ones. | `prp-dn-nondiagonal-scalar` (i), (ii), (iv) |
| scalar factors $A_\varsigma=A\rho_\varsigma g_\varsigma$ | **Verified** | Coordinator: formula as typed in the book text, both functional equations to $10^{-27}$, $n=4,5,6$ | `eq-dn-nondiagonal-scalar` |
| breather closed forms, $\Xi$ | **Verified** | Agent's check 5 (projection, $10^{-26}$); coordinator: $\Xi$ unitary and crossing-symmetric, dressed breather formula, $10^{-30}$ | `eq-dn-nondiagonal-breather` |
| continuation: $K^N_1=K_1(2-B)$, $K^{(-)}_1=K_1/K_1^{\rm cl}$; axioms, bootstrap | **Verified** (relations by hand from the block identities $(2h-x)=1/(x)$; axioms and bootstrap: agent's checks 6, 9) | | `eq-dn-neumann-amplitude` |
| Kim–Cho one loop | **Verified** | Coordinator read hep-th/9505138 TeX (eq. numbers (29), (40) by counting environments; arXiv page for the journal ref); own $O(\beta^2)$: $K^N_1$, $K^N_2$ agree to $10^{-17}$; $K^{(-)}_1$ only at $\theta+i\pi$; Kim–Cho's own self-dual $K_1$ and Kim's $d_5$–$d_7$ conjecture agree at one loop, differ by 8–16% at $B=0.6$ | `sec-dn-uniform`; Appendix D self-dual item extended; `kim_cho1996` added |
| book cited [@kim1995] for the $d_4$ one-loop amplitude | **Book error, corrected** | Kim 1995 has only the self-dual $d_n$ conjectures (fig. 7) | open problem 2 |
| $\Xi$ ambiguity ($x=k\omega/2$ vs $\frac12+k\omega/2$) | **Not integrated** (not rechecked); the text says only that $\Xi$ is fixed by one loop and not derived | | — |
| semiclassics: both $\varsigma$ give the Neumann time delay; $a_n$ $A^{(-)}$ too | **Partly integrated**: $d_4$ statement only (agent's check 8); the $a_n$ claim was not rechecked and is not in the book | | `sec-dn-uniform` |
| ch. 33 `prp-dn-r` nesting | **Verified** (markup error from the round-1 integration) | read the source | fixed |
| Task 3, $c_2$ | **Accepted, not integrated**: the folding argument is confirmed; the classical $c_2$ check only confirms the normalization in `sec-images-cn` | | — |

