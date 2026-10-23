# Chapter 42 drafting report: "Defects meet boundaries"

Files written (nothing else touched; no git commands, no quarto render):
- `part8-defects/42-defects-and-boundaries.qmd` (241 lines; skeleton replaced; all 73 cross-references and all citations resolve)
- `code/defect_doubling.py` (new; imports only `boundary_common` and `defect_common`)
- `briefs/draft_ch42_checkpoint.md`, `briefs/draft_ch42_report.md`
- Scratch scripts: session scratchpad `ch42/` (s1–s16, dd_engine.py); not needed for the book.

## 1. Chapter summary

**Intro.** The doubling trick turns a defect into a boundary; the chapter asks which symmetry the doubled boundary preserves, mirroring ch. 34. Three answers: invertible lines are diagonal-type QSPs; a jump defect keeps two Borel halves (one untwisted, one twisted by charge conjugation) but no twisted diagonal; a defect fused onto a boundary is a new boundary (exactly GZ for sine-Gordon solitons, scalar included).

**`sec-doubled-theory`.** The reflected copy carries $V'(u)=V(-u)$ with the action through $\Delta^{\rm op}$; doubled bulk algebra $A=U\otimes U^{\rm cop}$ with $\Delta_A$ (`eq-dbd-coproduct`); symmetries of the doubled boundary are left coideals of $A$ acting on $M\otimes\mathcal W$ by $(\sigma\otimes\chi)\Delta_A$. Via $\omega_q$ (Hopf iso $U^{\rm cop}\to U$) the reflected copy is the conjugate multiplet: parity maps solitons to antisolitons.
- `prp-reflected-copy` **[Established; numerically checked]**: (1) $P\check RP$ intertwines $\Delta^{\rm op}$; (2) doubled $\check R$ intertwines $\Delta_A$ of $U\otimes U^{\rm cop}$, not of $U\otimes U$; (3) $V(1/y)\circ\omega_q=\bar V(y)$, $P\check RP\propto\check R_{\bar V\bar V}$; (4) RE ⇔ triangle relations (cites `prp-doubling`, checked in `defect_axioms.py` d).

**`sec-diagonal-qsp`.** $\Delta(U)$ and $(1\otimes\psi)\Delta(U)$ ($\psi$ Hopf automorphism) are left coideals (`eq-dbd-diagonal`). Invertible lines $\mathbb T=g$ ↔ twisted diagonals (twist line: torus $t_1=1/t_0=e^{2i\alpha}$; charge conjugation: node swap $\tau$); particle-as-defect ↔ untwisted $\Delta(U)$. $\omega_q\otimes{\rm id}:U\otimes U^{\rm cop}\to U\otimes U$ is an algebra iso and coalgebra anti-iso mapping $\Delta(U)$ onto Bao–Wang's right coideal $\varphi(U)=(\omega\otimes1)\Delta(U)$. $q=1$: fixed points of $(x,y)\mapsto(\omega(y),\omega(x))$, second kind; Satake diagram = two copies, all white, $\tau$ exchanging (by hand). "What is proved, and where" paragraph (see §5).
- `prp-diagonal-coideal` **[Theorem (1)–(2); Established (3); Numerical (4)]**: (4) in a 21-monomial search space, the largest left coideal inside the commutant of the doubled K is exactly $(1\otimes\psi)\Delta(U)\cap X$ (5-dim), for transparent (sG spin-1/2, spin-1 probe, $a_2$), twist line, C line. Commutant itself is 17-dim (one particle), 7 or 12 (with two particles) — the coideal requirement does the work; parallel to Appel–Vlaar Thm 4.1.3.

**`sec-jump-coideal`.** Conjecture block kept verbatim-style, labelled **[Conjecture; refuted as stated by the test in `code/defect_doubling.py`]**, followed by the pre-registered expectation and the test description.
- `prp-diagonal-borel` **[Theorem (1), (4); Numerical (2), (3), (5)]**: (1) diagonal Borel $\Delta(U_q(\mathfrak b^\mp))$ is a symmetry with $\mathcal W$ the Weyl module; (2) NEW: $(1\otimes\tau)\Delta(U_q(\mathfrak b^\pm))$ is also a symmetry, on a second module structure $\tilde\rho$ of $\mathcal W$ (`eq-dbd-second-module`, explicit); (3) no twisted diagonal; (4) no coideal containing both acts on any $\mathcal W$; (5) limits $\theta_d\to\pm\infty$ restore the untwisted / $\tau$-twisted diagonal.
- Reading: the jump defect keeps the $\mathfrak b^-$ half of the twist line's diagonal and the $\mathfrak b^+$ half of the charge-conjugation line's diagonal (quantum form of `eq-dq-massless-defect`); $\tilde\rho$ is the $\mathfrak b^+$ symmetry of the same defect read as type I$_+$ after $\psi\to-\psi$ (`prp-two-borels` (3)) — labelled as an interpretation.
- "Beyond sine-Gordon": expectation for $a_r$ with the rotation $\pi$ in place of $\tau$ (untested; see obligations).
- RTT inside the doubled RE algebra: off-diagonal blocks $P\mathbb T=\mathbb L(\theta)$, $\bar{\mathbb T}P=\mathbb L(-\theta)^{-1}$; VV→V'V' block = RTT, V'V'→VV = RTT for $\mathbb L^{-1}$, mixed blocks follow; diagonal Borel coideal is again a Hopf algebra.
- $\mathcal W$ as a dynamical boundary module (cf. `sec-dynamical-boundary`; contrast Baseilhac–Koizumi's $q$-Onsager boundary module).
- New table `tbl-dbd-symmetries` (symmetries of doubled boundaries of sine-Gordon lines, with status).

**`sec-defect-boundary-fusion`.** `eq-dbd-fused-k`: $K_{\rm f}=\bar{\mathbb T}(1_{\mathcal W}\otimes K)\mathbb T$ (Sklyanin; BG (19); cz2012 (2.9)); Lee–Yang instance.
- `prp-defect-fusion-k` **[Theorem (1); Numerical (2)–(3); by hand (4)]**: (1) RE + unitarity; (2) sinh-Gordon identity $\bar{\mathbb T}K_D(E)\mathbb T=K_{\rm shG}(E,F=2i\theta_d/\pi)$; (3) sine-Gordon solitons: type-I$_-$ + Dirichlet $=$ GZ (`eq-bsg-k-full`) with $(\eta,\vartheta)=(\eta_D,\lambda\theta_d)$ sector by sector after a diagonal gauge, **scalar factor included** (no CDD); (4) with `cnj-bsg-uv-ir`'s classical limit and the classical fusion: $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$ = `cnj-defect-parameter` (sine-Gordon form), now from the soliton sector.
- Classical version `eq-dbd-dressed-k`: $K_\phi(\lambda)=\sigma(L_d(1/\lambda)^{-1})K_\psi(\lambda)L_d(\lambda)^{-1}$, a dynamical K-matrix (inverse defect matrices, as for the defect operator).
- `prp-fused-bdg` **[Theorem, by hand; numerically checked]**: (1) dressing satisfies the dynamical compatibility $\partial_tK=Ka_t(\lambda)+\sigma(a_t(1/\lambda))K$ for every $a_r$ boundary of `tbl-durham` incl. Neumann (checked type I$_-$, $a_1$–$a_3$) ⇒ new integrable dynamical $a_n$ boundaries (addresses Part VII problem 9; involution/quantization open); (2) sG defect potential = BDG $H_{\rm osc}$ at $(\hat p,\hat q)=(\beta\psi/2-i\theta_d,\beta\psi/2+i\theta_d)$; Dirichlet fusion = oscillator frozen at a complex point = static GZ; (3) Neumann fusion is NOT the BDG boundary (one variable conjugate to $\phi(0)$; extra $\psi_t$ in BC; BDG flow does not preserve $\hat p-\hat q$). **Answer to the brief's question: the fused system is not the BDG dynamical boundary.**

**`sec-tensor-k`.** Appel–Vlaar tensor K-matrices: $\mathbb K=\mathcal R^\psi_{21}(1\otimes\mathcal K)\mathcal R$ (Thm 5.6.1, (5.5)–(5.8)); Thm 6.4.2; Thm 6.5.2 trigonometric tensor K; Rem 6.5.3(a) (Ito–Terwilliger, $q$-Onsager) and (b). Physical reading: fusion of a particle-defect onto a boundary. One numerical check (sG coideal, $\mathcal W=V(\theta_d)$: intertwiner, 1-dim intertwiner space, RE). Jump defects lie outside AV's categories; link to AV §1.6.2 (Borel legs, boundary Q-operators).

**`sec-noninvertible-boundaries`.** Lines fused onto boundaries (CFT, Lee–Yang, Runkel 2010; $g^{-1}Kg$ for invertible lines); lines ending on boundaries (SY's weakly symmetric boundaries, Hom$(L\otimes B,B)\ne0$; example $\sum_c|L_c\rangle$); SY's modified boundary crossing `eq-dbd-sy-crossing` (their (3.29)), its origin (3.22)–(3.28); WT identities (4.2)–(4.3); $p\ge5$ (4.12), $p=4$ (4.26)–(4.27); minimal solutions (4.29), (4.35). Relation to ch. 36: the boundary modification is removed by the same rapidity-dependent kink renormalization as CCK's bulk one (by hand, `exr-dbd-boundary-crossing`). Relation to Part VIII defects; doubled boundaries of non-invertible lines open.

**Exercises (8):** `exr-dbd-cop`, `exr-dbd-diagonal-coideal`, `exr-dbd-invertible`, `exr-dbd-mixed-twist`, `exr-dbd-e-one`, `exr-dbd-fusion-gz`, `exr-dbd-classical-dressing`, `exr-dbd-boundary-crossing`.

**Notes** with history, a "stated here" list, and "Corrections and clarifications to the sources made in this chapter:" (§5 below).

**New ids** (all prefixed, no clashes): `eq-dbd-coproduct`, `eq-dbd-diagonal`, `eq-dbd-second-module`, `eq-dbd-fused-k`, `eq-dbd-dressed-k`, `eq-dbd-sy-crossing`, `tbl-dbd-symmetries`, `exr-dbd-cop`, `exr-dbd-diagonal-coideal`, `exr-dbd-invertible`, `exr-dbd-mixed-twist`, `exr-dbd-e-one`, `exr-dbd-fusion-gz`, `exr-dbd-classical-dressing`, `exr-dbd-boundary-crossing`. All 15 skeleton ids kept.

## 2. Script output (`code/defect_doubling.py`, 17 PASS, 0 FAIL, 30 s)

```
PASS doubled R-check (Rc on V V, P Rc P on V'V', flip on mixed pairs) intertwines Delta_A of U (x) U^cop for sine-Gordon and a_2 [8e-16]; with U (x) U instead it fails [1.66]
PASS reflected copy = conjugate multiplet: V(1/y) o omega_q = Vbar(y) exactly (a_2 [2e-16]; sine-Gordon: C V(y) C [7e-16]), and P Rc(1/y2, 1/y1) P = Rc on Vbar(y1) (x) Vbar(y2) up to a scalar [6e-16]
PASS largest left coideal inside the commutant of the doubled K (X: 21 monomials of degree <= 1 per node pair): transparent defect -> Delta(U) cap X (sine-Gordon spin 1/2 and spin 1 probes, a_2); twist line -> (1 (x) t)Delta(U) cap X, t_1 = 1/t_0 = e^{2 i alpha}; charge conjugation -> (1 (x) tau)Delta(U) cap X [max deviation 2e-15; commutant dims (one particle, with two particles, coideal) [(17, 7, 5), (17, 12, 5)]]
PASS prp-diagonal-borel: the doubled K of the type-I_-+ defect commutes with the diagonal Borel Delta(U_q(b^-+)) acting on M (x) W through Delta_A, W the Weyl module (windows N = 6, 8; s = 0, 0.3) [2e-15]
PASS cnj-twisted-diagonal as stated is refuted: for no twisted diagonal (1 (x) psi)Delta(U) (psi = torus or tau o torus) do both its e- and f-generators commute with the doubled jump-defect K for any theta-independent action on W; minimal residuals over windows N = 5, 7, 9: I- e/id 0.95-0.96, I- f/tau 0.96-0.96, I+ e/tau 0.93-0.93, I+ f/id 0.96-0.97
PASS what holds instead (same test): type I_- keeps f/id (diagonal b^-) and e/tau (tau-twisted b^+); type I_+ keeps e/id and f/tau [6e-15, all windows]
PASS second module rho~ of U_q(b^+) on W (type I_-): rho~(e_1)|n> = a q^-2n |n+1>, rho~(k_1)|n> = q^(-2n-2-s)|n+2>, rho~(e_0)|n> = b q^2n |n-1>, rho~(k_0)|n> = q^(2n-2+s)|n-2>, a = i y_d q^-s/(q^2-1), b = i y_d q^s/(q^2-1): k_0 k_1 = 1, k e k^-1 = q^(+-2) e, e_1 e_0 = e_0 e_1 central [1e-15]; and T intertwines (pi (x) rho~)Delta(x) with (rho~ (x) pi o tau)Delta(x) for x = e_j, k_j (s = 0.3) [2e-15]
PASS the two Borel symmetries of the I_- doubled boundary do not combine: rho(k_1) rho~(k_1) = q^4 rho~(k_1) rho(k_1) [1e-15; commutator 1.46], while k_1 (x) k_1 and k_1 (x) k_0 commute; any coideal containing both contains e_j (x) 1, which is not a symmetry [residual 0.96]
PASS limits (N = 6): as theta_d -> +inf the untwisted e-generator is restored (residual 4.4e-01, 6.3e-02, 8.0e-03 at theta_d = 1.5, 3, 4.5: the twist line), as theta_d -> -inf the tau-twisted f-generator (2.7e-01, 3.6e-02, 4.6e-03: the charge-conjugation line); decay ~ e^{-lambda |theta_d|}
PASS control: a particle frozen as a defect (W = V(theta_d), T = Rc) is invariant under the full diagonal Delta(U), W acting by pi_W [1e-15]
PASS defect + Dirichlet: the charge-sector blocks of K_f = Tbar (1 (x) K_D) T have label-independent gauge invariants (Corrigan-Zambon (2.18)) [2e-15] and equal Ghoshal-Zamolodchikov with (eta, vartheta) = (eta_D, lambda theta_d) [9e-16]; vartheta = -lambda theta_d or theta_d fail [0.03] (windows N = 5, 7, 9)
PASS the infinite-dimensional K_f solves the reflection equation on V (x) V (x) W, interior labels, windows N = 5, 7, 9 [8e-16, 7e-16, 1e-15]; a label-dependent perturbation fails [0.02]
PASS with the scalar factors (Phi of eq-sgd-scalar for T and Tbar; Dirichlet K_D = R_0 sigma(eta,u)/(2 cos eta) diag(y e^{i eta} + e^{-i eta}/y, ...), the vartheta -> inf limit of eq-bsg-k-full since sigma(i vartheta,u) -> 1 [2e-08 at vartheta = 24]) the diagonal entries of Tbar K_D T equal those of eq-bsg-k-full at (eta_D, lambda theta_d) [3e-09; the accuracy of the Richardson-summed Phi at complex rapidity is ~1e-8]
PASS sinh-Gordon (Bajnok-Simon): Tbar_shG K_D(E) T_shG = K_shG(E, F = 2 i theta_d/pi) of eq-shg-reflection [2e-16]
PASS classical fusion: K_phi(lam) = [L_d(1/lam)^-1]^T K_psi(lam) L_d(lam)^-1 satisfies d_t K = K a_t(lam) + a_t(1/lam)^T K on the fused defect + boundary system (a_1 with GZ C = (0.7, -1.3) and Neumann; a_2 with C = (1,-1,1) and Neumann; a_3 with C = (1,-1,1,-1)) [9e-10, finite differences]; L_d instead of L_d^-1, or K_psi undressed, fails [1.00]
PASS BDG: the sine-Gordon defect potential (type I_-) equals H_osc of eq-bsg-oscillator at (p^, q^) = (beta_SG psi/2 - i theta_d, beta_SG psi/2 + i theta_d); at psi = phi_D, C_1 = cos q^ = cosh(theta_d - i beta_SG phi_D/2); the BDG flow does not preserve p^ - q^ = -2 i theta_d [d(p^-q^)/dt = 0.34 at a sample point]
PASS tensor K-matrix: for the sine-Gordon coideal (generic eps-hat) and W = V(theta_d), Tbar (1 (x) K) T with T = Rc intertwines Delta(b_j) on V(theta) (x) W [3e-15], the intertwiner space is {1}-dimensional (Appel-Vlaar Rem. 6.5.3(b)), and the reflection equation holds [1e-15]
real 0m30.071s
```
Windows and stability: Weyl windows $|n|\le N$, equations on $|n|\le N-2$. Diagonal Borel: $N=6,8$ (2e-15). Twisted-diagonal tests: $N=5,7,9$; holding generators ≤6e-15 at every $N$; failing ones 0.93–0.97, stable under enlargement. Fusion: $N=5,7,9$ all at 1e-15. No untwisted trace is used anywhere. The only non-machine-precision number is the scalar identity (3e-9), limited by `defect_common.sg_transmission_scalar` at complex rapidity (K=2000 vs 8000 differ by 1–2e-8 there; 4e-14 at real rapidity).

Implementation notes: `boundary_coideal`, `qsp_checks`, `boundary_classical`, `boundary_sine_gordon` and `defect_axioms` run their checks at import (no `__main__` guard), so they are not imported. Small helpers are reimplemented inside the script: BPTT `R_0`, `sigma` (~25 lines, as in `boundary_sine_gordon.py`), a solver for `eq-boundary-k-equation` (~10 lines), the doubled R-check (cf. `defect_axioms.doubled_S`). Item (a) "doubled RE ⇔ triangle relations" is cited to `defect_axioms.py` check d, not duplicated.

## 3. Verification record

**By script (`defect_doubling.py`):** checks a1–a2 (`prp-reflected-copy` (2),(3)); b (`prp-diagonal-coideal` (4)); c1 (`prp-diagonal-borel` (1)); c2–c3 (`cnj-twisted-diagonal` refutation; `prp-diagonal-borel` (2),(3)); c4 ($\tilde\rho$ relations and exact intertwining, `eq-dbd-second-module`); c5 (`prp-diagonal-borel` (4): $q^4$ obstruction and $e_j\otimes1$ failure); c6 (limits, (5)); c7 (particle-as-defect control); d1–d4 (`prp-defect-fusion-k` (2),(3), incl. scalar identity and RE on windows); e1 (`prp-fused-bdg` (1)); e2 (`prp-fused-bdg` (2), and non-preservation in (3)); f (`sec-tensor-k` check).

**By hand (this agent):** proofs of `prp-reflected-copy` (1)–(3); `prp-diagonal-coideal` (1)–(2); the map $\omega_q\otimes{\rm id}$ onto Bao–Wang's $\varphi(U)$ and its anti-coalgebra property; the $q=1$ involution and second kind; Satake diagram of the diagonal type (stated "by hand"); `prp-diagonal-borel` (1) and the derivation of $\tilde\rho$ (matrix-element recursion, `exr-dbd-mixed-twist`) and (4) (the abstract argument that $e_1\otimes1$ lies in any coideal containing both, and its failure at order $y^2$); `prp-defect-fusion-k` (1) (Sklyanin-dressing argument) and (4) (classical limit of `cnj-bsg-uv-ir`); `prp-fused-bdg` (1) proof, (2) expansion, (3) comparison with BDG; the RTT-in-RE block discussion; the SY renormalization (`exr-dbd-boundary-crossing`: exponent bookkeeping done by hand twice; constant part gives exactly $\sqrt{d_d/d_b}$, $\theta$-parts cancel).

**Unchecked / quoted only:**
- Bao–Wang Rem. 4.10(2) and Wang's survey (content quoted, read in source text); Balagović–Kolb, Kolb, Appel–Vlaar 2022 indecomposability (read in source text, not re-derived).
- Appel–Vlaar 2024 Thms 4.1.3, 5.6.1, 6.4.2, 6.5.2, Rem. 6.5.3, §1.6.2 (read; quoted). Ito–Terwilliger quoted via AV.
- Shimamori–Yamaguchi results (read; quoted; not re-derived numerically). Their bulk factor (2.32) used as printed.
- Interpretation of $\tilde\rho$ via $\psi\to-\psi$ / BCZ two labels (heuristic, flagged in text); "quantum form of `eq-dq-massless-defect`" (interpretation).
- Boundary crossing-unitarity of $K_{\rm f}$ (only RE and unitarity claimed).
- Type I$_+$ classical dressing (only I$_-$ checked; the proof covers both); involution of the fused $a_n$ dynamical boundaries.
- $a_r$ expectation ("Beyond sine-Gordon": $\Delta(U_q(\mathfrak b^-))$ + $(1\otimes\pi)\Delta(U_q(\mathfrak b^+))$) — untested.
- CFT remarks: folded trivial defect = permutation-type boundary state preserving the diagonal chiral algebra; fusion of Verlinde lines onto Cardy boundaries by fusion coefficients (uncited general facts).
- Baseilhac–Koizumi boundary module statement (as in ch. 34).

**Sources.** `2504.08375` (Shimamori–Yamaguchi) and `2402.16676` (Appel–Vlaar tensor K) were not in `~/referee-copy/papers`; WebFetch saved their arXiv PDFs in the session tool-results directory and I converted them with pdftotext in the scratchpad. Recommend adding both to `~/referee-copy/papers` for the referee.

## 4. Outcome of `cnj-twisted-diagonal`

**Refuted as stated**, and the pre-registered expectation is confirmed, with one addition:
- No twisted diagonal $(1\otimes\psi)\Delta(U)$ ($\psi$ = torus or $\tau\circ$torus; torus factors absorbed by free $\kappa$) is a symmetry of the doubled boundary of either sine-Gordon type: for type I$_-$ the untwisted $e_j$ (0.95–0.96) and the $\tau$-twisted $f_j$ (0.96) fail; for type I$_+$ the $\tau$-twisted $e_j$ (0.93) and untwisted $f_j$ (0.96–0.97) fail — least squares over arbitrary $\theta$-independent weight-graded operators on $\mathcal W$, three rapidities, windows 5, 7, 9.
- The diagonal Borel $\Delta(U_q(\mathfrak b^\mp))$ holds, with $\mathcal W$ the Weyl module as a dynamical boundary module (2e-15) — as expected.
- **New:** the opposite Borel survives twisted by charge conjugation, $(1\otimes\tau)\Delta(U_q(\mathfrak b^\pm))$, on a second module structure $\tilde\rho$ (explicit; $\tilde\rho(k_1)$ shifts the defect label by 2). The two cannot be combined into one coideal acting on any module (any coideal containing both contains $e_1\otimes1$; also $\rho(k_1)$, $\tilde\rho(k_1)$ $q^4$-commute). Interpretation: one Borel half of each of the two invertible lines between which the defect interpolates; $\tilde\rho$ is the I$_+$ reading after $\psi\to-\psi$.
- In the limits $\theta_d\to+\infty$ / $-\infty$ the untwisted / $\tau$-twisted full diagonal is restored, residual decaying like $e^{-\lambda|\theta_d|}$.

## 5. Corrections and clarifications to the sources (in the Notes)

1. **Diagonal type (clarification).** Bao–Wang Rem. 4.10(2) is a finite-type remark ("one checks"); Kolb 2014 and Appel–Vlaar 2022 assume an indecomposable Cartan matrix, so the diagonal type of affine algebras ($\hat g\oplus\hat g$) is not covered as stated. Their right coideal $(\omega\otimes1)\Delta(U)$ = our left coideal $\Delta(U)\subset U\otimes U^{\rm cop}$ under $\omega_q\otimes{\rm id}$. Kolb–Yakimov's "diagonal type" = Nichols algebras of diagonal braiding (unrelated). Balagović–Kolb do not discuss the diagonal type.
2. **Corrigan–Zambon (2.18) (extension/clarification).** The fused type-I matrix is GZ with $(\eta,\vartheta)=(\eta_D,\lambda\theta_d)$ in book conventions, and the identity extends to the scalar factors (no CDD factor). Their $\hat T(\theta)=T(-\theta)^{-1}$ is the book's $\bar{\mathbb T}$.
3. **Shimamori–Yamaguchi (clarification).** Their modified boundary crossing (3.29) is equivalent to GZ's conventional crossing after the same rapidity-dependent kink renormalization that relates the bulk normalizations ((2.32)); WT identities, not crossing, distinguish normalizations (by hand).
4. **Book-internal.** `cnj-twisted-diagonal` refuted as stated; refines (does not contradict) ch. 40 `sec-perturbed-defects` ("no action of the $e_j$ extends $U_q(\mathfrak b^-)$ with the untwisted coproduct"): the $e_j$ survive $\tau$-twisted.

## 6. Requests (edits elsewhere) and obligations

**`tbl-boundary-defect-algebra` (ch. 40), third column "Defect, doubled (ch. 42)" — proposed row texts:**
- Symmetry: "left coideal of $U\otimes U^{\rm cop}$ (-@eq-dbd-coproduct): invertible lines and the transparent defect, a twisted diagonal $(1\otimes\psi)\Delta(U)$ (@prp-diagonal-coideal); type-I$_\mp$ jump defect, $\Delta(U_q(\mathfrak b^\mp))$ and $(1\otimes\tau)\Delta(U_q(\mathfrak b^\pm))$ on two module structures of $\mathcal W$, and no twisted diagonal (@cnj-twisted-diagonal, @prp-diagonal-borel)". Status column of that row: replace "doubled: to be decided" by "doubled: [Numerical] (sine-Gordon), with the module $\tilde\rho$ and the obstruction by hand".
- Module of the line: "$\mathcal W$ as a dynamical boundary module of $\Delta(U_q(\mathfrak b^\mp))$ (@sec-jump-coideal)".
- Matrix: unchanged.
- Linear equation: "intertwiner of the coideal on $(V\oplus V')\otimes\mathcal W$; for $\Delta(U_q(\mathfrak b^\mp))$ it is (-@eq-borel-intertwiner)".
- Quadratic equation: unchanged.
- Quadratic algebra: "reflection-equation algebra of $U\otimes U^{\rm cop}$; its off-diagonal part is the RTT algebra of $\mathbb L=P\mathbb T$ with $\mathbb L^{-1}$ adjoined (@sec-jump-coideal)".
- Universal object: replace "quasi-K-matrix of diagonal type = quasi-R-matrix [@bao_wang2018b]" by "diagonal type: quasi-K-matrix = Lusztig's quasi-R-matrix twisted by $\omega\otimes1$ [@bao_wang2018b] (finite type); for $\hat g$ the universal R-matrix plays this role (@sec-diagonal-qsp)".
- Operator on the circle: "boundary state of the doubled theory = defect operator [@bajnok_george2006]".
- Functional relations: "fusion onto boundaries, $K_{\rm f}=\bar{\mathbb T}(1\otimes K)\mathbb T$ (@prp-defect-fusion-k)".
- Classical limit: "boundary K-matrix dressed by $L_d^{-1}$ (-@eq-dbd-dressed-k), @prp-fused-bdg".
- In the sentence before the table: "where it is decided" → "where it is settled".

**Other pointers (ch. 40, 35, 29/30, 34):**
- ch. 40 `sec-perturbed-defects`, after "(the best candidate leaves a residual of 0.88 ...)": add "; twisted by charge conjugation they do survive, on a second module structure of $\mathcal W$ (@prp-diagonal-borel)".
- ch. 35 problem 9 (`sec-boundary-open-problems`, "Higher-rank dynamical boundaries"): add "Classically, a type-I defect fused onto a boundary of @tbl-durham gives integrable dynamical $a_n^{(1)}$ boundaries with $n$ boundary variables (@prp-fused-bdg); their quantization is open."
- ch. 29 `sec-bsg-dynamical` or ch. 30 `sec-dynamical-boundary`: "A defect fused onto a Dirichlet boundary is the oscillator boundary frozen at a complex point; fused onto Neumann it is a different dynamical boundary (@prp-fused-bdg)."
- ch. 34 `sec-universal-k` or Notes: pointer to `sec-diagonal-qsp` (diagonal type; defects as QSPs).
- ch. 37 `sec-sg-breather-transmission`: optional pointer to `prp-defect-fusion-k` (3)–(4) (soliton version gives the same $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$).

**`code/README.md`:** row for `defect_doubling.py` (ch. 42; 17 checks; ~30 s; imports `boundary_common`, `defect_common`).

**`briefs/part8_defects_status.md`:** ch. 42 row (drafted 2026-10-12, `defect_doubling.py`, 17 PASS (30 s), not refereed) and this verification record.

**Conventions brief, new §10 "Added from chapter 42":**
- Doubled theory: $V'(u)=V(-u)$ with $\Delta^{\rm op}$; $A=U\otimes U^{\rm cop}$, $\Delta_A(g\otimes h)=\sum(g_{(1)}\otimes h_{(2)})\otimes(g_{(2)}\otimes h_{(1)})$; doubled K ordered $(V\otimes\mathcal W)\oplus(V'\otimes\mathcal W)$, blocks $P\mathbb T=\mathbb L(\theta)$ and $\bar{\mathbb T}P=\mathbb L(-\theta)^{-1}$; symmetries are left coideals of $A$ acting by $(\sigma\otimes\chi)\Delta_A$.
- $\omega_q\otimes{\rm id}$ maps $\Delta(U)$ onto Bao–Wang's $(\omega\otimes1)\Delta(U)$; $V(1/y)\circ\omega_q=\bar V(y)$.
- Invertible lines ↔ $(1\otimes\psi)\Delta(U)$; sG twist line $t_1=1/t_0=e^{2i\alpha}$; C ↔ $\tau$ (node swap).
- Jump defect: $\Delta(U_q(\mathfrak b^\mp))$ on $\rho$ and $(1\otimes\tau)\Delta(U_q(\mathfrak b^\pm))$ on $\tilde\rho$ (formulas `eq-dbd-second-module`, sg normalization $A_n=q^{-n}$); not combinable.
- Fusion: $K_{\rm f}=\bar{\mathbb T}(1_{\mathcal W}\otimes K)\mathbb T$, $\bar{\mathbb T}(\theta)=\mathbb T(-\theta)^{-1}$; I$_-$ + Dirichlet($\eta_D$) = GZ($\eta_D$, $\vartheta=\lambda\theta_d$), scalar included.
- Classical dressing $K_\phi=\sigma(L_d(1/\lambda)^{-1})K_\psi L_d(\lambda)^{-1}$, dynamical compatibility $\partial_tK=Ka_t(\lambda)+\sigma(a_t(1/\lambda))K$.
- BDG dictionary $(\hat p,\hat q)=(\beta_{\rm SG}\psi/2-i\theta_d,\beta_{\rm SG}\psi/2+i\theta_d)$.
- Trap: `boundary_coideal`, `qsp_checks`, `boundary_classical`, `boundary_sine_gordon`, `defect_axioms` are not import-safe.

**Bib keys:** none missing. Optional additions (cited only via other papers now): Ito–Terwilliger (augmented tridiagonal algebra, Kyushu J. Math., quoted in AV Rem. 6.5.3(a)); Kolb, "Braided module categories via quantum symmetric pairs" (Proc. LMS 2020); Lusztig, *Introduction to quantum groups* (quasi-R-matrix).

**Obligations for ch. 43 (open problems):**
1. $a_r$ jump defects: test the expectation $\Delta(U_q(\mathfrak b^-))$ + $(1\otimes\pi)\Delta(U_q(\mathfrak b^+))$ ($\pi$ the diagram rotation) with the ch. 41 transmission matrices; identify the second module structure.
2. Type-II defects (generalized oscillator modules): symmetry of the doubled boundary; Corrigan–Zambon (2.19) fused matrices with irremovable charge dependence — which coideal and boundary module.
3. Quantization and involution of the classical dynamical $a_n$ boundaries obtained by fusing a defect onto a Durham/Neumann boundary (Part VII problem 9).
4. Doubled boundaries of non-invertible lines (CCK/SY kinks): their symmetry algebras.
5. Appel–Vlaar's tensor K-matrices with a Borel leg (§1.6.2) and boundary Q-operators, realized by defects in front of boundaries.
6. Mathematical gap: diagonal-type QSPs for affine $\hat g$ (decomposable Cartan matrix) are not covered by the published universal-K constructions.
7. The operator algebra generated by the two Borel symmetries of a jump defect (both present as operators, not as one coideal): structure and meaning.
8. Boundary crossing-unitarity of the fused infinite-dimensional $K_{\rm f}$ (not checked here).
