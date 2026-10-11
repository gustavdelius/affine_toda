# Referee report: chapters 32–33 (boundary quantum groups; soliton reflection matrices)

Report of the referee pass of 2026-10-11 on `referee_half_line_32_33.md`. Line numbers refer to the files as of commit 64ed4af. The referee's scratch scripts were in the session scratchpad (`referee_32_33/`); the ones used for integrated items were copied into `textbook_code/` during integration. The integration status is recorded at the end.

**Summary of real errors found**
- (a) The DM "(4.27) prints Δ(Q) twice" correction is wrong. The PDF has the bar; pdftotext dropped it.
- (b) The mechanism in the proof of `prp-dm-first-order` is wrong. DM's c is ∝ λ^{1/2}, not ∝ λ. The conclusion survives, and for sine-Gordon it can be derived exactly.
- (c) Delius–George (16)–(17) as printed is not an intertwiner.
- (d) `devega1994` is cited for a_{n−1} results it does not contain.
- (e) Doikou did not introduce SP/SNP names.
- (f) DM (4.21) uses a different q from their own (4.15).
- (g) There is a cross-chapter node-label inconsistency between ε̂_j and C_j (via ch. 29).

**Positive results:** the B-matrix expansion redone, the excited-state β_j, the general-n existence condition up to n=7, and the irreducibility hypothesis tested.

## 1. Neumann versus uniform boundary
**Verdict:** the conclusion is sound and [Conjecture] is the right label for n≥2. The proposed mechanism is **wrong**.

**Evidence**
- **Sine-Gordon (solid).** GZ fix the free boundary at ξ=0, η=π(λ+1)/2 from C-symmetry and an expected pole at iπ/2 (GZ §5b). Al. Zamolodchikov's relation (BPT (3.2), = `cnj-bsg-uv-ir`) reproduces this at μ̃→0, and BPT checked that relation against TCSA.
- **Exact map for sine-Gordon.** Combining `eq-bsg-eta-theta` with `cnj-bsg-uv-ir` gives ε̂_1 = −cos((λ+1)·arccos μ_+)/sin πλ, and ε̂_0 the same with μ_−, where μ_± = (μ̃/μ̃_crit)e^{±iβφ_0/2}. Consequences:
  - Each ε̂ depends on one boundary coupling only.
  - ε̂ = ε_* − (λ+1)/(2 sin(πλ/2))·μ_± + O(μ²). So the form of `prp-dm-first-order` holds exactly for sine-Gordon, with κ = −(λ+1)/(2 sin(πλ/2)) in units of μ̃/μ̃_crit.
  - ε̂=0 at φ_0=0 is at μ̃/μ̃_crit = cos(π/(2(λ+1))), which tends to C=1.
- **Caveat on classical limits.** The map ε̂(μ) oscillates like a Chebyshev function as λ→∞. Within O(β⁴) of C=1 it takes every value, including 0 and ±ε_*. So classical limits alone cannot choose between 0 and ±ε_*; only full amplitudes can.
- **a_2 at O(β²).** G98 §6 says what the book says: K^+ matches Kim's O(β²) Neumann result, and Gandenberger "therefore conjecture[s]" that it is Neumann; K^d matches Perkins–Bowcock's O(β²) result for b=+2 (C=+1). Both comparisons use minimal scalar factors.
- **Further support not cited by the book.** DG (hep-th/9904002, end of §3.4) compare their Neumann particle amplitude (3.20), obtained by duality, with Kim at O(β²) for general n, and find them identical.
- **Robustness of the semiclassical test.** The semiclassics and DG (3.16) are evidence for K=1 ↔ C=+1. The CDD factors (3.24) are ω-independent functions of θ, so they do not affect the leading semiclassical term, and the time-delay test is robust against them.
- **Prior art.** The a_2 identification is Gandenberger's own conjecture (1998). Say that the general-n conjecture extends Gandenberger's a_2 conjecture.

**Why the mechanism is wrong**
- DM (4.3) reads c = sqrt(λγ²(q_i²−1)/2πi) (checked against the PDF, p. 13), so c ∝ λ^{1/2}, not ∝ λ.
- In unnormalized charges the Neumann term is therefore ∝ c ε_* ∝ λ^{1/2}. This is non-analytic in the bulk coupling, so it is produced at no finite order of bulk perturbation theory.
- Dimensions agree: [λ^{1/2}] = [λ_b] = mass^{1−β̂²/2}, so the coefficient of q^T is a·λ^{1/2} + b·λ_b + ….
- With DM's image currents, the O(λ) bulk boundary contributions also cancel identically at x→0 (checked by hand).

**Corrections**
- ch. 32 line 252, replace the statement with: "**[Conjecture]** for $a_n$, $n\ge2$; for sine-Gordon it follows from @cnj-bsg-uv-ir. First-order boundary perturbation theory around the Neumann point determines $\partial\hat\epsilon_j/\partial\lambda_b$ at $\lambda_b=0$, not the value of $\hat\epsilon_j$ there: $\hat\epsilon_j=\epsilon_*+\kappa\,\lambda_b\epsilon_j+O(\lambda_b^2)$, with $\kappa$ the constant of [@delius_mackay2003, eq. (4.12)]. For sine-Gordon, (-@eq-bsg-eta-theta) and @cnj-bsg-uv-ir give exactly $\hat\epsilon=-\cos\big((\lambda+1)\arccos\mu_\pm\big)/\sin\pi\lambda$ with $\mu_\pm=(\tilde\mu/\tilde\mu_{\rm crit})e^{\pm i\beta_{\rm SG}\phi_0/2}$, so $\hat\epsilon=\epsilon_*-\frac{\lambda+1}{2\sin(\pi\lambda/2)}\mu_\pm+O(\mu^2)$."
- ch. 32 line 256, replace the proof's 3rd–4th sentences with: "In Delius and MacKay's normalization the constant is $c=\sqrt{\lambda\gamma^2(q^2-1)/2\pi i}\propto\lambda^{1/2}$ [@delius_mackay2003, eq. (4.3)]. A term $\epsilon_*q^{T_j}$ in the normalized Neumann charge is therefore a term $\propto\lambda^{1/2}$ in the unnormalized one, which has the same dimension as $\lambda_b$. Being non-analytic in the bulk coupling, it is invisible at every finite order of bulk perturbation theory, and the first-order computation in $\lambda_b$, done at $\lambda=0$, cannot see it." Keep the last sentence (the exact results force it).

## 2. DM corrections (ch. 32 notes, lines 297–303)
Script `main/dm_own.py` checks DM's own (4.15)–(4.16) equations in their own conventions.
- **q inverted: correct.** By (3.10) and footnote 3, q_DM = e^{8π²i/β_SG²} = 1/q_phys.
- **Factor 2: correct.** In (4.13), Q_i raises the T_i eigenvalue by 2, while (4.26) says α_i·α_i = 1.
- **"(4.27) prints Δ(Q) twice": WRONG.** The PDF (p. 16) prints Δ(Q̄_i) = Q̄_i⊗1 + q^{T_i}⊗Q̄_i, with an overline; the text extraction lost the bar. **Delete ch. 32 line 299.**
- **(4.19)–(4.20): correct.**
  - By hand and with sympy (`dm_own.py`), the (i,j) and (i+1,j) entries, j∉{i,i+1}, give two equations: ε̂_i(q−1)K_ij + x^{-1}K_{i+1,j} = 0 and xK_ij + ε̂_i(q^{-1}−1)K_{i+1,j} = 0.
  - Their determinant is 1 − ε̂_i²(1−q)(1−q^{-1}). This is the general-n origin of `eq-an-eps-condition`.
  - Suggested rewording of line 300: "their coefficients correspond to taking $q^{T_i}=q^{-1}e_{ii}+q\,e_{i+1,i+1}$, forgetting that $q^{T_i}$ is 1 on the other basis vectors".
- **Condition ε̂²=ε_*², a sign not a phase: correct.** In DM's own conventions, for n=2,3: all ±ε_* sign patterns give a unique K; |ε̂|=1 gives none; a phase times ε_* gives none.
- **"Arbitrary diagonal" contradicts (4.18): correct.** At ε̂=0, K=1 exactly.
- **(4.21) correct after the dictionary: correct, but sharpen.**
  - (4.21) does not solve DM's own corrected equations at the q of their (4.15).
  - It does solve them as printed (no basis reversal) for n=2…7 if the q in (4.21) is taken to be −(−q)^{1/2}, with q that of (4.15) (`main/dm_own2.py`, `main/dm_general_n.py`).
  - Suggested addition at ch. 33 line 139: "The $q$ in their (4.21) is not the $q$ of their representation (4.15): (4.21) solves their own (corrected) equations only after this substitution."

## 3. Other source corrections (ch. 33 notes)
- **DN (3.8): correct (script).** The re-run gives:
  - (3.8)–(3.13) solve (3.6)–(3.7) for j≤2 at real and complex η;
  - they fail at physical η for j≥1;
  - the square root taken over the second product fails;
  - the linear equations have a unique solution, which equals the book's K.
- **DG (3.9) prefactor: correct (by hand).** The two choices differ by (hω/4)_I². With a = hω/4, (a)_I(±μ+a) = ±cot(πμ/hω), so the square satisfies both conditions in (3.8). It is a special CDD factor.
- **DG (3.24) minus sign: correct, but a clarification rather than a correction.** DG print the minus. Without it, the product condition (3.23) acquires a factor (−1)^{n+1}, so it fails exactly for even n. Suggested rewording of ch. 33 line 261: "The overall minus sign that they print is necessary when $n$ is even."
- **Delius–George (22): correct.**
  - The omitted ε̂=0 solution and the torus conjugates are checked by script (PASS for d_4, d_5).
  - The η² = q^{−h} rescaling algebra of line 160 checks by hand. Their generator equals q(e + q^{-3}fk + q^{-1}ε̂k) − ε̂ with e = q^{-1/2}E and f = q^{1/2}F. This needs s² = q, which in the homogeneous gradation is a spectral shift by s^h.

## 4. `thm-k-reflection`
**Verdict:** the proof is complete given its hypotheses. Two hypotheses should be sharpened.
- **(ii).** Schur needs only dim End_B(V_μ(θ1)⊗V_ν(θ2)) = 1, together with the same for (μ̄,ν̄) or invertibility of one side; the target is then isomorphic to the source.
  - Computed numerically (`main/commutant.py`): dim End_B = 1 for a_2 and a_3, all four {V,V̄}² pairs, at ε̂=0, s=±1 and generic ε̂, and for a_1 spin ½⊗½ and ½⊗1.
  - At generic ε̂, End_B = C but no K exists, consistent with "B too large".
  - Whether generic irreducibility is proved in general is unclear. A general proof would come from the universal K-matrix of Balagović–Kolb / Appel–Vlaar, but only if B_ε̂ is a Letzter–Kolb QSP coideal. With s_j≠0 on adjacent a_n nodes (a_ij = −1) that is doubtful; Kolb's admissibility conditions should be checked.
- **(iii).** As stated it is vacuous, since one can always rescale an Ř. For physical normalizations the constant is 1 because C and P make the scalar factors on both sides equal: S0_{ν̄μ̄} = S0_{νμ} and S0_{μν̄} = S0_{νμ̄}.
- **Twisted case with four Ř's.** Domains, codomains and arguments (x1/x2 and x1x2) are consistent; no issue.
- **Suggested text for line 102 (ii):** "(ii) for generic $\theta_1,\theta_2$ the only $B_{\hat\epsilon}$-endomorphisms of $V_\mu(\theta_1)\otimes V_\nu(\theta_2)$ and of $V_{\bar\mu}(-\theta_1)\otimes V_{\bar\nu}(-\theta_2)$ are scalars (checked numerically for $a_1,a_2,a_3$)".

## 5. `eq-coideal-generators` (c_j = q_j^{-1}) and ε̂_j ↔ C_j
**c_j: correct.**
- Parity maps J ↔ J̄, and DM normalize Q and Q̄ with the same c, so the Q/Q̄ ratio is fixed.
- Write Q̄ = κ f k. In book conventions e(fk) − q^{-2}(fk)e = (k²−1)/(q−q^{-1}). Matching DM (4.26)'s form (q^{2T}−1)/(q_i²−1), after the T→T/2 fix, gives κ = q^{-1}.
- The parity argument "k^{1/2}fk^{1/2}" (line 65) is heuristic: parity is not an algebra automorphism of this form, and it does not square to 1. Better to cite (4.26).
- Independent confirmation: the B=L̄KL expansion of item 7 returns exactly q^{-1}f_jk_j.

**ε̂_j ↔ C_j: unclear, with a real inconsistency.**
- By ch. 29's own formulas (`eq-bsg-eta-theta`, `cnj-bsg-uv-ir`, and the C-dictionary at line 29, C_0 ∝ e^{+iβφ0/2}), ε̂_1 is a function of C_0 alone and ε̂_0 a function of C_1 alone.
- Ch. 32 lines 27 and 48 assert ε̂_j ∝ ε_j ∝ C_j, so one of the two places has a node-label or φ_0-sign convention flipped.
- The inconsistency is harmless at φ_0 = 0, where all the special points lie.
- **To settle:** trace Part VI's e_j ↔ J_j ∝ exp(…α_j·ϕ) against the soliton s of BPT, then fix either the sign of φ_0 in `eq-bsg-c-dictionary` or the sign of ϑ in `eq-bsg-eta-theta`.

## 6. General-n claims
**a_n: supported up to n=7.** `main/dm_general_n.py` (DM's corrected equations, which are the book's up to basis order) gives, for n=2…7:
- ε̂=0: a unique solution;
- random ± patterns of ε_*: a unique solution;
- mixtures, generic values and 1.1ε_*: none;
- the corrected (4.21) solves the equations.

Proof sketch for general n: the 2×2 determinant above; (4.18) then forces a constant diagonal, and chaining around the cycle kills mixtures. Suggest relabelling to "[Numerical] n≤7; [Proposition, sketch] in general" and adding the determinant argument.

**d_n: n=4,5 only.** Which classical d_n boundaries carry ε̂=0 and which carry s=±1 remains open.

## 7. B = L̄KL: redone in book conventions; correct, with one correction to the text
Scripts (exact, sympy): `item7/r_leading.py`, `item7/b_first_order_symbolic.py`.
- **Universal R.** For the Part VI coproduct, R = q^{G h⊗h}(1 + Σ_j (q−q^{-1}) f_j⊗e_j + …). The variants with a different Cartan placement or sign fail against the exact Ř (a_1–a_3, two values of q).
- **Consequence.** L and L̄ are power series in 1/y, not in x→0 as in DM.
- **Result.** B/κ = 1 + (q−q^{-1}) y^{-1} Σ_j π_V(f_j)|_{y=1} ⊗ (e_j + q^{-1}f_jk_j + ε̂_jk_j) + O(y^{-2}), exactly, for n=1…6, with K=1 (ε̂=0) and with `eq-an-k-matrix` (ε̂ = (sε_*, ε_*, …, ε_*)).
  - With R_21^{-1} the same b_j appear at O(y).
  - This confirms that `eq-boundary-charges` is exact, and confirms the parity-symmetric coefficient and the ε̂ values.
- **Replace ch. 32 line 218's last two sentences** with: "In book conventions the universal R-matrix is $q^{\sum G\,h\otimes h}(1+\sum_j(q_j-q_j^{-1})f_j\otimes e_j+\dots)$, so $\mathcal B^\mu_\theta$ is expanded in $1/y$: for the $a_n$ vector representation, with $K=1$ or (-@eq-an-k-matrix), $\mathcal B=\kappa\big[1+(q-q^{-1})y^{-1}\sum_j\pi_V(f_j)|_{y=1}\otimes b_j+O(y^{-2})\big]$ with exactly the $b_j$ of (-@eq-coideal-generators) and $\hat\epsilon_j=0$ resp. $\pm\epsilon_*$ (checked symbolically for $n\le6$)."

## 8. Entry-by-entry checks
- **DN (2.11) R^{(1/2,j)}: correct** (helper's run). Checks:
  - YBE holds for j≤2 at real, complex and physical η, and altered versions fail;
  - it is the unique intertwiner of DN's coproduct;
  - it equals the book's Ř^{(1/2,j)} after a reversed-and-diagonal basis change;
  - DN's mixed reflection equation holds with the spin-j K.
  - One FAIL line ("book K vs DN K after the same S") is attributed to applying the R-matrix basis change to a map V(y)→V(1/y). The book's own script already passes that comparison.
- **Delius–George (16)–(17): WRONG as printed.** This comes from a helper-written script whose output the referee read; the referee did not re-derive the formula.
  - For d_4, d_5 at generic and physical q, the printed matrix neither intertwines their coproduct on (12)–(14) nor satisfies YBE, and it lacks the TPG eigenvalues.
  - One consistent correction:
    - (1−q²) → (1−q^{−2});
    - x on the a<b sum instead of a>b;
    - the last sum over all a,b (their (17) defines b_aa, but the sum says a≠b);
    - the last term E_{ā,b}⊗E_{a,b̄} instead of E_{b̄,a}⊗E_{b,ā}.
  - The corrected matrix is their intertwiner and satisfies YBE. It equals the book's Ř with the same q and no gauge, and has eigenvalues 1, ⟨1⟩, ⟨1⟩⟨n−1⟩ with multiplicities n(2n+1)−1, n(2n−1), 1.
  - DG (12)–(14) coincide exactly with the book's `dn_vector`.
  - **Replace ch. 33 line 158** with: "Jimbo's explicit matrix as printed by Delius and George [@delius_george2002, eqs. (16)–(17)] contains misprints: as printed it is not an intertwiner. With $(1-q^2)\to(1-q^{-2})$, the factor $x$ on the $a<b$ terms, the last sum over all $a,b$ and its elementary matrices $E_{\bar a,b}\otimes E_{a,\bar b}$, it equals (-@eq-dn-vector-r) (checked for $n=4,5$)."
- **Gandenberger (3.3)–(3.4) in his conventions: unchecked.** To settle: transcribe, then verify his twisted reflection equation for n=2–4.
- **de Vega–González-Ruiz: quoted correctly, not checked numerically.** Their (29) in hep-th/9211114 is K_a^− = k sinh(ξ+θ)e^θ for a≤l_−, and k sinh(ξ−θ)e^{−θ} otherwise. But it is only in hep-th/9211114: `devega1994` is the XYZ/XXZ/XXX (a_1) paper. **Fix ch. 33 line 194** to cite `[@devega1993]` only.

## 9. `prp-boundary-module` vs BPTT: consistent (new check)
Script `main/excited.py`.
- At the boundary bound-state poles y = e^{−iλν_n}, n=0,1,2, π_y(b_0) and π_y(b_1) have a common eigenvector (residual ~1e-15).
- The resulting β_j are exactly the values that reproduce BPTT's excited-state matrix (s↔s̄, η→η̄ = π(λ+1)−η): β_{1,0} = −cos(η̄ ± iϑ)/sin πλ.
- β is the same for every n; the states differ only by the CDD factors a_n.
- The submodule appears at the reflected rapidity e^{−iλν_n}.
- **Replace ch. 32 line 134's "Which β_j each excited state carries is not checked here."** with: "At the poles $y=e^{-i\lambda\nu_n}$, $V(y)\otimes\mathbb C_{\hat\epsilon}$ has a one-dimensional $B$-submodule, with $\beta_{1,0}=-\cos(\bar\eta\pm i\vartheta)/\sin\pi\lambda$, $\bar\eta=\pi(\lambda+1)-\eta$, the same for all $n$; this is exactly the matrix part of [@bptt2002] (checked numerically)."
- **Open for a_n:** whether arbitrary β_j define one-dimensional B-modules for n≥2.

## 10. `tbl-k-types`
- **The assignment is correct.** Durham boundaries make solitons return as antisolitons, giving V→V̄ and the four-Ř twisted reflection equation. The Delius 1998b boundary gives V→V and the braid form.
- **The terminology is misattributed.** Doikou 2000 uses "soliton preserving" and "soliton non-preserving" in quotes, attributing the former setting to Delius (her ref. 4 = `delius1998b`, whose title already says "soliton-preserving"). The abbreviations SP/SNP do not occur in her paper.
  - **Fix ch. 33 line 194:** "…studied in the spin-chain setting by Doikou, who called it "soliton non-preserving", in contrast to the "soliton-preserving" case of Delius [@doikou2000; @delius1998b]".
  - **Fix table line 188:** relabel the row "Spin-chain name", or drop the abbreviations.

## 11. Semiclassical check
**Correct.** dδ/dθ = M sinhθ Δt and M_a = 2hm_a/β² give (4h/β²)ln(…). The limit ω→∞ is β→0. The leading term is O(ω) and the CDD and other corrections are O(1), so the normalized error is O(1/ω), as expected.

## 12. Exercises
All are well posed and solvable. Checked: exr-sg-neumann, exr-cdd-check (both conditions follow from tan(π/4+x) = 1/tan(π/4−x)), exr-dn-casimir (4n, 4n−4), exr-an-diagonal-crossing (v̄_i⊗v̄_i are top-component weight vectors, eigenvalue 1), exr-torus-gauge (n parameters = Gandenberger's E_i).

Caveats:
- **exr-boundary-unitarity-k:** the hint also needs End_B(V(y)) = C (uniqueness for π_y→π_y).
- **exr-semiclassical:** it asks about "attraction and repulsion" without stating a convention. With the usual convention (delay = repulsion), C=+1 has a positive delay, which clashes with ch. 31's "attractive C=+1 boundary binds". State the convention.

## 13. Bibliography
- All listed keys check out against INSPIRE and Crossref: authors, title, journal, volume, pages, year, DOI and eprint.
- `kim1995` (J. Phys. A 29 (1996) 2163) is correctly entered.
- `mezincescu_nepomechie1998` (hep-th/9709078, IJMPA 13 (1998) 2747) supports "first written down"; DM themselves say "first written down in [24]".
- **Citation-support problems:** ch. 32 line 256 (DM (4.3) does not say c ∝ λ); ch. 32 line 299 (delete); ch. 33 line 194 (`devega1994`, and the Doikou naming).

## 14. Prior art and later work
**Unclear/unchecked.** The INSPIRE citation lists were downloaded but not analysed. To settle:
- Read Kolb 1207.6036 §5, the admissibility conditions on s_i. For the Chevalley involution with adjacent a_ij=−1 these likely force s_i=0, which would mean the s=±1 coideals are not standard QSPs; this bears on items 4 and 9.
- Read Regelskis–Vlaar 1602.08471 on the a_n^{(1)} vector K-matrices.
- Scan the citing papers for the ε_* condition, the Neumann ↔ s=+1 identification, and the d_n ε̂=0 solution.

## Referee scripts (outputs summarized)
- `main/dm_own.py`: DM's equations in their own conventions, n=2,3. Gives 9 PASS (ε̂=0 → K=1; ±ε_* unique; |ε̂|=1 none; phase none). The printed (4.21) has residual ~1 at DM's own q.
- `main/dm_own2.py`: (4.21) matches the null vector to about 5e-16 with q_{4.21} = −(−q)^{1/2}.
- `main/dm_general_n.py`: n=2…7 all True, residuals ~1e-15; PASS.
- `main/commutant.py`: dim End_B = 1 in all 36 cases.
- `main/excited.py`: common eigenvector at y = e^{−iλν_n}, n=0,1,2; β = 0.7607 ± 0.3826i, equal to the BPTT prediction.
- Helpers: `item7/r_leading.py`, `item7/b_first_order_symbolic.py` (PASS); `item8a`, `item8b/check_A.py` (Delius–George (16)–(17)).

## Integration (2026-10-11)

**Integrated** (ch. 32 = `32-boundary-quantum-groups.qmd`, ch. 33 = `33-soliton-reflection.qmd`):
- Item 1: `prp-dm-first-order` restated as [Conjecture] for $a_n$, with the exact sine-Gordon map and its first-order expansion. The proof now uses $c\propto g^{1/2}$ from DM (4.3), writing $g$ for the bulk coupling to avoid a clash with sine-Gordon's $\lambda$. Gandenberger's $a_2$ conjecture and DG's general-$n$ duality comparison with Kim were added to the evidence for `cnj-neumann-coideal`.
- Item 2: the "(4.27) prints $\Delta(Q)$ twice" note was deleted; the (4.19)–(4.20) note was reworded; the sharpened (4.21) remark was added in ch. 33.
- Item 3: DG (3.24) reworded as a clarification.
- Item 4: `thm-k-reflection` (ii) is now stated as dim End$_B=1$, the proof uses invertibility plus (ii), and a remark on (iii) was added.
- Item 5 (c_j part only): the heuristic parity argument was replaced by the DM (4.26) argument, with the $\mathcal B$-expansion as independent confirmation.
- Item 6: `prp-an-k-existence` relabelled [Numerical] $n\le7$; [Proposition, sketch] in general. A proof sketch with the $2\times2$ determinant was added; the mixed-pattern step is stated as not written out.
- Item 7: the $\mathcal B=\bar LKL$ paragraph was replaced (universal R leading term, $1/y$ expansion, $n\le6$).
- Item 8: Delius–George (16)–(17) correction. This was confirmed independently: (16)–(18), (2) and (7) were read from the PDF/text of math/0208043 (the $a'$ of (18) included), the referee's `check_A.py` was rerun (72 PASS, 0 FAIL), and a compact port was added. Also added: with only $(1-q^2)\to(1-q^{-2})$ and the full last sum, the printed matrix is $P\check RP$.
- Item 9: the excited-state $\beta_j$ text in ch. 32.
- Item 10: the de Vega–González-Ruiz $a_{n-1}$ statement now cites `devega1993` only; the Doikou naming was fixed; the table row was relabelled "spin-chain name".
- Item 12: the `exr-boundary-unitarity-k` hint now asks for End$_B(V(y))=\mathbb C$; `exr-semiclassical` now states the convention (positive delay = repulsion) and its relation to ch. 31's "attractive $C=+1$". Ch. 31 is unchanged.
- From the 28–29 report, item 1c: the GZ (5.9) note in ch. 32 now says all four operators carry $-\theta$.

**New checks:**
- `textbook_code/boundary_coideal.py`:
  - the $a_n$ loop now runs over $n=2,\dots,7$;
  - §5: $2\times2$ systems (sympy, $a_3$); DM's corrected (4.15)–(4.16) are solved by (4.21) only with $q_{(4.21)}=-(-q)^{1/2}$, for $n\le7$;
  - §6: dim End$_B=1$, 36 cases;
  - §7: one-dimensional $B$-submodules at $y=e^{-i\lambda\nu_n}$, $n=0,1,2$, at $\lambda=4.3$, $\eta=8.1$ (all three $\nu_n$ in the strip). They give the BPTT $\beta_j$; the complementary values appear at $e^{+i\lambda\nu_n}$;
  - §8: universal-R leading term against the exact $\check R$ ($a_1$ spin $\frac12\otimes1$, $a_2$ all pairs, physical $q$ only), and the $\mathcal B$ first-order coefficients for $n\le6$ at the physical and one generic $q$. The referee's version was symbolic in $q$ and also checked $R_{21}^{-1}$; the port is numerical and omits $R_{21}^{-1}$;
  - §9: the sine-Gordon $\hat\epsilon(\mu_\pm)$ map.
  
  Result: 59 PASS, about 30–50 s.
- `textbook_code/reflection_matrices.py` §7: DG (16)–(17) for $d_4,d_5$ at the physical $q$ only (the referee also used a generic $q$). Result: 44 PASS, about 1–1.5 min.

**Not integrated:**
- **Item 5, the $\hat\epsilon_j\leftrightarrow C_j$ node-label inconsistency.** By ch. 29's formulas, $\hat\epsilon_1$ depends on $C_0$ alone; ch. 32 asserts $\hat\epsilon_j\propto\epsilon_j\propto C_j$. This is still open. The new sine-Gordon formula in `prp-dm-first-order` is written per node with $\mu_\pm$ and deliberately avoids assigning $\pm$ to $j$.
- Gandenberger (3.3)–(3.4) in his conventions (unchecked).
- The open question whether arbitrary $\beta_j$ give one-dimensional $B$-modules for $a_n$, $n\ge2$.
- Item 14 (prior art, Kolb admissibility, whether the $s=\pm1$ coideals are QSPs).
- The classical identification of the $d_n$ solutions.
- Minor points not on the integration list:
  - the Chebyshev caveat on classical limits (item 1);
  - the CDD robustness of the semiclassical test (item 1);
  - DN (2.11) $R^{(1/2,j)}$, verified by the referee's helper; ch. 33 still says "not checked here entry by entry".
- `textbook_code/README.md` was not updated; its runtime and ids columns for the two scripts are now stale.
