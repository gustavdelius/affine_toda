# Brief: fixes and gaps from the Part V reading (Chapters 16–19)

Status: done 2026-10-10 (uncommitted, branch `reading-fixes`; see Resolution at the end). Written 2026-10-10 after a reading of Part V by a reader with graduate-level mathematical physics and QFT, who knows Parts I–IV as written. Chapters 18 and 19 are long and technical; this brief records what was checked and what was not.

Files: `part5-semiclassics/part.qmd`, `16-dhn.qmd`, `17-thimbles.qmd`, `18-jordan-chains.qmd`, `19-one-loop-masses.qmd`.

Scope: fix the items below in place. Items A3, A6 and A7 need the literature or a calculation before any change, and the report must say what was checked. Do not edit anything else without checking with the author. The Part V results checked by hand or by reference are in section D. Numerical claims were not rerun in this review.

## A. Errors and likely errors

**A1. Channel form of the transmission-factor formula (`16-dhn.qmd`, Section `sec-dhn-an`, and `19-one-loop-masses.qmd`, Proposition `prp-one-loop-mass`).** Equation `eq-dhn-psi` is stated for a single field with mass $m$, $\Delta M=\frac m2\sum_w\operatorname{ord}_w(X)\psi(\arcsin w)$. Section `sec-dhn-an` then applies it "channel by channel" to $a_n^{(1)}$. In the channel sum the prefactor must be the channel mass $m_b$, not $m$. I checked this for $a_2^{(1)}$, $a=1$: the channel sum with $m_b=\sqrt3m$ gives $\Delta M_1=m\big(\frac14-\frac{3\sqrt3}{2\pi}\big)$, which equals $-m_1[\frac h{2\pi}-\frac14\cot\frac\pi h]$ (Hollowood's formula with the $\frac14$ coefficient). Fix: write the channel form with $m_b$ and state the reduction.

**A2. Attribution to Hollowood's equations (`16-dhn.qmd`, paragraph after `eq-hollowood-mass`).** The text says Hollowood's paper "prints $\frac12\cot\frac\pi h$, because its bound-state sum and continuum integral (his (3.2) and (3.26)) are evaluated as twice their actual values". The coefficient check in A1 supports $\frac14$, but the specific equation numbers and the factor-of-2 explanation were not checked against Hollowood (1993). Check the paper, and either give the exact location of the error or remove the sentence and cite the correction in MacKay–Watts (1995) only.

**A3. Proof of the thimble Gaussian (`17-thimbles.qmd`, Proposition `prp-thimble-gaussian`).** The proof begins "Triangularize $B$ (Schur form)". For a complex symmetric $B$ the quadratic form $\frac12 q^{\mathsf T}Bq$ is preserved only by complex orthogonal changes of variable, and a unitary Schur form does not preserve it. The correct tool is the complex-orthogonal normal form of complex symmetric matrices (every complex symmetric matrix is complex-orthogonally similar to a direct sum of symmetric Jordan blocks; see Horn–Johnson, *Matrix Analysis*, the section on complex symmetric matrices). Rewrite the proof with that normal form, and state the branch convention for $\prod\lambda^{-m_\lambda/2}$ (principal branch, since the spectrum avoids $(-\infty,0]$). Also check the sign of the upward flow in the thimble definition against the claim $\frac{d}{dt}S=|S'|^2$.

**A4. "Verified" in place of "checked numerically" (`18-jordan-chains.qmd`, Section `sec-counting-rule-transmission-factors` and Section `sec-hirota-form-dn1-en1`; also `19-one-loop-masses.qmd`, the paragraph after `cnj-counting-rule`).** The text says hypothesis (D) "is verified for every case used in this book", "verified for $d_n^{(1)}$ with $n\le8$", and "$(D)$ is verified for ... the 19 foldings". The evidence is numerical: evaluation at 96 random complex $k$ per pair, to $10^{-25}$. That is strong evidence but not a verification, and the theorems of Part V are conditional on it. Replace "verified" by "checked numerically at 96 random points per pair" wherever it occurs, and state at each theorem whether it depends on (D).

**A5. The $a_{2n}^{(2)}$ folding (`18-jordan-chains.qmd`, Section `sec-hirota-form-dn1-en1`, Foldings paragraph; `19-one-loop-masses.qmd`, table row $a_{2n}^{(2)}$).** The text says $a_{2n}^{(2)}$ is obtained from $d_{2n+2}^{(1)}$ "by an order-4 automorphism that moves the affine node". Twisted algebras are defined by automorphisms of order $k$ (Part II, `sec-twisted`), and $a_{2n}^{(2)}$ has $k=2$. This is the same issue as item A1 in the Part II brief: the parent in `tbl-foldings` (Part II) is $d_{2n+2}^{(1)}$, the node count is inconsistent with an order-2 fold for $n=1$ ($d_4^{(1)}$ has 5 nodes, $a_2^{(2)}$ has 2), and Part V now justifies it by an order-4 automorphism that is not a twist of the stated kind. Resolve this together with Part II A1, and make Parts II and V agree.

**A6. Sign of the one-loop floating Coxeter shift for $b_n^{(1)}$ and $c_n^{(1)}$ (`19-one-loop-masses.qmd`, table and the paragraph before it, versus `14-non-simply-laced.qmd`, `tbl-floating-h`).** Chapter 14 gives, for the particle masses, $H=2n-\tfrac12B$ for $b_n^{(1)}$ (decreasing to $h^\vee=2n-1$) and $H=2n+\tfrac12B$ for $c_n^{(1)}$. Chapter 19 gives $\delta H=+\frac1{4\pi}$ for $b_n^{(1)}$ and $\delta H=-\frac1{4\pi}$ for $c_n^{(1)}$, with $H_s=h+\delta H\,\beta^2$. The signs are opposite for both algebras, while $g_2^{(1)}$ and $f_4^{(1)}$ agree between the two chapters. Chapter 19 writes "$H=2n+\beta^2/4\pi$" for the spinor of $b_n^{(1)}$ and compares it to the exact [@gmw1996]. Decide whether the chapter 14 and chapter 19 values refer to particles or solitons, state it in both places, and correct the sign that is wrong for the quantity compared.

**A7. The coupling map $1/\omega$ and its relation to $B$ (`19-one-loop-masses.qmd`, Section `sec-one-loop-masses-folded`, paragraph "Comparisons use the normalization ...").** The text says the gradation formula gives $\omega=4\pi/\beta^2-1$, so $1/\omega=\beta^2/4\pi+O(\beta^4)$, and then "$H_s=6+2/\omega\approx6+\beta^2/2\pi$". Part IV gives $B=\beta^2/2\pi+O(\beta^4)$ and $2/\omega=\beta^2/2\pi$ at leading order. These agree, but the relation $B=2/\omega+O(\beta^4)$ is used without being stated. State it once, with its origin.

**A8. Symbol $u$ for the breather rapidity and the fusing angle (`17-thimbles.qmd`, Section `sec-which-saddles-contribute`, proposition `prp-static-saddles`).** The text says "bound states of two species-1 solitons at rapidities $\pm iu$ ... tends to $M_2$ at the fusing angle $u=\pi/4$". Here $u$ is half the fusing angle $U=2u=\pi/2$, which is the convention of Part III, chapter 10 (breathers) and not of Part I (fusing angles). Write "$u=\pi/4$, so that the fusing angle is $U=\pi/2$".

**A9. Reduced-model numbers (`17-thimbles.qmd`, Section `sec-which-saddles-contribute`, the bullet "[Numerical]" on $a_2^{(1)}$ with $L=1$).** The text gives "54 complex critical points", "exactly four contribute", and a residual "$6.6\times10^{-9}$ of the total at $\hbar=0.15$". These come from `foundations_code/thimbles/`. Rerun the script and report the numbers; they were not rerun in this review.

**A10. Phase-shift claims (`19-one-loop-masses.qmd`, Remark "Phase shifts").** The text says the naive energy-weighted phase-shift sum "for the sine-Gordon kink gives $0$ instead of $-m_1/\pi$, and for the $\phi^4$ kink it is too large by $3m/2\pi$". These follow from the cut-off prescription described, but no calculation is shown. Give the calculation for the sine-Gordon case (it is short) or mark the claim as numerical.

**A11. Numerical "closed form" agreement (`19-one-loop-masses.qmd`, Section `sec-one-loop-masses-folded`, "Calibration").** The text says the closed form, unpaired form, direct phase-shift form and heat-trace route "agree to $10^{-12}$", and gives other agreement figures. These come from `foundations_code/oneloop/heat_trace_checks.py` and `phi4_check.py`. Rerun them and state the date of the run.

## B. Gaps for a reader following Parts I–IV

**B1. Thimbles and intersection numbers.** `17-thimbles.qmd` gives the thimble, dual thimble and $n_\sigma$ in one variable. The reader needs the definition of the intersection number $\langle\mathbb R,\mathcal K_\sigma\rangle$ in several variables, with orientation, and one example of a non-trivial case beyond the zero-dimensional $i\phi^3$. Add a sentence and a pointer to Pham (1983) and Witten (2010).

**B2. Operator-theoretic terms used without definition (`18-jordan-chains.qmd`).** Each needs one line: Riesz projection (defined by $-\frac1{2\pi i}\oint R$, state it); root space $\mathcal R(\lambda_0)$ (defined); Birman–Schwinger operator $K(\lambda)$ (defined in `prp-net-count`); Konno–Kuroda formula (named; state the resolvent identity it expresses); m-sectorial; Jost–Pais formula; Fredholm determinant of trace class operators; Evans function (named only); "threshold resonance" and "half-bound state" (defined only by example).

**B3. Hypothesis roadmap.** Chapter 18 has hypotheses (R), (H), (D), (X1)–(X3), and the grouping into mass levels. These are the core of the chapter, and the reader must track them through `thm-jordan-chains`, `prp-net-count`, `thm-hirota-factorization`, `prp-hirota-form` and `prp-one-loop-mass`. Add a table at the start of chapter 18: each hypothesis, where it is used, and whether it is proved, checked numerically, or assumed for a family.

**B4. Trace-class and resolvent assumptions.** The proof of `prp-net-count` part 2 uses Hilbert–Schmidt and trace-class properties of $V_jR_0$ with a citation to Simon (2005). Say in one line why the exponential decay of $M(x)-M_0$ gives the required $\ell^1$ bounds.

**B5. Counterterm conventions.** `16-dhn.qmd` introduces the counterterm $\delta M_{\rm ct}$ and `17-thimbles.qmd` introduces the term $-\frac14\operatorname{Tr}[(\mathcal A-\mathcal A_0)\mathcal A_0^{-1/2}]$. The reader needs one sentence linking them to the normal ordering of Part IV, chapter 12, and noting that the two are the same quantity.

**B6. Numerical evidence.** Many claims are labelled "[Numerical]" in Part V, and the scripts are in `foundations_code/`. A reader who does not run them cannot check the claims. State in the chapter introduction which claims are numerical, and give the script for each in one place (a table would do).

**B7. Conjecture and theorem labels.** `cnj-counting-rule` has a Part 1 that is a theorem given (D), and a Part 2 that is a conjecture. This is clear in the text but easy to misread at a glance. Add "Part 1 is a theorem given (D); Part 2 is a conjecture" to the statement heading, as it is partly done.

## C. Notation clashes

**C1. $\psi$.** The fluctuation solutions $\psi_c$, $\psi_b^\mp$ (`18-jordan-chains.qmd`) and the function $\psi(t)=(t\cos t-\sin t)/\pi$ of `eq-dhn-psi` and `prp-one-loop-mass`. Rename the function, for example $\Psi(t)$ or $\varpi(t)$, everywhere in Part V.

**C2. $C$.** The set of closed channels (`18`, Section `sec-counting-rule-transmission-factors`), the affine Cartan matrix $C$ (`18`, Hirota form section, "$NC$"), the cubic couplings $C_{abc}$ (Part III), the Casimir $C(\vec\lambda)$ (Part II), and the constants $c_b$ (`19`). Rename the Cartan matrix to $\mathsf C$ or $A$ in the Hirota section, and the set of closed channels to $\mathcal C$.

**C3. $\mathcal A$ and $A$.** The fluctuation operator $\mathcal A$ and the angle $A=\pi a/h$ in the example of `sec-example-an1-zeros-alone` and `sec-example-cn1-jordan-blocks`. Also $A_{ab}$ (interaction coefficient, Part III). Rename the angle to $\alpha$ or $\vartheta_a$.

**C4. $B$.** The Gaussian matrix $B$ (`17`), the angle $B=\pi b/h$ (`18`, examples), the parameter $B(\beta)$ (Part IV), and $B$ in the comparison with breathers (`19`). Rename the angle to $\gamma$ or $\vartheta_b$.

**C5. $K$.** The dual thimble $\mathcal K_\sigma$ (`17`) and the Birman–Schwinger operator $K(\lambda)$ (`18`). Rename the Birman–Schwinger operator to $\mathcal K(\lambda)$ or $Q(\lambda)$.

**C6. $S$.** The action $S(x)$, $S_E$ (`17`), and the S-matrix $S_{ab}$ (Parts I, IV). Note the clash at the first use in `17`, or rename the action to $\mathcal S$.

**C7. $R$.** The resolvent $R(\lambda)$, $R_0$ (`18`), the functions $R_b^\mp$ in hypothesis (H) (`18`), and the R-matrix (Part II). Rename $R_b^\mp$ to $Q_b^\mp$.

**C8. $\omega$.** The mode frequency $\omega_n$ (`16`, `19`), the gradation parameter $\omega=4\pi/\beta^2-1$ (`19`, comparison section), and the $h$-th root of unity $\omega=e^{2\pi i/h}$ (Part III, `eq-an-soliton`). Rename the gradation parameter to $\Omega$ or $\nu$.

**C9. $M$.** The potential matrix $M(x)$ (`16`, `18`) and $M_0$, against the masses $M_a$ and classical soliton mass $M$ of Part III. Rename the matrix to $\mathsf M$ or $U(x)$ in Part V, as in the Part III brief (C7).

**C10. $n$.** The rank in $\hat g$, the channel index $n$ in $c_n^{(1)}$ (`18`, `19`), the Kac labels $n_j$, the multiplicities $n_{\mu,w}$ and $n_b$ (`18`, `19`), and the intersection numbers $n_\sigma$ (`17`). Keep $n_j$ for Kac labels, use $\mu_w$ for multiplicities, and $\nu$ for the intersection numbers.

**C11. $V$.** The potential $V$ (`16`), the difference $V=M-M_0$ (`19`, Born lemma), and the factors $V_1,V_2$ with $\mathcal A-\mathcal A_0=V_1V_2$ (`18`). Use $\delta M$ for $M-M_0$.

**C12. $\kappa$ and $z$.** $\kappa_b=\sqrt{m_b^2-\lambda}$ (`16`, `18`) and $\kappa$ as a parameter of the Pöschl–Teller well (`18`, `sec-counting-rule`, part 3). Also $z=\kappa/m$ (`16`, `18`, `19`), which is not the spectral ratio $z$ of Part II. Note both in the first use of each section.

**C13. $\chi$.** The Jordan chains $\chi_{c,j}$ (`18`) and the characters $\chi_g$ of the centre (`18`, Hirota section). Rename the characters to $\varepsilon_g$ or $\theta_g$, and check the $\varepsilon$ used for the perturbation parameter in `18`.

**C14. $E$.** The Hirota exponential $E=e^{m_ax+\xi}$ (`16`, `18`, `19`) and the energy $E$ of Part III, and the sum $E_i$ in (H). Consistent with Part III except for the energy; see the Part III brief (C5).

**C15. $T$.** The heat-trace parameter $t$ (`18`) and the gradation parameter $T$ (Part VI, cited in `19`). Only a notational point, but worth a sentence in `19` where $T$ appears.

## D. Items checked and found correct (do not re-check unless you change them)

- Chapter 16: the mode-sum formula (`eq-dhn-mode-sum`); the Pöschl–Teller transmission factor $X(z)=\prod_j(z-z_j)/(z+z_j)$ with $z_j=\sin t_j$ and $\omega_j=m\cos t_j$; the Cahill–Comtet–Glauber formula (`eq-ccg`) and its pole–zero form (`eq-dhn-psi`), including the statement that a zero at $z_j$ and a pole at $-z_j$ each give half of (`eq-ccg`).
- Sine-Gordon $\Delta M=-m/\pi$; $\phi^4$ kink $\Delta M=m(\frac1{4\sqrt3}-\frac3{2\pi})$, with $\ell=2$, $\mu=m/2$ in the Pöschl–Teller form (exercise `pt-transmission`), which matches $V''=m^2(1-\frac32\operatorname{sech}^2\frac{mx}2)$.
- Chapter 16 `eq-an-transmission`: $X_{ab}$ is the two-soliton interaction coefficient (Part III, `eq-interaction-coefficient`) at $z=\cos$-form; the $a_4^{(1)}$ values $\lambda=1.25\,m^2$, $m_1^2=1.382\,m^2$ and $m_2=2m\sin\frac{2\pi}5$ are reproduced, with $X_{12}$ zero and $X_{11}$ pole at the same energy.
- Chapter 17: the one-variable example $S=\frac{x^2}2+\frac{ix^3}3$: critical points $0$ and $i$ with $S(i)=-\frac16$; the coefficient $-\frac5{6\lambda}$ from $\langle x^6\rangle=15/\lambda^3$; the antilinear map $x\mapsto-\bar x$.
- Chapter 17 Mellin representation of $\sqrt\lambda$ and $\lambda^{-1/2}$, and the combined formula (`eq-one-loop-mellin`).
- Chapter 17 `eq-one-loop-mass` has the normal-ordering counterterm with coefficient $-\frac14$.
- Chapter 18, Lemma `lem-wronskian`: $W=-2ik$ at $-\infty$ and $-2ikX_b(k)X_{\bar b}(-k)$ at $+\infty$, so $X_bX_{\bar b}(-k)=1$.
- Chapter 18 `thm-jordan-chains`, Step 4: self-orthogonality from bilinear symmetry.
- Chapter 18 `prp-net-count` part 1: the Riesz projection and $\operatorname{ord}_{\lambda_0}d=\operatorname{tr}P$.
- Chapter 18 $a_3^{(1)}$ closed form: $(1+e^{-2t})\operatorname{erf}\sqrt{2t}$ from $X_{11}X_{13}$ at $\mu=\sqrt2$ and $X_{12}$ at $\mu=2$.
- Chapter 18 `prp-cn-threshold-blocks`: $X_n(z)=(\frac{z-\sin A}{z+\sin A})^2$ has a double zero at $z=\sin A$, and $m_n^2\cos^2A=m_{n-a}^2$; $\tau_j=1+2\cos(\pi ja/n)E+\cos^2A\,E^2$ is the coincident pair with $A_{a,2n-a}(0)=\cos^2A$.
- Chapter 19 Born lemma: $X_b(k)=1-\frac{ic_b}{2k}$; for sine-Gordon $c=-4$ gives $X\approx1+2i/k$, which matches $(k+i)/(k-i)$.
- Chapter 19 counterterm identity $\sum_{m_b=\mu}c_b=-4\mu\sum_wn_{\mu,w}w$.
- Chapter 19 `prp-one-loop-mass`, step 3–5: $F_w$, the Frullani integral giving $\frac w{2\pi}\int_0^1\log(1-w^2+w^2u^2)du$, and $\int_0^1\log(\cos^2\theta+u^2\sin^2\theta)du=2\theta\cot\theta-2$, giving $\psi(\theta)$.
- Chapter 19 phase-shift integral: $\int_0^\infty\frac{dk}{\sqrt{1+k^2}(k^2+w^2)}=\frac{\arccos w}{w\sqrt{1-w^2}}$ (check at $w\to1$ gives $1$).
- Chapter 19 table: $\phi^4$ and sine-Gordon calibration; $a_2^{(1)}$ channel check (A1).
- All 57 `@` labels and all 33 citation keys used in Part V resolve. The gap ids H7 and H8 match the claims in chapters 17 and 19. The code files cited in chapters 17 and 19 exist (`textbook_code/semiclassics.py`, `foundations_code/thimbles/`, `foundations_code/oneloop/`, `foundations_code/hirota/`).

## E. Verification checklist for the agent who makes the fixes

1. Re-run the label check on `part5-semiclassics/` (exclude `.claude/` and `_site/`).
2. For A1, redo the $a_2^{(1)}$ channel sum (see D) and check the $a_3^{(1)}$ closed form in code before editing.
3. For A3, check the complex-orthogonal normal form statement in Horn–Johnson before rewriting the proof.
4. For A6, check the gradation and floating-Coxeter conventions in the cited papers before changing either chapter.
5. For A9–A11, rerun the named scripts and record the run date and output in the report.
6. For each notation rename (C1–C15), grep the whole book (excluding `.claude/` and `_site/`) for the old symbol in the same meaning, and change only those occurrences. Several symbols (C, A, B, S, R, K) are used in Parts II–IV with other meanings, so check each use.
7. Report what changed, what was left, and anything you could not verify.

## Resolution (2026-10-10)

- **A1.** Done. `sec-dhn-an` now displays the channel form $\Delta M_a=\frac12\sum_bm_b\sum_w\operatorname{ord}_w(X_{ab})\varpi(\arcsin w)$, says that it reduces to `eq-dhn-psi` for one channel and that `prp-one-loop-mass` proves it. `prp-one-loop-mass` already had $m_b$. Rechecked the $a_2^{(1)}$, $a=1$ sum by hand ($m(\frac14-\frac{3\sqrt3}{2\pi})$, equal to `eq-hollowood-mass`) and the $a_3^{(1)}$ closed form $(1+e^{-2t})\operatorname{erf}\sqrt{2t}$ (levels $\sqrt2$, $w=1$ and $2$, $w=1/\sqrt2$).
- **A2.** Checked against Hollowood, arXiv:hep-th/9209024 (PDF read). His (3.2) is $\sum\frac12\nu$, correct as written, but his evaluation (3.3) is exactly twice it (checked numerically for $h=3,\dots,7$). His (3.27) is twice the value of the integral (3.26) (numerical quadrature, $h=3,5,7$; for even $h$ after including the threshold delta term). With (3.2) and (3.26) evaluated correctly the total is $-m_a[\frac h{2\pi}-\frac14\cot\frac\pi h]$. The text now names (3.3), (3.27) and (3.28). MacKay–Watts (hep-th/9411169, Table 3) have $1-\frac{\beta^2}{8h}\cot\frac\pi h$ times $M^{\rm cl}$, that is the $\frac14$ coefficient; the sign convention of their $h/2\pi$ term was not compared.
- **A3.** Done. The proof of `prp-thimble-gaussian` is rewritten with Takagi's factorization $B=U\Sigma U^{\mathsf T}$ (a congruence): the thimble is $\bar U\mathbb R^n$, the integral is $(2\pi)^{n/2}\det\bar U\prod\sigma_j^{-1/2}$, its square is $(2\pi)^n/\det B$, and the sign is fixed by transporting the orientation along $(1-s)I+sB$, which gives the principal branch. The statement now gives the branch. A remark adds the complex-orthogonal normal form (symmetric Jordan blocks), and the exercise asks for the Takagi form instead of the Schur form. Checked numerically (orientation transport, the matrix of `exr-jordan-determinant` and random $3\times3$ cases). New bib entry `horn2012` (DOI 10.1017/CBO9781139020411, verified via doi.org); cited as Sec. 4.4 without theorem numbers, which were not checked. The upward flow $\dot z=\overline{S'}$ gives $\frac d{dt}S=\lvert S'\rvert^2$; no change needed.
- **A4.** Done. "Verified" is removed throughout ch. 18–19. Correction to the brief: the scripts (`hirota/fluct_scan.py`, `fold.py`) use $3(h+2)$ random complex $z$ per pair, from 24 ($d_4^{(1)}$) up to 96 ($e_8^{(1)}$), not 96 per pair. The text now says this. The dependence on (D) is stated at `prp-hirota-form`, in the "Which solitons" bullet, after `cnj-counting-rule` and in the H8 paragraph, and in the new hypothesis table.
- **A5.** Done, consistent with the Part II resolution (`sec-affine-folding`, "The parent of $a_{2n}^{(2)}$"). Ch. 18 now refers to that section, and says that both foldings are foldings of the Toda theory, not the order-2 twist of `sec-twisted`. It also says that the one-loop masses agree (checked in `oneloop/analysis_twisted.log`: $a_{4,6}^{(2)}$ from $d_{2n+2}^{(1)}$ and from $a_{2n}^{(1)}$ give the same $\Delta M/M_{\rm cl}$). Table row in ch. 19: "$d_{2n+2}^{(1)}$ or $a_{2n}^{(1)}$".
- **A6.** No sign error. The brief's claim is wrong: ch. 14 gives the particle floating Coxeter numbers at real coupling, ch. 19 the soliton ones at imaginary coupling. They are related by the Lie-duality rule of ch. 19: soliton $b_n^{(1)}$ ($+1/4\pi$) = particle $c_n^{(1)}$ ($+1/4\pi$), soliton $c_n^{(1)}$ ($-1/4\pi$) = particle $b_n^{(1)}$ ($-1/4\pi$), and $g_2^{(1)}$, $f_4^{(1)}$ are self-dual. Both places now say which quantity they mean (ch. 19 before the table and in the $b_n^{(1)}$ sentence of the resolution section; one sentence in the `tbl-floating-h` caption in ch. 14).
- **A7.** Done. Ch. 19 now states $2/\omega=B(\beta)+O(\beta^4)$, and exactly $2/\omega=-B(\beta_r)$ with $\beta_r^2=-\beta^2$. It also gives the $g_2^{(1)}$ example $3B(\beta/\sqrt3)=2/\omega+O(\beta^4)$. Note that the relation the brief proposed, $B=2/\omega$, holds for $B(\beta)$ but not for the $B(\beta_*)$ of `tbl-floating-h`.
- **A8.** Done: "$u=\pi/4$, so that the fusing angle is $U=2u=\pi/2$".
- **A9.** Rerun 2026-10-10 in a scratch copy. `thimbles/check_cg.py` (broad search) finds 54 critical points, while `run_a2L1_final.py` finds 23 with its narrower search, so the text now says "a broad search finds 54". `run_a2L1_final.py` reproduces its stored output: four saddles contribute with $\lvert\nu_\sigma\rvert=1$. At $\arg\hbar=+0.02$, however, the signs of the complex pair are reversed relative to $\arg\hbar=0$, so "the $n_\sigma$ are unchanged" was too strong. It now says that the same four saddles contribute, with $\lvert\nu_\sigma\rvert=1$. `check_th4.py` cannot run as committed: it loads `cps_a2L1_k1.00.npy` (not in the repo) and hard-codes $n_c,n_d=1,-1$ from a missing `check_th2.py`. After regenerating the file with the search of `run_a2L1_final.py`, $\hbar=0.15$ gives a residual of $-3.0977\times10^{-7}$ (stored $-3.0976\times10^{-7}$) and a complex-pair prediction of magnitude $3.0315\times10^{-7}$ but of the opposite sign. The magnitudes reproduce the stated $6.6\times10^{-9}$ (2%); the sign depends on the tangent-frame orientation of the regenerated critical points and was not settled independently. The $\hbar=0.1$ rows were not rerun. Suggest committing the `.npy` file and `check_th2.py`.
- **A10.** Done. The sine-Gordon calculation is now given (the continuum integrand vanishes identically, $c=-4m_1$, so the whole answer is $c/4\pi=-m_1/\pi$), and the $\phi^4$ statement follows from $c=-6m$.
- **A11.** Rerun 2026-10-10 in a scratch copy. `heat_trace_checks.py` output is identical to the stored log (Born $\le2.6\times10^{-11}$, Mellin $\le4.3\times10^{-10}$, phase-shift form $2.0\times10^{-13}$). `phi4_check.py` gives $1.6\times10^{-15}$. The run date is now in the text. The `--grid 768` run was not repeated.
- **B1.** Done: a paragraph "Several variables" in `sec-picard-lefschetz` with the orientation convention, relative homology, the Airy example ($a>0$: $t=i\sqrt a$ contributes, $-i\sqrt a$ does not) and a pointer to the two-variable reduced model.
- **B2.** Done: one-line definitions of the Riesz projection, resolvent, Fredholm determinant (`simon2005`, Ch. 3), Birman–Schwinger operator, Konno–Kuroda formula (via the second resolvent identity), m-sectorial, Jost–Pais formula, Evans function, half-bound state and threshold resonance.
- **B3.** Done: `tbl-jordan-hypotheses` at the start of ch. 18.
- **B4.** Done: one line in ch. 18 (Hilbert–Schmidt property of $f(x)g(p)$ with $f,g\in L^2$) and in ch. 19 (definition of $\ell^1(L^2)$ and why exponential decay gives it).
- **B5.** Done: ch. 16 derives $\delta M_{\rm ct}$ from normal ordering and identifies it with the $-\frac14\operatorname{Tr}$ term; ch. 17 refers back.
- **B6.** Done: new section `sec-dhn-numerics` with `tbl-part5-numerics` at the end of ch. 16 (script names taken from the docstrings in `foundations_code/`).
- **B7.** Done: the heading of `cnj-counting-rule` now reads "[Part 1 is a theorem given (D); Part 2 is a conjecture]".
- **C1.** $\psi(t)\to\varpi(t)$ throughout Part V, with a note at `eq-dhn-psi`.
- **C2.** Closed channels $C\to\mathcal C$. Cartan matrix $C\to\mathsf C=(a_{ij})$, and the matrix inverted in the centre characters is now written $\mathsf C_{\rm fin}$ (the finite Cartan matrix, as in `hirota/lie.py`).
- **C3, C4.** The angles $A,B,C=\pi a/h,\pi b/h,\pi c/h$ are now $\vartheta_a,\vartheta_b,\vartheta_c$ in ch. 16 and 18. The coefficients $A_1',A_2'$ and the Gaussian matrix $B$ of ch. 17 are unchanged.
- **C5, C7.** Birman–Schwinger $K(\lambda)\to Q(\lambda)$, and $R_b^\mp\to Y_b^\mp$, not $Q_b^\mp$: the brief suggested $Q$ for both.
- **C6.** A note at the first use of $S$ in ch. 17.
- **C8.** A note, not a rename, because Part VI uses $\omega$ for the gradation throughout.
- **C9, C11.** $M(x)\to\mathsf M^2(x)$ and $M_0\to\mathsf M^2$, the vacuum mass matrix of Parts III–IV; $M-M_0$ (and $V$ in the Born lemma) $\to\delta\mathsf M^2$. $V_1V_2$ is unchanged.
- **C10.** Multiplicities $n_{\mu,w}\to o_{\mu,w}$ ("order"); $\mu_w$ would clash with the mass level $\mu$. Intersection numbers $n_\sigma\to\nu_\sigma$ in ch. 17, and in ch. 33 (one occurrence). Levinson's $n_b\to N$ in ch. 16.
- **C12.** The Pöschl–Teller $\kappa$ was dropped ("depth parameter $\mu w$"), and there is a note on $z$ in ch. 16.
- **C13.** $\chi_g\to\gamma_g$. $\varepsilon_g$ was not used, because $\epsilon$ is the perturbation parameter in the same section.
- **C14.** A note at the first $E$ and $\omega$ in ch. 18, following Part III. The frequency $\omega$ in $e^{ikx-i\omega t}$ is replaced by $\sqrt\lambda$, because it clashed with the root of unity in the same sentence.
- **C15.** A sentence at $H_b=6T/(T+1)$.
- **Files outside Part V:** `references.bib` (new `horn2012`); `part4-real-coupling/14-non-simply-laced.qmd` (one sentence in the `tbl-floating-h` caption, A6); `part7-foundations/33-outlook.qmd` ($n_\sigma\to\nu_\sigma$).
- **Verification.** Label check (script over all `.qmd` anchors and `references.bib` keys) on the seven touched `.qmd` files: 0 problems. Scratch render (`quarto render --to html` in a copy): exit 0, no warnings or errors in the log, no `?@` in the rendered Part V, ch. 14 or ch. 33 pages. The render predates the last one-sentence A9 wording edit in ch. 17, which has no references.
