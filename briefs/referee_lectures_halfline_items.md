# Phase 2 of the referee pass on the half-line lecture notes: the coordinator's items

This continues `briefs/referee_lectures_halfline.md`. Your phase-1 report, `briefs/referee_lectures_halfline_phase1_report.md`, is complete and committed; do not change it.

The rules, notation, sources and standard of evidence of the phase-1 brief still apply, with one change: the blinding rule is lifted. You may now read the earlier referee material on Part VII of the book (`briefs/referee_half_line_*_report.md`, `briefs/referee_bsg_dirichlet_sign_report.md`, `briefs/reader_review_part7.md`). Do not take their verdicts as evidence.

Write this phase to `briefs/referee_lectures_halfline_report.md`.

## For each item below

- **Overlap.** Was the item, or part of it, in your phase-1 risk list? Give the entry number and rank, and whether you checked it.
- **Check.** If your phase-1 check already covers the item fully, say so and point to the section; do not repeat it. Otherwise check it now: verdict (correct / wrong / unclear), your own computation and its output, and, for errors, replacement text ready to paste, with file and line.

The items were written before your phase 1 started. Their SHA-256 was recorded in `briefs/README.md` at that time. Line numbers refer to the lectures as on `main` at 035351d.

## Items (most important first)

1. **What the lectures present as settled, against the book's labels.** Several results the lectures build on carry weaker labels in the book:
   - `cnj-bsg-uv-ir` (the UV–IR relation of boundary sine-Gordon), `cnj-neumann-coideal` (the Neumann boundary has $\hat\epsilon_j=\epsilon_*$, the uniform one $\hat\epsilon=0$) and `cnj-neumann-amplitude` (the Neumann particle amplitudes) are **[Conjecture]**;
   - `prp-solitonic-vacua`, `prp-a2-excited`, `prp-solitonic-amplitude` and `prp-ly-reflection` are **[Numerical]**;
   - Appendix E says that the Neumann particle amplitudes were derived "given @cnj-neumann-coideal, the continuation in the coupling, and a CDD factor chosen for exact duality with the uniform amplitude, which one loop supports but does not fix".

   Against this, the lectures write, among other things: "we will resolve both puzzles … uncovering a breathtaking weak–strong boundary duality" (`index.qmd`, idea 5); "Why is that assumption wrong?", "confirming both the intercept $\epsilon_*$ at Neumann" and "the triumphant confirmation of this dictionary" (Lecture 5, section "Resolving the Neumann vs. Uniform Boundary Puzzle"); "we discover that weak–strong duality … swaps the uniform boundary … and the Neumann boundary" (Lecture 6, callout); "produces the exact real-coupling particle reflection amplitudes … for the Neumann boundary" (Lecture 6, "What You Will Learn"). Also: "the purely Dirichlet boundary condition ($\phi(0,t)=0$) is *never* integrable for $n\ge2$" (`index.qmd`, idea 2) and "full Dirichlet boundary conditions are never integrable in Sklyanin's framework" (Lecture 2), against the scope of `prp-no-dirichlet` (field-independent $K$, finite-dimensional evaluation representations without trivial summands).
   - Go through the six lectures and the index. List every sentence that states as established, derived or discovered what the book labels [Conjecture] or [Numerical], or that states a result more broadly than the book. For each, give the book statement and its label, and replacement text.
   - Is the duality $K_a^N(\theta;B)=K_a(\theta;2-B)$ a consequence of the soliton calculation of Lecture 6, or of the CDD choice that Appendix E describes? Reconstruct the step from $K^{(+)}_{B_1}$ to $K_1^N$ in Lecture 6, and say where a CDD factor enters ("with the CDD factor $P_{\omega/2}=\dots$", Blackboard Calculation 2, step 2).
   - Do the $B\to0$ limits ($K_a\to$ the classical uniform factor, $K_a^N\to1$) and the one-loop agreement identify the boundaries, or only check consistency?

2. **Lecture 1, the boundary Lee–Yang calculation** (Blackboard Calculation 3, `eq-lec1-s2theta` to `eq-lec1-ly-shift`).
   - Is `eq-lec1-s2theta`, $S_{\rm LY}(2\theta)=-(2)(4)$ with blocks of argument $\theta$, right? Derive $S_{\rm LY}(2\theta)$ as a product of blocks of argument $\theta$.
   - The text between `eq-lec1-s2theta` and the conclusion $\epsilon=\epsilon^2=+1$ contains steps that are abandoned or left open ("Wait—let us be extra careful!", "$(x)=-(12-x)$? No, …", "$=\cdots$", "Better yet"). Does the text as printed derive $\epsilon=+1$? Write a complete derivation that a student can follow, as replacement text.
   - The claim that Ghoshal and Zamolodchikov printed $-(1)(3)(4)^{-1}$ in their (3.51): quote the source, and compare with Appendix D.
   - Check the residue $-2i(2\sqrt3-3)$, the coupling $g=2i\sqrt{2\sqrt3-3}$, and $\Gamma^2=-2\sqrt3\,(2+\sqrt3)$ (`eq-lec1-ly-g` and the sentence after it) against the residue conventions of `eq-lec1-one-point-residue`, of the book (`eq-one-point-residue`) and of Ghoshal–Zamolodchikov (3.45). In which normalization is $\Gamma$ defined, and is it the one in which $S_{\rm LY}$ has residue $i\Gamma^2$?
   - Check "the identity cardy boundary", the dimension $x_B=-\frac15$ of $\phi_B$, and $c=-\frac{22}5$.
   - Check `eq-lec1-ly-shift` ($K_b^{(1)}=K_{b+2}$) against Dorey–Pocklington–Tateo–Watts and the book.
   - `exr-lec1-crossing-forms`, part 4: is $f_x(\theta)=(2x)(6-2x)$ in the blocks of `eq-lec1-block-def`? Is $K_b/K_{(\mathbb 1)}=f_{(1+b)/6}/f_{(1-b)/6}$?

3. **Lecture 1, conservation laws, reflection equations and the rotation.**
   - With $x^\pm=\frac12(t\pm x)$ and $\partial_\pm=\partial_t\pm\partial_x$, rewrite $\partial_-T_{s+1}=\partial_+\Theta_{s-1}$ in $(t,x)$. Is the conserved density $T_{s+1}+\Theta_{s-1}$ or $T_{s+1}-\Theta_{s-1}$? Then check `eq-lec1-flux-balance`, `eq-lec1-surface-condition`, "$Q_1+Q_{-1}=E$" next to "$P^++P^-=2E$", and "$\Sigma_1$ is the boundary potential energy". Compare with the book's section on half-line charges.
   - `eq-lec1-reflection-braid`: the text says that multiplying by $P$ "puts every reflection matrix in the **first** tensor slot", and the display has $\mathbf 1\otimes K$. Check against the book's `eq-reflection-braid`.
   - `eq-lec1-reflection-twisted`: trace the spaces through both sides with the stated rule $\check R_{A,B}:V_A\otimes V_B\to V_B\otimes V_A$ (this is `exr-lec1-twisted`). Do both sides typecheck? Compare with the book's `eq-reflection-twisted` and with the version in Lecture 4, Step 2 of Blackboard Calculation 1.
   - The rotation: Lecture 1 and its figure caption write $(x,t_E)\mapsto(-t_E,-x)$, `index.qmd` writes $(x,t_E)\mapsto(-t_E,x)$, and the book's figure `fig-boundary-state-rotation` uses $(x,y)\mapsto(y,-x)$. Which of these are rotations? Check the signs in `eq-lec1-rapidity-rotation`.
   - Compare `eq-lec1-boundary-state`, `eq-lec1-pair-amplitude` and `eq-lec1-crossing-unitarity` with the book's `eq-boundary-state` and `eq-boundary-crossing-unitarity` (index order, $S$ against $\check S$, the measure).
   - `exr-lec1-crossing-forms` calls `eq-lec1-cu-cd` the "Corrigan–Delius form"; the text calls it the Corrigan–Dorey–Rietdijk–Sasaki form. Which source has it in that form? Which has `eq-lec1-cu-dg`?
   - Compare the bootstrap `eq-lec1-bulk-bootstrap` (with $S_{b\bar a}$) with the one in Lecture 6, Blackboard Calculation 1, Step 3 (with $S_{ab}$) and with the book's `eq-boundary-bootstrap-bulk`.
   - `exr-lec1-elastic`: is the argument (dividing by $\frac12e^{s\theta'_{\max}}$, $s\to\infty$) correct and complete? What about the particles with $\theta'_j=\theta'_{\max}$ and the requirement "same mass"?

4. **Lecture 3, mirror pairs and vacua.**
   - Step 1 of Blackboard Calculation 1 says $k=m_a\cosh\theta=2ms$; `exr-lec3-mirror-square` says $k=m_a\cosh\theta=2m\frac{\sin(\pi a/h)}s$. Check both against $s=\sin\frac{\pi a}h/\cosh\theta$, and check the exercise's conclusion with the correct $k$.
   - The factorization "$\zeta^{ja}E_1\big(1+\zeta^{-2ja}AE_2\big)$" before `eq-lec3-virtual-collision`: check, and check $x_{\rm mirror}$, $d$ and `eq-lec3-time-delay`, with the signs of the interpretations ("delayed", "early", "behind the wall"), against the book's `eq-boundary-time-delay`, `eq-virtual-collision` and `fig-boundary-time-delay`.
   - `eq-lec3-cj-from-chi` and the three intervals of $\chi$ for $a_2^{(1)}$ with their sign patterns $(C_0,C_1,C_2)$. Recompute the pattern in each interval directly, by evaluating the real-coupling boundary condition on `eq-lec3-vacuum-pair`. Is the index rule $C_j=\epsilon_{j+1}\epsilon_{j-1}$ well defined, given that $\operatorname{sgn}\cos(\frac\chi2+\frac{\pi ka}h)$ is not periodic in $k$ with period $h$ when $a$ is odd? Compare with `prp-solitonic-vacua` ($C_j=(-1)^{\vec\alpha_j\cdot\vec\lambda}$) and with the caption of `fig-lec3-mirror-images` ($\chi=0.25$, pattern $(-1,+1,-1)$). The drafter's script `code/lectures/halfline_l3_f1_mirror_images.py` picks, for each test value of $\chi$, whichever of the three patterns fits best; it does not test the stated assignment.
   - Is the positivity argument for $\tau_j(x)$ on $x<0$ (item 1 of Blackboard Calculation 2: "both roots … lie at $u\ge1$") complete?
   - The parity argument: is the constant $\frac i\beta\sum_j\vec\alpha_j\ln\zeta^{ja}$ a shift in $\frac{2\pi}\beta\Lambda_W^\vee$ for every branch of the logarithm, and does the conclusion "(a, θ, ξ) ↦ (ā, −θ, −ξ)" agree with `prp-antisoliton-return`?
   - `eq-lec3-linear-reflection`: check "$=-\big[(a)(h-a)\big]^{-C}$" at $C=0$ and $C=-1$, and the linearized boundary condition before it.
   - `exr-lec3-breathers` and `exr-lec3-self-conjugate`: are they true? In the second, check "exact energy $E=0$" and the zero mode against `prp-self-conjugate-vacuum`.

5. **Lecture 4, the coideal and boundary sine-Gordon.**
   - The derivation of the charges from boundary perturbation theory: the parity-image field (the text switches between $\tilde\Phi=\Phi+\bar\Phi(-x)$ and $\tilde\phi=\phi(x)+\phi(-x)$), the contact identity, the double-pole argument, the normalization of $T_j$, and the dictionary $\bar Q_j\mapsto q_j^{-1}f_jk_j$. Check the chiralities against `AGENTS.md` and `sec-boundary-charges`.
   - The display in "Why Boundary Crossing-Unitarity Fixes the Sign of $q$", $K_{ab}(\theta)=\sum S_{cd}^{ba}(2\theta)K_{\bar d\bar c}(i\pi-\theta)$: is it consistent with `eq-lec1-crossing-unitarity` and with the book's `eq-boundary-crossing-unitarity`?
   - "Under $q\mapsto-q$, the six-vertex bulk $\check R$-matrix is conjugated by $\Pi\otimes1$ … $K(\hat\epsilon;-q)=K(-\hat\epsilon;q)$. Therefore, the reflection equation is completely blind to the sign of $q$": check, given that $q\to-q$ is a Klein-factor twist. Is "holds if and only if $q=-e^{-i\pi\lambda}$" shown, or only that one sign fails? Check the claim about the Klein-twisted conjugation $\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)$. Compare with `sec-sg-k-from-symmetry`.
   - The one-dimensional submodule: the text contains an abandoned step ("$\det\pi_y(b_1)=\hat\epsilon_1^2-1$ wait: …"). Check $y_n=i(-1)^ne^{-i\eta}$, the simultaneous eigenvector, $\beta_{1,0}$, $\bar\eta=\pi(\lambda+1)-\eta$, and "$\vartheta\mapsto-\vartheta$ (which swaps $P^+\leftrightarrow P^-$)", against `eq-bsg-excited-k` and `prp-bsg-excited`. Give replacement text for the paragraph.
   - The 4-line proof that intertwiners solve the reflection equation: are the hypotheses of `thm-k-reflection` stated, and does "matching highest weights fixes the constant to 1" hold for $V_\mu\otimes V_\nu\to V_{\bar\mu}\otimes V_{\bar\nu}$?
   - `exr-lec4-neumann-eps`: check $E_n-E_0=\frac12m_{2n}$ and the mirror-image hint. What is the $n=0$ "pole" at $\theta=\frac{i\pi}2$?
   - The paragraph on node labels was rewritten in ed0585a after its own referee pass. Check only that it agrees with the book after that commit.

6. **Lecture 5, the $a_2^{(1)}$ intertwiners and the Serre relation.**
   - Are the two branches exhaustive? The text treats "$\hat\epsilon_j^2\ne\epsilon_*^2$ for all $j$" and "$\hat\epsilon_j^2=\epsilon_*^2$ for all $j$"; what about mixed cases? Compare with `prp-an-k-existence` (labelled "[Numerical for $n\le7$; Proposition, sketch in general]" in Appendix E). Check "only the product of signs $s$ is gauge-invariant".
   - Check the general-$N$ Toeplitz formula against `eq-an-k-matrix`, and the explicit $3\times3$ $K_s$. In `exr-lec5-toeplitz-solve`, part 1, do the node-0 equations alone give $K_{31}$?
   - Boundary unitarity in the form $K_s(-\theta)^TK_s(\theta)=U_s\mathbf 1$: is the transpose the reflection matrix $\bar V\to V$ in the book's conventions? Check $U_s$, and $\kappa(\zeta)$ in `exr-lec5-cyclic-unitarity`.
   - The Serre computation: check every case of $C(X,Y)$ and the result against `prp-serre-obstruction`, and `exr-lec5-serre-ej`. Is the argument that $k_i^2k_j\in B_{\hat\epsilon}$ forbids an intertwiner correct (what does conjugation by powers of $k_i^2k_j$ separate, and for which $q$)?
   - The callout "Quantum Symmetric Pairs vs. Generalized $q$-Onsager Algebras":
     - "the degree-4 quantum Serre polynomial is antisymmetric under word reversal, so its weight-zero terms cancel in pairs for any $(\hat\epsilon_0,\hat\epsilon_1)$" (compare `prp-q-dolan-grady`);
     - Kolb's parameter set (quote Kolb, arXiv:1207.6036, the definition of the parameter sets);
     - the universal K-matrix citations (do Balagović–Kolb and Appel–Vlaar cover this affine coideal?);
     - "$F_{ij}(b_i,b_j)=b_j$ … Baseilhac and Belliard's generalized $q$-Onsager algebra … obtained from $B_0$ by a one-dimensional character twist $(\mathrm{id}\otimes\chi)\circ\Delta$" (compare `prp-one-dim-modules` and Baseilhac–Belliard).
   - "Delius and MacKay made that very assumption", and "the overall normalization constant scales as $c\propto g^{1/2}$": compare with Delius–MacKay and `prp-dm-first-order`.
   - The two meanings of $s$: classical $C_0C_1C_2$ (Lecture 2) and quantum $\prod_j\hat\epsilon_j/\epsilon_*$ (Lectures 5–6). In "Resolving the Neumann vs. Uniform Boundary Puzzle", which classical class does each quantum solution belong to, according to the book, and is the lecture consistent with it? What does the book say about $s=-1$?

7. **Lecture 6, from solitons to particles.**
   - Check $A(\mu)$, its two equations, $K_{B_1}$, the continuation, $K_1$ and $K_a$ against `eq-an-diagonal-scalar-eqs`, `eq-an-diagonal-scalar`, `eq-breather-k-imaginary`, `eq-breather-continuation` and `eq-uniform-amplitude`.
   - The singlet argument: "On identical-colour states … $\check R$ acts as the identity", the order $(1\otimes K(\mu_a))\check R(2\mu)(1\otimes K(\mu_b))$, and "$K_{\bar V}=K_V^T$".
   - Check $c_+(\mu)$ (`prp-an-nondiagonal-scalar`), the equations for $A^{(+)}$, $K^{(+)}_{B_1}$ (`eq-an-neumann-breather`), $K_a^N$ (`eq-neumann-amplitude`), and the $s=-1$ formula with "the amplitude conjectured by Corrigan, Dorey, Rietdijk, and Sasaki for $C\equiv-1$".
   - The solitonic factor $M_a^{(b)}$ (`eq-solitonic-factor`), and the description "trapping a stationary imaginary-coupling soliton … in front of a uniform ($C\equiv+1$) wall", against Lecture 3 and `sec-solitonic-boundaries`.
   - The $a_2^{(1)}$ spectrum: $K_1^{b_{0,0}}$ and "poles inside the physical strip … (block labels $0<x<\frac32$)"; `eq-a2-excited` as quoted; $F(k)$; $e_{b_{n,-n}}-e_{b_{0,0}}=m_1\sin^2(n\delta)/\sin\delta$; the limit and $I=\frac6{\beta^2m_1}$. What do negative indices in $b_{n,m}$ mean, when "$n$ and $m$ count bound particles"?
   - The rotor callout: "Computing the classical kinetic moment of inertia $I=\int|\partial_\chi\vec\phi_\chi|^2dx$ of Lecture 3's exact vacuum solution gives $I=\frac6{\beta^2m_1}$ ([Section 36.6])". Compute the integral. Does the book contain it? Is "$\chi$ runs around the full $2\pi$ circle connecting all three quantum-equivalent solitonic patterns" consistent with Lecture 3, where different intervals of $\chi$ satisfy different boundary conditions?
   - The Coleman–Thun paragraph: the formula "order $=P-2L-Z$" (quote Dorey–Tateo–Watts); the window $\frac1{2n-1}<B<\frac2{2n-1}$; the pole "$p=\frac12+\frac{3B}2=\frac{47}{40}$" for $B=\frac9{20}$, and whether it lies in the physical strip in the block labels used; "Below $B=\frac13$, this pole creates the bound state $b_{-1,2}$" against the state $b_{-2,2}$ discussed; the residue $-1.1875$; $\theta_R$; and the count $P=5$, $L=2$, $Z=1$. Compare with `sec-boundary-spectrum-a2`.
   - The opening paragraphs: "only later were they recovered from the imaginary-coupling soliton S-matrices"; "In the early 1990s, several authors [@sasaki1993; @fring_koberle1994a; @kim1995] tried to fix the CDD ambiguity by assuming … self-dual"; and the claim about $\frac12B(2-B)$ and one loop. Quote the sources.

8. **Lecture 2, Sklyanin's construction.**
   - Derive $d\tau/dt$ from `eq-lec2-dt-t` with the stated transport convention. Is `eq-lec2-compatibility` the right condition, and are "$x_-\to-\infty$ … take $K_-=\mathbf 1$" and the claim about the classical reflection equation right?
   - Branch I and Branch II: do they exclude mixed patterns (some $C_j^2=1$, others not), and so prove that $a_2^{(1)}$ allows only $C\equiv0$ or all $C_j=\pm1$? Compare with `prp-durham`. Which one-parameter slice does the right panel of `fig-lec2-sklyanin-rigidity` show, and does the caption say so?
   - "for $d_4^{(1)}$ … *none* of the 32 integrable sign choices $C_j=\pm1$ satisfies $C_j\propto\sqrt{n_j}$": is $C_j=\pm1$ the integrable set for $d_4^{(1)}$ with this normalization (`tbl-durham`)? Check `exr-lec2-normalization`.
   - The Dirichlet argument ("a diagonal matrix … cannot intertwine … unless $\lambda^{2h}=1$") against `prp-no-dirichlet`.
   - The soliton-preserving section: `eq-lec2-sp-sum`, `eq-lec2-soliton-preserving`, the coweight extension, the image relation $\vec\phi(x,t)=2\vec v-\overline{\vec\phi(-x,t)}$, and the attribution of the generalization to arbitrary $\sigma$ to Delius (hep-th/9809140). Compare with `prp-soliton-preserving` and `prp-soliton-preserving-half`.
   - `exr-lec2-coweight-symmetries`: does conjugation by $D$ flip "the signs of $C_j$ and $\lambda$"? Does $K\mapsto K^{-1}$ preserve or flip $C_0C_1C_2$? Is the claimed partition into two classes compatible with the title ("generated from the single uniform solution")? Is "$s=+1$ (four solitonic sign patterns)" consistent with the book's definition of solitonic patterns?
   - `exr-lec2-a1-dirichlet`: check the dictionary against `tbl-boundary-dictionary` and the sign of $\phi_0$ in `eq-bsg-dirichlet`.

9. **Prior art and attributions.** The lectures present several results of the book's Appendix E ("Half-line") as discoveries of named authors, or as their own, and say in places who found what. For each attribution below, quote the source and say whether it is right, overstated or misattributed:
   - "uncovering a breathtaking weak–strong boundary duality [@corrigan1998; @delius_gandenberger1999]" (`index.qmd`) and "It was Corrigan, Gandenberger, and Delius and Gandenberger who discovered the true path" (Lecture 6), against the Appendix E items on the Neumann amplitudes and on the $s=-1$ duality;
   - "When Corrigan, Dorey, Rietdijk, and Sasaki and Bowcock, Corrigan, Dorey, and Rietdijk … tested this in 1994–1995, they discovered a shock" (`index.qmd`);
   - "The resolution, discovered by Mezincescu and Nepomechie [@mezincescu_nepomechie1998] for sine-Gordon at the free-fermion point and developed into a general theory by Delius and MacKay" (Lecture 4);
   - "(generalized to arbitrary involutive anti-automorphisms σ by Delius [@delius1998b])" (Lecture 2);
   - the value $\epsilon_*$ and the generalized $q$-Onsager relations (Lecture 5), against the Appendix E statement that the value "was first published by Baseilhac and Belliard";
   - "Following Delius [@delius1998a] and Bowcock [@bowcock1998]" (Lecture 3) and "[@ssw1995; @delius_gandenberger1999]" for mirror images (`index.qmd`);
   - "(and at $O(\beta^2)$ reproduces Kim's one-loop Neumann amplitudes [@kim1995])" (Lecture 6).
   Search for prior art on any result the lectures present without attribution.

10. **The other citations.** Check that each remaining citing sentence says what the source says, in particular: `gz1994` (eqs. (3.45), (3.51), and as the source of the 90° argument); `dptw1998`; `cdr1995` for the linearized reflection factor; `devega1994`; `gandenberger1999a` for $K^{(\pm)}$; `delius_gandenberger1999` for the shifted crossing form and the solitonic patterns; `bptt2002` for $(\eta,\vartheta)$; `dtw1999`; `bernard1991` for $U_q(\hat{\mathfrak g}^\vee)$; `cherednik1984` and `sklyanin1988`; `letzter1999` and `kolb2014` for "left coideal subalgebra"; `baseilhac_koizumi2005`.

11. **Exercises.** Work all 18 in outline. For each, say whether it is true, well posed and solvable with what the lectures provide, and whether the book object in its heading link is the corresponding one:
    - Lecture 1: `exr-lec1-elastic`, `exr-lec1-twisted`, `exr-lec1-crossing-forms`;
    - Lecture 2: `exr-lec2-normalization`, `exr-lec2-a1-dirichlet`, `exr-lec2-coweight-symmetries`;
    - Lecture 3: `exr-lec3-mirror-square`, `exr-lec3-breathers`, `exr-lec3-self-conjugate`;
    - Lecture 4: `exr-lec4-general-cj`, `exr-lec4-neumann-eps`, `exr-lec4-1d-submodule`;
    - Lecture 5: `exr-lec5-toeplitz-solve`, `exr-lec5-serre-ej`, `exr-lec5-cyclic-unitarity`;
    - Lecture 6: `exr-lec6-bootstrap-a2` (including the "doubling identity" and $S_{11}=\{1\}_B$), `exr-lec6-kim-oneloop` (is part 3 a well-posed question? is the expansion in part 2 right?), `exr-lec6-rotor-spectrum` (the binding rapidities, "$m$ particles of species 2 at the same rapidities", and the degeneracy at $B=\frac12$).

12. **Claims that a script checks something.** For each sentence in the lectures that says a script verifies a claim ("checked numerically in", "verified in", "verified to $10^{-14}$", "to machine precision", "Verified and plotted by"), read the script and say whether it checks that claim, in the conventions of the text, or something weaker. Examples: the braid form in Lecture 1; the panels of `fig-lec2-sklyanin-rigidity`; the vacuum patterns in Lecture 3 (see item 4); crossing-unitarity at both signs of $q$ in Lecture 4; $U_s$ and matrix crossing-unitarity in Lecture 5; $c_+(\mu)$ and the Coleman–Thun residue in Lecture 6.

## Report format (`briefs/referee_lectures_halfline_report.md`)

1. A header with the date, your model and the commit reviewed.
2. One section per item above, with the overlap note and the check.
3. **Overlap table.** One row per item above and per entry of your phase-1 risk list, with these columns:
   - the claim;
   - in the coordinator's list (item number, or no);
   - in your phase-1 list (entry number, or no);
   - verdict.
4. **What each list found.** List the errors (verdict "wrong") found only through your phase-1 list, only through this list, and through both. Say which kinds of error each list was better at finding.
5. The claims you could not check, and why.
6. Other observations, labelled as unverified.
