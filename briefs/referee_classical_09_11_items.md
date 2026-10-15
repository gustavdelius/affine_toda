# Phase 2 of the referee pass on Part III: the coordinator's items

This continues `briefs/referee_classical_09_11.md`. Your phase-1 report, `briefs/referee_classical_09_11_phase1_report.md`, is complete and committed; do not change it.

The rules, notation, sources and standard of evidence of the phase-1 brief still apply, with one change: the blinding rule is lifted. You may now read the earlier pass on these chapters: `briefs/referee_textbook_report.md` (section "Item 13", and "Item 12" for the exercises) and `briefs/part3_reading_fixes.md`. Do not take their verdicts as evidence.

Write this phase to `briefs/referee_classical_09_11_report.md`.

## For each item below

- **Overlap.** Was the item, or part of it, in your phase-1 risk list? Give the entry number and rank, and whether you checked it.
- **Check.** If your phase-1 check already covers the item fully, say so and point to the section; do not repeat it. Otherwise check it now: verdict (correct / wrong / unclear), your own computation and its output, and, for errors, replacement text ready to paste, with file and line.

The items were written before your phase 1 started. Their SHA-256 was recorded in `briefs/README.md` at that time.

## Items (most important first)

1. **The $\Theta$ selection rule and its entry in Appendix E** (`sec-reality-energy-momentum`, bullet "Selection rule"; first item of `appendices/E-original-results.qmd`).
   - Appendix E labels the rule **[Theorem, given additivity of the higher charges]**, but in chapter 11 it is an unlabelled bullet. Given additivity, is there a complete proof in the text? Write the proof out. It must show:
     - that a charge depends on the data only through $\{(a_i,\theta_i)\}$, not on the $\xi_i$;
     - that every local conserved current has a definite parity $\epsilon_Q=\pm1$ under $\vec\phi\mapsto-\bar{\vec\phi}$.
     
     Test the second point on `eq-spin3-current` and on spin-2 and spin-4 currents of $a_2^{(1)}$.
   - Should the rule be a labelled statement, so that `check-refs.py` can match the label in Appendix E against it? At present the item cites a section.
   - Check "This closes gap H6 of @sec-critique, apart from the question of singularities" against H6 in `part9-foundations/44-twelve-gaps.qmd`.
   - Check "The complex-energy solutions of Khastgir and Sasaki have asymptotic data that are not $\Theta$-invariant". Construct one of their solutions (hep-th/9507001) and read off its data.
   - **Prior art.** The notes say "The formulation as a selection rule from $\Theta$ is new here". These papers use PT symmetry to explain the real energies of complex solitons, and none of them is cited in the book:
     - Cen and Fring, "Complex solitons with real energies", arXiv:1602.05465;
     - Cen, Correa and Fring, "Time-delay and reality conditions for complex solitons", arXiv:1608.01691;
     - Fring and Taira, "Complex BPS solitons with real energies from duality", arXiv:2007.15425;
     - Correa, Fring and Taira, "Linearly stable and unstable complex soliton solutions with real energies in the Bullough–Dodd model", arXiv:2110.06825. The Bullough–Dodd model is $a_2^{(2)}$ Toda theory.
     
     Read them, and search for further work. Say exactly what they prove, and what, if anything, of the book's rule is new. Give replacement text for the notes of chapter 11 and for the Appendix E item.

2. **`prp-hirota-energy` and the energy–momentum tensor.**
   - Claim: $E=\sum_iM_{a_i}\cosh\theta_i$ and $P=\sum_iM_{a_i}\sinh\theta_i$, for complex rapidities, cited as "[@tw1999, eq. 4.6]". Quote the source.
   - Check numerically on an $a_2^{(1)}$ breather and on a single soliton with complex rapidity. Integrate the energy density of the actual field; do not use the asymptotic formula.
   - Is the label **[Established; reformulated]** right? What exactly is established, and by whom?
   - The notes (`sec-classical-literature`, paragraph "Reality") give $T_{\mu\nu}=(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$, with $C=-\frac2{\beta^2}\sum_j\frac2{\vec\alpha_j^2}\ln\tau_j$, in signature $(+,-)$. Check the sign and the normalization against the canonical $T_{\mu\nu}$ on a Hirota solution.
   - Check "$C$ involves $\ln\tau_0$, which is not a local function of $\vec\phi$", and "for the trace, $T^\mu{}_\mu\propto V$, it follows directly from the bilinear equations". Quote Olive–Turok–Underwood or Zhu–Caldi for the formula.

3. **The additivity assumption, tested** (`sec-reality-energy-momentum`, bullet "A weaker condition for $E$ and $P$").
   - Claim: the spin-2 charge of a single $a_2^{(1)}$ soliton is $q_ae^{2\theta}$, where $(q_1,q_2)$ is the eigenvector of the $A_2$ Cartan matrix with eigenvalue 3 (Freeman), so $q_2=-q_1$. Additivity then gives $2iq_1e^{2\theta_0}\sin2u$ on a breather.
   - Construct the spin-2 conserved current of $a_2^{(1)}$ yourself. Integrate it on an explicit breather solution, without assuming additivity. Does the result equal $2iq_1e^{2\theta_0}\sin2u$? The label in Appendix E rests on this assumption.
   - Quote Freeman (hep-th/9408092) on single solitons.

4. **Time delays and the semiclassical phase** (`sec-multi-solitons`, `fig-two-soliton-shift`).
   - Check the signs from the asymptotics of `eq-n-soliton`: $\Delta x>0$ (forward), $\Delta t<0$, and $E_1\Delta x=-p_1\Delta t=-(2h/\beta^2)\ln A$. Do this for both solitons; the caption says each is displaced forward along its own direction of motion.
   - Claim: Wigner's relation gives $d\delta_{ab}/d\theta=(2h/\beta^2)\ln A_{ab}(\theta)$ at leading order, for every $a_n^{(1)}$. The earlier pass verified this only for sine-Gordon and left $n\ge2$ "unclear".
     - Expand the $a_n^{(1)}$ soliton S-matrix to leading order in $1/\beta^2$ and compare, for $a_2$ and $a_3$, with $a=b$ and $a\ne b$. Use Hollowood's (`hollowood1993b`, "Quantizing $SL(N)$ solitons and the Hecke algebra"), or the book's worked example `sec-an-example` (ch. 21).
     - Check the derivation of $d\delta/d\theta=p_1\Delta t$ from $\Delta t=d\delta/dE_1$.
   - Sine-Gordon: is the leading term of $-i\,d\ln S_0/d\theta$, from `eq-sg-s0`, equal to $(16/\beta_{\rm SG}^2)\ln\tanh\frac\theta2$?

5. **The Lax pair, boosts and the spins of the charges** (`sec-lax-pair`, `sec-conserved-charges`).
   - `eq-lax-pair` with $c_j^2=n_j\vec\alpha_j^2/2$ and `eq-zero-curvature` should be equivalent to `eq-toda-eom`. Check this with explicit matrices for one twisted algebra, for example $a_2^{(2)}$.
   - Check the boost claim $\mathcal A_\pm[\vec\phi_\vartheta](x;\lambda)=e^{\mp\vartheta}\mathcal A_\pm[\vec\phi](x';e^\vartheta\lambda)$, with "no gauge transformation". Check also that another gradation needs a constant gauge transformation.
   - `eq-spin3-current`:
     - Derive it.
     - Is $Q_3=\int(T_4+\mathcal X_2)$ the charge it conserves?
     - Is it consistent with the form "$\partial_-T_{s+1}=\partial_+\mathcal X_{s-1}$" stated after `eq-toda-eom`?
     - Check its value on a one-particle state.
   - Check the spin sets by computing the centralizer of $E_1$:
     - "$a_{2n}^{(2)}$: the odd integers not divisible by $2n+1$" (do $a_4^{(2)}$);
     - "the two members of a dual pair have the same spins" (do $(c_2^{(1)},d_3^{(2)})$ and $(g_2^{(1)},d_4^{(3)})$);
     - "$a_n^{(1)}$, $n\ge2$: every integer not divisible by $h$".
   - Check the Drinfeld–Sokolov description, and the citations for quantum conservation (`feigin1993`, `dgz1992`, `niedermaier1994`).

6. **Masses, couplings and Dorey's rule** (`sec-masses-couplings`, `sec-dorey-rule`, `fig-dorey-d4`).
   - Check the closed forms for the $e_8^{(1)}$ masses, the $d_n^{(1)}$ and $c_n^{(1)}$ masses, and the sum rule.
   - Check that every cubic coupling of $V_{\rm imag}$ is $i$ times the real-coupling one.
   - Check the charge-basis prescription for degenerate masses, and the area rule $\lvert C_{abc}\rvert=\frac{4\beta}{\sqrt h}\Delta_{abc}$. Is the normalization of $\Delta_{abc}$ ("with $m=1$") consistent with the mass table? Do $a_4^{(1)}$ and $e_6^{(1)}$ in the charge basis, and show that the rule fails in a real orthonormal basis.
   - Check the derivation of $\vec\gamma_a$ for white and black nodes.
   - Check the caption of `fig-dorey-d4`: which particles fuse, with which angles?
   - Is the label **[Established]** of `prp-dorey-rule` right?
   - Check the attributions to `braden1990`, `freeman1991`, `dorey1991`, `dorey1992` and `fring1991`.

7. **Solitons, charges and Lie duality** (`sec-hirota`, `sec-vertex-solitons`, `sec-topological-charges`, `sec-soliton-masses`).
   - Check `eq-hirota-bilinear` for non-simply-laced and twisted algebras (the factor $n_j\vec\alpha_j^2/2$), and the remark on $a_1^{(1)}$.
   - Check the "embedded sine-Gordon kink": the species-$n$ soliton of $a_{2n-1}^{(1)}$ at $\operatorname{Im}\xi=\pi/2$, with $\tau_{2k+1}=\bar\tau_{2k}$. Is it real? Do $n=2$ and $n=3$.
   - Check the $c_2^{(1)}$ $\tau$-functions printed in `sec-vertex-solitons`.
   - Check the degree statements: $n_j^\vee$ for folded theories, and a coincident orbit having degree $\lvert O\rvert n_j^\vee$ (test $c_2^{(1)}$ and $g_2^{(1)}$). Check the claim that in Part III only these statements depend on @hyp-D.
   - Check the formula for $\vec\phi(+\infty)$, the count $h/\gcd(a,h)$ and the caption of `fig-topological-charges`. Does McGhee (hep-th/9307035) give the count $h/\gcd(a,h)$?
   - Check `eq-soliton-mass`, and the hint of `exr-soliton-mass`: the first integral, and $V=\frac1{\beta^2}\sum_j\frac2{\vec\alpha_j^2}\partial_x^2\ln\tau_j$.
   - Lie duality:
     - Check the $c_n^{(1)}$ folding argument.
     - Check the claim that for the other untwisted non-simply-laced algebras Lie duality "is quoted from [@kneipp1996]". Kneipp and Olive treat twisted algebras; do they state it?
   - Check the two one-sentence pointers to Part VIII: `prp-stored-charge` ("moves a soliton's charge along the same cycle") and `prp-bianchi` ("$A_{11}$ is the square of the inverse of its coefficient").

8. **Chapter 11, other statements.**
   - "A $\Theta$-invariant field … is purely imaginary, and $i\vec\phi$ solves the real-coupling equation." The earlier pass proposed $-i\vec\phi$. Which is right?
   - Check `eq-theta-data` and `exr-theta-interaction`.
   - Singular solutions:
     - Check the double pole with zero residue.
     - Check "the Hirota data continue analytically through the singularity".
     - The text reads Takács–Watts as "the classical image of a quantum sector in which $\Theta$ is broken". Is that supported, or is it interpretation stated as fact?
   - Check the missing-charge section against `sec-which-saddles-contribute` and `cnj-missing-charge` (ch. 17).
   - Check the claim of `exr-breather-real` that the real sine-Gordon breather has $\xi_2=\bar\xi_1+i\pi$.

9. **Citations and exercises.**
   - Check every citing sentence. In particular:
     - `tw1999` eq. 4.6;
     - `mcghee1994` and `harder1995`;
     - `fring1994` and `hollowood1992` for $A_{ab}$;
     - `khastgir1996`: "solutions with complex energy that become singular in finite time, and solutions with negative energy";
     - `zhu1995` and `freeman1995`;
     - `mikhailov1981`;
     - `olive1985` for the spins;
     - `kac1990` for the exponents of $X_N^{(k)}$ modulo $kh$.
   - Work all 17 exercises in outline, and say whether each is true and solvable with what the chapters provide:
     - chapter 9: `exr-lax`, `exr-spin3`, `exr-dn-masses`, `exr-cn-masses`, `exr-shg-cubic`, `exr-dorey-a2`;
     - chapter 10: `exr-bilinear`, `exr-moving-soliton`, `exr-soliton-mass`, `exr-charges-a3`, `exr-lie-duality`, `exr-fusion`, `exr-breather-sg`;
     - chapter 11: `exr-theta-solutions`, `exr-theta-interaction`, `exr-breather-real`, `exr-singular-points`.

## Report format (`briefs/referee_classical_09_11_report.md`)

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
