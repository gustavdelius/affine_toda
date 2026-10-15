# Part III (chapters 9–11): the two phase-1 referee reports, combined

- **Date:** 2026-10-15.
- **Compiled by:** Claude Opus 5.5, the coordinating session. This is not a referee report.
- **Sources.**
  - Gemini 4 Argon (high): `briefs/referee_classical_09_11_phase1_report.md` in the main checkout. Reviewed `a03edd4`. Scripts in `referee_scratch/classical_09_11/` (main checkout), no saved outputs.
  - Codex GPT-6 Astra (high): `briefs/referee_classical_09_11_codex_phase1_report.md` in the worktree `.claude/worktrees/referee-briefs-part-iii`. Reviewed `bcd82b9`. Scripts and outputs in `briefs/referee_classical_09_11_codex_scripts/`, with a frozen copy of its risk list and SHA-256 snapshots of the chapters.
- **Commits.** Between `bcd82b9` and `a03edd4` Part III changed only in script citations and one qualifier ("for simply-laced $\hat g$", 09:98, added in f7739bc). The substance reviewed is the same. Line numbers below are those of `main` (`a03edd4`).
- **Why this file is here.** It sits in `~/referee-held/`, not in `briefs/`, so that neither referee can read the other's findings before phase 2. Move it into `briefs/` after phase 2.
- **Status.** Nothing here has been integrated. Rule 1 still applies at integration. I spot-checked every point where the two reports disagree and every finding I mark "verified" below. My script is `~/referee-held/referee_classical_09_11_combined_spotchecks.py`.

Notation: **G3** is Gemini's rank 3, **C4** is Codex's item 4.

## 1. Findings, by severity

### Major

**M1. The Θ selection rule and its Appendix E entry** (11:43–49, 11:26, 11:96; `appendices/E-original-results.qmd:8`). G2, C3, C18.

- "Every conserved charge" is false. Codex's counterexample is the sine-Gordon kink. Its data $\{(1,\theta)\}$ are Θ-invariant, since $\bar a=a$, and its topological charge has parity $\epsilon=-1$. The rule would make the charge imaginary, but it is 1. *Verified by reasoning.* The book already says, at 11:26, that Θ flips the kink's topological charge. Codex only.
- The missing hypothesis is the same in both reports. The charge's value must depend only on $\{(a_i,\theta_i)\}$, not on the $\xi_i$. Additivity alone does not give this. Olive–Turok–Underwood state the property for $E$ and $P$: the contribution is "independent of both the complex number $Q_i$ and the sign $\epsilon_i$" (hep-th/9305160, sec. 8). Neither referee cites this sentence.
- The charges need a basis of definite Θ-parity. Both reports say this. Gemini adds that it matters where an exponent is degenerate, as in $d_{2n}^{(1)}$.
- In chapter 11 the rule is an unlabelled bullet, but Appendix E gives it **[Theorem, given additivity of the higher charges]** and cites a section. Gemini proposes a `prp-theta-selection-rule` div with a proof. That structure is right.
- At 11:26 the text conflates a Θ-invariant *field* with Θ-invariant *data*. Codex only.
- Θ-invariance gives reality, not positivity. The finite negative-energy solutions of Khastgir–Sasaki (hep-th/9507001, eq. 5.13) start from $\vec\phi=0$, $\partial_t\vec\phi=i\vec f$, with $E=-\frac12\int\vec f^2<0$. These data are Θ-invariant. *Source verified.* Gemini only. It bears on 10:180, 11:47 and the "physical sector" sentence at 11:49.
- **Prior art.** Codex found, unprompted, Cen–Fring (1602.05465), Cen–Correa–Fring (1608.01691) and Fring–Taira (2007.15425). It also distinguished Zhu–Caldi's "pairing condition", which is about the density. Gemini did not search the literature. One more item, not in either report: Takács–Watts (hep-th/9810006, p. 11, just before eq. 4.8) already state the reality condition for $a_2^{(1)}$ two-solitons. Real energy needs rapidities that are "both real, or … a complex conjugate pair".
- **Recommendation.** Make the rule a labelled statement. Take its structure and proof from Gemini (G2). Take its hypotheses from Codex (C3): definite parity; a value fixed by $\{(a_i,\theta_i)\}$; the regular domain of M2; topological charges excluded. Rewrite 11:96 and the Appendix E item once the phase-2 prior-art reading is done. Choose the label after that reading. Gemini's **[Theorem, …]** is defensible as a conditional statement only if the novelty survives.

**M2. `prp-hirota-energy` at complex rapidity** (11:32–41). G1, C2.

- The formula $E=\sum_iM_{a_i}\cosh\theta_i$ needs $\operatorname{Re}\cosh\theta_i>0$. It also needs a finite, nonzero leading product $\prod A$ and no zero of any $\tau_j$ on the slice.
- Codex's counterexample: sine-Gordon at $\theta=i\pi$, $\xi=i\pi/2$. This is a real, regular static antikink. Quadrature gives $E=8$; the formula gives $-8$. *Verified* (my quadrature: 8.0).
- Gemini's general rule is that each term carries a factor $\operatorname{sgn}(\operatorname{Re}\cosh\theta_i)$. At the endpoint $A=0$ the energy halves.
- The source states the condition. Olive–Turok–Underwood: "holds good even if $e^{\eta_i}$ is complex, providing that its real part is positive" (sec. 8.3). *Quote verified.*
- **Fix.** Add the hypotheses. Gemini's replacement is usable, but drop its "so $\operatorname{Re}E$ remains positive" aside.

**M3. Singular solutions: "integrated charges are unaffected"** (11:60, `exr-singular-points` 11:90–92). G7, C1.

- A double pole with zero residue does not make the real-line integral converge. For static sine-Gordon with $\tau_{0,1}=1\pm e^{2x}$, the density is $-8\operatorname{csch}^2 2x$. With $(-\epsilon,\epsilon)$ removed, the integral is $8(1-\coth2\epsilon)\approx-4/\epsilon$. *Verified* (−392.05 at ε=0.01, −3992.005 at 0.001). Only a contour or finite-part prescription gives a finite value.
- When a zero crosses the real line, the topological charge jumps by a coroot. Gemini only; *verified by reasoning*.
- **Fix.** Merge the two replacements. Use Codex's wording on divergence and the extra definition it needs. Use Gemini's sentence on the coroot jump. For the exercise, use Codex's explicit sine-Gordon version.

**M4. The breather band** (10:128–136, 11:23, caption 11:28, figure `ch11_f16`). G3, C4.

- A regular neutral (type B) breather needs $A_{a,h-a}(2iu)=1-\sin^2(\pi a/h)/\cos^2u<0$. That is $|\frac\pi2-\frac{\pi a}h|<u<\frac\pi2$, so the mass lies in $(0,2M_a\sin\frac{\pi a}h)$, not $(0,2M_a)$. *Verified:* the formula, the a₂ values, and Harder–Iskandar–McGhee §3.2, which says it is "not possible to have type B breathers with positive interaction coefficient".
- Codex also proves that no choice of phases removes the zero when $0<A<1$.
- The figure and its script must be redone.
- **Do not paste Gemini's type-A sentence.** Its formula $A_{aa}(2iu)=-\tan^2(\pi a/h)/\tan^2u<0$ for all $u$ is wrong. From `eq-interaction-coefficient`, $A_{aa}(2iu)=-\sin^2u/(\sin^2\frac{\pi a}h-\sin^2u)$, which changes sign at $u=\pi a/h$. For $a=1$, $h=4$, $u=1.2$ the book's formula gives $+2.36$ and Gemini's gives $-0.15$. HIM (3.5) also says the type-A coefficient changes sign at a critical velocity. The remark that the type-A family is omitted is still fair.

### Medium

**M5. Spin of `eq-spin3-current`** (09:96, `exr-spin3`). G8, C9.

- With $\partial_\pm=\partial_t\pm\partial_x$, $\partial_+f(x\cosh\theta-t\sinh\theta)=e^{-\theta}f'$. So the displayed current gives $Q\propto e^{-3\theta}$, i.e. spin −3. Both referees independently find $Q=-\frac{16m^3}{3\beta^2}e^{-3\theta}$ on the kink. *Verified* (the scaling).
- "Replace $\beta^2$ by $-\beta^2$": Gemini calls this false. That is overstated. $\cosh\beta\phi$ is even in β, so the substitution is consistent when read as acting on β². It is still worth making explicit, as both replacements do.
- Both note that the example's sinh-Gordon $m,\beta$ are not the root-normalized ones of `exr-lax`.

**M6. Fusion** (10:126, `exr-fusion`). G19, C6.

- The text needs $a+b<h$. At $a+b=h$ the limit is a vacuum and the pole is double (Codex). For $a+b>h$ the channel is species $a+b-h$. I checked that both referees' rapidity choices give $m_{a+b-h}$; Codex's lies in the physical strip.
- Gemini adds that the limit must be joint in $\theta_1-\theta_2$ and $\xi_{1,2}$.

**M7. The vertex-operator formula has no spacetime dependence** (10:72). C16 only; Gemini (G17) judged it correct.

- The displayed $\tau_j$ contains no $x$ or $t$. Chapter 6 (06:168) has the factor $W_k(x,t)$. Olive–Turok–Underwood obtain $g(t)$ by multiplying each $\hat F^i(z)$ by $W$ (eqs. 2.28–2.29). *Verified.*
- Gemini's defence says OTU eq. (6.3) absorbs the spacetime dependence into $F^a$. That misreads the source: (6.3) is the level-$x$ tensor-power identity used for nilpotency.

**M8. The area rule needs "nonzero"** (09:132–138). C7 only; Gemini (G11) judged it correct.

- In $a_2^{(1)}$, $C_{112}=0$ by $\mathbb Z_3$ charge, although the triangle exists. The displayed equality is therefore false as written. The book's own sentence at 09:142 ("Many triangles exist that carry no coupling") contradicts it. Dorey's sequel (hep-th/9110058) says "which, if non-vanishing, obey the 'area rule'".
- Codex also notes the "(with $m=1$)" normalization.
- Gemini's separate point is valid. For $d_{2m}^{(1)}$ the charge basis is real, so "in a real orthonormal basis … neither holds" needs "an arbitrary" or a restriction to conjugate pairs. The phrase "the complex eigenvectors of the centre" is also wrong there. *Verified for $d_4$:* the Z₂×Z₂ characters are ±1. The nonzero couplings are $HHH$, $HL_iL_i$ (three) and $L_1L_2L_3$, as Codex lists.
- Gemini's printed $d_4$ couplings are mislabelled. It lists $C_{1,1,s'}$ and $C_{s,s,s'}$ as nonzero, which Z₂×Z₂ forbids.

**M9. Lax-pair generators for twisted algebras** (09:57). C5 only.

- "For $j\ge1$ … Chevalley generators of $g$; $e_0,f_0$ are $E^{\mp\vec\theta}$ times the loop variable" describes only the untwisted case. In $a_2^{(2)}$, $\vec\alpha_0=-2\vec\alpha_1$ is not a root of the rank-1 finite algebra. *Verified by reasoning.* The curvature computation itself is fine; both referees checked it.

**M10. Degrees of folded τ-functions** (10:77). G6 only.

- `tbl-foldings` includes $a_{2n}^{(2)}$, and the caption defines $n_j^\vee=k\,n_j\vec\alpha_j^2/2$. With $n=(1,2,\dots,2)$ and $\vec\alpha^2=(2,1,\dots,1,\frac12)$ this gives $n^\vee=(2,\dots,2,1)$.
- A coincident pair from $a_{2n}^{(1)}$ has degree 2 at every node, but $|O|n_j^\vee=4$ at nodes $0,\dots,n-1$. *Verified from chapter 6's definitions.*
- **Fix.** Exclude $a_{2n}^{(2)}$ (chapter 6 already calls it "an exception to both statements"), or state the degree in parent labels.
- A smaller point: $\tau_1=\tau_3$ uses parent indexing.

**M11. Time advance sign** (10:122; caption of `fig-two-soliton-shift`). C8; Gemini (G13) judged it correct.

- "$\Delta t=-\Delta x/v_1<0$" and "shifted forward" need $v_1>0$. So does the caption's "each soliton is displaced forward along its own direction of motion", which needs $v_1>0>v_2$. The figure's ±0.6 satisfies this. *Verified by reasoning.*
- Gemini notes the clash of $E_1$ as an exponential and as an energy in the same paragraph.

### Minor (one referee each; plausible; not re-derived by me unless stated)

- 09:30: "makes $\vec\phi=0$ a minimum" for a complex potential. Say "the real part, on real fields" (C17).
- 10:59: the position at $t=0$ is $-\operatorname{Re}\xi/(m_a\cosh\theta)$, not $-\operatorname{Re}\xi/m_a$ (C, unranked). *Verified by reasoning.*
- 10:174: "and the energies are real" needs "for regular real-rapidity solutions" (C13).
- 10:180, 11:47: describe Khastgir–Sasaki correctly. Their negative-energy solutions are not Hirota solitons (G1).
- 11:49: "closes gap H6". Codex proposes "a sufficient reality criterion … does not select a quantum sector".
- 11:72: "no static saddle" should read "no regular static solution approaching vacua at both ends" (C19). Codex proved this independently, with a symbolic static three-species family for $a_3^{(1)}$.
- 09:88, the twisted spin lists. Gemini (G10) reports a centralizer computation; Codex (C15) leaves it unclear. Not re-run.

## 2. Claims both referees checked and found correct

- Topological charges and the count $h/\gcd(a,h)$ (G9, C10; Codex enumerated $h\le14$).
- The $d_n$ and $e_8$ mass tables and the sum rule (G14, C20).
- Dorey's rule, $\vec\gamma_a$ and the $d_4$ caption (G12, C14).
- The stress tensor as a total derivative, pointwise (G4, C13).
- Lie duality for $c_n$ (G5, C11).
- The $c_2$ τ-functions (G6, C16).
- The zero-curvature computation and the boost (G16, C5).
- The semiclassical sine-Gordon phase (G13, C8; Codex derives the cosine transform).

Single-referee "correct" verdicts:

- Gemini: the embedded sine-Gordon kink (G15) and the sign of $q$ (G18). The Takács–Watts citations (G20) also check out; I confirmed in the source that they cover complex soliton–breather time delays and missing excited solitons.
- Gemini also built the $a_2^{(1)}$ spin-2 current and integrated it on an explicit breather. The result matches $2iq_1e^{2\theta_0}\sin2u$ to $5\times10^{-14}$ without assuming additivity. This is the check on which the Appendix E label rests.
- Codex: the Hirota equations (C12), the vacua (C17) and the $a_3$ static obstruction (C19).

## 3. Relation to the held phase-2 list

Compared with `referee_classical_09_11_items.md` (SHA-256 verified).

- **Found blind:**
  - the Fring-group prior art (item 1), by Codex;
  - the spin-2 additivity test on a breather (item 3), by Gemini;
  - the Khastgir–Sasaki question (item 1), by Gemini;
  - the complex-rapidity energy (item 2), by both;
  - the spin-3 value (item 5), by both;
  - the area rule and degenerate basis (item 6), by Codex and Gemini respectively;
  - the twisted Lax pair (item 5), by Codex.
- **Found blind but not targeted by the held list:**
  - the breather band;
  - the divergence of the singular integrals (the list asks only about the residue);
  - the topological-charge counterexample;
  - the missing $W_k$;
  - the $a_{2n}^{(2)}$ degree;
  - the fusion range (reached only through the exercise list);
  - the frame dependence of the time advance.
- **Still open for phase 2:**
  - the Wigner relation for $a_n^{(1)}$, $n\ge2$, against Hollowood's S-matrix (item 4);
  - $i\vec\phi$ against $-i\vec\phi$ (item 8);
  - H6;
  - "$C$ is not local";
  - Drinfeld–Sokolov;
  - the attributions to `braden1990`, `freeman1991` and `fring1991`;
  - `mikhailov1981` and `olive1985`;
  - all 17 exercises.

## 4. Judgement: Codex did the better job

**Recall is similar. Precision and evidence are not.**

- **Shared core.** Both found the four major errors (M1–M4) and the spin and fusion errors (M5, M6). These are the most consequential findings, and the overlap is reassuring.
- **Codex's extra finds are confirmed errors.** The missing $W_k$, the area-rule qualifier, the twisted Lax generators, the topological-charge counterexample and the frame dependence of the time advance are all real. Gemini marked four of these five "correct".
- **Gemini's extra finds are real but smaller.** The Θ-invariant negative-energy data, the coroot jump, the $a_{2n}^{(2)}$ degree and the $d_{2m}$ basis are all valid. So are the call for a labelled statement and the spin-2 breather integral. The spin-2 integral is the most useful single computation in either report.
- **The prior-art test.** The brief named originality claims as fair game. Codex searched, found three of the four papers on the held list, and declined to certify novelty. Gemini did not put originality on its risk list at all. This was the experiment's main test.
- **Errors made by the referee.** I found none in Codex's report. In Gemini's I found three: the type-A breather formula, which it would have pasted into the book; the misreading of OTU eq. (6.3) behind its "correct" verdict on the vertex formula; and the mislabelled $d_4$ coupling table. It also overstated the β² → −β² point and proposed **[Theorem]** without weighing prior art.
- **Evidence and honesty.** Codex quotes exact source lines, keeps a frozen risk list and chapter hashes, saves outputs, and marks two items "unclear" with reasons. Gemini mostly cites line ranges without quotes. It saved no outputs. It wrote "None. All 20 entries were checked in full", although several checks (G5's twisted cases, G17, G18) rest on reading rather than computation.
- **Replacement text.** Gemini's is more paste-ready and closer to the book's markup. Codex's is terse, uses backtick math and needs rewriting in the book's style. This is Gemini's clear advantage.
- **Process.** Codex worked in the worktree named in the brief. Gemini ran in the main checkout, which is harmless here.

**Verdict.** Codex is the more reliable referee. Its findings can be integrated after the usual Rule 1 checks. Gemini found real problems that Codex missed, and its report is worth integrating as a second source. Treat its derivations and replacement texts as pointers, and check each before use. This matches the Part II calibration.

## 5. Phase-1 reports frozen before phase 2

SHA-256 of each phase-1 report, taken before the phase-2 items were copied into `briefs/` (2026-10-15):

- Gemini, `briefs/referee_classical_09_11_phase1_report.md` (main checkout): `1a4da4425316699dfd04ee1505fe60cd35125af49b59cb8ffaa6fbd46d0396b3`
- Codex, `briefs/referee_classical_09_11_codex_phase1_report.md` (worktree): `26c44d6352239cde46bb1ae60c9f66567a41fe7b79b752a4993b68bdc4faa24b`

The items file was then copied into `briefs/` in both checkouts and its hash checked against `briefs/README.md`. Commit the phase-1 reports, the items file and these hashes on the integration branch.
