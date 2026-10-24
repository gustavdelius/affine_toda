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

## Integration

- **Date:** 2026-10-24 (started 2026-10-15; both integration agents were stopped by usage limits and restarted from their checkpoints, with every check redone).
- **Integrated by:** Claude Opus 5.5. The coordinating session ran two integration agents one after the other, A (chapters 9–10, branch `part3-int-ch9-10`, commit 81a93cd) and B (chapter 11, Appendix E, `references.bib`, H6, branch `part3-int-ch11`, commit 89a2cf8), then merged both into `part3-referee-integration` and did the checks below. Procedure and prompts: `briefs/part3_integration_coordinator.md`, `briefs/part3_integration_A_prompt.md`, `briefs/part3_integration_B_prompt.md`. Checkpoints: `briefs/part3_integration_A_checkpoint.md`, `briefs/part3_integration_B_checkpoint.md`.
- **Base:** `origin/main` b7087eb, merged into the branch before the agents started (64f2cdb). Main had already fixed the spin of the sinh-Gordon example and the sign of the conservation law in chapter 9 (0c1e204, 2fa0691), independently of these reports.
- **Agents' own checks (Rule 1):** `briefs/referee_classical_09_11_integration_scripts/integration_A_*.py` and `integration_B_*.py`, with outputs. They import no referee script. The book's new checks are in `code/toda_classical.py` (extended) and `code/classical_reality.py` (new).
- **Status labels.** No label was upgraded. Changes, none of them yet seen by a referee in their new form:
  - `prp-hirota-energy`: kept **[Established; reformulated]**, with a narrower statement (hypotheses (i)–(iii)).
  - `prp-theta-reality`: **new, [Established; reformulated]**, with a proof. It replaces the unlabelled "selection rule" bullet.
  - `prp-breather-charges`: **new, [Numerical]**.
  - Appendix E: the item **[Theorem, given additivity of the higher charges]** is removed; a **[Numerical]** item citing `prp-breather-charges` replaces it.

### The coordinator's checks of the agents' work

- **New and changed statements, re-derived by hand.**
  - `prp-hirota-energy`: the boundary evaluation under (i)–(iii); the relabelling $\tau_j(a,\theta+i\pi,\xi)=\omega^{ja}e^{-X}\tau_j(\bar a,\theta,-\xi)$, where the constant factor $\omega^{ja}$ shifts $\vec\phi$ by a period ($\vec\alpha_k\cdot\sum_jj\vec\alpha_j\in\{0,\pm h\}$); $A_{ab}(\theta+i\pi)A_{\bar ab}(\theta)=1$.
  - `prp-theta-reality`: the one-line proof, the parity of the $E$, $P$ densities, and the existence of a basis of definite parity (the map $Q\mapsto\overline{Q[\Theta\,\cdot\,]}$ is an antilinear involution on local charges). The label **[Established; reformulated]** fits: the proof is one line, and the mechanism is prior art.
  - Fields fixed by $\Theta$: with $\vec\phi=-i\vec\chi$, $\vec\chi$ solves the real-coupling equation and $E_{\rm imag}[\vec\phi]=-E_{\rm real}[\vec\chi]\le0$, since $e^y-1\ge y$ and $\sum_jn_j\vec\alpha_j=0$.
  - `prp-breather-charges`: recomputed with my own script, `integration_coord_breather_charges.py` (in the same directory), which builds the $a_2^{(1)}$ currents in the variables $y_k=\vec\alpha_k\cdot\vec\psi$ and evaluates them through $\partial_\pm E_i=\lambda^\pm_iE_i$ on the $\tau$-functions. Result: **confirmed.** With $q_1$ the charge of a species-1 soliton, $q_2/q_1=-1$ at spins 2 and 4, a soliton at complex rapidity gives $q_1e^{s\theta}$, and three breathers ($(\theta_0,u)=(0.3,0.8),(-0.5,1.1),(0.1,1.4)$, i.e. $A=-0.55,-2.6,-25$, with different $\xi_1$ and times) give $Q_s/q_1=2ie^{s\theta_0}\sin su$ to $10^{-11}$ for $s=2,4$. At spin 3 every solution of the current equations has zero charge (total derivatives). The field-equation residual of the $\tau$-functions is below $10^{-14}$. Output in `integration_coord_breather_charges_output.txt`; about 15 minutes, most of it the exact spin-4 nullspace.
- **Key edits, spot-checked by hand.**
  - Chapter 9: the kink value $Q_{-3}=-\frac{16m_0^3}{3\beta_{\rm SG}^2}e^{-3\theta}$ (the two integrals $\int(f''^2-\frac{\beta^2}4f'^4)ds$ and $\int\cos(\beta f)f'^2ds$ both equal $-\frac8{3\beta^2}$, which is what makes the charge proportional to $e^{-3\theta}$); the dictionary $m_0=2m$, $\beta_{\rm SG}=\sqrt2\beta$; the Feigin–Frenkel scope (hep-th/9310022, line 154 of the text: "are odd, and the Coxeter number is even").
  - Chapter 10: $A_{a\bar a}(2iu)=1-\sin^2\frac{\pi a}h/\cos^2u$; the regular family $\xi_2=\bar\xi_1+i\pi$; the ellipse proof of the zeros for $0<A<1$; the three fusion cases; the signs of $\Delta x_1$ and $\Delta x_2$ from the asymptotics; the gauge shift of $C$ by $-2hf/\beta^2$; Takács–Watts on $\sqrt3M$ against $2M$ (hep-th/9810006, line 718 of the text).
  - The breather band in chapters 10 and 11 is the same statement, in the same symbols. The domain of `prp-hirota-energy` covers its use in chapter 10 ($\theta_0\pm iu$, $0<u<\frac\pi2$). The Appendix E item carries the label of the statement it cites (`check-refs.py`).
- **Merge.** Both branches merged without conflict (e7e8d9c, 06b609d). `code/README.md` got the agents' rows.
- **Final checks on the merged branch:** `code/classical_reality.py` 20 PASS; `code/toda_classical.py` 24 PASS; `code/algebra_data.py` 6 PASS; `code/smatrix_simply_laced.py` 11 PASS; the figure scripts `ch9_f13_dorey_d4.py`, `ch10_f14_topological_charges.py`, `ch10_f15_two_soliton_shift.py`, `ch11_f16_breather_theta.py` and `ch44_f33_gap_map.py` all PASS (96 PASS, 0 FAIL in all; regenerated SVGs whose diff was only ids and dates restored); `python3 check-refs.py` passes (one older WARN on `citation.qmd`); `./render-scratch.sh` renders the HTML with no unresolved cross-reference (one older TeX warning, from a title in `references.bib`).

### Decisions for the editor

1. **The label of `prp-theta-reality`.** [Established; reformulated] (the coordinator's default), not [Theorem]. *Decided by the editor on 2026-10-24: kept.* The antilinear-symmetry mechanism is in Cen–Fring, Cen–Correa–Fring and Fring–Taira, and the reality of $E$, $P$ for conjugate pairs is in Olive–Turok–Underwood §8.3 and Takács–Watts §4.2.
2. **Appendix E.** The [Theorem, …] item is replaced by a [Numerical] item for `prp-breather-charges`, the evaluation of the $\Theta$-odd charges of spins 2 and 4 on $a_2^{(1)}$ breathers, which none of the papers read contains. The alternative was to drop the item. *Decided by the editor on 2026-10-24: the [Numerical] item is kept.*
3. **The title of chapter 11** is now "The classical theory: reality from $\Theta$" (was "a selection rule from $\Theta$").
4. **Chapter 21 (not edited).** Agent A observed that the product $c\,f$ of `sec-scalar-factor` builds in self-conjugate crossing, $f(T-t)=f(t)$, so for $a_n^{(1)}$, $n\ge2$, it does not reproduce Hollowood's scalar factor (scaled leading phase derivative $-0.85$ against $\ln A_{11}=-1.93$ at $h=3$). `sec-an-example` defers to Hollowood and is consistent with chapter 10; only the sentence "the scalar factor can place its zero at $t=\omega$ itself" invites the wrong reading. Not verified by the coordinator; a separate task.
5. **Optional.** A sentence in chapter 10 saying that Hollowood's $a_2^{(1)}$, $a_3^{(1)}$ amplitudes have the leading term $(2h/\beta^2)\ln A_{11}$ (evidence `integration_A_wigner.py`; a book check would need a new script in `code/`).
6. **Not recommended:** Appendix D entries for Takács–Watts eq. (4.6) ("θ ∈ C" without the restriction, which they impose in §4.2) or for the sign in Hollowood's §5, which disagrees with his eq. (4.5); an Appendix E line for the negative energy of fields fixed by $\Theta$ (a two-line computation).
7. **Small leftovers.** Chapter 38, line 250 writes "$\beta^2\to-\beta^2$" for the spin-3 current; it is correct, and "$\beta\to i\beta$" would match chapter 9. Chapter 39, line 145 ("respects the physical sector") was left unchanged.

### The two-phase experiment: which errors each phase and each referee found

Rows are the errors and overclaims that the integration verified and corrected (Rule 1). G*n* is Gemini's phase-1 rank *n*, C*n* Codex's phase-1 item *n*, GP2-*n* Gemini's phase-2 item *n*, C2§*n* Codex's phase-2 section *n*. "Held list" says whether the coordinator's phase-2 items (`briefs/referee_classical_09_11_items.md`) targeted the error. A dash means not found.

| # | Error (as corrected) | Ch. | Gemini P1 | Codex P1 | Gemini P2 | Codex P2 | Held list |
|---|---|---|---|---|---|---|---|
| 1 | Θ rule stated for "every conserved charge": fails for the topological charge (kink) and for charges without definite parity | 11 | G2 (parity) | C3 (kink, parity) | GP2-1 | C2§1 | item 1 |
| 2 | Θ rule: missing hypothesis that the charge depends only on the reduced data | 11 | G2 | C3 | GP2-1 | C2§1 | items 1, 3 |
| 3 | Prior art for the Θ mechanism; "new here" and the [Theorem] item in Appendix E | 11, E | — | C18 (3 of 4 papers, unprompted) | GP2-1 (papers supplied) | C2§1 | item 1 |
| 4 | Θ gives reality, not positivity: Khastgir–Sasaki's negative-energy data are Θ-invariant | 10, 11 | G1 | — | GP2-1 | C2§1, §9 | item 1 |
| 5 | Domain of `prp-hirota-energy` ($\operatorname{Re}\cosh\theta_i>0$, $\prod A\ne0,\infty$, no zero) | 11 | G1 | C2 | GP2-2 | C2§2 | item 2 |
| 6 | Singular solutions: zero residue does not make the integral converge | 11 | G7 | C1 | GP2-8 | C2§8 | item 8 (residue only) |
| 7 | Topological charge jumps by a coroot when a zero crosses the real line | 11 | G7 | — | — | — | — |
| 8 | Zeros of real-rapidity solitons are lines, not isolated points | 11 | G7 | — | — | — | — |
| 9 | Neutral-breather band $(0,2M_a\sin\frac{\pi a}h)$, not $(0,2M_a)$ | 10, 11 | G3 | C4 | GP2-8 | C2§8 | — |
| 10 | Spin of the sinh-Gordon example is $-3$ (fixed on main before integration) | 9 | G8 | C9 | GP2-5 | C2§5 | item 5 |
| 11 | Fusion needs $a+b<h$; channels $a+b>h$ and $a+b=h$; joint limit | 10 | G19 | C6 | — | C2§9 (exercise) | — (exercise list only) |
| 12 | Vertex-operator formula without $W_k(x,t)$ | 10 | — (G17 judged it correct) | C16 | — | C2§7 | — |
| 13 | Area rule only for nonzero couplings; "(with $m=1$)" | 9 | — (G11 judged it correct) | C7 | — | C2§6 | item 6 |
| 14 | $d_{2m}$: the charge basis is real | 9 | G11 | — | GP2-6 | — | item 6 (partly) |
| 15 | Lax generators described for untwisted algebras only | 9 | — | C5 | — | C2§5 | item 5 |
| 16 | Folded degrees: $a_{2n}^{(2)}$ is an exception; $\tau_1=\tau_3$ in parent labels | 10 | G6 | — | — | C2§7 ("unclear") | item 7 |
| 17 | Sign of the time advance depends on the frame | 10 | — (G13 judged it correct) | C8 | — | C2§4 | item 4 |
| 18 | "Minimum" of a complex potential | 9 | — | C17 | — | — | — |
| 19 | Soliton position at $\theta\ne0$ is $-\operatorname{Re}\xi/(m_a\cosh\theta)$ | 10 | — | C (unranked) | — | C2§13 (unverified) | — |
| 20 | "Energies are real" needs the normalized τ gauge and convergence | 10 | — | C13 | — | C2§2 | item 2 |
| 21 | "Closes gap H6" | 11, 44 | — | C (minor) | — (GP2-1 claimed it closes H6) | C2§1 | item 1 |
| 22 | Missing charges: "no static saddle" → no regular static vacuum-to-vacuum solution | 11, 17 | — | C19 | — | C2§8 | item 8 |
| 23 | The sinh-Gordon example's $m,\beta$ are $m_0,\beta_{\rm SG}$; the β→iβ substitution made explicit | 9 | G8 | C9 | GP2 obs. 2 | — | item 5 |
| 24 | Generic conservation law against the displayed current (fixed on main before integration) | 9 | — | — | — | C2§5 | item 5 |
| 25 | Scopes of the quantum-conservation citations | 9 | — | — | — | C2§5, §9 | item 9 |
| 26 | Takács–Watts's singular final states as a Θ-broken sector: an interpretation | 11 | — | — | — | C2§8 | item 8 |
| 27 | Symbol $E_1$ used for the exponential and the energy | 10 | — | — | GP2 obs. 1 | — | — |
| 28 | Learning goals and introduction of ch. 11 overclaim | 11 | — | — | — | C2§1 | item 1 |
| 29 | Chapter 9 claimed area-rule checks for $a_n$, $d_n$, $e_6$ that `toda_classical.py` did not contain | 9 | — | — | — | — | — (found by agent A) |

**Contested, not an error.** The Wigner sentence (10:122) for higher rank. Codex (C2§4) found the opposite sign in Hollowood's §5 and correctly left it "unclear"; Gemini (GP2-4) "confirmed" it with an integral representation of its own. Agent A showed that Hollowood's eq. (4.5) equals the book's $S_T$ at $n=2$ and gives the book's sign for $h=2,3,4$, so the sentence stands, and Hollowood's §5 is inconsistent with his (4.5).

**Referee claims that were false** (rejected at integration):
- Gemini: the type-A breather coefficient $-\tan^2(\pi a/h)/\tan^2u$ (P1); "correct" for the vertex formula, from a misreading of OTU eq. (6.3) (P1); the mislabelled $d_4$ coupling table (P1); "replace β² by −β² is false" (P1, overstated); the label [Theorem, given additivity] (P1, P2); "real and positive" as part of the general rule; "real negative energy" for a soliton at $\operatorname{Im}\theta=\pi$ (the $i\pi$ antikink has $E=+8$); "half energy at $A=0$" (P1); "closes H6 for Hirota soliton solutions" (P2); real $\tau_j$ stay positive for small $u$ (GP2-8); its summary of Correa–Fring–Taira (two-solitons and breathers, *Front. Phys.*; in fact one-solitons, *Nucl. Phys. B* 979) (P2); "None [could not be checked]" (P2).
- Codex: none found.

**Counts.** Of the 29 corrected errors, 28 were found by at least one referee, and 23 of these already in phase 1, including every major one (rows 1–6, 9). Phase 2 alone found five (rows 24–28), all minor: citation scope, interpretations, notation, one sign that main had already fixed. By referee, over both phases: both found 10 (rows 1–6, 9–11, 23; Gemini found row 3 only when phase 2 supplied the papers), Codex alone 13, Gemini alone 5 (rows 7, 8, 14, 16, 27). Codex made no false claim that the integration found; Gemini made about a dozen.

**Outcome.**
- **The blind phase did most of the work.** It found 23 of the 29 errors, and Codex found the Fring-group prior art without being pointed to it. Gemini did not search the literature until phase 2 gave it the papers, and then misread one of them.
- **The directed phase found few new errors** but settled the open questions with evidence: the spin-4 additivity on breathers, Freeman's antisoliton parity read in the source, the Hollowood comparison and the citation scopes. Its new findings are of one kind, scope and attribution, which a blind referee rarely targets.
- **Calibration.** As in Part II, Codex is reliable and Gemini finds real problems but needs every claim checked. Gemini's phase 2 added new false claims (on H6, on real $\tau_j$ at small $u$, on Correa–Fring–Taira) as well as confirmations.
- **For future briefs.** Keep the blind phase first. Make the directed phase short and aim it at what blind referees miss: prior art, citation scope, interpretive sentences, and claims that a book script checks something.

### Integration record A: Part III referee reports, chapters 9 and 10

- **Date:** 2026-10-24 (started 2026-10-15; resumed after a usage-limit stop, with every earlier check redone).
- **Integrated by:** Claude Opus 5.5 (integration agent A), worktree `.claude/worktrees/part3-int-ch9-10`, branch `part3-int-ch9-10`, base 30b61dd.
- **Reports:** Gemini phase 1 and 2 (`briefs/referee_classical_09_11_phase1_report.md`, `briefs/referee_classical_09_11_report.md`), Codex phase 1 and 2 (`briefs/referee_classical_09_11_codex_phase1_report.md`, `briefs/referee_classical_09_11_codex_report.md`), and the coordinator's merge `briefs/referee_classical_09_11_phase1_combined.md` (M1–M11). G*n* is Gemini's phase-1 rank *n*, C*n* Codex's phase-1 item *n*, C2§*n* Codex's phase-2 section *n*.
- **Own checks (Rule 1):** written from the chapters' definitions, not importing the referees' scripts, in `briefs/referee_classical_09_11_integration_scripts/`:
  - `integration_A_ch9.py` (items 1, 2, 5, 6, 7), 7 s, all PASS;
  - `integration_A_ch10.py` (items 9, 10, 11, 13, 14), 1 min, all PASS;
  - `integration_A_wigner.py` (item 12), 4 min, all PASS; output in `integration_A_wigner_output.txt`.
- **Sources read for this record:** Feigin–Frenkel hep-th/9310022 (pp. 4–5, §4.7); Delius–Grisaru–Zanon hep-th/9201067 (§2); Niedermaier hep-th/9401078 (abstract, §1); Olive–Turok–Underwood hep-th/9305160 (eqs. 2.27–2.30); Hollowood hep-th/9203076 (eqs. 3.21, 4.4–4.7 and §5, from the PDF); Harder–Iskandar–McGhee hep-th/9409035 (§3); Takács–Watts hep-th/9810006 (§4.2); Khastgir–Sasaki hep-th/9507001 (§5); Kneipp–Olive hep-th/9404030 (§6); McGhee hep-th/9307035 (p. 2).
- **Labels:** no status label changed in chapters 9–10.

#### Items

| # | Item and location | Sources in the reports | Verdict | How checked | Files changed |
|---|---|---|---|---|---|
| 1 | Spin of `eq-spin3-current`, sine-Gordon substitution, dictionary, `exr-spin3` | G8, C9, M5; C2§5 | **Verified; modified.** Main (0c1e204, 2fa0691) had already made the charge $Q_{-3}$ of spin $-3$ with $T_4$, $\mathcal X_2$ in the convention of ch. 2; that text is right. Added: the example's $m,\beta$ are $m_0,\beta_{\rm SG}$ of `sec-sg-dictionary`; the sine-Gordon current follows by $\beta\to i\beta$ ($\cosh\to\cos$, $\frac{\beta^2}4\to-\frac{\beta^2}4$, $V\to$ the potential of `eq-sg-lagrangian`); `exr-spin3` now states the kink value $Q_{-3}=-\frac{16m_0^3}{3\beta_{\rm SG}^2}e^{-3\theta}$. Gemini's "replace $\beta^2$ by $-\beta^2$ is false" rejected (cosh is even in $\beta$); made explicit instead. | Symbolic: sine-Gordon current conserved after $\beta\to i\beta$; potential maps to `eq-sg-lagrangian`. Quadrature of $\int(T_4-\mathcal X_2)dx$ on the boosted kink at $\theta=0,0.6,-0.45$ equals $-16m^3e^{-3\theta}/(3\beta^2)$ to $10^{-15}$. | ch. 9; `code/toda_classical.py` |
| 2 | Conservation-law sign at 09:41 versus the displayed current | C2§5 | **Verified; no edit.** After 2fa0691 the generic law $\partial_-T_{s+1}=\partial_+\mathcal X_{s-1}$ and the example's $\mathcal X_2=-m^2\cosh(\beta\phi)(\partial_+\phi)^2$ agree, and the density $T_4-\mathcal X_2$ is the conserved one (as in `eq-charge-eigenvalue`'s section). | By hand: $\partial_-T=\partial_+\mathcal X\Rightarrow\partial_t(T-\mathcal X)=\partial_x(T+\mathcal X)$; symbolic check already in `toda_classical.py`. | none |
| 3 | Quantum conservation citations (09:32) | C2§5, C2§9 | **Verified; modified.** The three sources have different scopes. Feigin–Frenkel prove Theorem 4.7.3 at generic $\beta$ when "all the exponents of g are odd and the Coxeter number is even" and say "the general proof … will be published separately" (also p. 5: "In the remaining cases we can prove these results by a more refined argument"). Delius–Grisaru–Zanon quote the spin-3 current of $a_3^{(2)}$ from their ref. [9] and write "Assuming that higher-spin quantum-conserved currents … exist". Niedermaier computes eigenvalues from a free-field realization of the ZF algebra, with first-order perturbative checks. The sentence now states the three scopes. | Read the three sources. | ch. 9 |
| 4 | "Minimum" of $V_{\rm imag}$ (09:30) | C17 | **Verified; modified.** On real fields $\operatorname{Re}V_{\rm imag}=\frac{m^2}{\beta^2}\sum_jn_j(1-\cos\beta\vec\alpha_j\cdot\vec\phi)\ge0$; the quadratic term is $\frac12\vec\phi^{\mathsf T}\mathsf M^2\vec\phi$. Text now says the quadratic term is positive definite and the real part is minimal at 0. | By hand. | ch. 9 |
| 5 | Lax-pair generators for twisted algebras (09:57); the gradation remark (09:77) | C5, M9; C2§5 | **Verified; modified.** For $a_2^{(2)}$ ($\vec\alpha_0=-2\vec\alpha_1$, labels $(1,2)$), $H=\frac1{\sqrt2}\operatorname{diag}(1,0,-1)$, $e_1=\sqrt2(E_{12}+E_{23})$, $e_0=E_{31}$ satisfy the Chevalley relations; `eq-lax-pair` with $c_j=\sqrt{n_j\vec\alpha_j^2/2}=(1,1/\sqrt2)$ has zero curvature exactly when the field equation holds; $e_0$ has weight $-2\vec\alpha_1$, not a root of the $A_1$ of $e_1,f_1$. Wording now separates untwisted and twisted. The gradation remark is right: $\lambda^R(e_1+\lambda^3e_0)\lambda^{-R}=\lambda(e_1+e_0)$, so in the homogeneous parameter $z=\lambda^3$ a boost needs the constant conjugation $e^{\vartheta R}$; left unchanged. | Symbolic 3×3 computation. | ch. 9; `toda_classical.py` |
| 6 | Area rule: nonzero couplings, normalization, degenerate bases (09:130–138, 09:160) | C7, G11, M8; C2§6 | **Verified; modified.** In the charge basis of $a_2,a_4,d_4,e_6$ (joint eigenvectors of $\mathsf M^2$ and the centre, bilinear couplings) every nonzero coupling obeys $\lvert C\rvert=\frac{4\beta}{\sqrt h}\Delta$ with $\Delta$ the area of the triangle of the actual masses (checked at $m=1.7$, so "(with $m=1$)" was wrong and is removed); the mass triples of the nonzero couplings are exactly those allowed by Dorey's rule (2, 4, 5, 18 triples); in $a_2$, $C_{112}=0$ although the triangle exists; the $d_4$ charge basis is real (characters $\pm1$) and a generic orthonormal basis of its light triplet violates the rule. Gemini's point on $d_{2m}$ adopted ("complex" and "real orthonormal basis" reworded). The sentence at 09:160 claimed checks for $a_n,d_n,e_6$ that `toda_classical.py` did not contain; the script now does $a_2,a_4,d_4,e_6$ and the sentence says exactly that. Gemini's printed $d_4$ coupling table (mislabelled, M8) not used. | numpy, own construction of the centre maps from permutations of the affine nodes. | ch. 9; `toda_classical.py` |
| 7 | Twisted spin lists (09:88) and simply-laced qualifier (09:98) | C15, G10; C2§5 | **Verified; no edit.** Spot check $a_4^{(2)}$: in $\mathfrak{sl}_5[z,z^{-1}]$ with $E_1$ the cyclic matrix and $\sigma(X)(z)=-KX(-z)^{\mathsf T}K^{-1}$ ($K$ antidiagonal, found by search to fix $E_1$), the centralizer of $E_1$ has dimension 1 in principal degrees 1, 3, 7, 9 and 0 in the others up to 10: odd degrees not divisible by 5, as printed. The qualifier at 09:98 is already there. | numpy nullspaces. | none |
| 8 | Vertex-operator formula without $x,t$ (10:72) | C16, M7; G17 (its defence via OTU (6.3) rejected) | **Verified; modified.** The displayed matrix element had no spacetime dependence. OTU (2.28)–(2.29) multiply each $\hat F^i(z)$ by $W(i,z,x^\mu)$; ch. 6 has $W_k(x,t)$. The formula now has $Q_kW_k(x,t)F^{a_k}(z_k)$, with $W_k$ referred to `sec-vertex-operators`. | Read OTU pp. 8–9 and ch. 6. | ch. 10 |
| 9 | Degrees of folded $\tau$-functions (10:77) | G6, M10; C2§7 | **Verified; modified.** With `tbl-foldings`' $n_j^\vee=kn_j\vec\alpha_j^2/2$, $a_{2n}^{(2)}$ has $n^\vee=(2,\dots,2,1)$, while a coincident pair $(a,2n+1-a)$ of $a_{2n}^{(1)}$ has $E^2$ coefficient $A_{a,\bar a}(0)\ne0$ at every node, so degree 2, not $\lvert O\rvert n_j^\vee=4$. The statement is now restricted to the foldings other than $a_{2n}^{(2)}$, with $a_{2n}^{(2)}$ named as the exception. The $c_2$ $\tau$-functions are right and are now said to use the parent labels. Nothing depending on @hyp-D was upgraded. Codex's $g_2$ degrees $(3,3,6)$ not recomputed (consistent with the restricted statement). | Exact $\tau$-functions from `eq-n-soliton`. | ch. 10; `toda_classical.py` |
| 10 | Position $-\operatorname{Re}\xi/m_a$ (10:59) | C P1 §4 obs., M (minor) | **Verified; modified.** $\lvert E\rvert=1$ on $x=t\tanh\theta-\operatorname{Re}\xi/(m_a\cosh\theta)$. | By hand and numerically. | ch. 10; `toda_classical.py` |
| 11 | Time-delay signs (10:122), caption of `fig-two-soliton-shift`, the symbol $E_1$ | C8, M11; G13 and Gemini phase-2 obs. 1 | **Verified; modified.** On the two-soliton $\tau$ ($a_5$, $\theta=\pm0.6$), soliton 1 emerges with $\xi_1+\ln A$; $\Delta x_1=-\ln A/(m_{a_1}\cosh\theta_1)>0$ and $\Delta x_2=+\ln A/(m_{a_2}\cosh\theta_2)<0$. In $a_2$ at $\theta=(-0.4,-1)$, $\Delta x_1>0$ but $\Delta t_1=+3.10>0$. The paragraph now gives both displacements, states the frame, and says that in the rest frame of soliton 2 $\Delta t_1<0$. The energy is now $\mathcal E_1$ (the note at 10:61 reserves $E$). The caption states the frame and speaks of a distance. The figure script was not changed: the caption's checked claims are unchanged. | Asymptotics of the exact $\tau$-functions. | ch. 10; `toda_classical.py` |
| 12 | Wigner relation versus the quantum S-matrix (10:122) | C2§4, Gemini phase-2 item 4 | **Verified: the sentence is right (option (a)); no discrepancy with ch. 21.** Section `sec-an-example` uses Hollowood's $a_n^{(1)}$ amplitude. (i) Hollowood's (4.5) at $n=2$ equals the book's $S_T$ of `eq-sg-smatrix` exactly, including the imaginary part of its log-derivative (λ = √3, e), so his rapidity and in/out conventions are the book's. (ii) The same formula gives $\frac{\beta^2}{2h}\big(-i\,d\ln S_{ik\to ki}/d\theta\big)\to\ln A_{11}(\theta)<0$ for $h=2,3,4$, with error falling like $1/\lambda$ (e.g. $h=3$, $\theta=0.7$: $-1.878$, $-1.910$ at λ = 28.3, 56.6; $\ln A_{11}=-1.928$). This is the book's sign. (iii) The ch. 21 scalar factor $c\,f$ at $h=2$ reproduces $S_0$ of `eq-sg-s0`, whose leading term is $(16/\beta_{\rm SG}^2)\ln\tanh\frac\theta2$. Codex's opposite sign comes from Hollowood's own §5, eqs. (5.2)–(5.5), where $\delta'=-(n/\beta^2)\ln A_{11}$ and $\Delta t>0$; that section disagrees with his (4.5) as printed. Gemini's "confirmation" used an integral representation it wrote itself and is not evidence. No correction to a published paper is proposed (it would need the full AGENTS.md procedure). | Digamma series for (4.5) (read from the PDF, p. 23), cosine transform for (-@eq-sg-s0), the $c\,f$ product of ch. 21. | none |
| 13 | Fusion (10:126, `exr-fusion`) | G19, C6, M6 | **Verified; modified.** For $h=3,\dots,8$: $a+b<h$ gives species $a+b$ at a simple pole $\theta_{12}=i\pi(a+b)/h$; $a+b>h$ gives species $a+b-h$ with $\theta_1=\theta+i\pi(h-b)/h$, $\theta_2=\theta-i\pi(h-a)/h$, simple pole at $i\pi(2h-a-b)/h$; $a+b=h$ gives a vanishing exponent (vacuum) and a double pole. The joint limit ($\theta_{12}\to$ pole, $\xi_{1,2}\to-\infty$, $A e^{\xi_1+\xi_2}$ fixed) reproduces the single soliton with error $\propto\sqrt\epsilon$ ($a_2$: 1+1→2, 2+2→1; $a_4$: 1+2→3). Text and exercise rewritten. | Exact $\tau$-functions. | ch. 10; `toda_classical.py` |
| 14 | Breathers (`sec-breathers`) | G3, C4, M4; HIM §3.2 | **Verified; modified.** $A_{a\bar a}(2iu)=1-\sin^2(\pi a/h)/\cos^2u$ and $A_{aa}(2iu)=-\sin^2u/(\sin^2\frac{\pi a}h-\sin^2u)$ for $h\le8$. For $A<0$, $\xi_2=\bar\xi_1+i\pi$ gives $\operatorname{Re}\tau_j=1+\lvert A\rvert r^2>0$ in the rest frame. For $0<A<1$, Codex's ellipse argument checked by hand ($f(A^{-1/2})=4A\cos^2\delta/(s+s^{-1})^2<1$; the segment case $s=1$) and numerically (zeros found for three phase choices). Band: $0<2M_a\cos u<2M_a\sin\frac{\pi a}h$, all of $(0,2M_a)$ only for $a=h/2$ — the statement the coordinator gave for both chapters. Added: type A breathers (HIM §3.1) with the correct $A_{aa}$; Takács–Watts's remark that the $a_2$ band ends at $\sqrt3M_a$ while quantum breathers reach $2M_a$. Gemini's type-A formula $-\tan^2(\pi a/h)/\tan^2u$ rejected (wrong; e.g. $a=1$, $h=4$, $u=1.2$). The quantum sentence ("the lowest breathers are the elementary particles") stays consistent: they lie at the bottom of the band. | By hand; exact $\tau$-functions. | ch. 10; `toda_classical.py` |
| 15a | Notes, "Reality" (10:174) | C13, C2§2 | **Verified; modified.** For $a_n^{(1)}$, $\ln\tau_j\to\ln\tau_j+f$ with $\partial_+\partial_-f=0$ leaves $\vec\phi$ and the bilinear equations unchanged (columns of the affine Cartan matrix sum to 0) but changes $C$ by $-2hf/\beta^2$; with $f=x^2+t^2$ it changes $T_{\mu\nu}$. Text now says "with the $\tau_j$ normalized as in `eq-n-soliton`", and replaces "the energies are real" by boundary terms where the integrals converge and real energies for regular real-rapidity solutions. Freeman's sentence made explicit. | By hand. | ch. 10 |
| 15b | Notes, Khastgir–Sasaki (10:180, 10:183) | G1, M1; C2§9 | **Verified; modified.** Khastgir–Sasaki §5: c-solitons (complex rapidity) "develop singularities after a finite time"; (5.13) $\phi=0$, $\partial_t\phi=if$ gives $E=-\frac12F<0$, and "We do not know if the solutions … remain finite or not". These data are Θ-invariant. Bullet rewritten; the closing sentence now says ch. 11 gives a sufficient condition for real charges of regular Hirota solutions and does not single out the physical solutions. | Read the source. | ch. 10 |
| 16 | Other ch. 9–10 points | C2§7, C2§9, G5, G9 | `kneipp1996` for untwisted Lie duality **verified** (Kneipp–Olive §6, p. 21: $M_{\langle i\rangle}=\lvert\langle i\rangle\rvert M_i$ for untwisted non-simply-laced folding, via coincident constituents). McGhee's count **verified** (hep-th/9307035, p. 2: $\tilde h_a=h/\gcd(a,h)$ "turns out to be equal to the number of topological charges"). Exercises (C2§9, Gemini item 9): only `exr-spin3` and `exr-fusion` needed changes (items 1, 13); the others were judged correct by both referees and nothing in this pass contradicts them. | Read the sources. | none |

#### What was not integrated, and why

- **`exr-breather-sg` and the remaining exercises.** No change; both referees found them solvable as stated.
- **Ch. 38, line 250** writes "(-@eq-spin3-current) with $\beta^2\to-\beta^2$". This is correct (cosh is even in $\beta$), and ch. 38 is touched by the open PR #42, so it was left alone. If wanted, "$\beta\to i\beta$" would match ch. 9.
- **The figure `fig-two-soliton-shift`.** Its script checks that $0<A_{ab}\le1$ "for every species of $a_5^{(1)}$"; at $\theta=0$ and $a=b$, $A_{aa}=0$. The caption's plotted claim is about distinct species and rapidities in practice; not changed here.
- **Codex's statement that Part III relies on @hyp-D only for degree claims** is reasonable; no edit.
- **A sentence on higher-rank amplitudes in ch. 10.** The item-12 check shows that Hollowood's amplitude has the leading term $(2h/\beta^2)\ln A_{11}$ for $h=3,4$. The book's sentence compares only with sine-Gordon, which is what it claims; a sentence on $a_2^{(1)}$, $a_3^{(1)}$ could be added, but its check takes about a minute per point and would need a script in `code/` other than `toda_classical.py`, which this task may not create. Left for the editor.

#### Pending edits in files this agent may not touch

- **Ch. 11 (other agent).** State the same neutral-breather band as ch. 10: regular iff $A_{a\bar a}(2iu)<0$, i.e. $\lvert\frac\pi2-\frac{\pi a}h\rvert<u<\frac\pi2$, so $2M_a\cos u\in(0,2M_a\sin\frac{\pi a}h)$, all of $(0,2M_a)$ only for $a=h/2$; the caption of `fig-breather-theta` ("fills a continuous band") and its figure. The Khastgir–Sasaki text in ch. 11 should match ch. 10's bullet (Θ-invariant data of negative energy; regularity undetermined), and ch. 11 should present itself as a sufficient condition for real charges of regular Hirota solutions. Ch. 10's breather section cites `prp-hirota-energy` for $\theta_0\pm iu$, $0<u<\frac\pi2$, where $\operatorname{Re}\cosh(\theta_0\pm iu)>0$, so the hypotheses that ch. 11 is adding hold.
- **Ch. 21 (observation, no error claimed).** `sec-an-example` defers the $a_n^{(1)}$ scalar factor to Hollowood, and that factor agrees with ch. 10 (item 12). But the sentence "the scalar factor can place its zero at $t=\omega$ itself" invites the reader to use the $c\,f$ of `sec-scalar-factor` with $a_1=\omega$. For $h\ge3$ that product does not reproduce Hollowood's factor: its leading semiclassical phase derivative is not $(2h/\beta^2)\ln A_{11}$ (scaled values $-0.85$ against $-1.93$ for $h=3$, and tending to 0 for $h=4$; scratch script `wigner.py`, not copied). The reason is that $f(T-t)=f(t)$ builds in self-conjugate crossing, while $V_1^*=V_n$ for $a_n^{(1)}$, $n\ge2$. The coordinator may want ch. 21 to say that the $c\,f$ form applies to self-conjugate multiplets.
- **Appendix E, references.bib, Part IX:** nothing needed from chapters 9–10. No new bibliography entry.
- **Appendix D:** no correction to a published paper. Hollowood's §5 sign is inconsistent with his own (4.5), but it is not recorded as a correction without the full procedure.

#### `code/README.md` row

The row of `toda_classical.py` should read (ids and description):

| `toda_classical.py` | `eq-lax-pair`, `eq-toda-eom`, `eq-spin3-current`, `exr-spin3`, `sec-masses-couplings`, `prp-dorey-rule`, `eq-area-rule`, `eq-hirota-bilinear`, `eq-an-soliton`, `sec-vertex-solitons`, `eq-n-soliton`, `eq-interaction-coefficient`, `sec-multi-solitons`, `sec-breathers`, `eq-soliton-mass`, `eq-theta-data`, `exr-breather-real` | Lax pair ⇔ field equation for $a_2^{(1)}$ and $a_2^{(2)}$; sinh-Gordon $\partial_+$ current and the spin $-3$ of its charge; sine-Gordon current and $Q_{-3}$ on the kink; masses = Perron–Frobenius vector and sum rule; Dorey's rule and $\lvert C_{abc}\rvert=\frac{4\beta}{\sqrt h}\Delta_{abc}$ for $e_7,e_8$, and in the charge basis for $a_2,a_4,d_4,e_6$ (nonzero couplings only; a generic $d_4$ basis fails); one- and two-soliton $a_n^{(1)}$ solutions; soliton position; folded degrees ($c_2$, $a_4^{(2)}$); time-delay signs and their frame dependence; fusion for $a+b<h$, $>h$, $=h$; neutral-breather band ($A<0$ regular, zeros for $0<A<1$); soliton energy $2hm_a/\beta^2$ at complex $\xi$; $\Theta$ on Hirota data; sine-Gordon breather reality | about 4 minutes |

#### Final checks (2026-10-24)

- `code/toda_classical.py`: 24 PASS, 0 FAIL (about 4 minutes).
- `code/algebra_data.py`: 6 PASS; `code/smatrix_simply_laced.py`: 11 PASS. Both exec only the root-data part of `toda_classical.py`, which is unchanged.
- The figure scripts were not run or changed: no caption now claims a check that its script does not make.
- `python3 check-refs.py`: all PASS (one WARN on `citation.qmd`, present before this change).

### Integration record B: Part III referee reports, chapter 11

- **Date:** 2026-10-24 (started 2026-10-15; resumed twice after usage-limit stops; every earlier check redone from scratch, since the first scratch directory was lost).
- **Integrated by:** Claude Opus 5.5 (integration agent B), worktree `.claude/worktrees/part3-int-ch11`, branch `part3-int-ch11`, base 30b61dd.
- **Reports:** Gemini phase 1 and 2 (`briefs/referee_classical_09_11_phase1_report.md`, `briefs/referee_classical_09_11_report.md`), Codex phase 1 and 2 (`briefs/referee_classical_09_11_codex_phase1_report.md`, `briefs/referee_classical_09_11_codex_report.md`), the coordinator's merge `briefs/referee_classical_09_11_phase1_combined.md` (M1–M4 are this agent's). G*n* is Gemini's phase-1 rank *n*, GP2-*n* Gemini's phase-2 item *n*, C*n* Codex's phase-1 item *n*, C2§*n* Codex's phase-2 section *n*.
- **Own checks (Rule 1).** Written from the chapters' definitions; no referee script imported.
  - `code/classical_reality.py` (new, cited in ch. 11): 20 PASS, 8 s.
  - `code/figures/ch11_f16_breather_theta.py` (rewritten): 8 PASS, 2 s.
  - Scratch verifications, with outputs, in `briefs/referee_classical_09_11_integration_scripts/`: `integration_B_currents.py` (a₂ currents of spin 1–4 solved from the field equation), `integration_B_evalq.py` (canonical $E$, $P$ and charges on Hirota solutions), `integration_B_run1_charges.py`, `integration_B_run2_energy_domain.py`, `integration_B_run3_relabel_density.py`.
- **Sources read.** Olive–Turok–Underwood hep-th/9305160 (= `olive1993b`), §8.3; Takács–Watts hep-th/9810006, eq. (4.5)–(4.13), §4.2; Freeman hep-th/9408092, abstract, (3.44)–(3.49), Conclusions; Khastgir–Sasaki hep-th/9507001, §5 and (5.12)–(5.13); Harder–Iskandar–McGhee hep-th/9409035, §3.2; Zhu–Caldi hep-th/9404097, introduction and (2.5); Cen–Fring 1602.05465, Cen–Correa–Fring 1608.01691, Fring–Taira 2007.15425, Correa–Fring–Taira 2110.06825 (quotations checked on my own `pdftotext` of the PDFs); Cen–Fring 1812.02111 (fetched from arXiv). Metadata of the five new entries from the arXiv API and INSPIRE.

#### Items

| # | Item and location | Sources in the reports | Verdict | How checked | Files changed |
|---|---|---|---|---|---|
| 1 | `prp-hirota-energy` (11:32–41) | M2; G1, C2; GP2-2, C2§2 | **Verified; modified.** The formula needs (i) $\operatorname{Re}\cosh\theta_i>0$, (ii) finite nonzero $\prod A$, (iii) no zero of $\tau_j$ on the slice. Statement restricted to (-@eq-n-soliton) with (i)–(iii), boundary derivation added, OTU §8.3 quoted ("holds good even if $e^{\eta_i}$ is complex, providing that its real part is positive"; independence of $Q_i$, $\epsilon_i$), TW eq. 4.6 kept for $a_2$. Counterexample in the text: sine-Gordon $\theta=i\pi$, $\xi=i\pi/2$, $E=+8$ against $-8$. Relabelling $(a,\theta+i\pi,\xi)\simeq(\bar a,\theta,-\xi)$ stated, with $A_{ab}(\theta+i\pi)A_{\bar ab}(\theta)=1$ for multi-solitons. Label kept **[Established; reformulated]**. Gemini's "half energy at $A=0$" and its "real negative energy at $\operatorname{Im}\theta=\pi$" not used (the latter is false: the $i\pi$ antikink has $E=+8$). | Canonical quadrature of $\frac12(\phi_t^2+\phi_x^2)+V$ and $-\phi_t\phi_x$ from exact $\tau$ derivatives: complex single solitons of both species, a non-conjugate complex pair, breathers ($10^{-14}$); $i\pi$ antikink $E=8.0000000000$; $a_2$ soliton and a two-soliton with one $\operatorname{Re}\cosh<0$ constituent give $-M\cosh\theta$ for it; $\tau_j$ identity checked pointwise; density $=\frac2{\beta^2}\sum_j\frac2{\alpha_j^2}\partial_x^2\ln\tau_j$ pointwise ($10^{-14}$). | ch. 11 |
| 1a | Other uses of `prp-hirota-energy` | prompt | **Checked.** Ch. 39:19 already restricts to $\lvert\operatorname{Im}\theta\rvert<\frac\pi2$ "as in @prp-hirota-energy": consistent. Ch. 49:112, 126 use it for real masses: consistent. Ch. 17 does not cite it. Ch. 10 (agent A) cites it at $\theta_0\pm iu$, $0<u<\frac\pi2$: (i) holds there; I confirm A's remark. Ch. 30:319 cited it for the reality of all charges of $\Theta$-invariant data, which is the other statement: now cites `prp-theta-reality` and says "real energy, and higher charges real or imaginary according to their parity". | grep of ids. | ch. 30 |
| 2 | The $\Theta$ selection rule (11:43–49) | M1; G2, C3, C18; GP2-1, C2§1 | **Verified; modified.** New statement `prp-theta-reality`, **[Established; reformulated]**, with proof. Hypotheses: definite parity $Q[\Theta\phi]=\epsilon_Q\overline{Q[\phi]}$; a class of regular solutions closed under $\Theta$ on which $Q$ depends only on $D=\{(a_i,\theta_i)\}$; $\Theta$-invariant $D$. Remarks: $(1+i)Q$ has no parity (Codex); topological charge excluded (sine-Gordon kink, Codex); phase independence argued for single solitons (complex translation) and separated real-rapidity solitons, not proved for complex rapidities; the weaker $\{(M,\theta)\}$ condition kept. Field versus data distinguished (11:26): a field fixed by $\Theta$ is $-i\vec\chi$ with $\vec\chi$ a real-coupling solution, and its energy is $-E_{\rm real}[\vec\chi]<0$ (new two-line computation, not claimed as new). Freeman (3.49) verified in the source: "Q−N(antisoliton) = (−1)^{N−1} Q−N(soliton)", with real normalization; so $\epsilon_{Q_s}=(-1)^{s-1}$, and for $a_n$ $\sin\frac{s\pi(h-a)}h=(-1)^{s-1}\sin\frac{s\pi a}h$. Gemini's **[Theorem, given additivity]** and "real and positive" rejected; positivity on the class of `prp-hirota-energy` comes from (i), not from $\Theta$ (now said in the text). | Proof by hand. Kink: topological charge 1 from the $\tau$s, reversed by $\Theta$. $\Theta$-fixed field: functional identity $E_{\rm imag}[-i\chi]=-E_{\rm real}[\chi]$ on a random $a_2$ profile. Freeman parity: my currents give $q_2/q_1=-1$ at spins 2 and 4, raw $q_1$ imaginary for real-coefficient currents. | ch. 11 |
| 2a | Prior art (11:96) | C18, C2§1, GP2-1 | **Verified; modified.** What each paper does: Cen–Fring 2016, complex KdV/mKdV/sine-Gordon one- and two-solitons with real dispersion parameters; real energies explained by PT restored by space-time shifts; N-soliton case conjectured ("we conjecture that one might be able …"; "Here we have not allowed complex dispersion relations"). Cen–Correa–Fring 2017, KdV: all Gardner charges real on complex multi-solitons by PT plus asymptotic separation ("Thus PT-symmetry and integrability guarantee the reality of all charges"). Fring–Taira 2020: complex BPS solitons; self-dual and anti-self-dual solutions related by a modified CPT symmetry ("the energy is the same for either choice"). Correa–Fring–Taira 2022: complex **one**-solitons of Bullough–Dodd (types I and II), energies $-6\lvert k\rvert$, linear stability; a symmetry mapping a solution "to itself or possibly to a second degenerate solution". Cen–Fring 2020: multicomplex KdV; antilinear symmetries force components of the charges to vanish. OTU §8.3 already has real $E,P$ for equal-mass conjugate pairs (breathers, any affine Toda), TW §4.2 the $a_2$ two-soliton condition. Gemini's summary of Correa–Fring–Taira (2-solitons and breathers; *Front. Phys.*) is wrong: one-solitons only; *Nucl. Phys. B* 979 (2022) 115783. Zhu–Caldi select by reality of the energy *density* (pairing condition): different and stronger. | Read the sources; quotations on my own PDF extraction. | ch. 11, `references.bib` |
| 3 | Appendix E, Part III item | prompt; C3, C2§1, GP2-1 | **Modified.** The [Theorem, given additivity…] item removed. Replaced by **[Numerical]** for `prp-breather-charges`: the $\Theta$-odd charges of spins 2 and 4 of $a_2^{(1)}$ are imaginary on breathers and additive at complex rapidity. Reason: the general argument is prior art (item 2a); the breather evaluation of higher charges appears in none of the papers read (OTU §8.3 evaluates only $E,P$; HIM, Gandenberger hep-th/9501136 and an arXiv search for affine Toda breathers found no classical higher-charge evaluation). The editor should confirm. | `check-refs.py`: original-results items match their statements (39 items). | `appendices/E-original-results.qmd` |
| 4 | Additivity test (11:46) | GP2-3, C2§3 | **Verified; stated as checked.** New `prp-breather-charges` **[Numerical]**: with $q_1=1$, $Q_s=e^{s(\theta_0+iu)}-e^{s(\theta_0-iu)}=2ie^{s\theta_0}\sin su$, $s=2,4$. "Assumed here" removed. | Currents solved from the field equation alone in $\vec\psi=i\beta\vec\phi$ (one primitive current at spins 1, 2, 4; none at 3; exact rational nullspace); integrated on single solitons (rest, complex rapidity: $q_ae^{s\theta}$, $q_2=-q_1$) and on breathers $\theta_0=0.23$, $u=\pi/3$ with three phase choices ($\xi_2=\bar\xi_1+i\pi$, shifted real parts, and an unrelated pair), each at a different time: errors $6\times10^{-15}$ (spin 2), $3\times10^{-14}$ (spin 4). Same numbers as both referees (up to the current's normalization). | ch. 11, `code/classical_reality.py` |
| 5 | Khastgir–Sasaki (11:47) | M1, G1; GP2-1; C2§1, C2§9 | **Verified; modified.** KS §5: "all c-soliton solutions develop singularities after a finite time"; (5.13) $\phi=0$, $\partial_t\phi=if$, $E=-\frac12F<0$; "We do not know if the solutions … remain finite or not". These data are fixed by $\Theta$, so the solution is a $\Theta$-fixed field, hence of negative energy by item 2. Text now: complex-rapidity solitons have $E,P$ not both real and non-invariant data (consistent with a sufficient condition); negative-energy data are $\Theta$-fixed; $\Theta$ gives reality, not positivity; the parallel with complex eigenvalues is a heuristic. Matches agent A's ch. 10 bullet. | Read the source; uniqueness of the Cauchy problem gives $\Theta\phi=\phi$. | ch. 11 |
| 6 | H6 (11:49; 44:60, 44:81) | C2§1, GP2-1 | **Verified; modified.** "closes gap H6" → partly addresses it, with what remains open (necessity, singular solutions, choice of physical solutions). Table cell now "partly open: …"; it must start with one of the words that `code/figures/ch44_f33_gap_map.py` classifies, so "partly addressed" would have failed its check. The gap map was regenerated: H6 is now orange (6 PASS). 44:60 rewritten as a sufficient condition. 44:33 (the gap's definition) unchanged. Gemini's "closes H6 for Hirota soliton solutions" rejected. | Gap-map script. | ch. 11, `44-twelve-gaps.qmd`, `ch44_f33_gap_map.svg` |
| 7 | Singular solutions (11:58–62), `exr-singular-points` | M3; G7, C1; GP2-8, C2§8 | **Verified; modified.** Double pole, zero residue, divergent integral: static sine-Gordon density $-\frac{8m^2}{\beta^2}\operatorname{csch}^22mx$, cut integral $\frac{8m}{\beta^2}(1-\coth2m\varepsilon)$; only a contour or finite-part prescription gives a finite value. Added: slices that avoid the zeros give the boundary value (checked on a singular-band $a_2$ breather); a zero crossing the line changes the winding of $\tau_j$ by $2\pi$ and the topological charge by a coroot (Gemini); isolated zeros need independent gradients, real-rapidity single solitons have lines (Gemini); "the Hirota data continue analytically" replaced by: the $\tau$ do, the field and the integrals do not (Codex). The Takács–Watts sentence now marks the $\Theta$-broken-sector reading as an interpretation. Exercise rewritten after Codex. | Density at three points ($10^{-10}$); cut integrals at $\varepsilon=10^{-2},10^{-3}$ (geometric grid, $10^{-6}$ relative); $u=0.2$ breather at $t=0,1.9,5.1$; complex-rapidity soliton: zero located, Jacobian $\ne0$, winding jump $2\pi$. | ch. 11 |
| 8 | Breathers and the figure (11:23, 11:28, `fig-breather-theta`) | M4; G3, C4; GP2-8 | **Verified; modified.** Same statement as ch. 10 (agent A): regular iff $A_{a\bar a}(2iu)=1-\sin^2\frac{\pi a}h/\cos^2u<0$, mass in $(0,2M_a\sin\frac{\pi a}h)$, full band only for $a=h/2$; zeros of $\tau$ for $0<A<1$ whatever the phases. Right panel redrawn for $a_5^{(1)}$ (bars per species, excluded part hatched); caption, alt text and the example bullet updated. Gemini's type-A formula not used. Gemini's claim (GP2-8) that for $\xi_2=\bar\xi_1$ and small $u$ the real $\tau_j$ stays positive ("negative discriminant $4-4A<0$") rejected: $A<1$ always, so $1-2r+Ar^2$ is negative for some $r$; `exr-breather-real` holds for every $u$, as Codex says. | Figure script: $A$ formula and sign on a grid, regular family ($\operatorname{Re}\tau_j\ge1$), a zero of $\tau_0$ at $a=1$, $u=0.524$ ($A=0.667$), band tops. `classical_reality.py` repeats these for three phase choices. | ch. 11, figure script and SVG |
| 9 | Missing charges (11:72) | C19 | **Verified (by reading Codex's argument and ch. 17); modified.** "no static saddle of the right energy" → "no regular static solution of the right energy that approaches vacua at both ends". `prp-static-saddles` (ch. 17) said "no static saddle of energy $M_2$"; made the same minimal change there. Codex's unstable-manifold argument not recomputed. | Reading. | ch. 11, ch. 17 |
| 10 | Notes (11:96), learning goals (11:5–6), `index.qmd:43`, 44:60 | C18, C2§1, GP2-1 | **Modified.** "The formulation as a selection rule from $\Theta$ is new here" replaced by the prior art of item 2a and what the chapter adds (the statement on reduced data; the checks of `prp-breather-charges` and of the negative energy of $\Theta$-fixed fields). Learning goals: sufficient condition; reality versus regularity, positivity and physical selection. Intro paragraph and chapter title ("reality from $\Theta$" instead of "a selection rule from $\Theta$") adjusted. `index.qmd:43` no longer says "guarantees". | — | ch. 11, `index.qmd`, 44 |
| 11 | Other ch. 11 points | GP2-9, C2§8, C2§9 | **Verified.** $+i\phi$ is right ($\vec\chi=i\vec\phi$ solves the real-coupling equation; derived in item 2). `eq-theta-data` and `exr-theta-interaction`: cos is even and $2\pi$-periodic; `code/toda_classical.py` already checks $\bar\tau_j(a,\theta,\xi)=\tau_j(h-a,\bar\theta,\bar\xi)$; no change. `exr-theta-solutions`: correct, no change. `exr-breather-real`: correct; rewritten in two parts (the $\xi_2=\bar\xi_1$ singular case with $1-2r+Ar^2$; the regular case $\xi_2=\bar\xi_1+i\pi$, $A<0$) to match the band. | `classical_reality.py` (real $\tau_j$ take negative values; regular family). | ch. 11 |

#### Status labels changed or created

- `prp-hirota-energy`: kept **[Established; reformulated]**; the statement is narrower (hypotheses (i)–(iii), $a_n^{(1)}$ $N$-solitons). Not yet seen by a referee in this form.
- `prp-theta-reality`: **new, [Established; reformulated]**, with a proof (one line plus the $E,P$ parity). The coordinator's default label. Its structure follows Codex's conditional derivation (C3, C2§1); the text itself has not been refereed.
- `prp-breather-charges`: **new, [Numerical]**. Both referees obtained the same values independently (spin 2: both; spin 4: Codex); the statement text has not been refereed.
- Appendix E: the item **[Theorem, given additivity of the higher charges]** is removed; a **[Numerical]** item for `prp-breather-charges` replaces it.
- No label was upgraded.

#### Not integrated, and why

- **TW eq. (4.6)** is stated for "θi ∈ C" without the restriction to $\lvert\operatorname{Im}\theta_i\rvert<\frac\pi2$; with $\theta=i\pi$ it fails as in item 1. I did not add an Appendix D correction: TW themselves restrict $\alpha$ to the "fundamental region" $[-\pi/2,\pi/2]$ in §4.2, and the slip is one of parametrization. The coordinator may decide otherwise; the AGENTS.md procedure would be needed.
- **Gemini's "half energy at $A=0$"** (G1): not used. At $A=0$ condition (ii) fails and the solution is no longer a two-soliton; `prp-hirota-energy` now simply excludes it.
- **Gemini's extension of `eq-theta-data` to $d_n$, $e_n$** (via MacKay–McGhee): not checked; the chapter keeps the $a_n^{(1)}$ statement.
- **Ch. 39:145** speaks of a defect that "respects the physical sector" for $\Theta$-invariant data. That is the book's terminology there, not a repetition of the corrected claim; left unchanged (Part VIII, near the open PR #42).
- **A complete proof of phase independence for complex rapidities** (needed to turn `prp-breather-charges` into a theorem) is not attempted.

#### Pending edits in files this agent may not touch

- **Ch. 10 (agent A, branch `part3-int-ch9-10`).** Nothing further needed: A's breather band, Khastgir–Sasaki bullet and closing sentence agree with ch. 11 as written here. A's ch. 10 Notes provide the formula $T_{\mu\nu}=(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$ that ch. 11 now cites via `@sec-classical-literature`; on `main` the same formula is in ch. 10's Notes, so the reference resolves either way.
- **`code/README.md`** (rows below).

#### `code/README.md` rows

New row (after `toda_classical.py`):

| `classical_reality.py` | `prp-hirota-energy`, `prp-theta-reality`, `prp-breather-charges`, `sec-theta-classical`, `sec-singular-solutions`, `exr-breather-real`, `exr-singular-points` | canonical $E$, $P$ by quadrature: additive for complex rapidities with $\operatorname{Re}\cosh>0$; the $i\pi$ sine-Gordon antikink has $E=+8$; a constituent with $\operatorname{Re}\cosh<0$ contributes $-M\cosh\theta$ (relabelling, $A_{ab}(\theta+i\pi)A_{\bar ab}(\theta)=1$); density $=\frac2{\beta^2}\sum\frac2{\alpha_j^2}\partial_x^2\ln\tau_j$; breather $E,P$ real; the kink's topological charge has parity $-1$; a $\Theta$-fixed field has $E<0$; $a_2^{(1)}$ currents of spin 2 and 4 built from the field equation, $q_2=-q_1$, breather charges $2ie^{s\theta_0}\sin su$; divergence $\frac{8m}{\beta^2}(1-\coth2m\varepsilon)$ at a $\tau$ zero; regular slices of a singular breather; winding jump at a zero crossing; $a_5^{(1)}$ breather band, regular family, $\tau$ zeros for $0<A<1$ | 8 s |

Changed figure row:

| `figures/ch11_f16_breather_theta.py` | `fig-breather-theta` | `sec-theta-classical`, `sec-breathers` | breather $E$, $P$ real with $E=2M_a\cosh\theta_0\cos u$, $E^2-P^2=(2M_a\cos u)^2$; $\Theta$ maps the pair to itself; $M_a=M_{h-a}$; for $a_5^{(1)}$: $A_{a,h-a}(2iu)=1-\sin^2\frac{\pi a}h/\cos^2u$, negative exactly in the band; regular family in the band; a zero of $\tau_0$ below it; band tops $2M_a\sin\frac{\pi a}h$ | seconds |

The row of `figures/ch44_f33_gap_map.py` is unchanged (its checks are the same; only the H6 colour changed).

#### For the coordinator to decide

1. The label of `prp-theta-reality` (**[Established; reformulated]**, the default) and the Appendix E choice (**[Numerical]** item for the breather charges, rather than no item).
2. Whether the observation that every non-vacuum field fixed by $\Theta$ has negative energy deserves an Appendix E line (I think it is too elementary; it is a two-line computation).
3. Whether TW (4.6)'s "θi ∈ C" warrants an Appendix D entry (I recommend not).
4. The chapter title change, "The classical theory: reality from $\Theta$".

#### Final checks (2026-10-24)

- `code/classical_reality.py`: 20 PASS, 0 FAIL (8 s).
- `code/figures/ch11_f16_breather_theta.py`: 8 PASS; SVG regenerated and viewed on light and dark backgrounds.
- `code/figures/ch44_f33_gap_map.py`: 6 PASS; SVG regenerated (H6 orange).
- `python3 check-refs.py`: all PASS (one WARN on `citation.qmd`, present before).
- No Quarto render (as instructed).
