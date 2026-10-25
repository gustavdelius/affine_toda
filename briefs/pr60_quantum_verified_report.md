# PR #60: independently verified quantum-defect and lecture recommendations

- Date: 2026-10-25.
- Referee: Codex / GPT-6, quantum subagent; no subagents spawned.
- Reviewed snapshot: `eecc48fdf06bc079af2175fac53ea2d5d2d92abf`, the merge of PR #60. Original PR head: `9f44dedfcdf88b4a84950c1188899e9082f749bc`. Sources were read in `referee_scratch/pr60_snapshot/`; line numbers below refer to that snapshot, not the original reports' `9b912ee`.
- Reports reviewed: every item 1–7c of `briefs/referee_defects_38_41_pr52_report.md`; sections 3, 4 and 5.1 of `briefs/referee_lectures_defects_report.md`. The coordinator handles exercise counts. The Gemini lecture report is excluded at the editor's request because it has already been integrated.
- No existing book, lecture, script, or report file was edited.

## Editor-authorized integration exception

**For the verified conclusions and the specific modified replacements in this report, Claude may integrate directly without repeating Rule 1 of `AGENTS.md`, step 3 of `.agents/skills/integrate-referee-report/SKILL.md` (also exposed through `.claude/skills/`), or that skill's requirement to re-derive label changes. Gustav Delius explicitly requested this exception for this verification round. It is not a global change to the repository rules.**

The independent computations are complete for the dispositions below. “Unclear” means that the proposed stronger assertion is **not cleared**; its specified disposition is to retain a weaker label, omit the proposed addition, or narrow the claim. Claude need not undertake fresh mathematical research to carry out those conservative dispositions. Checking that the target passage still matches, running regression scripts, checking references, and rendering the resulting text remain ordinary integration work. The exception does not turn the original Gemini reports or their proposed replacement text into certified evidence.

## Essential corrections to the incoming reports

1. **Do not paste the incoming Item 2 proof.** Its integral `cos(ku)/(k cosh(pi k/2))` diverges at zero. Its differentiation and its asserted value at zero are also wrong. The existing book uses the correct sine integral; preserve it.
2. **Do not paste incoming Item 7b.** It describes a different, incorrect boundary matrix and says off-diagonal shifts cancel although its own displayed product has nonzero off-diagonal entries. The actual snapshot already contains a complete proof of the correct Dirichlet-to-general-Ghoshal–Zamolodchikov fusion statement.
3. **Do not paste the incoming Nemes Appendix D bullet.** Its purported quotations and equation numbers do not match the primary source. The four extra misprints are not already in the snapshot's Notes, contrary to the report.
4. **Do not promote mutual commutation of two Weyl defect operators to a theorem on the supplied argument.** The train proof with one finite auxiliary soliton and the local fusion proof of TQ are valid. The extra sentence about “the same power” on charge-neutral paths is not a proof of the two-defect commutator.
5. **Do not copy incoming Item 5's Weyl-pair definitions.** They differ from the book, omit essential factors of $i$, give the wrong state eigenvalues, and assert a nonzero commutator between two multiples of the same field combination. The book's actual formulas pass the independent BCH calculation.

## Evidence and execution

The following are new computations from the displayed definitions. None imports a book or previous referee helper.

| Independent file | What it establishes |
|---|---|
| `briefs/pr60_quantum_verified_scripts/check_quantum.py` | Six exact BCH brackets; localized-mode determinants; scalar crossing, unitarity, Gamma-cutoff convergence and B1 ratio; pole scaling; scalar and matrix boundary fusion; finite-chain TQ, t–D and shift covariance; numerical D–D; Wigner and classical-breather limits |
| `briefs/pr60_quantum_verified_scripts/check_local_blocks.py` | Exact kernel and quotient fusion blocks with independent formal variables $t=q^n,r=q^s$; exact obstruction to the homogeneous intertwiner solution |
| `briefs/pr60_quantum_verified_scripts/check_triangle.py` | Exact RTT relation on all four two-soliton input spins at arbitrary (q,t,r,Y,Z); all four crossing entries and the wrong-sign control |

Run each with `python3 -u <path>`. Corresponding `.out` files contain the complete outputs. Representative results:

```text
PASS BCH U,W residual = 0
PASS BCH Us,X residual = 0
PASS BCH U,X residual = 0
PASS BCH Us,W residual = 0
PASS BCH W,X residual = 0
PASS BCH U,Us residual = 0
PASS odd determinant residual = 0
PASS even determinant residual = 0
PASS zero-rapidity kernel residual = 0
PASS scalar crossing 0.31 residual = (-5.67407592e-28 + 4.19188247e-37j)
PASS B1 ratio 0.31 residual = (-3.76158192e-37 - 7.52316385e-37j)
PASS gamma cutoff errors 0.00045306436 0.000226533136 0.000113266687
PASS GZ scalar identity residual = (1.75328045e-28 + 3.56079433e-28j)
PASS Dirichlet matrix dressing residual = 0
PASS TQ, tD, shift N=5 residual = 2.987304636342383e-15
PASS DD numerical only N=5 residual = 8.673617379884035e-18
PASS fusion block e input 1 all residuals zero
PASS fusion block e input -1 all residuals zero
PASS fusion block f input 1 all residuals zero
PASS fusion block f input -1 all residuals zero
uniqueness homogeneous obstruction = c*hom*t*(q - 1)*(q + 1)*(q**2 + 1)/q**3
PASS exact triangle input 1 1
PASS exact triangle input 1 -1
PASS exact triangle input -1 1
PASS exact triangle input -1 -1
PASS all four crossing entries with Yprime=-qY
```

Numerical representation tests use $\lambda=\sqrt7$, not a root-of-unity decimal, and use exact paths on the infinite integer lattice, not a truncated Weyl module. The cutoffs in the Gamma calculation are explicitly finite: the errors halve with each doubling; the computation does **not** claim 30-digit Gamma-product convergence. Symbolic identities and the derivations below provide the all-parameter evidence.

As a secondary comparison only, selected functions of `code/defect_sg_scalar_poles.py` were run in the frozen snapshot: `check_pole_origin`, `check_widths`, `check_odd_jump`, `check_emission_tail`, `check_inverted_oscillator`, `check_two_weyl_pairs`, and `check_faddeev_normalization`. All seven printed PASS; output is in `book_comparison.out`. This is not the independent evidence.

## Item-by-item verification

### Q1 / lecture report 3.3: the scalar integral

**Verdict: correct book formula; wrong intermediate signs and limit argument in the incoming report.** Targets: `part8-defects/38-defect-sine-gordon.qmd:225` (`prp-sgd-scalar-integral`); `lecture-notes/defects/07-quantizing.qmd:415` (`eq-l7-Phi`). Keep the displayed formula and its existing proof.

My derivation starts with
$$
\log\Gamma(a-z)-\log\Gamma(a+z)
=\int_0^\infty[-2ze^{-t}+2e^{-at}\sinh(zt)/(1-e^{-t})]\,dt/t.
$$
For $z=i\lambda u/(2\pi)$ the sine term has **positive** sign. The signed exponential sum for the four Gamma arguments is
$$
\frac{(1-e^{-\lambda t/2})(e^{-t/4}-e^{-\lambda t/2-3t/4})}{1-e^{-\lambda t}}.
$$
Thus the phase of $r/\bar r$ is $+\int H(k)\sin(ku)\,dk/k$, not the negative integral printed in incoming Item 1. Splitting
$$
H(k)=\frac{\tanh(\pi k/2)}{2\sinh(\pi k/(2\lambda))}
+\frac1{2\cosh(\pi k/(2\lambda))}
$$
and adding the prefactor gives exactly the book's integral. The $\log\cosh$ and Gudermannian Fourier identities follow by differentiating and integrating their elementary cosine transforms, with zero value at the origin.

At zero the integrand is $\lambda(i\pi u-u^2)/(2\pi)+O(k^2)$. At infinity its three relevant decay rates are $\pi/(2\lambda)+\operatorname{Im}u$, $\pi+\pi/(2\lambda)-\operatorname{Im}u$, and $\pi/(2\lambda)$. This proves local uniform convergence in the stated open strip and zero-freeness after exponentiation. The combined (k)-summands in the Malmstén calculation have fixed sign before the common sine and are (O(t)) in total at (t=0), so the exchange of their sum with the integral is justified.

For the normalization at $u\to-\infty$, the phase kernel is (h(k)/k), with $h(0)=\lambda/2\ne0$. It is not an (L^1) amplitude to which Riemann–Lebesgue applies. Subtract (h(0)e^{-k}): the remainder divided by (k) is integrable, and $\int_0^\infty e^{-k}\sin(ku)\,dk/k=\arctan u\to-\pi/2$. Hence the phase tends to $\pi\lambda/4-\pi\lambda/4=0$. The book already uses this valid Dirichlet-integral limit. Crossing follows by adding the integrands at (u) and $u+i\pi$, and unitarity from the real part. These are mathematical identities for the specified scalar, not a derivation of the Lagrangian quantum theory.

**Additional lecture precision:** the blanket uniqueness/minimality assertion at Lecture 7 line 429 exceeds this proof: $-\Phi$ alone is a counterexample to “any other solution” having a new pole, zero or faster growth. Replace that sentence with:

```markdown
We choose the normalized minimal scalar factor of Konik and LeClair. Multiplication by a defect CDD factor preserves the scattering equations; the integral fixes that choice and its normalization, rather than deriving it from the defect Lagrangian.
```

This does not change the chosen amplitude.

### Q2 / lecture report 3.4 and 5.1: breather transmission

**Verdict: correct formula and existing B1 proof; wrong incoming replacement proof; label modification warranted.** Targets: book `38-defect-sine-gordon.qmd:335`, `prp-sg-breather-transmission`; Lecture 7 `thm-l7-breather`; Lecture 9 `eq-l9-TBn`.

Acting with the two matrix parts on the singlet at $Y_1=qY_2$ cancels the charge-changing channels and leaves (q^s$1+qY_1Y_2$). At $s=0$, crossing reduces the scalar to
$$
\Phi(u+i\pi/2-i\pi/(2\lambda))/\Phi(u+i\pi/2+i\pi/(2\lambda)).
$$
Both arguments are in the convergence strip when $|\operatorname{Im}u|<\pi/2$. In the difference of logarithms the numerator is
$2i\sin(ku)\sinh(\pi k/(2\lambda))$, not a cosine. Therefore
$$
\log T_{B_1}=i\pi/2+i\int_0^\infty\frac{\sin(ku)}{k\cosh(\pi k/2)}\,dk,
\qquad T_{B_1}=\frac{i-\sinh u}{\cosh u}.
$$
The integral vanishes at zero and its derivative is $1/\cosh u$. The incoming cosine integral diverges logarithmically at zero, and its proposed final exponential gives (-1), rather than (i), at (u=0). These are fatal errors in that proposed proof, not in the book proof.

For $1\le n<\lambda$, the bootstrap recurrence
$$
T_{B_n}(u)=T_{B_{n-1}}(u+i\pi/(2\lambda))T_{B_1}(u-i(n-1)\pi/(2\lambda))
$$
inserts exactly the last factor into the shifted (n-1) product. This proves the displayed (n)-fold expression by induction, conditional on the usual breather bootstrap. Pairing the conjugate shifts proves real-axis unitarity; applying $T_{B_1}(w)T_{B_1}(w+i\pi)=1$ factor by factor proves crossing. A physical (B_n) requires $n<\lambda$; outside that range the same product is only a meromorphic formula.

Primary source: Bowcock–Corrigan–Zambon, [hep-th/0506169](https://arxiv.org/abs/hep-th/0506169), pp. 20–21, (5.10), (5.12), (5.13). The exact source line preceding the last two formulas is: “In summary, repeatedly using (5.3) one finds,”. The odd formula has ((-1)^{s+1}i/(1)) times paired reciprocal blocks; the even formula has ((-1)^s) times paired reciprocal blocks, exactly as obtained by pairing this product. **The incoming citation to (5.8) as the breather bootstrap is false: (5.8) is $S_{11}$; (5.9) defines the blocks.**

**Integration:** preserve the current B1 proof, change the book label to `**[Established; reformulated, with an analytic proof]**`, add `For each physical breather, $1\le n<\lambda$,` before the higher formula, and append to its proof:

```markdown
For the higher breathers, the fusion of $B_{n-1}$ at $\theta+i\pi/(2\lambda)$ with $B_1$ at $\theta-i(n-1)\pi/(2\lambda)$ gives the product by induction, for $1\le n<\lambda$. Pairing the shifts gives real-axis unitarity, and applying the crossing identity of the lightest breather to each factor gives the stated crossing identity. The resulting products are those of [@bcz2005, eqs. (5.10), (5.12)–(5.13)].
```

Do not replace the first-word label with Theorem merely because an identity between bootstrap amplitudes now has a proof. In Lecture 7 use the same Established first word for `thm-l7-breather`, retaining the stable id. Appendix E must describe the analytic reformulation/proof, not claim the known breather factors as new, and must cite the proposition id rather than its parent section.

### Q3 / lecture report 3.3: poles and widths

**Verdict: correct.** Targets: book `38-defect-sine-gordon.qmd:278`, `prp-sgd-pole-ladder`, `eq-sgd-widths`; Lecture 7 `thm-l7-poles`. Keep Theorem for the pole statement.

The zero-free strip from Q1 and
$\Phi(u)=[(1+e^{2\lambda u})\Phi(u+i\pi)]^{-1}$
give the continuation down to $-\pi-\pi/(2\lambda)$. The new denominator has exactly the simple zeros $u_n=-(2n+1)i\pi/(2\lambda)$, $0\le n<\lambda$, with derivative $-2\lambda$. Hence the residue is the stated $-[2\lambda\Phi(u_n+i\pi)]^{-1}$. The fixed-charge matrix has nonzero diagonal entries and determinant (q^{1+s}(1+Y^2)), so it has rank one, irrespective of whether its trace happens to vanish. (The incoming report's `trace != 0` is unnecessary and not generally true.)

At even (n), the Gamma argument $1/4-i\lambda u_n/(2\pi)=-n/2$ is a nonpositive integer. At odd (n), (Y_n=i) and (1+iY_n=0). All other lower-half-plane Gamma singularities start at or below the excluded lower boundary. This separately proves the origin and completeness of the two families.

The independent residue errors at $\varepsilon=10^{-4},10^{-6}$ are, for (n=0,1,2), respectively $9.86732\cdot10^{-5},9.86769\cdot10^{-7}$; $1.21481\cdot10^{-4},1.21482\cdot10^{-6}$; $7.65705\cdot10^{-5},7.65755\cdot10^{-7}$. Thus the singularity was tested through scaling, not a threshold.

Taking the imaginary part of $M\cosh(\theta_d-ia_n)$ gives $\Gamma_n=2M\sinh\theta_d\sin a_n$. For a physical breather $2n+1<\lambda$, this is $m_{2n+1}\sinh\theta_d$. At fixed (n), $m_{2n+1}/m_1\to2n+1$. The interpretation as decaying resonances needs positive real energy and width, as the book explicitly says; do not generalize it to every pole or to negative $\theta_d$.

### Q4: odd-jump modes and inverted oscillator

**Verdict: correct for nonzero rapidity, with an omitted zero-rapidity qualification in the proof and surrounding prose.** Target: book `38-defect-sine-gordon.qmd:309`, `prp-sgd-odd-jump`.

My independent matching matrix is
$$
\begin{pmatrix}\gamma+\kappa-m\sigma&-\gamma+\kappa-m\sigma\\
\gamma-\kappa+m/\sigma&\gamma+\kappa-m/\sigma\end{pmatrix}.
$$
On $\gamma^2=\kappa^2-m^2$, its determinant is $4\kappa[\kappa-m(\sigma+\sigma^{-1})/2]$. Reversing the vacuum sign replaces the minus by plus and excludes localized even-vacuum modes. At $\theta_d=0$, $\gamma=0,\kappa=m$, the **whole odd matrix vanishes**, giving two independent static amplitudes. There is no positive growth exponent there. At nonzero rapidity the positive and negative branches force respectively one of (A,B) to vanish, as stated for $\theta_d>0$.

For the nonlinear tail, $d(4\arctan e^z)/dz=2\operatorname{sech}z$ and $\sin(2\arctan e^z)=\operatorname{sech}z$; substitution into the two defect equations with $z=m(x\cosh\theta_d-t\sinh\theta_d)+c$ verifies emission. Expanding $4\arctan e^z=2\pi-4e^{-z}+O(e^{-3z})$ yields its growing tail. The oscillator algebra gives $[X,P]=i\hbar$ and $H=\gamma(XP+PX)/2$, so $H X^n=-i\hbar\gamma(n+1/2)X^n$. These generalized resonant functions do not prove a reduction of the field theory to one oscillator. The book correctly calls the match heuristic.

**Replacement:** in the proof add `for $\theta_d\ne0$` to both amplitude conclusions and append:

```markdown
At $\theta_d=0$ the matching matrix vanishes for $\gamma=0$ and $\kappa=m_0$, so $A$ and $B$ are independent static zero modes.
```

Replace the first sentence after the proof by:

```markdown
For $\theta_d\ne0$ the odd defect has one growing and one decaying localized mode, with rates $\pm m_0|\sinh\theta_d|$. At $\theta_d=0$ these are static zero modes. For the decaying-resonance discussion below we take $\theta_d>0$.
```

Keep Theorem for the linearized-mode statement; do not promote the oscillator interpretation.

### Q5 / lecture report 3.3: the two Weyl pairs

**Verdict: correct book algebra; wrong incoming derivation.** Target: book `41-defect-quantum-groups.qmd:191`, `prp-dq-two-weyl-pairs`; Lecture 7 `sec-l7-weyl`.

Write an exponent as its coefficient vector in $(\phi,\psi)$. The actual book has
$$
U=(-b/(2\sqrt2),b/(2\sqrt2)),\quad W_+=(b/\sqrt2,0),\quad
U_s=(-\sqrt2\pi/(\hbar b),\sqrt2\pi/(\hbar b)),\quad
X_+=(2\sqrt2\pi/(\hbar b),0).
$$
The determinant pairing $-i\hbar(a_1c_2-a_2c_1)$ gives $i\hbar b^2/4, 4i\pi^2/(\hbar b^2), i\pi, i\pi, 0, 0$ for ((U,W),(U_s,X),(U,X),(U_s,W),(W,X),(U,U_s)). This is the complete independent symbolic calculation in `check_quantum.py`.

At $b=i\beta_{SG}/\sqrt2$, $\hbar=1$, $q_s=e^{-8\pi^2i/\beta_{SG}^2}=q$, **not (-q)**. In a jump eigenstate, the exponent of (U_s) is $-4\pi^2i(2n+s)/\beta_{SG}^2$, giving (q^{n+s/2}), not its inverse. Further, $[\phi-\psi,X_+]=-4\pi X_+/\beta_{SG}$, so the jump increases by $4\pi/\beta_{SG}$. Squaring the first pair makes both cross-phases unity. Its multiplier is $e^{i\hbar b^2}=e^{2\pi i\tau}$, and that of ((U_s,X_-)) is $e^{-4\pi^2i/(\hbar b^2)}=e^{-2\pi i/\tau}$, $\tau=\hbar b^2/(2\pi)$.

The report substituted different exponent vectors, including real rather than imaginary sine-Gordon exponents, and asserted nonzero brackets between proportional vectors. Its printed evidence cannot support its verdict. The replacement evidence above does support the book's actual identities.

**Integration:** preserve formulas; use `**[Established; reformulated]**` for this canonical-commutator/BCH construction, consistently with the repository rule against theorem labels for convention derivations. Add at the start `As identities of formal exponentials with a central commutator,` so that no spectral-domain or Hilbert-space representation theorem is implied. The physical identification with the Lagrangian quantum theory remains the explicitly incomplete correspondence discussed after the proposition. Update Appendix E's first word if retaining the defect-specific reformulation there; do not claim the general modular-double construction as new.

### Q6: source corrections and Appendix D/E

**Verdict: mixed; existing three corrections supported, proposed Nemes expansion rejected as written.**

Primary sources read directly: Nemes [arXiv:0909.3268v1](https://arxiv.org/pdf/0909.3268), downloaded and text-extracted in this referee's scratch directory; Faddeev [math/9912078v1](https://arxiv.org/pdf/math/9912078), pp. 4–5; Bowcock–Corrigan–Zambon [hep-th/0506169](https://arxiv.org/abs/hep-th/0506169), local primary text. The following are source quotations, not quotations copied from the incoming report.

| Subclaim | Direct evidence and independent conclusion | Disposition |
|---|---|---|
| Nemes resonance count | His (67) is $n_p\sim\gamma/4$. Q3 independently counts both even and odd families, giving $#\{n\ge0:2n+1<\gamma\}\sim\gamma/2$. | **Correct.** Keep existing book Notes and Appendix D:165. |
| Nemes phase variable (Q) | The phase $Q^\alpha$ occurs in (33), while (34) is the scalar (f). A different (Q) was defined earlier as topological charge. | **Wrong as quoted/numbered.** Do not add the proposed `(34)` bullet. |
| Nemes (54) | The printed right side is $\ln((\eta-\theta)/2)$, and the left side is $\int_0^\infty\cos((\theta-\eta)t)\tanh(\pi t/2)\,dt/t$. There is no printed $-\gamma$. With the usual Abel prescription the integral is $\ln\coth(|\theta-\eta|/2)$, not the report's $-\gamma\ln\tanh((\eta-\theta)/2)$. | **Wrong quotation and proposed text.** Omit proposed expansion. |
| Nemes alleged (56)/(61) phase errors | Actual (56) is an odd-defect phase integral; actual (61) is $\theta<\eta\Longleftrightarrow\theta-\eta<0$. The phase manipulations occur in (58)/(63). For example (63) prints $e^{i\pi(1+\gamma)/4}/(1+i e^{\gamma(\theta-\eta)})=-e^{i\pi\gamma/4}$. Its elementary prefactor limit is $e^{i\pi/4}e^{i\pi\gamma/4}$, so that equality is erroneous, but it is not the report's quotation. | **Wrong as supplied.** Omit proposed four-equation bullet rather than introduce an unrequested source-correction expansion. |
| Claimed four misprints already in chapter Notes | Snapshot Notes, book lines 444–454, contain the pole-count correction but not these four extra corrections. | **Wrong.** There is no book/Appendix D synchronization debt for those absent claims. |
| Nemes narrow integral strip | His (49) has a positive (i) exponential with the kernel equivalent to Q1's (H(k)); its tail is $e^{-\pi k/(2\gamma)}$. | **Correct strip**, $|\operatorname{Im}(\theta-\eta)|<\pi/(2\gamma)$; incoming negative-sign quotation is wrong. |
| BCZ absolute width | Exact source sentence: “In the classical limit the width goes to zero so the only process is absorption.” Since $M\sim(\gamma+1)m_0/\pi$, $-2\operatorname{Im}\Delta E=2M\sinh\eta\sin(\pi/(2\gamma))\to m_0\sinh\eta$. | **Correct existing clarification**, distinguishing absolute units from units of (M); keep book Notes:449 and Appendix D's BCZ entry. |
| Faddeev normalization | Exact formulas (14), (18): $[Q,P]=2\pi i\hbar I$, $\alpha\tilde\beta=\tilde\alpha\beta=\hbar$. Independent BCH gives cross-phase $e^{-2\pi i\hbar\alpha\tilde\beta}$, so the printed choice leaves $e^{-2\pi i\hbar^2}$. Integer multiples of $1/\hbar$ are required. Opposite unit integers give the book's inverse modular parameter. | **Correct existing clarification**, book:350 and Appendix D:85. |
| Appendix E labels “already match” | The B1 item refers to a section and the incoming report treats partial labels inconsistently. Independent formulas do not excuse label/reference mismatch. | **Wrong as blanket audit.** Coordinator's Appendix E audit governs; use the statement ids and final labels above. |

The Nemes quotation errors are not merely line-number drift: their mathematical content is different. No new Appendix D entry is authorized by this report for those four purported misprints. This is a completed rejection, not a request for Claude to investigate them again.

### Q7a / lecture report 5.1(1),(3): chirality and the second Borel

**Verdict: incoming correction superseded; current text correct.** Target: book `43-defects-and-boundaries.qmd:110`, `prp-diagonal-borel`; Lecture 9 `sec-l9-borel`, lines 112–136.

With $\phi=L(t+x)+R(t-x)$, $\psi=l(t+x)+r(t-x)$, the massless equations $\phi_x=\psi_t, \psi_x=\phi_t$ give (l'=L'), (r'=-R'). Thus $\phi+\psi=2L$, $\phi-\psi=2R$, up to constants. This is an independent two-equation elimination. The current book already says the passed left-moving half has untwisted (f_j), and the conjugated right-moving half has twisted (e_j). The limits also agree: $\theta_d\to+\infty$ means $Y\to0$, whereas $\theta_d\to-\infty$ means $Y\to\infty$. There is no sign flip to make.

For the second module, direct application of (B^{-1}) then (C) gives $CB^{-1}|n\rangle=q^{-2n-2-s}|n+2\rangle$. Its (e_1) coefficient (a q^{-2n}), (a=iy_dq^{-s}/(q^2-1)), has difference (-iy_dq^{-2n-2-s}), as required. The (e_0) partner gives (BC^{-1}), their product is one, conjugation weights are $q^{\pm2}$, and (e_1e_0=e_0e_1=abq^2); hence the Serre relation cancels term by term. This verifies the substantive “second Borel” claim without interpreting it as one combined coideal. Indeed $\rho(k_1)\tilde\rho(k_1)=q^4\tilde\rho(k_1)\rho(k_1)$, so they cannot be images of the commuting ambient Cartans in one representation.

The incoming lecture report names `thm-l9-doubled` and `prp-l9-general`, which are not ids in this snapshot. Do not create statements under those invented/stale ids. Match the actual discussion and the actual `thm-l9-crossing`, `thm-l9-fusion`, `thm-l9-dirichlet`.

### Q7b / lecture report 5.1(2): fusion onto a Dirichlet boundary

**Verdict: correct current statement and proof; wrong and superseded incoming replacement.** Target: book `43-defects-and-boundaries.qmd:143`, `prp-defect-fusion-k`; Lecture 9:285 `thm-l9-dirichlet`.

Independently invert (A(Y')) with $A(v)=\begin{pmatrix}1&-iv\\-iv&1\end{pmatrix}$, and multiply $A(Y')^{-1}\operatorname{diag}(k_+,k_-)A(Y)$, (Y=y/y_d), (Y'=1/(yy_d)), $k_\pm=ye^{\pm i\eta}+y^{-1}e^{\mp i\eta}$. The resulting entries are exactly
$$
\frac1{1+y^{-2}y_d^{-2}}
\begin{pmatrix}
k_++y_d^{-2}k_-&iy_d^{-1}e^{i\eta}(y^{-2}-y^2)\\
iy_d^{-1}e^{-i\eta}(y^{-2}-y^2)&k_-+y_d^{-2}k_+
\end{pmatrix}.
$$
Comparing $\cosh\vartheta=(y_d+y_d^{-1})/2$ and $\sinh\vartheta=(y_d-y_d^{-1})/2$ proves $\vartheta=\lambda\theta_d$, with $\eta$ unchanged, up to the stated diagonal conjugation. This is a sectorwise general GZ matrix; it is not a scalar on all of $\mathcal W$.

In particular, for $K=1$, the off-diagonal coefficient is $i(Y'-Y)/(1+Y'^2)$, independently computed as `-I*d*y*(y-1)*(y+1)/(d**2*y**2+1)`. It is generically nonzero. A diagonal change of basis cannot make a nonzero entry zero. This directly refutes the incoming paste-ready proof.

For the scalar, sum the four Gamma exponentials in the GZ product: the (l)-th sum is $-4e^{-t/2}e^{-(2l-1)\lambda t}\sinh(\lambda t/2)\sinh(xt/\pi+\lambda t/2)$. The geometric series gives the integral `eq-dbd-sigma-integral`. At $x=i\lambda\theta_d,u=-i\theta$, subtract the two Phi logarithms. The difference of the two sine integrals is
$$
-2\int_0^\infty\frac{\sin(k\theta)\sin(k\theta_d)}{k\sinh(\pi k/(2\lambda))}\,dk
=\log\cosh\lambda(\theta-\theta_d)-\log\cosh\lambda(\theta+\theta_d).
$$
The elementary prefactors then cancel, proving `eq-dbd-scalar-identity`. Independent quadrature agrees to $4\cdot10^{-28}$. The limit $\sigma(i\vartheta,-i\theta)\to1$ follows only **after** evaluating the singular sine term; the book's current proof does that correctly.

The current book already labels (1),(3) Theorem, leaves (2) Numerical and makes (4) conditional on UV–IR matching. **Keep this scoped label, not the incoming blanket Theorem.** The lecture's matrix/scalar theorem can stand. Do not promote its finite-coupling UV–IR conjecture or interpret fusion as deriving the physical defect from its Lagrangian.

### Q7c / lecture report 4.2–4.3: traces, commutation and TQ

**Verdict: correct finite-chain t–D and TQ proofs after a small continuation qualification; unclear proposed D–D proof; wrong blanket upgrade of the full transfer proposition.** Targets: book `41-defect-quantum-groups.qmd:228–262`; Lecture 8:181,239,284.

**Certified theorem hypotheses.** Take a finite chain with $N\ge1$, real irrational $\lambda>0$, and the fixed physical phase $q=-e^{-i\pi\lambda}$, so $q$ is not a root of unity. Fix the logarithm $\log q=i\pi(1-\lambda)$ and use this same logarithm for every power of $q$, including $q^{1/2}$, $q^s$ and the sector twist. Let the offset $s$ be fixed and real. All spectral parameters $y_1,\ldots,y_N,y_d$ are nonzero complex numbers. Set $Y_j=y_j/y_d$. For t–D commutation assume $1+Y_j^2\ne0$ for every site; the auxiliary spectral parameter is also nonzero and the transfer matrix is evaluated away from its poles. For the displayed TQ identity additionally assume $Y_j^2\notin\{-q^2,-q^{-2}\}$ for every site, so both shifted inverse-L products and the fused transfer matrix are defined. In charge sector $m$, the commutation twists are $\omega=\pm q^{-m/2}$; the written TQ formula uses $\omega=q^{-m/2}$. These are the hypotheses to include explicitly in the upgraded book statements. The algebraic identities extend meromorphically where removable limits exist, but such extensions are not needed for the certified theorem as stated.

For two banded scalar auxiliary matrices (A,B) with bounds $C_A,C_B$, all terms of a truncated trace commutator cancel except paths crossing a cutoff edge. Their number is bounded independently of the cutoff; division by (2K+1) makes the difference zero. This proves trace-per-state cyclicity when the limit exists. For physical-space matrix entries cyclicity is only used to move an auxiliary-only matrix such as $X$, not arbitrarily to commute two operator-valued factors. $X,L_j^{-1}$, and the fusion basis transformations are banded and bounded for phase $q$, fixed $s$, and nonzero spectral denominators. Thus no unbounded two-defect intertwiner is hidden in the t–D or TQ proof.

The explicit inverse entries used in the independent path computation are
$$
L^{-1}|+,n\rangle=\frac{q^n|+,n\rangle+iYq^{-n-1-s}|-,n+1\rangle}{1+Y^2},
\quad
L^{-1}|-,n\rangle=\frac{q^{-n-s}|-,n\rangle+iYq^{n-1}|+,n-1\rangle}{1+Y^2}.
$$
Shifting the auxiliary origin multiplies each row by $q^{\text{outgoing spin}}$, proving $\langle n|M|n\rangle=q^{mn}D$. The Cesàro trace selects $\omega^2q^m=1$, and ordinary bilateral trace convergence fails for every nonzero twist because one of the two tails fails the term test. For (X) invertible, the exact RTT relation proved in `check_triangle.py` gives (X T_aM=MT_aX). The auxiliary twist commutes with (X), so bounded cyclicity proves $[t_\omega,D]=0$. At $1+e^{2\lambda(\theta-\theta_d)}=0$, continuity extends the identity wherever $t_\omega$ is finite. This is the minor completion identified by the incoming report.

At fusion (Y_a=i), write (t=q^n,r=q^s). The kernel vector is $e_n=|+,n\rangle-q^{-2n-1-s}|-,n+1\rangle$. In the bases (q^ne_n) and $q^{-n}[|-,n\rangle]$, independent symbolic action of (F_j=R_{aj}L_j^{-1}) gives, **for all four input/block cases**, precisely
$$
F_j|_S=\frac{1+q^2Y_j^2}{1+Y_j^2}L^{-1}_{s+1}(qY_j),\qquad
F_j|_{\mathrm{quotient}}=q^{-1}L^{-1}_{s-1}(Y_j/q).
$$
The change of basis and its inverse each shift (n) by at most one and have bounded entries, so the normalized trace splits into the two diagonal blocks. Their twists are $\omega^{2n+s\pm1}$. Since $D_{s\pm1}=q^{\mp(N-m)/2}D_s$, division by $\omega^s$ gives (q^{-N/2}) for both branches. This is a complete all-$N$ proof of the displayed TQ relation, under the explicit certified hypotheses above, including $Y_j^2\notin\{-1,-q^2,-q^{-2}\}$. Exact local identities, not the numerical N≤5 check, support the theorem label.

The statement that neutral paths have “the same power” in two defect orders does not exhibit a bijection of intermediate physical spins and two independent defect paths. Neither bounded cyclicity nor one-soliton fusion supplies that missing proof. Independent numerical D–D tests pass for N=2–5, but do not establish all N. Keep `prp-l8-dd` Numerical and retain the corresponding book clause as Numerical. Leave the Fock statement at Proposition/sketch with its trace/intertwiner qualification.

**Precise integration actions:**

- In Lecture 8, replace the proof of `thm-l8-commute` with the existing proof plus the incoming report's one-sentence invertibility/continuation qualification; use the complete condition $e^{2\lambda(\theta-\theta_d)}\ne-1$, not just two selected rapidities.
- At book `prp-defect-transfer`, split the label and clauses explicitly: Fock remains `**[Proposition, sketch; numerically checked]**`; bilateral trace factorization and t–D commutation become `**[Theorem]**` under the explicit certified hypotheses above; D–D becomes `**[Numerical]**`. Keep the id for the combined proposition if convenient, but do not put a single Theorem label over all claims.
- Replace the bilateral formal-trace argument with the bounded trace-per-state proof above, retaining the distinction from the divergent ordinary trace. For the positive Borel counterpart the shift exponent is (-mn), so the selected twist is $\omega^2=q^m$; inversion of all phase exponents gives the same bounded argument.
- Change `cnj-defect-tq` to `**[Theorem]**` with the explicit assumptions above. Preserve its id to avoid reference churn. Copy the current Lecture 8 `thm-l8-tq` proof, including the two verified block formulas and trace-basis justification, into the book. It is a finite-chain algebraic theorem, not a theorem about existence of the continuum QFT.
- In the book status table distinguish t–D/TQ (Theorem), D–D (Numerical), and Fock (sketch). Do not add a new TQ originality claim to Appendix E. Do not use incoming Item 7c's blanket replacement of Appendix E.
- The spectral identities $q^{\pm1}y_d=-e^{\lambda(\theta_d\mp i\pi)}$ follow directly from the physical q and (D(-y_d)=D(y_d)) follows because closed paths have equally many raising and lowering conversions. They identify **scalar-free** shifted modules. Restoring scalar factors and physical resonance data is additional information. Preserve the book's cautious qualification at line 262. In Lecture 8 replace its emphatic physical identification after `thm-l8-tq` with:

```markdown
The kernel and quotient therefore have the spectral parameters of the two odd-sector shifts. This proves the scalar-free operator relation. Identifying its terms with the full physical defect bootstrap also requires the scalar transmission factors and their pole normalization.
```

### Lecture report 3.1–3.2: Borel obstruction, uniqueness and triangle

**Verdict: correct, independently proved.** Targets: Lecture 7 `thm-l7-half`, `thm-l7-T`, `prp-l7-triangle`.

For non-root-of-unity q, the (k_1)-eigenspaces are one-dimensional. A putative (e_1) has $e_1|n\rangle=a_n|n+1\rangle$. The relation ([e_1,f_0]=0) makes (a_n) constant, and then ([e_1,f_1]=0), contradicting ((k_1-k_1^{-1})/(q-q^{-1})) on all but at most one state. This is a complete algebraic contradiction; retain Theorem.

For uniqueness, the weight equation gives four channels. Solving the elementary shift equations gives $A_n=a q^{-n}+Hq^n$, $D_n=q^{2n+s}A_n$; the remaining equation has obstruction $cHq^n(q^4-1)/q^3$. Therefore $H=0$ for generic q. The remaining amplitudes are $q^{-n},q^{n+s},-iYq^{-n},-iYq^{n+s}$, up to a scalar. `check_local_blocks.py` computes this obstruction symbolically. `check_triangle.py` then checks RTT exactly on every two-soliton input with unconstrained $q,Y,Z,q^n,q^s$, so the triangle conclusion need not rely on an unproved uniqueness assertion for the larger tensor product. Retain the current labels.

### Lecture report 3.3: crossing, Weyl phases, Wigner asymptotics

**Verdict: correct actual formulas; wrong phase shorthand in the incoming report.**

The inverse-L formulas above prove all four crossing equations after setting $Y'=e^{-i\pi\lambda}Y=-qY$: multiplication by (q^s(1+q^2Y^2)) matches the original four amplitudes. Changing to (Y'=+qY) leaves the diagonal entries and reverses the converting ones, so a common scalar cannot repair the sign with symmetric charge conjugation. This is `check_triangle.py`'s exact calculation. Retain `thm-l7-crossing` as a statement in the fixed conjugation conventions.

The lecture's actual canonical pair is $U=e^{2\pi i\Delta/(\hbar\beta^2)}$, $X_\pm=e^{\mp4\pi i\Phi/(\hbar\beta^2)}$. Its multiplier is $e^{-8\pi^2i/(\hbar\beta^2)}$. The incoming report's `-exp(-4*pi^2*i/(hbar*beta^2))` does not use this lecture's $\beta=\beta_{SG}$ and must not be copied. Q5 independently proves the correctly normalized brackets.

For the Wigner slope, put $L=\lambda$:
$$
\delta'(0)=\frac L\pi I(L),\quad I(L)=\int_0^\infty\tanh(Lt)/\sinh t\,dt.
$$
Subtract the sharp cutoff at (t=1/L). The hyperbolic-denominator constant is $\log2$; the smoothing constant is $\gamma_E+\log(4/\pi)$. One direct derivation of the latter uses
$\int_0^\infty x^s\operatorname{sech}^2x\,dx=2^{1-s}\Gamma(s+1)\eta(s)$, where $\eta(0)=1/2$, $\eta'(0)=\frac12\log(\pi/2)$; integration by parts gives the constant as minus its derivative at zero. Expanding $1/\sinh t=1/t-t/6+O(t^3)$ in the smoothing region makes the remainder (O(L^{-2})). Thus
$\delta'(0)=L[\log(8L/\pi)+\gamma_E]/\pi+O(L^{-1})$.
Independent errors at L=80,160,320 are 0.00027269911, 0.00013635277, 0.00006817679. This verifies the asymptotic including its constant, not just the leading logarithm. Keep the lecture formula.

### Lecture report 4.1: turned defect

**Verdict: correct lecture; incoming verbal sign explanation should not be copied.** Target `08-q-operators.qmd:56`, `thm-l8-turned`.

The rotated equations require (G_a=iB). Comparing coefficients of $\cos(\Sigma/2)$ and $\cos(\Delta/2)$ gives respectively $a=-i\sigma$ and $a^{-1}=i/\sigma$, the same condition. The incoming phrase “replaces sigma by i sigma, i.e. B_{-i sigma}” conflates the Lax-matrix sigma with the Bäcklund parameter; the current lecture distinguishes them. The generating function $F=\beta^{-2}\int(\Phi\Psi_x-iB)dx$ gives $\delta F/\delta\Phi=(\Psi_x-iB_\Phi)/\beta^2$ and $-\delta F/\delta\Psi=(\Phi_x+iB_\Psi)/\beta^2$, exactly the rotated matching conditions. Finally differentiation of $G(x)=T^\Phi(x,-\ell)L(-\ell)-L(x)T^\Psi(x,-\ell)$ gives $G'=-a_x(\Phi)G$ and $G(-\ell)=0$; periodicity of L yields conjugate monodromies away from its isolated determinant zeros, and their traces extend by spectral continuity. These steps independently prove the stated classical result. No change required.

### Lecture report 4.3: Baxter polynomial and Bethe roots

**Verdict: correct polynomial/eigenvalue reduction; root equation needs regularity hypotheses.** Target `08-q-operators.qmd:284`, `prp-l8-baxter`.

Each conversion in $\prod_j(1+Y_j^2)D$ contributes one (y_d^{-1}). Closed auxiliary paths have equal raising and lowering counts, at most the number of minority spins. This proves degree at most ((N-|m|)/2) in (y_d^{-2}). A simple eigenvalue of one transfer matrix defines a one-dimensional invariant subspace for every commuting matrix, proving the eigenvalue TQ relation. But evaluating its rational coefficients at a zero requires that both that coefficient and the transfer eigenvalue are finite, and division by (Q(qy_k)) requires it nonzero. The incoming blanket endorsement omitted this issue.

Replace the root sentence of the proposition with:

```markdown
At a nonzero zero $y_k$ of this eigenvalue, assume that the transfer-matrix eigenvalue and the coefficient $\prod_j(1+z_j^2)/(1+q^{-2}z_j^2)$ are finite there, and that $\mathcal Q(qy_k)\ne0$. Then
```

Keep its displayed Bethe equation and proof, appending `These assumptions justify evaluation at the zero and division by $\mathcal Q(qy_k)$.` Do not conclude that every polynomial has the **maximum** degree merely from an upper bound: replace the subsequent sentence “Their number is the number of solitons of the minority charge” by `Their number is at most the number of solitons of the minority charge; the generic examples below attain that bound.` Keep Theorem for the qualified polynomial statement.

### Lecture report 5.1: the actual Lecture 9 statements

**Verdict: substantive formulas correct, but the original report did not review the actual Lecture 9 file and cannot certify its stated ids.** I checked the actual statements as follows.

- `thm-l9-crossing` (line 80): in the doubled boundary only mixed-copy pair amplitudes are nonzero. The mixed-copy bulk operator is a flip, so the boundary crossing component is $T_{c\alpha}^{\beta b}(\theta)=\bar T_{\alpha\bar b}^{\bar c\beta}(i\pi-\theta)$. It leaves both defect labels in place. The same-copy components vanish; the two diagonal blocks of $K(\theta)K(-\theta)$ are precisely the two defect-unitarity equations. This is a complete component proof. Retain the statement and its complete component proof; label it `**[Established; reformulated]**` as a reformulation of the known crossing conventions.
- The second-Borel discussion is verified in Q7a; it does not assert a full diagonal quantum-group symmetry.
- `thm-l9-fusion` (line 152): in the scalar case substitute $T(i\pi/2-\theta)=\bar T(i\pi/2+\theta)$ and its partner into $\bar TKT$; scalar commutativity extracts $S(2\theta)$, giving the boundary crossing equation. The product at opposite rapidities cancels the two defect-unitarity pairs and K-unitarity. This is an exact algebraic consequence of the assumed axioms; the current proof is complete. Label this known dressing identity `**[Established; reformulated]**`, retaining the complete proof and stable id.
- `thm-l9-dirichlet` (line 285) is verified in Q7b, including its scalar, with its theorem label scoped to the specified matrices. Its UV–IR sequel remains Conjecture.
- `thm-l9-breather` (line 413): put $g(x)=\log[-\coth((u+ix)/2)]$, on a continuous branch for $\alpha\le x\le\pi-\alpha$. There are no zeros or poles on this closed interval for real u and $0<\alpha<\pi/2$. The midpoint rule proves $(\pi/\lambda)\sum g\to\int g$. Applying it to $\partial_u g=-i\partial_x g$, and using $z(u+i(\pi-\alpha))=1/\overline{z(u+i\alpha)}$, yields the phase-derivative limit $2\log|z(u+i\alpha)|$. Comparing the endpoint trapezoidal sum for (n+1) with the midpoint sum for (n) leaves $[g(\alpha)+g(\pi-\alpha)]/2+O(1/\lambda)$, whose exponential is the asserted clock factor. Independent clock errors at n=40,80,160 are 0.0082004112, 0.004101471, 0.0020508938. Dividing by $m_n\sinh\theta$, with $m_n=2(\lambda+1)m_0\cos\alpha/\pi$, gives the classical delay wherever $\sinh\theta\ne0$. Retain Theorem as the mathematical product-asymptotic statement. Add `with a continuous branch of the logarithm along the string` to its first limit, and `for $\sinh\theta\ne0$` to its delay conclusion. Do not read it as a proof of the full quantum-to-classical field-theory correspondence.

## Exact Appendix E replacements

These refer to `appendices/E-original-results.qmd` in the reviewed snapshot. They replace only the named bullets; the other bullets remain outside this packet unless separately verified by the coordinator.

**Lines 43–44:** replace the scalar/breather combined bullet and the odd-saddle bullet with these three bullets. The first two separating labels prevent the existing section reference from concealing a mismatch.

```markdown
- **[Theorem]** An integral representation of the chosen minimal sine-Gordon defect scalar factor in a strip of width $\pi+\pi/\lambda$, and the complete ladder of its nearby poles at $\theta_d-(2n+1)i\pi/(2\lambda)$ for $0\le n<\lambda$. In the decaying-resonance range their widths are $2M\sinh\theta_d\sin((2n+1)\pi/(2\lambda))$. Bowcock, Corrigan and Zambon discussed the first pole, and Nemes gave the integral near the real axis (@prp-sgd-scalar-integral; @prp-sgd-pole-ladder).
- **[Established; reformulated]** An analytic proof, using the scalar integral, that the lightest sine-Gordon breather has the classical transmission factor at every coupling where it exists. Breather fusion gives the higher factors for $1\le n<\lambda$. The transmission factors themselves were obtained by Bowcock, Corrigan and Zambon (@prp-sg-breather-transmission).
- **[Theorem]** For $\theta_d\ne0$, the odd sine-Gordon defect has one growing and one decaying localized linear mode, with rates $\pm m_0|\sinh\theta_d|$. At $\theta_d=0$ the two localized modes are static. For $\theta_d>0$ the fixed-$n$ resonance widths tend to $(2n+1)m_0\sinh\theta_d$; identifying this whole ladder with an inverted oscillator remains heuristic (@prp-sgd-odd-jump).
```

**Line 47**, the two-Weyl-pair bullet: replace with:

```markdown
- **[Established; reformulated]** The formal canonical algebra of the type-I$_-$ defect values realizes two Weyl pairs, with multipliers $\hat q=e^{i\hbar b^2/4}$ and $q_{\rm s}=e^{4\pi^2i/(\hbar b^2)}$. Their cross-relations are signs, removed by squaring the first pair. This gives the defect-specific relation to Faddeev's modular-double construction; it does not establish a canonical quantization of the full defect field theory (@prp-dq-two-weyl-pairs).
```

**The fusion bullet at line 50:** preserve its existing label and its statement of the proved scalar identity. It already states the correct $\vartheta=\lambda\theta_d$ and attributes the matrix construction. Do not replace it with PR60's nonexistent scalar-on-$\mathcal W$ Dirichlet claim. Its separate fused-spin-three clause is not newly certified by this quantum packet and is not changed here.

**TQ:** there is no old transfer/TQ bullet in this snapshot, despite the original report's purported quotation. **Do not add a TQ item to Appendix E in this integration.** The independent proof certifies the mathematical theorem, not originality of this particular proof or normalization. The chapter itself reserves novelty pending a detailed comparison with Rossi and Weston. Preserve its attribution of the bilateral factorization to Rossi and Weston and of the TQ mechanism to Bazhanov, Lukyanov and Zamolodchikov and Weston. The completed conservative disposition is to upgrade the finite-chain theorem and its chapter status table, with the proof and assumptions in Q7c, while making no new original-result claim. No further literature investigation is required for that disposition. Do not add original-result claims for D–D commutation, Fock TQ, or the physical bootstrap interpretation either.

**Ch41 status table:** line 294's final cell becomes `Fock: [Proposition, sketch]; Weyl t–D commutation: [Theorem]; Weyl D–D commutation: [Numerical]`. Line 295's final cell becomes `Fock: [Numerical]; Weyl finite-chain TQ: [Theorem]`. Keep all qualifiers in prose outside status-label brackets if cross-references are added.

## Integration ledger and checks for Claude

| Incoming item | Final disposition | Necessary destinations |
|---|---|---|
| Quantum 1 | Verified formula; reject sign and Riemann–Lebesgue mistakes | No book formula change; qualify Lecture 7 minimality sentence |
| Quantum 2 | Modified: preserve existing proof, add induction and physical n-range; Established label | Ch38, Lecture7 breather label, Appendix E |
| Quantum 3 | Verified | No formula change; preserve physical-resonance qualifications |
| Quantum 4 | Modified: zero-rapidity qualification | Ch38 proof and post-proof prose |
| Quantum 5 | Modified: actual book formulas verified, report derivation rejected; formal BCH/Established label | Ch41, Appendix E |
| Quantum 6 | Existing three source corrections verified; extra Nemes bullet rejected | Preserve existing D entries; no new Nemes expansion |
| Quantum 7a | Superseded/no-op | Preserve current Ch43 and Lecture9 chirality |
| Quantum 7b | Superseded/no-op; incoming replacement rejected | Preserve current Ch43 partial label and proof, Appendix E scope |
| Quantum 7c | Modified: prove t–D/TQ; D–D remains Numerical; Fock remains sketch | Ch41, its status table (no new Appendix E item); Lecture8 continuation and bootstrap qualification |
| Lecture 3.1–3.2 | Verified | No formula change |
| Lecture 3.3 | Verified actual formulas; qualify minimality | Lecture7 |
| Lecture 3.4 | Modified label; correct existing proof retained | Lecture7 and Ch38 in sync |
| Lecture 4.1 | Verified | No change |
| Lecture 4.2 | Modified continuation qualification; D–D Numerical retained | Lecture8 and Ch41 in sync |
| Lecture 4.3 | Modified root regularity and degree wording | Lecture8 |
| Lecture 5.1 | Replace stale-id endorsement with actual statement audit above | Lecture9 small qualifications; no phantom ids |

The only all-N upgrade cleared here is the **finite-chain bilateral Weyl t–D/TQ statement**, with its explicit hypotheses and supplied full proof. It does not clear Fock TQ, mutual D–D commutation, the defect UV–IR relation, physical scalar bootstrap identification, or a canonical quantization of the full field theory.

When integrating the corrected formulas, add regression checks for the zero-rapidity matching matrix and nonzero Dirichlet-fusion off-diagonal entry. The three new independent scripts provide directly reusable check logic, but repository check scripts are to be written/adapted by Claude and indexed in `code/README.md` or `code/lectures/README.md`. Book prose must cite those production paths, not `briefs/` or `referee_scratch/`. For the unchanged proven identities, the existing check scripts already cover the formulas; do not duplicate their whole implementations.

Run the affected `code/defect_sg_scalar_poles.py`, `code/defect_transfer.py`, `code/defect_fused_boundary.py`, and relevant Lecture 7–9 check scripts after editing; update their README rows only if checks change. Run `python3 check-refs.py` and `python3 check-refs.py lecture-notes/defects`, then scratch renders of book and the Claude lecture project. These are integration/regression checks, not a request to re-referee the mathematical conclusions above. No figure regeneration is required by this report.

**Independent quantum verification pass complete.**
