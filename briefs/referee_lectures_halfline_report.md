# Half-line lectures: directed phase-2 referee report

- Date: 2026-10-25.
- Referee: Codex (GPT-6; exact runtime model identifier not exposed).
- Snapshot: the brief identifies `035351d`. `git rev-parse HEAD` again fails because this export has no `.git`; the hash is supplied, not independently certified.
- Phase 1: [frozen report](referee_lectures_halfline_phase1_report.md), hereafter **P1**. Its SHA-256 on beginning this phase's written report was `de6f2ffaa8f8543bca62d67e1b2c9e8b535c834dfce5edb2738710525158de4b`.
- No lecture, book, existing script or figure was edited. No Git mutation, Quarto render, author figure-script execution, or subagent was used. Earlier referee verdicts were not used as evidence. This report contains recommendations, not integrated changes.
- L1–L6 mean the numbered files in `lecture-notes/half-line/`; Index means `index.qmd`. Locations refer to this export. P1 entry numbers are also their ranks.
- New executable evidence: `referee_scratch/lectures_halfline/phase2_checks.py` and `phase2_additional.py`, with corresponding `.log` files. They import no author helper. Algebraic derivations below are additional independent evidence. A checked formula does not certify its physical identification, all-rank validity or historical attribution.

## 1. Status, scope and the CDD-dependent duality

**Overlap:** P1 risks 5, 17, 18, 2, 24 and 25; all checked, with explicit limits. P1 §§2.5, 2.17–18 establish the central scope and status problems. This phase expands the sentence inventory and independently tests the CDD ambiguity. **Verdict: wrong as an account of established physical identifications; many displayed identities are correct.**

### 1.1 Sentence-level status inventory

The following inventory covers the affected assertions in the index and all six lectures, including captions and teaching summaries. Multiple quotations in a row identify repeated assertions requiring the same qualification. “Numerical” is the book's label, not a claim that an elementary identity can never be proved. In particular, a repaired proof can legitimately establish the Lee–Yang sign or an explicit matrix identity; it cannot establish a Lagrangian identification merely by calling the formula exact.

| Location and assertion | Book statement / status | Replacement text |
|---|---|---|
| Index:12, rigidity “governs all higher-rank affine Toda theories” | `prp-durham` **[Numerical]**, with an analytic simply-laced necessity result; `tbl-durham` has continuous non-simply-laced branches | “The $a_2^{(1)}$ example exhibits the rigidity of the simply-laced $a_n^{(1)}$ boundaries. Other affine types can retain continuous parameters.” |
| Index:30, couplings “cannot be tuned continuously at all” and Dirichlet is “never integrable” | `prp-durham`; `prp-no-dirichlet` **[Theorem]** with specified representations and field-independent $K$ | “Within the stated field-independent Lax construction, the $a_n^{(1)}$, $n\ge2$, boundary potentials have discrete nonzero couplings. The Dirichlet obstruction applies in finite-dimensional evaluation representations without trivial summands.” |
| Index:34, “a standard $a_n^{(1)}$ boundary reflects ... into its antisoliton”; stationary pair becomes “the exact ... background ... needed” | `prp-antisoliton-return`, `prp-solitonic-vacua` **[Numerical]**; the latter does not classify all ground states | “Parity conjugates the Hirota species. The mirror construction gives antisoliton reflection for the boundaries considered here and produces explicit regular stationary backgrounds on the allowed modulus intervals.” |
| Index:46, “we will resolve both puzzles”, WKB checks, “literally exchange places” | `cnj-neumann-coideal`, `cnj-neumann-amplitude`, `prp-dm-first-order` **[Conjecture]**; `prp-semiclassical-uniform` **[Numerical]** | “We will examine the proposed boundary/coideal dictionary and the amplitudes exchanged by $B\leftrightarrow2-B$. Their physical interpretation depends on the continuation and CDD choice. The book supplies numerical semiclassical evidence; these lectures do not perform a WKB calculation.” |
| Index:46, “only half ... admit real non-singular mirror-soliton vacua” | `prp-solitonic-vacua` plus `prp-self-conjugate-vacuum`; missing patterns remain | “The mirror-pair construction yields particular solitonic sign patterns. It is not a classification of all regular vacua.” |
| L1:353,404, “minimal solution” and sign conclusion; L1:464–485, family and excited shift | `prp-ly-reflection` **[Numerical]** | “The following factors solve the three scalar consistency equations, as the block calculation shows. Their assignments to the identity and perturbed conformal boundaries use the cited physical evidence.” Retain the sign and shift after repairing the proof; P1 §§2.3–4. |
| L2:4,176,313, the calculation proves higher-rank rigidity / “completely ... classically integrable” | `prp-durham` **[Numerical]** | “The calculation proves the parameter restriction in the defining $a_2$ representation after the mixed-branch step is included. The broader classification is summarized in the book with its stated evidence.” |
| L2:73, a “true ground state ... must be” a nonconstant background | `prp-solitonic-vacua` **[Numerical]**, not an existence/minimality theorem for every pattern | “If a ground state exists for such a pattern, it cannot be the zero field. Lecture 3 constructs stationary backgrounds for some patterns.” |
| L2:325, “never integrable in Sklyanin's framework” for all $\hat g\ne a_1^{(1)}$ | `prp-no-dirichlet` **[Theorem]** | Use P1 §2.5 replacement, including untwisted algebras, field-independent matrices and representation restrictions. |
| L3:4, “we proved ... an $a_n^{(1)}$ ... is classically integrable”; generic wall is “secretly a whole-line two-soliton solution” | `prp-durham`, `prp-images` **[Numerical]** | “Lecture 2 treated $a_1$ and $a_2$. We now study the explicit mirror-pair solutions of the uniform and Neumann $a_n$ boundaries.” |
| L3:55, “every standard boundary condition in Table 31.1” returns the antispecies | `prp-antisoliton-return` **[Numerical]**, checks $a_2,a_3$ and specified $C$ | “For the mirror construction considered here, parity sends $(a,\theta,\xi)$ to $(h-a,-\theta,-\xi)$ up to a vacuum shift.” This is an analytic statement actually derived here; it is not a theorem about every entry of the table. |
| L3:118–127, identity “for every $j$ and every $t$”, three-line proof and verification | `prp-images` **[Numerical]** | “With the corrected $k$ and a consistent square-root/vacuum branch, the following algebra checks the mirror boundary condition.” P1 §2.10 supplies the corrected steps. |
| L3:191–243, “remains completely free”, properties of the vacuum and sign pattern | `prp-solitonic-vacua` **[Numerical]** | “The displayed tau functions are positive on $x<0$. A regular wall excludes the endpoint phases, and a fixed boundary potential selects one open interval of $\chi$.” Use the corrected intervals in P1 §2.1. |
| L3:238, product-one patterns called “the solitonic sign patterns” without excluding the uniform pattern | `prp-dg-weights` **[Proposition]**, definition in `sec-real-coupling-vacua` | “Solitonic patterns have product $+1$ and are not all $+1$; this construction need not realize every such pattern.” |
| L4:8, “complete ... reflection matrix” from eight linear equations; :253 scalar “fixed” by unitarity/crossing | `prp-sg-k-unique` **[Theorem]** for matrix part; scalar has CDD freedom | “The linear equations determine the matrix part up to a scalar. Unitarity and crossing constrain that scalar, leaving CDD freedom.” |
| L4:267–280, UV–IR relation “states”, “what this implies”, “at the Neumann point ... equals”, zero “is ... uniform” | `cnj-bsg-uv-ir`, `cnj-neumann-coideal` **[Conjecture]** | Begin “Assuming the conjectured UV–IR relation...” and qualify its consequences accordingly. At finite coupling, say that the zero-coideal point **approaches** the classical uniform boundary. |
| L4:10, “generate the entire tower”; :310 “every excited boundary state $\lvert n\rangle$ ... same ... $\bar\eta$” | `prp-bsg-excited` **[Numerical]**, first-level states $\lvert n\rangle$; full bootstrap has multi-index states | “The invariant lines explain the matrix structure of the first-level states $\lvert n\rangle$. Further levels require the boundary bootstrap and scalar factors.” |
| L5:6,58,97–108, existence dichotomy for all $a_n$ | `prp-an-k-existence` **[Numerical for $n\le7$; Proposition, sketch in general]** | “For $a_2$ the mixed-branch argument below completes the restriction. The book gives checks through $n=7$ and a general sketch.” |
| L5:7,114–131, three matrices “valid for any $N$” with all associated scattering properties | `prp-an-nondiagonal` **[Numerical]**, $n=2,3,4$ | “The displayed candidate Toeplitz formula extends to general $N$. The $a_2$ calculation verifies it directly; the book records the ranges of its other checks.” |
| L5:10,232–251, “resolve”, “Why ... wrong?”, both bold boundary identifications, missed intercept, “confirming both”, “triumphant confirmation”, “exact particle reflection factor” | `cnj-neumann-coideal`, `prp-dm-first-order`, `cnj-neumann-amplitude` **[Conjecture]** | “Under the proposed dictionary, $\hat\epsilon=0$ corresponds to the uniform boundary and $s=+1$ to Neumann. The nonanalytic normalization explains how the proposed intercept could evade the perturbative argument; it does not prove the intercept. The SG expansion follows conditionally from UV–IR. Fusion with the stated CDD choice is consistent with this interpretation.” |
| L5:212 caption, “exact ... curve”, uniform point $\mu\to1$ | `cnj-bsg-uv-ir` **[Conjecture]** | “The curve implied by the UV–IR conjecture has its zero at $\mu_0=\cos[\pi/(2(\lambda+1))]$, tending to 1 classically.” |
| L6:6–9, scalar factors “determine”, continuation “produces the exact ... amplitudes”, swaps boundaries, quantizes modulus into rotor | `prp-an-diagonal-scalar`, `prp-an-nondiagonal-scalar`, `prp-uniform-amplitude`, `prp-solitonic-amplitude` **[Numerical]**; Neumann **[Conjecture]**; rotor is interpretation | “We study scalar solutions and their CDD ambiguity, continue them to proposed particle amplitudes, and compare their classical limits. The duality and rotor interpretation have the qualifications explained below.” |
| L6:22,34, “discovered the true path”, “we discover ... swaps” | same, plus Appendix E's explicit CDD condition | “These works motivate soliton fusion and the conjectured boundary duality. The chosen amplitudes obey exact algebraic duality; identifying the corresponding boundary theories needs additional input.” |
| L6:53,81,101, Gamma products “yields”, “every single ... telescopes”, bootstrap “produces ... every particle” | `prp-an-diagonal-scalar` **[Numerical]** ($n=2,3,4$); `prp-breather-amplitude` **[Established]**; `prp-uniform-amplitude` **[Numerical]** | “Use the stated scalar solution. Its breather projection gives the three-block expression; a complete all-rank telescoping proof is not supplied here. The particle product is the corresponding bootstrap candidate.” |
| L6:111, “landed squarely on the uniform boundary” | `prp-uniform-amplitude` **[Numerical]** and continuation assumption | “Its classical limit agrees with the uniform-boundary linear reflection factor.” |
| L6:119–134, general $c_+$ and telescoping, “full multiplet ... Neumann” | `prp-an-nondiagonal-scalar` **[Numerical]**, $n=2,3,4$; `cnj-neumann-amplitude` **[Conjecture]** | “With the chosen scalar factor, including $P_{\omega/2}$, fusion yields the following proposed Neumann amplitudes. The script checks the matrix projection for $a_2$.” |
| L6:151, one loop taken as identification | `cnj-neumann-amplitude` **[Conjecture]** | “The classical limit and checked one-loop coefficients support, but do not uniquely determine, the Neumann assignment.” |
| L6:165–181, “classical ground state”, dressing “simply” gives amplitudes and entire lattice | `prp-solitonic-vacua`, `prp-solitonic-amplitude`, `prp-a2-excited` **[Numerical]** | “For the solitonic backgrounds considered here, the stationary-soliton construction motivates these reflection factors and their bootstrap spectrum. State labels obey the physical ranges below.” |
| L6:206–210, rotor “must”, full circle, exact prefactor establishes quantization | book `sec-boundary-spectrum-a2` explicitly calls this an interpretation | Use P1 §2.2 replacement. The independently verified inertia does not establish the global Hilbert space. |
| L6:229 caption, physical boundaries swap and “quantized rotor tower” | same conjectures / interpretation | “Algebraic duality of the proposed amplitudes; boundary-state energies with a rotor-like weak-coupling limit.” |
| L6:254, “true Neumann amplitude” and “failed at two loops” | `cnj-neumann-amplitude` **[Conjecture]**; no two-loop computation in the lecture | “Compare the proposed Neumann factor with a self-dual candidate that agrees at one loop but differs at $O(\beta^4)$.” |

L2's explicit $a_1,a_2$ matrices and L3's local Robin calculation do not require a conjectural status merely because a larger classification in the book is numerical. Likewise L4's coproduct is an exact algebra identity. Conversely a link whose text says “Conjecture” does not cure an unqualified assertion elsewhere.

### 1.2 Where the CDD choice enters

Write $A^{(+)}=A^{(+)}_{\rm base}P_{\omega/2}$. Multiplying by a soliton CDD factor $P_x=(T/2-x)_I(T/2+x)_I$ leaves the two scalar functional equations unchanged. Fusion multiplies the breather factor by $P_x(\mu_a)P_x(\mu_b)$, not by a constant. After continuation let $\chi=2x/\omega$ and $v=\pi(1/2+B/(2h))$. Direct multiplication of the shifted hyperbolic blocks gives
\[
 D(\chi)=(h+\chi+B/2)(\chi-B/2)(h-\chi+B/2)(-\chi-B/2).
\]
Thus the specific $P_{\omega/2}$, with $\chi=1$, is already built into L6:124 before its three-block formula is displayed. Continuation then gives
\[
 (-n)(n+B/2)(-B/2)=\frac{(h-1+B/2)}{(h-1)(B/2)}=K_1^N.
\]
Bootstrap yields the displayed $K_a^N$, and substitution $B\mapsto2-B$ proves its equality to $K_a(2-B)$. This proves an identity **after a scalar choice**, not uniqueness of that choice from the coideal or fusion.

Independently, $D(\chi)|_{B=0}=1$ for every $\chi$, so replacing $\chi=1$ by $1-B/2$ changes $D$ only at $O(B^2)$. In `phase2_checks.py`, at $h=3,\theta=0.71$, the absolute change divided by $B^2$ is $0.7494831,0.7473872,0.7471793$ for $B=10^{-2},10^{-3},10^{-4}$. Hence neither the classical limit nor one loop fixes this freedom. This agrees with Appendix E, but the computation does not rely on its claim.

**L6:124 insertion:** “The functional equations leave CDD freedom. We choose $P_{\omega/2}$ because it gives exact duality with the uniform candidate after continuation. Other admissible choices can share the same classical and one-loop limits. The Neumann identification remains conjectural.”

## 2. Boundary Lee–Yang calculation

**Overlap:** P1 risks 3, 4, 22 and 25, all checked. P1 §§2.3–4 supply the complete replacement bootstrap proof, the correct doubling identity and the independently verified excited shift; §2.22 supplies the exact source line for GZ (3.51) and bulk residue. Those checks need not be repeated. **Verdict: mixed; additional one-point-residue sign error.**

- $S_{\rm LY}(2\theta)=-(1)(2)/[(4)(5)]$ in the lecture's $h=6$ blocks. The abandoned steps do not prove the conclusion. Use the ready-to-paste proof in P1 §2.4. It establishes $\epsilon=\epsilon^2=+1$.
- GZ (3.51) really prints the minus sign. Appendix D's GZ entry says the same; it is not a misquotation. P1 quotes the TeX source and checks the bootstrap independently.
- The residue at $i\pi/2$ and $g=2i\sqrt{2\sqrt3-3}$ are correct in the GZ convention, up to the simultaneous sign freedom of couplings. GZ's source explicitly prints `$g_{(1)} = 2i\sqrt{2\sqrt{3}-3}$` in (3.54). Their (3.45) is a formula for the **pair amplitude**, `$K^{c_1 c_2}(\theta) \simeq -{i\over 2}\,{{g^{c_1} g^{c_2}}\over \theta}$`; changing variable to the reflection rapidity gives $+ig^2/2$ at $i\pi/2$.
- **New:** exactly the same Jacobian is missing in L1:291's `eq-lec1-one-point-residue`. The book writes $k(t)\sim i\Gamma g/[2(t-iu/2)]$. Since $t=i\pi/2-\theta$, the reflection pole has **$-i\Gamma g/2$**, not $+i\Gamma g/2$. With $\Gamma=i\sqrt{2\sqrt3}$ and the above $g$, independent block evaluation gives $\operatorname{Res}_{i\pi/6}K=i(3-\sqrt3)=1.26794919243112i$, exactly $-i\Gamma g/2$. The positive sign fails. `phase2_checks.py` checks the limit at displacement $10^{-25}$.

**Replace L1:291 by**
\[
K_a^b(\theta)\sim-\frac i2\frac{\Gamma_{a\bar b}^{c}g_c}{\theta-i(\pi/2-u_{a\bar b}^{c}/2)}.
\]
**Replace L1:452's coupling value by** `$\Gamma^2=-2\sqrt3$ in the convention $\operatorname{Res}_{2\pi i/3}S=i\Gamma^2$`. The printed extra $(2+\sqrt3)$ is not this normalization, and no different normalization is defined in the lecture. These are separate issues: changing the sign of $\Gamma$ cannot repair a wrong square.

The CFT data are **correct**. Independently the minimal-model formulas give $c=1-6(5-2)^2/(2\cdot5)=-22/5$ and $h_{1,2}=((5-4)^2-9)/40=-1/5$. A boundary primary has dimension $h$, not $h+\bar h$. The two Cardy labels are $\mathbb1,\phi$; the identity boundary's self-channel contains only the identity, whereas the $\phi$ boundary admits the perturbing boundary field of dimension $-1/5$. DPTW's source `Paper4.tex`, discussion of the conformal boundaries preceding `eq:rfb`, makes this distinction. Use “Cardy” with a capital C. Equality of the massive factors at $b=0$ does not identify the two boundary conditions (P1 §2.25).

$K_b^{(1)}=K_{b+2}=K_{4-b}$ is **correct**; P1 §2.3 gives the block cancellation and source comparison. For Exercise 1.3(4), the block identity is **wrong**, but the stated ratio is **correct**:
\[
 f_x=-(6x)(6-6x),\qquad
 \frac{f_{(1+b)/6}}{f_{(1-b)/6}}=
 \frac{(1+b)(5-b)}{(1-b)(5+b)}=K_b/K_{(\mathbb1)}.
\]
This follows by factoring $\sinh\theta\pm i\sin\pi x$ and needs no numerical assumption.

## 3. Charges, reflection equations and Euclidean crossing

**Overlap:** P1 risks 11, 20, 21 and 25, checked. The density signs, braid convention, twisted domain failure and omission of zero-momentum particles are already established in P1 §§2.11,2.20–21. **Verdict: wrong at several additional normalization/index/channel steps.**

### Charges and normalization

The correct laws are $\partial_t(T-\Theta)=\partial_x(T+\Theta)$ and $\partial_t(\bar T-\bar\Theta)=-\partial_x(\bar T+\bar\Theta)$. Consequently the bulk density is $T-\Theta+\bar T-\bar\Theta$ and the surface flux is $T+\Theta-\bar T-\bar\Theta$. Set that flux to $-\dot\Sigma_s$. P1 §2.20 has the replacement displays.

In the book $Q_{\pm1}=E\pm P$, so $Q_1+Q_{-1}=2E_{\rm bulk}$ and $\Sigma_1=2\mathcal B$. Independently the canonical energy flux is $\dot E_{\rm bulk}=\dot\phi\phi'|_0=-\dot{\mathcal B}$, proving the factor two. **Replace L1:57:** “With the book's normalization, $\widehat Q_1=2(E_{\rm bulk}+\mathcal B)$. Dividing by two gives the total energy.” One may instead normalize both chiral charges by $1/2$, but must change their eigenvalues consistently.

The boundary contribution on asymptotic states is a constant $H^{(s)}$, not literally an operator that does nothing. It cancels for elastic scattering on the same boundary state; it need not cancel when the boundary changes state. Add that restriction at L1:63.

### Braid and twisted equations

The display has $1\otimes K$, the **second** slot. P1 §2.11 derives the needed $\check S=S_{12}P=PS_{21}$ convention for that displayed fixed-slot equation. The printed $PS_{12}$ conversion is not general. The L1 twisted RHS needs $\check R_{\nu,\bar\mu}$, as in the book and L4:146. Its intermediate spaces are
\[
 V_\mu V_\nu\to V_\nu V_\mu\to V_\nu V_{\bar\mu}
 \to V_{\bar\mu}V_\nu\to V_{\bar\mu}V_{\bar\nu}.
\]
L4's typed equation is correct; its scalar-normalization proof needs the repair in §5 below.

### Rotation, pair measure and a new crossing index error

The matrices of $(-t_E,-x)$, $(-t_E,x)$ and $(t_E,-x)$ have determinants $-1,+1,+1$. Only the last two are rotations. Index's rotation maps $x\le0$ to **negative** new Euclidean time; the book's clockwise rotation maps it to positive new time. Use $(\tilde x,\tilde t_E)=(t_E,-x)$ in Index:26, L1:188 and the caption. P1 §2.20 independently checks the inconsistent rapidity implication; use its replacement and derive the crossed-particle continuation separately.

The measure $\frac12\,d\theta/(2\pi)$ is correct for states normalized to $2\pi\delta(\theta-\theta')$. Rescaling each creation operator by $\sqrt{2\pi}$ relative to the book's $\delta$ normalization requires exactly this measure. Add that normalization and the possible one-particle term. The pair coefficient $\mathcal K^{ab}=K_{\bar a}^{b}(i\pi/2-\theta)$ agrees with the book.

**New wrong index order:** the lecture's own ZF rule gives
\[
 A_c(+t)A_d(-t)=S_{cd}^{ef}(2t)A_f(-t)A_e(+t).
\]
Therefore the coefficient of $A_a(-t)A_b(+t)$ is $S_{cd}^{ba}\mathcal K^{cd}(-t)$, not $S_{dc}^{ba}\mathcal K^{cd}(-t)$. **Replace the lower pair $dc$ by $cd$ in L1:232,238.** This agrees with the book's vector equation $k=\check S k(-t)$. In `phase2_additional.py` an independently solved normalized six-vertex matrix gives proportionality residual $1.8642\times10^{-16}$ for the corrected equation and $0.8520064$ for the printed one. A scalar factor cannot repair that failure.

### Sources of scalar crossing forms and bootstrap channels

The three diagonal identities are algebraically equivalent with bulk unitarity, crossing, C-invariance and periodicity as used in the book. However CDRS prints the **shifted** equation, source `9404108/source.tex`, `\eqn\crossing{K_1^0(\theta )K_2^0(\theta -i\pi )= S_{11}(2\theta ).}`. The centered product form is used by Corrigan–Delius (1999), (2.7), and the shifted form by Delius–Gandenberger, (3.17). **Replace L1:255's source/name by Corrigan–Delius `@corrigan_delius1999`**; the exercise's name is the appropriate one. CDRS supports an equivalent relation, not the particular source attribution as printed.

**New wrong channels:** for diagonal particle reflection the species is unchanged at the wall. The internal collision therefore uses $S_{ab}$ in the bulk-fusion bootstrap and $S_{ba}$ in **both** factors of the excited-boundary bootstrap. **Replace $S_{b\bar a}$ at L1:282 and :299 by $S_{ba}$.** L6 and the book have the right channels. At $h=3,B=\sqrt2/4,\theta=0.71$, the correct $S_{11}$ bootstrap error is $3.93\times10^{-51}$; replacing it by $S_{12}$ gives $0.147034953$. This is a direct counterexample from the stated blocks (`phase2_checks.py`). If a boundary conjugates species, one must retain the full typed matrix bootstrap instead of inserting a bar into the diagonal formula.

### Elasticity exercise and scope

For the exercise's stipulated equal charge coefficients, let $r$ particles attain $M=\theta'_{\max}$. Dividing by $e^{sM}/2$ makes the RHS tend to $r$, independently of whether $r=1$. Equality forces $\theta_1=M$, after which the LHS tends to 1, hence $r=1$. Subtraction leaves a sum of strictly positive coshes, so there are no other particles. **Replace Exercise 1.1(2)'s “finite ... iff ... only one” by these two successive steps.** Equal mass alone does not imply equal higher-spin coefficients for distinct species. Write “one species, or species with identical normalized higher-charge eigenvalues”. The surrounding L1:77 assertion for arbitrary $\gamma_b^{(s)}$ is not proved by linear independence of functions: its coefficients depend on $s$. Replace it by “Under the charge-separation hypotheses of the elasticity argument, reflection preserves the rapidity and mass multiplet; Exercise 1.1 proves the equal-coefficient case.”

## 4. Mirror pairs, parity and vacua

**Overlap:** P1 risks 1, 10, 15, 23, 2 and 25, all checked. P1 §§2.1,2.10,2.15,2.23 fully cover the wrong $k$ identities, tau derivative, vacuum sign assignments, root-positivity proof, Robin block expression, and the two vacuum/breather exercises. **Verdict: mixed; parity shift and time-delay signs survive independent checks.**

Use $k=2m\sin^2(\pi a/h)/s$, the derivative and algebraic square-root branch in P1 §2.10. The factor in the asymptotic tau polynomial is $1+\zeta^{-ja}AE_2$, not $1+\zeta^{-2ja}AE_2$.

The trajectory calculation is independently short. Put $k=m_a\cosh\theta$, $v=\tanh\theta$. Incoming and outgoing intercepts give $\Delta t=-(\xi_1+\xi_2)/(kv)$; the incoming mirror has intercept $-(\xi_2+\ln A)/k$. Their intersection is $d=-(\xi_1+\xi_2+\ln A)/(2k)$. Substituting the three shifts gives precisely L3's printed delays and $d$: $C=+1$ delays and has $d>0$; $C=-1$ advances and has $d<0$; Neumann advances with $d=0$. The book agrees. Thus “mirror behind the wall” must be restricted to the appropriate boundary/asymptotic trajectory, not used as a universal geometric claim in Index:34 or L3:4.

**Parity and logarithm branches are correct for $a_n$.** In units $2\pi/\beta$, the constant shift is
\[
 v=-\frac ah\sum_{j=0}^{h-1}j\alpha_j-\sum_j\ell_j\alpha_j,\qquad\ell_j\in\mathbb Z.
\]
The affine $a_n$ Cartan matrix applied to the first coefficient vector has integer entries: interior second differences vanish and the two wraparound differences are integer multiples of $a$. Every branch change adds a root; in this simply-laced system its pairing with each simple root is integral. Hence $v\in\Lambda_W^\vee$ for every branch. `phase2_additional.py` checks all species for $h=3,\ldots,9$ using explicit roots. The displayed transformation $(a,\theta,\xi)\mapsto(h-a,-\theta,-\xi)$ is correct. Parity of a solution alone is not a proof that every possible wall returns this species; that is the scope distinction in §1.

The $\epsilon_j$ used in $C_j=\epsilon_{j+1}\epsilon_{j-1}$ are extended integer-index signs, with $\epsilon_{j+h}=(-1)^a\epsilon_j$. Their product formula is nevertheless periodic and has product $+1$. Both wraparound signs must be retained. This agrees with $C_j=(-1)^{\alpha_j\cdot\lambda}$ after choosing the continued field's vacuum branch; it does not justify reducing each $\epsilon$ index modulo $h$. P1 §2.1 recomputes the wall condition directly, rather than choosing the best sign as the author script does.

The linearized condition is $\eta_a'(0)=-C m_a^2\eta_a(0)/(2m)$ and
\[
 K_a^{\rm cl}=\frac{ip+C m_a^2/(2m)}{ip-C m_a^2/(2m)}
 =[-(a)(h-a)]^{-C}.
\]
This gives 1 at $C=0$ and $-(a)(h-a)$ at $C=-1$. The printed placement of the minus sign fails at $C=0$.

For the self-conjugate family, differentiating its energy with respect to $\xi$ gives only the boundary variation $[\phi'+\nabla\mathcal B]\partial_\xi\phi=0$. Its limit at $\xi\to-\infty$ is the zero field of energy zero, so the entire regular $\xi<0$ family has $E=0$. Its derivative is the zero mode; the pole at $i\pi/2$ is at the strip **edge**, not inside. The open-half-line positivity in the breather exercise is true, but the $\rho_-=0$ solution has wall zeros at recurring times. Its stronger regular-boundary interpretation and the claimed $C=-1$ endpoint sign change are wrong; the factorization and replacement are in P1 §2.23.

## 5. Coideal charges and boundary sine-Gordon

**Overlap:** P1 risks 12,16,17,21,25 were checked; the CPT normalization, theorem hypotheses and Neumann exercise were incompletely covered. **Verdict: exact coproduct and selected matrix identities correct; CPT presented too strongly, crossing display wrong, physical-spectrum scope incomplete.**

### CPT: what is and is not established

The two image formulas are compatible if $\bar\Phi$ means the **parity transform of the operator** $\Phi$, not complex conjugation. For a scalar, the full-field image is $\phi(x)+\phi(-x)$; for a chiral operator its image has the opposite chirality. Differentiating this full-field sum gives zero normal derivative at the wall. L4 does not define the uppercase bar, so **insert after L4:39:** “Here the bar on $\Phi$ denotes parity transformation; it is not complex conjugation. For the scalar field this reduces to $\tilde\phi(x)=\phi(x)+\phi(-x)$.”

The contact identity is **correct**. From $(u-i0)^{-1}-(u+i0)^{-1}=2\pi i\delta(u)$, differentiating $n-1$ times with respect to $t'$ gives the printed formula, including its sign. Upon integration, $\int dt'\,\partial_{t'}\delta(t-t')f(t')=-f'(t)$. The independent regulated $n=2$ integral against $e^{-t'^2}(1+0.3t')$ tends to $-0.6\pi i$ (values $-1.84297656i,-1.88070735i$ at regulators $10^{-2},10^{-3}$; target $-1.88495559i$).

“Spin-1 non-local currents” at L4:58 is **wrong**. In SG, from the canonically normalized chiral correlator $\langle\phi_R(z)\phi_R(0)\rangle=-\log z/(4\pi)$, the vertex with exponent $a=8\pi/\beta_{\rm SG}$ has weight $a^2/(8\pi)=8\pi/\beta_{\rm SG}^2=\lambda+1$; its integrated charge has spin $\lambda$. The double pole follows from the product of the two vertex charges and the image contraction, not from integer spin one. **Replace by:** “The appropriately normalized non-local vertex currents have a double-pole OPE with the same-node boundary perturbation. Their weights depend on the coupling; in SG the current weight is $\lambda+1$.” The sentence at L4:24 about a left-mover hitting the right wall reverses the direction: replace it by “A right-moving excitation hits the wall and returns left-moving.”

The displayed $T_1$ specializes correctly to $\beta_{\rm SG}[\phi(0)-\phi(-\infty)]/(2\pi)$, the integer soliton charge in SG. The simply-laced coweight normalization is likewise consistent. **However the claimed literal matching of the contact vertex with $q^{T_j}$ in all conventions is unclear.** The lecture omits the free-field normalization, the image/physical-field identification at the wall and the branch of $\log q$. These matter because $T_j$ at a boundary is not integer-valued. In Delius–MacKay's source, (3.25)–(3.28) use $e^{i(1/\hat\beta-\hat\beta)\phi}$ and $q=\exp[2\pi i(1-\hat\beta^2)/\hat\beta^2]$; matching these is direct in their normalization. Replacing the complex number $q$ by an equivalent exponential with a $2\pi i$-shift changes its noninteger power. I have not completed a convention-faithful OPE translation sufficient to certify the lecture's coefficient. No correction to the published source is inferred.

There is also a definite perturbative-order omission. The source says exactly “is conserved to first order in perturbation theory” after its displayed $\hat Q_+$. The book explicitly keeps $O(\lambda_b^2)$ before the charge display; the lecture drops it. **Replace L4:68–75 by:** “To first order in the boundary coupling, the contact term can be written in the conventional charge normalization as $\partial_t(Q_j+\bar Q_j)=-\hat\epsilon_j\partial_tq^{T_j}+O(\lambda_b^2)$. This motivates the conserved-charge ansatz below; the first-order computation alone does not prove its exact conservation.” Restrict the field-theory discussion to the simply-laced case unless the dual-root factors are specified.

The algebraic dictionary and coproduct are **correct**: $\Delta(fk)=fk\otimes1+k\otimes fk$, and $k^{1/2}fk^{1/2}=q^{-1}fk$. Thus the resulting $b_j$ generate a left coideal in the stated coproduct convention, independently of the unresolved field-theory normalization.

### Crossing and sign of q

P1 §2.16 constructs the six-vertex braid independently: $R(-q)=(\Pi\otimes1)R(q)(\Pi\otimes1)$ for the normalized matrix, and $K(\hat\epsilon;-q)=K(-\hat\epsilon;q)$ by direct substitution. These are valid **matrix identities at fixed spectral convention**, not an ordinary gauge identification of the physical theories. The reflection equation admits both sign conventions after corresponding matrix transformations. The generic matrix crossing test distinguishes them when the ordinary $C$ is fixed; replacing $C$ by the antisymmetric version restores matrix proportionality. Its exceptional-parameter cases and full scalar crossing are not established by the script, so “if and only if” for all data is too strong. Use P1 §2.16's qualification.

**New:** L4:261 is not the reflection-rapidity version of L1's centered crossing formula. Substituting $t=i\pi/2-\theta$ into the corrected equation of §3 gives
\[
 K_{\bar a}^{b}(\theta)=\sum_{c,d}S_{cd}^{ba}(i\pi-2\theta)K_{\bar c}^{d}(i\pi-\theta).
\]
**Replace L4:261 with this display**, or simply use the centered vector equation $k(t)=\check S(2t)k(-t)$ with $k$ explicitly defined. The printed $S(2\theta)$ has the wrong argument for reflection rapidity; its $K_{ab}$ indices are also undefined.

### Invariant line, excited parameters and reflection-equation theorem

P1 §2.12 fully solves the simultaneous eigenvector problem and supplies replacement text for L4:299–304: $v_n=(1,(-1)^ne^{-\vartheta})^T$ with the printed $\beta_j$ and $\bar\eta$. The determinant must be $\det(b_j-\beta_j1)$; there are other algebraic reducibility branches. Substituting $\vartheta\mapsto-\vartheta$ in $P_0^\pm$ exchanges them and leaves $Q_0$ unchanged. This is **correct at the shifted value $\bar\eta$**; it does not say that changing $\eta$ alone exchanges the original diagonal entries. First-level states $|n\rangle$ have this matrix structure; a common eigenvector by itself does not select a physical pole or exhaust all bootstrap levels.

The coideal proves that the two reflection histories are intertwiners. To conclude proportionality requires generic nonzero/invertible intertwiners and a scalar endomorphism algebra (or generic irreducibility). It does **not** fix the scalar function of the rapidities. A simple independent test of the logical gap is to rescale one of the mixed bulk intertwiners by a nonconstant scalar: all linear intertwining equations survive, while one side of the claimed reflection equation is rescaled. Agreement at a single numerical rapidity cannot fix that function either.

**Replace L4:149's last sentence:** “Under generic invertibility and the scalar-endomorphism hypothesis, the histories are proportional at each pair of rapidities. Equality additionally requires compatible normalization of the four bulk intertwiners; this must be verified for the mixed representations. A highest-weight normalization does not automatically supply that comparison when the modules change.” The book's `thm-k-reflection` includes a normalization hypothesis, but its phrase “agree at one point” should explicitly mean a specified nonzero matrix element **for every generic rapidity pair**, not one parameter value. Otherwise that theorem's proof has the same logical gap.

### Neumann poles and corrected node labels

With $\eta=\pi(\lambda+1)/2$, direct substitution gives $\nu_n=\pi/2-n\pi/\lambda$ and $E_n-E_0=M\sin(n\pi/\lambda)=m_n\cos(n\pi/(2\lambda))$. The formal mass expression $\frac12m_{2n}$ agrees. A bulk breather $B_{2n}$ exists only for $2n<\lambda$; a positive physical boundary pole requires $0<n<\lambda/2$. **Add those bounds in Exercise 4.2.** At $n=0$, $\nu_0=\pi/2$ is a strip-edge, zero-energy pole, not a new positive-energy excited state and not a bulk $B_0$. A mirror-image interpretation is a classical heuristic; it is not a proof of the exact quantum identity or of every state's existence. For $\epsilon_*$ choose the analytic square-root branch $1/[2\cos(\pi\lambda/2)]$, rather than a positive-real square root of its square.

L4:281 agrees with the post-ed0585a book: $\hat\epsilon_1$ depends on $\mu_-\propto C_1$, $\hat\epsilon_0$ on $\mu_+\propto C_0$, and $\eta=-4\pi\phi_0/\beta_{\rm SG}$. The node assignment also follows directly by matching the two exponential coefficients in the boundary cosine. No new sign correction is warranted here.

## 6. Higher-rank intertwiners, Serre relations and q-Onsager scope

**Overlap:** P1 risks 8,9,17,18,19,25; the first two were fully checked for their advertised targets. Mixed quantum branches, the Toeplitz exercise and source scope are new checks. **Verdict: formulas largely correct, proof and attribution qualifications needed.**

### Mixed branches and Toeplitz matrix

The two cases as printed do **not** exhaust patterns. For $a_2$ the omitted argument is elementary. If node 1 has nonzero spectator determinant, $K_{13}=K_{23}=K_{31}=K_{32}=0$. In the node-2 row/column systems these zeros force $K_{12}=K_{21}=0$ regardless of node 2's determinant. Thus $K$ is diagonal; the active off-diagonal equations make all diagonals equal, and the active diagonal equations force every $\hat\epsilon_j=0$ for generic $q$. Relabeling covers either other initial node. **Insert this paragraph after L5:108.** This proves the $a_2$ dichotomy without assuming all determinants simultaneously nonzero. A general-rank propagation proof is not supplied here.

Independent numerical checks enumerate all 64 $a_2$ patterns with node values $0,\epsilon_*,-\epsilon_*,0.37\epsilon_*$ at $\omega=\sqrt2$. Exactly the all-zero and eight all-special patterns have nullspaces; largest allowed smallest-singular-value $8.56\times10^{-16}$, smallest disallowed $0.19525$. This is a test of the proof, not a numerical proof at arbitrary parameters.

Under Cartan sign changes the rescalings on a cyclic root system have product 1. They preserve $\prod_j\hat\epsilon_j/\epsilon_*$, and any even set of node signs is obtainable by assigning successive diagonal signs around the cycle. Hence the claimed **single discrete invariant** is correct within the fixed parity-symmetric normalization. It is not a claim about arbitrary continuous rescaling of the $c_j$.

The general-$N$ formula and explicit $3\times3$ matrix match the book. Direct multiplication in `phase2_checks.py` for $N=3,\ldots,8$, both signs, gives maximum intertwining-entry residual $1.47\times10^{-15}$; P1 §2.9 checks unitarity and its $U_s$. The transpose really is the reverse-multiplet matrix in these conventions; it is not complex conjugation. The cyclic eigenvalue formula is checked in P1 §2.9; the error is $Z^T=sZ^{-1}$, which must be $Z^T=Z^{-1}$.

**New exercise defect:** node-0 **spectator** equations involve $K_{32},K_{12},K_{21},K_{23}$, not $K_{31}$. They give $K_{32}=sXK_{12}$ and $K_{21}=sXK_{23}$. Then a node-1 row equation gives $K_{31}=X^{-1}K_{32}$. **Replace L5:263's instruction by:** “Use node 0's spectator equations to find $K_{21}$ and $K_{32}$, then node 1's spectator-row equation to find $K_{31}$.” Parts 2–3 then give the stated $\kappa_s$.

### Serre computation and obstruction

P1 §2.8 independently checks the quadratic polynomial, including the Cartan coefficients and tensor-product tests. Exercise 5.2 is valid: $\varphi_i e_j=q^{-1}e_j\varphi_i$, the $e_i^2$ terms vanish by the bulk Serre relation, the $\varphi_i^2,k_i^2,k_i$ cross contributions cancel after commuting them, and the remainder is
\[
-q^2D_ie_j+e_jD_i
=\frac{c_i}{q-q^{-1}}e_j[-(k_i^2-q^2)+(k_i^2-1)]
=q c_i e_j.
\]
All nine coefficients of the Cartan calculation are consistent. The missing scalar in $\kappa=(F-b_j)/d$ and unauthorized universal inversion remain the errors identified in P1. In a finite-dimensional represented algebra, Cayley–Hamilton supplies the inverse of an invertible Cartan matrix; conjugation then separates the three components because the Vandermonde eigenvalues $1,q^3,q^{-3}$ are distinct. This requires generic $q$ (in particular $q^3\ne1$, $q^6\ne1$), with further genericity for propagation. It does not put inverses into the abstract coideal without proof.

### Rank one and mathematical sources

The quartic relation is **correct** for arbitrary boundary parameters. Independent 8-dimensional tensor-product checks at generic $c_0,c_1,\hat\epsilon_0,\hat\epsilon_1$ give errors $1.57\times10^{-15}$ and $1.12\times10^{-15}$ in
\[
 b_i^3b_j-[3]_qb_i^2b_jb_i+[3]_qb_ib_jb_i^2-b_jb_i^3
 =q c_i(q+q^{-1})^2[b_i,b_j].
\]
Word-reversal antisymmetry **alone** is not a proof of the cancellation: a weight-zero commutator such as $[e_i,f_i]k_i$ need not vanish. **Replace L5:207's explanatory clause:** “Expanding the degree-four polynomial and using the quantum-group commutation relations cancels the Cartan obstruction and gives the q-Dolan–Grady relation. This holds for arbitrary $(\hat\epsilon_0,\hat\epsilon_1)$.” Kolb Example 7.6 supplies the precise relation; the lecture should display it rather than relying on reversal alone.

From [Kolb](https://arxiv.org/abs/1207.6036), equations (5.9)–(5.11), the parameter sets are, transcribed with modern typography,
\[
\mathcal C=\{c\in(K(q)^\times)^{I\setminus X}:c_i=c_{\tau(i)}\text{ if }\tau(i)\ne i\text{ and }(\alpha_i,\Theta\alpha_i)=0\},
\]
\[
 I_{\rm ns}=\{i\in I\setminus X:\tau(i)=i,\ \alpha_i(h_j)=0\ (j\in X)\},\quad
\mathcal S=\{s:s_i\ne0\Rightarrow i\in I_{\rm ns},\ a_{ij}\in-2\mathbb N_0\ (j\in I_{\rm ns}\setminus\{i\})\}.
\]
The exact source definition calls $B_{c,s}$ a quantum symmetric pair coideal for “$c\in C$ and $s\in S$”. For the split diagram $X=\varnothing,\tau=\mathrm{id}$, $I_{\rm ns}=I$, and the lecture's odd-neighbor exclusion follows. It is **overbroad as a statement of Kolb's general theory** without these conditions. Add “For the split diagram used here” to L5:208.

[Balagović–Kolb](https://arxiv.org/abs/1507.06276) says in its abstract “a universal K-matrix if g is of finite type”; it does not by itself cover this affine case. [Appel–Vlaar](https://arxiv.org/abs/2007.09218) extends the construction to symmetrizable Kac–Moody type and a twisted/generalized reflection equation, using completions and additional diagram data. It does not by citation alone identify every physical scalar factor or ordinary numerical evaluation matrix here. **Replace the parenthesis at L5:208:** “Universal K-matrices were constructed in finite type by Balagović–Kolb and extended to the Kac–Moody setting, with a generalized reflection equation, by Appel–Vlaar.”

The character twist is **correct as an image construction**. Evaluate the scalar Serre relation: $(2-q-q^{-1})\beta_i^2\beta_j=\beta_j$. Connectedness gives $\beta=0$ or every $\beta_i^2=\epsilon_*^2$. Then
\[
 ({\rm id}\otimes\chi)\Delta(e_i+q^{-1}f_ik_i)
 =e_i+q^{-1}f_ik_i+\beta_i k_i.
\]
This is an algebra map with the desired image; no injectivity or automatic equality with a universal presentation follows. **Change “defining ... algebra ... obtained” in L5:209 to “realizing the generalized q-Onsager relations; this embedded coideal is the image of $B_0$ under the character twist”.**

[Baseilhac–Belliard](https://arxiv.org/abs/0906.1215), Proposition 2.1 and its simply-laced application, explicitly give `$\widehat{\epsilon}_j=\pm\frac{i}{q^{1/2}-q^{-1/2}}$`. Squaring it gives $1/[(1-q)(1-q^{-1})]$, independently matching $\epsilon_*^2$. Thus the value has prior publication and must not be implied to originate in these notes. The relation's cubic expansion is still useful teaching.

Delius–MacKay's printed `$c=\sqrt{\lambda\gamma^2(q_i^2-1)/2\pi i}$` directly establishes the $g^{1/2}$ scaling. Their first-order formula has no intercept and their Neumann discussion uses $Q_j+\bar Q_j$. Thus the description of what their calculation assumes is fair, but the proposed nonzero intercept remains the conjecture discussed in §1, not a disproved published theorem.

Finally the two signs called $s$ are **not the same dictionary**. Proposed assignments are quantum $\hat\epsilon=0\leftrightarrow$ classical product-$+1$ class (including uniform and its solitonic backgrounds), and quantum $s=+1\leftrightarrow$ Neumann. The book leaves quantum $s=-1$ and classical product-$-1$ without a proven pairing. Add that sentence at L5:228; the numerical similarity of “three classes” does not pair the third entries.

## 7. Scalars, particles, spectra and Coleman–Thun

**Overlap:** P1 risks 2,6,7,13,14,17,18,24,25, all checked with limitations. **Verdict: most formulas survive; interpretation, matrix-space bookkeeping and a source attribution need correction.**

### Scalar and fusion formulas

P1 §§2.7,2.18 compare and independently test the diagonal scalar product, singlet fusion and continuation. The displayed $A$, its two equations, $K_{B_1}$, $K_1$ and $K_a$ agree with the specified book objects. The apparent middle-block sign in continuation is resolved by period $2h$, not by changing the final factor. The diagonal infinite-product check had controlled cutoff convergence; it was not an all-rank proof.

**New space error at L6:48:** for soliton-conjugating $K:V\to\bar V$, the crossed boundary pair lies in $\bar V\otimes\bar V$. Replace “$\sum_i|\bar v_i\rangle\otimes|v_i\rangle$” by “$\sum_i|\bar v_i\rangle\otimes|\bar v_i\rangle$”. The printed mixed-multiplet vector cannot be acted on by an identical-soliton matrix. This is distinct from the incoming breather singlet, which does lie in $V\otimes\bar V$.

The subsequent fusion ordering is **correct**: $V\bar V\to VV\to VV\to V\bar V$, with the right constituent reflected first. The normalized braid acts as 1 on $v_i\otimes v_i$, and $K_{\bar V}=K_V^T$ is the correct reverse map. P1 §2.18 explicitly constructs the affine singlet at the stated rapidity separation and checks the $a_2$ matrix projection, including $c_+$. L6's $A^{(+)}$ equations and three-block $K_{B_1}^{(+)}$ agree with the book, but the lecture never provides the actual $A^{(+)}$ product or $\Psi_+$, so “telescopes” is a cited result, not a reproducible blackboard step. The full non-diagonal Gamma-product computation remains **unclear** in this pass; its CDD dependence is checked in §1.

The algebraic identity $K_a^{(-)}=-(a)(h-a)K_a$ is consistent with the continued three-block minus solution and subsequent bootstrap. Its classical limit is 1, yet this does not identify it as Neumann. The attribution of its **dual** to the CDRS $C=-1$ candidate is directly verified for $a_2$: CDRS writes `$-{(1)(2+B/2)\over(B/2)}$`, which is exactly $-(1)(2)K_1^N$. Gandenberger's source explicitly notes this relation in the footnote in §6. The general attribution is also correct: CDRS, equation `\general`, gives $\prod_{c=1}^a\langle c-1/2\rangle/\langle h-c+1/2\rangle$. Dividing this by the lecture’s $K_a^N$ telescopes to $\prod_{c=1}^a(c)(h-c)/[(c-1)(h-c+1)]=-(a)(h-a)$. Thus the dual of the minus candidate is the general CDRS candidate. This is an algebraic comparison of candidates, not an identification of the minus coideal with a classical boundary.

### Solitonic backgrounds and state labels

$M_a^{(b)}$ agrees with the book. At $a=b=1,h=3$, its four factors reduce using $(-x)=(x)^{-1}$ and period 6 to the four blocks displayed for $K_1^{b_{0,0}}$. This independent cancellation also verifies the ground-state specialization of the eight-block formula. The two positive block labels are in $(0,3/2)$ for $0<B<1$. Whether a given full amplitude retains a pole requires checking cancellations; the complete product was used in P1's residue tests.

The “stationary soliton in front of a uniform wall” is the **imaginary-coupling construction with its vacuum shift**, not a real stationary soliton simply added to the zero-field uniform vacuum. On continuation the branch of the field changes the real boundary signs, precisely as L3's wall calculation shows. **Replace L6:165–167's introductory wording:** “Start with the stationary imaginary-coupling soliton and its mirror at the uniform wall, including the vacuum shift. Continuing this solution gives a real stationary background for the associated solitonic sign pattern. Factorized passage through that soliton motivates the dressing factor below.” The resulting amplitude is still subject to the book's stated numerical/continuation evidence.

The eight-block spectrum, $F(k)$, diagonal energy and inertia are already checked independently in P1 §§2.2,2.24. The book **does contain** the inertia integral and its value at chapter 36:282. The lecture's assertion that the global coordinate must be a $2\pi$ circle is not established; use P1's replacement. Negative indices are bootstrap lattice coordinates, not negative counts of bound particles. **Replace L6:181:** “For nonnegative indices, $n,m$ count the two types of successive binding. Negative labels describe states reached along diagonal bootstrap paths; the energy formula extends to them, subject to the physical existence bounds.” Those bounds are given in P1 §2.24 and must be included.

### Coleman–Thun check and source scope

The window, label $p=47/40$, $\theta_R=7i\pi/120$, negative residue and residue identity were independently checked in P1 §2.6. $47/40<3/2$, so the pole is **inside** the strip in these labels; it is not outside merely because $47/40>1$. Directly,
\[
 F(-1)-F(-2)=\cos[\pi(1/2+3B/2)/3].
\]
Thus the transition from $b_{-2,2}$ to $b_{-1,2}$ has exactly the pole energy. The sentence naming $b_{-1,2}$ is **correct**, not a state-index typo. Its existence bound $2<1/(2B)+1/2$ gives $B<1/3$; its threshold relative to $b_{-1,1}$ has binding rapidity $\pi(1/2-3B/2)/3$, which reaches zero at $B=1/3$. `phase2_checks.py` checks the energy difference independently at $B=0.23,0.45$.

The arithmetic $5-2(2)-1=1$ is **wrong**. P1 §2.6 carefully leaves the separate $P,L$ counts unclear, rather than inventing a sixth line. Its replacement is still recommended. The source [Dorey–Tateo–Watts](https://arxiv.org/abs/hep-th/9810098), discussion of Figure 8, says “Thus the prediction is demoted to first order, as required.” It discusses a naive second-order pole multiplied by a vanishing reflection factor; I found **no literal formula $P-2L-Z$ in that paper**. Thus I cannot quote that formula as requested: the attribution of that exact formula is unsupported. Ordinary two-dimensional power counting gives one singular denominator per on-shell propagator and two integration directions per loop under the usual nondegeneracy assumptions; zeros then reduce the order. Application to a boundary graph needs its actual integrations and geometry. **Replace L6:214:** “Dorey, Tateo and Watts explain how a zero of an internal reflection factor can lower a boundary Coleman–Thun pole's order. For a specified graph, ordinary on-shell power counting must include both propagators and integration variables.” This is a source-scope correction, not a rejection of the verified residue identity.

### One loop and historical opening

The change $B\mapsto B(2-B)/2$ preserves the first derivative at zero because its derivative is 1 there. It does not establish that all historical self-dual proposals matched perturbation theory. Kim's source states of the earlier minimal $a_3$ factors: “None of these agrees with the perturbative result.” His own ansatz assumes self-duality. Therefore **L6:20's group attribution and blanket one-loop claim are wrong**. Replace by: “Early bootstrap constructions left a boundary-identification and CDD problem. Kim used one-loop Neumann calculations together with a self-duality assumption. Later non-self-dual candidates shared the one-loop coefficients, showing that this test does not determine the higher-order terms.”

Independent differentiation of the proposed $a_3$ amplitudes reproduces Kim's explicit light and heavy one-loop expressions in his Appendix, errors $4.38\times10^{-47},6.57\times10^{-47}$ at $\theta=0.71$ (`phase2_additional.py`). The $a_2$ coefficient was checked in P1 §2.13 and is printed explicitly by Gandenberger in §6. The actual Kim $a_2$ ansatz is $[1/2][3/2]$, whose duality follows by exchanging its denominator blocks. It is **not** the purported self-dual product in Exercise 6.2(3). P1 supplies a valid replacement exercise. Neither a different $O(\beta^4)$ coefficient nor a conjectured alternative is itself a completed two-loop field-theory calculation. The broad priority sentence about when all bulk particle amplitudes were first recovered by soliton fusion remains **unclear**; it needs an algebra-by-algebra historical scope, not a universal chronology.

## 8. Sklyanin construction and classical exercises

**Overlap:** P1 risks 5,19,25; all checked for $a_1,a_2$ matrices, Dirichlet obstruction and figure. **Verdict: main transport calculation correct; genericity, coweight exercise and one book explanation need repair.**

Differentiating the two transports gives exactly
\[
\dot\tau=-\operatorname{tr}[(Ka_++\sigma(a'_+)K)TK_-\sigma(T')]
+\operatorname{tr}[KT(a_-K_-+K_-\sigma(a'_-))\sigma(T')].
\]
This follows by the product rule and reversing products under $\sigma$, so the stated compatibility equation has the correct plus sign. For transpose and the zero vacuum, $a_\infty(\lambda)+a_\infty(1/\lambda)^T=0$, and $K_-=1$ is appropriate. For an arbitrary involution or other vacuum that cancellation must be checked, not assumed. The Poisson-commuting conclusion also assumes the bulk ultralocal $r$-matrix bracket and the compatible equation at the left end. **Add to L2:126–130:** “For the transpose involution and vacuum asymptotics used below the left condition is satisfied by $K_-=1$. Involution of the charges additionally uses the bulk Poisson bracket and both boundary reflection equations.” The displayed finite-interval differentiation itself is sound.

P1 §2.19 proves the missing mixed-branch propagation for classical $a_2$ and checks the matrices symbolically. The right figure has two different slices: all-equal couplings, and the slice with two fixed opposite signs. The latter has no Neumann zero at $C_0=0$. Use P1 §2.25's caption correction. “Invertible for arbitrary parameters” at L2:238 must mean **for generic spectral parameter**; isolated determinant zeros are allowed. At $C_0=C_1=0,\lambda=1$ the printed matrix is even zero, so the unqualified statement is false.

The $d_4$ normalization claim is **correct**, not a trap to be reversed. With roots $-e_1-e_2,e_1-e_2,e_2-e_3,e_3-e_4,e_3+e_4$, labels $(1,1,2,1,1)$, the zero-field boundary condition requires $C_j\sqrt{n_j}\propto n_j$. All $|C_j|=1$ fail at the central node. Independently constructed $8\times8$ vector Chevalley matrices give invertible solutions for **all 32 normalized sign patterns** at $\lambda=1.317$; maximal smallest singular value $9.68\times10^{-16}$, minimal normalized determinant modulus $3.49\times10^{-10}$. Omitting $\sqrt{n_j}$ replaces the central coefficient by $C_2/\sqrt2$; all 32 unit-sign patterns then have no kernel, minimum smallest singular value $0.184869$. This verifies the requested concrete exercise conclusion, with a faithful representation stated explicitly.

**New book/exercise wording issue:** the altered equations still decouple **node by node**; they simply contain $C_j/\sqrt{n_j}$. Thus the book's `exr-node-local` and nearby explanation that the factors “couple the nodes” are misleading. **Replace that clause in chapter 31:145 and the exercise:** “The equations remain node-local, but the allowed values of the newly normalized parameters acquire $n_j$ dependence.” In Lecture Exercise 2.1, make the vector representation and its spectator subspaces explicit; a determinant derived only for one defining $a_2$ module cannot be applied to an unnamed $d_4$ representation without explanation. The resulting nontrivial branch has $C_2^2=2$ in the altered normalization.

For Dirichlet, use P1 §2.5: the transpose Cartan condition is an **anticommutator**, not a diagonal commutation condition. The book's general dual-representation argument and its $\lambda^{2h}$ obstruction have a broader, explicitly stated scope; the lecture has not reproduced that argument.

The soliton-preserving sum is **correct** for real $\lambda$. Writing $\phi=u+iv$ makes the Cartan coefficient $-\beta v'$ and each step coefficient proportional to $e^{-\beta\alpha\cdot v/2}\sin(\beta\alpha\cdot u/2)$. Vanishing gives $v'=0$ and $u\in(2\pi/\beta)\Lambda_W^\vee$. Conjugation by $e^{i\pi\nu\cdot H}$ changes the latter to $u\in(\pi/\beta)\nu+(2\pi/\beta)\Lambda_W^\vee$, since its root eigenvalue is $(-1)^{\alpha\cdot\nu}$. The image relation imposes exactly these Dirichlet/Neumann components. Parity and $\phi\mapsto-\bar\phi$ each conjugate the Hirota species, so their product preserves it. Delius's source Proposition 1 explicitly starts “Let $\sigma$ be a linear or anti-linear involutive anti-automorphism”; the historical generalization is correctly attributed, with the Poisson/antilinear implementation qualifications stated in the book.

**Exercise 2.3 is wrong as written.** Conjugating the node equation by $D$ sends its step combination $A_j$ to $t_jA_j$, where $t_j=d_ad_b=\pm1$. Multiplying that node equation by $t_j$ gives $C_j\mapsto t_jC_j$ with the **same global $\lambda$**. It does not create independently flippable spectral parameters. For $a_2$, $\prod_jt_j=1$, so these flips alone have two four-element orbits. Inversion sends every $C_j$ to $-C_j$, hence sends their product $s$ to $-s$. Including inversion therefore gives **one eight-element orbit**, consistent with generation from the uniform solution, not two. In the product-$+1$ orbit only three patterns are solitonic; $(+++)$ is uniform. **Replace part 2's conclusion:** “Coweight flips alone preserve $s$ and give two four-element classes. Inversion exchanges these classes, so together the operations generate all eight patterns from $(+++)$. The $s=+1$ class consists of the uniform pattern and three solitonic patterns.”

Exercise 2.2's dictionary and Dirichlet limit are **correct**. Put $z=e^{i\beta_{\rm SG}\phi_0/2}$. With the division specified there,
\[
 K_{\rm Dir}=-\begin{pmatrix}0&z\lambda-z^{-1}\lambda^{-1}\\z^{-1}\lambda-z\lambda^{-1}&0\end{pmatrix}.
\]
It anticommutes with $H$. At fixed $\phi_0$ its off-diagonal part is proportional to $a_t$, while $a_t(1/\lambda)^T$ reverses the step part, proving compatibility for arbitrary normal derivative. `phase2_additional.py` verifies this symbolically. Neumann is $M_0=0$; uniform is $M_0=4m_0/\beta_{\rm SG}^2,\phi_0=0$ modulo the appropriate cosine period. This agrees with $C_0/C_1=e^{i\beta_{\rm SG}\phi_0}$ and the book's negative sign in the quantum Dirichlet relation.

## 9. Prior art and attribution audit

**Overlap:** P1 risks 17,18,22,25 included attribution concerns; only the Lee–Yang sources were fully checked there. The following fills that gap. **Verdict: mixed; the principal fault is converting conjectures and special cases into settled general discoveries.** Source files below are primary arXiv TeX under `referee_scratch/lectures_halfline/sources/`; article links identify the versions/sources independently of the scratch paths. Exact short quotations identify the relevant source claims; the calculations in §§1–8 test their application to the lecture conventions.

1. **Boundary duality — overstated, with substantial genuine prior art.** [Corrigan](https://arxiv.org/abs/hep-th/9707235), equations `neumannK`, `neumannKdual`, prints $K_N=(1+B/2)/[(B/2)(1)]$ and its dual $(2-B/2)/[(1-B/2)(1)]$, saying “Clearly, this is not self-dual”. This is the rank-one case of the lecture's products, as substitution $h=2,a=1$ independently verifies. [Delius–Gandenberger](https://arxiv.org/abs/hep-th/9904002), around `KNeum`, says “Therefore, we conjecture that the dual reflection amplitude”. Their general product agrees term by term with $K_a^N$. [Gandenberger](https://arxiv.org/abs/hep-th/9806003), §6, likewise says “we therefore conjecture that” before identifying the plus particle factor as Neumann. Thus neither the general product nor the idea of Neumann/uniform duality originated in these lectures or solely in the book. The book's added coideal/CDD interpretation is not the same priority claim. Replace Index:49 and L6:22 by: “Corrigan, Gandenberger, and Delius–Gandenberger proposed these boundary identifications and duality relations in the cases described below. We examine their general product formulas and the book's coideal interpretation, retaining the CDD and identification assumptions.” The $s=-1$ relation to the CDRS candidate is checked for general rank at §7 above; its coideal interpretation requires separate credit and qualification.

2. **Classical rigidity — correct for the simply-laced examples, overstated for all higher rank.** [CDRS](https://arxiv.org/abs/hep-th/9404108), `anboundary`, prints $|A_i|=2$ for every node or $A_i=0$ for every node. Substituting $A_i=2C_i$ gives precisely the lecture's $a_n$ restriction. [BCDR](https://arxiv.org/abs/hep-th/9501098), abstract, describes “a complete set of conjectures for permissible boundary conditions”. Its Lax calculation supports the normalized node equations; §8 independently checks the small-rank examples and $d_4$. Neither citation warrants the unrestricted no-Dirichlet claim. Use the scope replacements in §1, not a historical claim that the two papers prove every higher-rank assertion made in the index.

3. **Boundary symmetry — correct lineage, overly compressed account of the coideal.** [Mezincescu–Nepomechie](https://arxiv.org/abs/hep-th/9709078), abstract: “We construct integrals of motion (IM) for the sine-Gordon model with boundary at the free Fermion point”. The paper studies a boundary quantum group and also proposes extensions away from that point. [Delius–MacKay](https://arxiv.org/abs/hep-th/0112023) gives the boundary charges and explicit coideal construction used here. The independent coproduct calculation in §5 verifies the latter. Replace L4's “resolution, discovered ...” sentence by: “Mezincescu and Nepomechie constructed boundary integrals of motion at the sine-Gordon free-fermion point. Delius and MacKay developed the boundary quantum-group/coideal method used here.” Do not attribute the complete arbitrary-affine coideal result to the free-fermion calculation.

4. **Arbitrary involutive anti-automorphism — correct attribution.** Delius [hep-th/9809140](https://arxiv.org/abs/hep-th/9809140), Proposition 1, explicitly permits a linear or antilinear involutive anti-automorphism; the exact introductory phrase is quoted in §8. The independently differentiated double row and the $\dagger$ example there support the lecture's use. This is not a finding against Delius's generalization.

5. **Special coideal value and generalized Onsager relations — prior publication confirmed.** [Baseilhac–Belliard](https://arxiv.org/abs/0906.1215), Proposition 2.1 and subsequent simply-laced specialization, prints $\widehat\epsilon_j=\pm i/(q^{1/2}-q^{-1/2})$. Squaring gives $1/[(1-q)(1-q^{-1})]$, exactly the lecture value. Their relations specialize, with the normalization of §6, to $F_{ij}(b_i,b_j)=b_j$. The character-image construction is independently proved in §6; the additional assertion of an isomorphism is not needed. **Add at L5:96 and 209:** “This special value and the generalized $q$-Onsager relations were published by Baseilhac and Belliard; here they arise from the intertwiner/Serre calculation.” No claim of first publication by the book is justified.

6. **Classical images and static vacua — mostly correct, with incomplete index citation.** [Saleur–Skorik–Warner](https://arxiv.org/abs/hep-th/9408004), abstract: “method of images still works, despite the non-linearity of the problem.” Its subject is sine-Gordon. [Delius](https://arxiv.org/abs/hep-th/9807189), abstract: “mirror must be placed either at, in front, or behind the boundary.” This is the general affine-Toda construction used in L3. [Bowcock](https://arxiv.org/abs/hep-th/9609233), abstract: “using solitons which are analytically continued from imaginary coupling theories”. It constructs classical backgrounds and reflection data. DG builds on these backgrounds and relates their signs to weights. The direct tau-function tests in P1 and §4 confirm the lecture's construction after its branch/index repairs. **Replace Index's image citation by** `[@ssw1995; @delius1998a; @bowcock1998; @delius_gandenberger1999]`, and distinguish SG images from general Toda images. L3's “Following Delius ... and Bowcock” is supported; it does not excuse the incorrect interval assignment.

7. **Kim's one-loop amplitudes — checked for the explicit $a_2,a_3$ examples, not a universal all-rank calculation.** [Kim](https://arxiv.org/abs/hep-th/9506031) actually assumes self-duality and prints perturbative results and candidate factors for specified examples. The explicit $a_3$ coefficients independently match the lecture product derivatives (§7); the $a_2$ coefficient agrees with the expression printed in Gandenberger §6 and with the actual Kim product $[1/2][3/2]$. Kim's sentence about the earlier minimal factors is quoted in §7 and contradicts the lecture's blanket historical story. **Replace L6:151's parenthesis:** “For the explicitly checked $a_2$ and $a_3$ cases, its $O(\beta^2)$ coefficients agree with the one-loop expressions discussed by Kim and Gandenberger. This does not fix higher orders or the CDD factor.” Do not describe a two-loop failure without a two-loop calculation.

**Additional prior-art search:** targeted searches for affine-Toda boundary $K$-matrices, general $a_n$ conjugating solutions, and the centered crossing form found two relevant primary sources beyond the supplied list: [Gandenberger, hep-th/9911178](https://arxiv.org/abs/hep-th/9911178), whose abstract explicitly promises “explicit solutions” for the $a_n^{(1)}$ boundary Yang–Baxter equation, and [Corrigan–Delius, hep-th/9909145](https://arxiv.org/abs/hep-th/9909145), equation (2.7), for the centered crossing form. Their TeX was downloaded. The latter equation is compared independently in §3. The former should be cited alongside the general-rank matrix discussion at L5:114; I have **not** completed a convention-by-convention gauge equivalence between all of its matrices and the lecture's Toeplitz matrix, so its exact formula-level identification remains **unclear**. Discovery of prior art is not itself proof of such equivalence. No unqualified novelty recommendation is made for the all-rank matrix, character-twist interpretation or rotor quantization.

## 10. Remaining citation scope

**Overlap:** P1 risks 3,17,18,22,25; GZ/DPTW and the LY amplitude were checked there, most remaining citations were not. **Verdict: the following table separates source-supported statements from incomplete or overstated ones.** A “correct” source match does not overrule a mathematical correction elsewhere in this report.

| Source and citing claim | Exact source anchor and independent comparison | Verdict / replacement |
|---|---|---|
| GZ: space/time interchange, boundary state and crossing | Their §2 says “one could take $x$ to be the ‘euclidean time’”; §§3.2–3.3 supply the pair-state and crossing construction. The determinant, ZF index, integration-measure and Jacobian calculations in §§2–3 check its translation. Equations (3.45),(3.51),(3.54) are explicitly quoted/compared in P1 §2.22 and this report §2. | Attribution correct; printed lecture formulas still require the corrections already given. GZ's opposite LY sign is a real source statement, not an invented correction. |
| DPTW: identity Cardy boundary, $R_b$, shifted excited amplitude | Their `eq:rfb` is quoted in P1 §2.3. In their blocks $R_{(1)}=(1/2)(3/2)(2)^{-1}$; halving the lecture block labels gives exactly that product. Fusion multiplication in P1 proves $R_b^{(1)}=R_{b+2}=R_{4-b}$. The conformal weights are checked in §2. | Correct attribution. Equal reflection factors at $b=0$ do not identify the two conformal boundary conditions; repair Figure 1 as in P1 §2.25. |
| CDR: linearized reflection | [Source](https://arxiv.org/abs/hep-th/9407148), `\reflection`: $K_a=(ip_a+A m_a^2/(4m))/(ip_a-A m_a^2/(4m))$. Put $A=2C$, $p_a=m_a\sinh\theta$: this is exactly the lecture's rational factor. | Correct. The lecture's unified block formula fails at $C=0$ independently of the citation. |
| de Vega–González-Ruiz: SG matrix | [Source](https://arxiv.org/abs/hep-th/9306089), abstract, says “reflection matrices” for “eight vertex and six vertex models”; its six-vertex solution supplies the matrix part. The SG intertwiner calculation in P1 and §5 verifies the corresponding two-parameter matrix structure. | L4:158 should say “GZ found the physical SG reflection amplitude; de Vega and González-Ruiz solved the corresponding six-vertex reflection-matrix problem.” Their lattice matrix alone is not a physical SG scalar factor. |
| Gandenberger: $K^{(\pm)}$ and continued factors | hep-th/9806003 gives the two non-diagonal $a_2$ branches and their scalar/breather analysis. Its §6 explicitly identifies the plus branch only conjecturally (quote in §9). P1 §§2.9,2.18 and §6 independently verify the $3\times3$ matrices and singlet factor in the book's gauge. | Correct for $a_2$ matrices and candidate factors; not a source proving the all-rank coideal/Lagrangian dictionary. |
| DG: shifted crossing and solitonic signs | Equation (3.17) is the shifted diagonal crossing relation used in §3. Their weight formula gives $C_j=(-1)^{\alpha_j\cdot\lambda}$, hence $\prod C_j=1$ for $a_n$ because $\sum\alpha_j=0$. P1 §2.1 evaluates the actual boundary condition, not just this necessary product. | Correct; the lecture's erroneous modulus-to-pattern table is its own addition. Uniform must be excluded from the word “solitonic”. |
| BPTT: $(\eta,\vartheta)$ | [Source](https://arxiv.org/abs/hep-th/0106070), reflection-matrix and UV–IR displays use $\cos(\eta/(\lambda+1))\cosh(\vartheta/(\lambda+1))$ and the corresponding sine/sinh product. Matching the SG matrix to its $P^\pm$ gives the cosines in the lecture's $\hat\epsilon_j$ dictionary; this is a conversion, not BPTT's original coideal notation. | Correct parameter attribution, subject to the book's spatial-orientation sign and UV–IR status. Replace “Following ... we parametrize” by “Using BPTT's physical parameters, and the book's orientation and node conventions, we write ...”. |
| DTW: zeros reducing Coleman–Thun pole order | Their source says “Thus the prediction is demoted to first order, as required.” The local order subtraction and residue test are independently checked in §7/P1 §2.6. | Correct qualitative mechanism. The source does not print the claimed general $P-2L-Z$ formula or prove this particular $a_2$ topology. Use §7's qualified replacement. |
| Bernard–LeClair: non-local charges and $U_q(\hat g^\vee)$ | Full cited 1991 paper not obtained; abstracts and later citations cannot certify its precise dual algebra or chirality conventions. The coproduct and SG weight calculation were checked independently in §5. | **Unclear citation scope**. Retain an explicitly qualified attribution pending a source read; do not claim its chirality labels match the book. |
| Cherednik (1984), Sklyanin (1988): factorization and double row | Bibliographic/abstract records found, but no complete cited primary text was reliably inspected in this pass. The operator histories and transport derivative were independently checked in §§3,8. | **Unclear source-level attribution**, not a claim that standard historical credit is false. The derivations are independently assessable. |
| Letzter (1999), Kolb (2014): coideal direction | Letzter's cited paper was not obtained in full. Kolb's definition is explicitly quoted in §6 and uses $\Delta(B)\subset B\otimes U$ in its coproduct convention. Direct expansion here gives $U\otimes B$ (the spatial right boundary). | Letzter citation **unclear**; Kolb scope correct only after specifying opposite coproduct. “Left” is a convention, not a proof that the two displayed coproducts coincide. |
| Baseilhac–Koizumi (2005) | [Source](https://arxiv.org/abs/hep-th/0507053), abstract: “XXZ open spin chain with general integrable boundary conditions” and “q-deformed analogue” of Onsager's algebra. The rank-one Dolan–Grady relation is independently tested in §6. | Correct rank-one antecedent. It is not the source of the general simply-laced special-parameter constraint; cite Baseilhac–Belliard for that. |
| Kolb, Balagović–Kolb, Appel–Vlaar: parameter sets and universal $K$ | Full primary texts were inspected; exact definitions and scope qualifications are in §6. The split-case reduction and character construction are derived there. | Repair the generality of the odd-neighbor condition; separate finite-type universal $K$ from the later Kac–Moody/completed twisted construction. Neither alone fixes a physical scalar factor. |
| Sasaki, Fring–Köberle, Kim: early bootstrap/self-duality story | Primary TeX was read. Sasaki and Fring–Köberle supply bootstrap constructions; Kim explicitly distinguishes their minimal factors from his perturbative result (quoted §7). | The collective claim of matched Neumann one-loop amplitudes is **wrong**. Whether each earlier author consciously imposed the same self-duality assumption is **unclear** where no explicit source statement was located; do not manufacture such a statement. |

These qualifications are source-specific. Reading later papers that cite an unavailable original does not count as reading that original.

## 11. All eighteen exercises: solution outlines and link audit

**Overlap:** P1 risks 4,8–15,19,20,23,24 covered most mathematical subparts. Exercises 2.1,2.2,2.3 and 4.1 receive the additional derivations in this phase. **Verdict: several are sound after small repairs; a number are impossible as printed.** “Solvable” below means with the displayed definitions and the indicated book link, at generic parameters. No link in a heading is treated as a proof.

### Lecture 1

**1.1 `exr-lec1-elastic` — wrong intermediate limit; repairable.** The outline is §3: let $r$ be the multiplicity of the maximal outgoing rapidity, divide, force equality of maximal rapidities, then $r=1$, then subtract the surviving cosh. All remaining terms are positive. The printed claim that a finite limit requires only one maximizer is false. “Same mass” is insufficient for identical higher-charge coefficients if there are multiple species. With the equal-coefficient assumption made explicit this is a short, solvable exercise. Heading target `exr-boundary-elastic` is the corresponding book exercise.

**1.2 `exr-lec1-twisted` — wrong as printed.** Starting at $V_\mu\otimes V_\nu$, apply the rightmost operator first, exchange the two labels at each braid, and conjugate the label only at reflection. One printed mixed braid has the wrong domain. P1 §2.11 supplies the full four-step space trace and the corrected equation. After that repair the exercise is true and entirely solvable from the definitions. Heading target `exr-twisted-domains` is appropriate.

**1.3 `exr-lec1-crossing-forms` — parts 1–3 correct algebraically; part 4 wrong.** The first two transitions use $\theta\mapsto\theta\pm i\pi/2$, bulk crossing and $S(\theta)S(-\theta)=1$, and the converse also uses boundary unitarity. The appropriate attribution of the centered form is §3. In the $h=6$ convention, direct sinh factorization gives $f_x=-(6x)(6-6x)$, not the printed $(2x)(6-2x)$. Therefore $K_b/K_{(1)}=f_{(1+b)/6}/f_{(1-b)/6}$ is correct even though the preceding block identity is false: the two minus signs cancel. P1 §2.4 gives the ready replacement. The heading's `exr-crossing-forms` target is appropriate for the crossing part; its text purports to link two book exercises while linking only one. Replace the label by the title of the actual linked exercise, and add a separate CDD exercise link if desired.

### Lecture 2

**2.1 `exr-lec2-normalization` — correct numerical conclusion; incomplete derivation setup.** Equating root-exponential coefficients cancels $c_j=\sqrt{n_j}$ in the normalized simply-laced equation. Omitting it gives $C_j/\sqrt{n_j}$ instead. In the $8$-dimensional $d_4$ vector representation, each root acts on two doublets and leaves spectators; the same spectator determinant enforces the nontrivial branch $(C_j/\sqrt{n_j})^2=1$, hence $C_2^2=2$ at the central node. Independent matrices and all 32 patterns are tested in §8. State this representation so a student knows where the spectator calculation applies. Heading `exr-node-local` corresponds, but the book's claim that omission destroys node locality is wrong; §8 supplies its replacement.

**2.2 `exr-lec2-a1-dirichlet` — correct.** Substitute $\beta_{\rm SG}^2=2\beta^2$, $m_0=2m$ into the two exponentials of the boundary potential and combine them into $-M_0\cos[\beta_{\rm SG}(\phi-\phi_0)/2]$ with the specified potential sign. The ratios are $C_0/C_1=e^{i\beta_{\rm SG}\phi_0}$ and the magnitudes give Neumann/uniform values in §8. Dividing the explicit classical matrix and taking $M_0\to\infty$ gives the displayed $K_{\rm Dir}$ there; direct multiplication verifies compatibility. Its sign agrees with the book's quantum orientation convention, $\eta=-4\pi\phi_0/\beta_{\rm SG}$. The heading target `exr-gz-dictionary` covers the dictionary, but the text “Exercises ... & ...” supplies only that one link. Link the Dirichlet object separately; the calculation itself is well posed.

**2.3 `exr-lec2-coweight-symmetries` — wrong class-count wording and spectral-parameter claim.** Let $D=\operatorname{diag}(d_1,d_2,d_3)$, $d_i=\pm1$. The node-step signs are $t_j=d_ad_b$, with $\prod t_j=1$. Conjugation changes $C_j$ to $t_jC_j$ while leaving $\lambda$ fixed; inversion changes all $C_j$ to $-C_j$. Hence coweights alone give two classes and coweights plus inversion one class. There are three solitonic patterns and one uniform pattern in the $s=+1$ class. The replacement in §8 makes the exercise true and solvable. There is no heading link to audit.

### Lecture 3

**3.1 `exr-lec3-mirror-square` — wrong $k$ and derivative, correct target identity after repair.** Put $z=\zeta^{ja}e^{\xi_1-m_at\sinh\theta}$. The boundary tau polynomial is a perfect square after using the prescribed shift. Differentiate the *bulk* tau first, then evaluate at zero; differentiating the simplified boundary square with respect to $x$ loses information. Use $k=2m\sin^2(\pi a/h)/s$ and the corrected derivative in P1 §2.10. Substitution proves the boundary condition with the prescribed square-root/vacuum branch. The heading points to `exr-mirror-neumann`, which is the **Neumann** exercise, not this $C=\pm1$ calculation. Replace the heading link by `[Mirror construction](book:sec-images-boundary)` (the section exists), or the actual uniform-boundary proposition. As printed, a student cannot derive the target using the supplied false intermediate identities.

**3.2 `exr-lec3-breathers` — partially true, regularity conclusion too broad.** Summing the two exponentials gives $2e^{y+S/2}\cosh(\rho_-+i\psi)$. At $\rho_-=0$, use $\cos\psi\ge-1$ to obtain the positive lower bound on $y<0$. This does not exclude a wall zero when $\psi=\pi$. For $C=-1$, $s<1$, factor the special phase as $(1-e^y)(1-(1+s)e^y/(1-s))$, giving an interior zero; for $s>1$ the indicated wall value is negative. At $s=1$ the normalization used in the question is undefined. P1 §2.23 gives the corrected cases and endpoint warning. The computation excludes regularity for this real $\rho_-=0$ slice, not every complex phase/modulus choice without additional argument. Heading `prp-boundary-breathers` is the correct subject. Replace the title/conclusion by “Regularity tests for the displayed mirror-pair family” unless the omitted parameter-range proof is supplied.

**3.3 `exr-lec3-self-conjugate` — correct, with strip-edge wording.** $u'_j=2mu_j$ and neighboring $u$'s are $-u_j$, giving the printed logarithmic derivative and square-root ratio. The boundary combination is constant in $j$, so its root sum vanishes. The energy variation in $\xi$ is zero by the bulk equation and boundary condition, and the vacuum limit gives $E=0$. The pole has $p=2mi$, frequency zero, at $\theta=i\pi/2$; $\partial_\xi\phi$ is its normalizable zero mode. Say “at the edge of the physical strip”, not inside it. Heading `prp-self-conjugate-vacuum` corresponds. The energy argument requires the standard variational identity, now supplied in §4.

### Lecture 4

**4.1 `exr-lec4-general-cj` — product constraint correct generically; uniqueness of the sign is overstated.** For $\hat\epsilon=0$,
\[
 A_1(y)=\begin{pmatrix}0&y\\c_1q/y&0\end{pmatrix},\qquad
 A_0(y)=\begin{pmatrix}0&c_0q/y\\y&0\end{pmatrix}.
\]
For $K=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)$, the equations give $d=y^2a$, $a=y^2d$, $c=c_1qb$, $b=c_0qc$. If $y^4\ne1$, $a=d=0$ and a nonzero solution exists exactly when $c_0c_1=q^{-2}$; it is invertible. At exceptional $y$, extra solutions exist, so state genericity. The coproduct follows by multiplying $\Delta(f)\Delta(k)$ and needs no coefficient constraint. Parity gives $c_0=c_1$, hence **both** signs $c_0=c_1=\pm q^{-1}$; selecting plus is an additional normalization choice. **Replace the last sentence:** “For generic $y$, existence fixes the product. We choose the parity-symmetric normalization $c_0=c_1=q^{-1}$; parity alone would also allow their simultaneous negatives.” Independent generic sample residuals are $4.11\times10^{-16}$ on the constraint and $0.0741$ off it. There is no heading link; the body link to `prp-sg-k-unique` is appropriate. With these repairs the exercise is solvable.

**4.2 `exr-lec4-neumann-eps` — correct conditional calculation with missing bounds/branch.** Insert $\eta=\pi(\lambda+1)/2$ into $-\cos\eta/\sin\pi\lambda$ to get $1/[2\cos(\pi\lambda/2)]$. Its square is $\epsilon_*^2$; the displayed square root means its analytic branch, not necessarily the positive root. Then $\nu_n=\pi/2-n\pi/\lambda$ and the elementary double-angle identity gives $M\sin(n\pi/\lambda)=m_n\cos(n\pi/(2\lambda))=m_{2n}/2$. Require $0<n<\lambda/2$ for an interior pole and the indicated breather. The $n=0$ edge singularity is not a new zero-energy bound state or a particle $B_0$. Images explain the half-energy kinematics but do not prove quantum state existence from the classical solution. There is no heading link. Mark its input as the UV–IR conjecture and apply §5's replacements.

**4.3 `exr-lec4-1d-submodule` — correct for the specified branch.** Solving the two row equations for $v_n=(1,\alpha_n)^T$ gives the common eigenvector and $\beta_{1,0}$ written in P1 §2.12. Substituting $y_n=i(-1)^ne^{-i\eta}$ and $q=-e^{-i\pi\lambda}$ turns these eigenvalues into $-\cos(\bar\eta\pm i\vartheta)/\sin\pi\lambda$. Unlike the surrounding prose, this exercise does not claim these are the only reducibility points. It is solvable at generic parameters. No heading link; the body link `eq-sg-coideal-rep` is correct. Repair the abandoned determinant paragraph in the body using §5.

### Lecture 5

**5.1 `exr-lec5-toeplitz-solve` — final formula correct; node attribution in part 1 wrong.** The node-0 spectator gives $K_{32}=sXK_{12}$ and $K_{21}=sXK_{23}$. One also needs the node-1 spectator equation $K_{31}=X^{-1}K_{32}$ to get $K_{31}=sX^{-1/2}$. The printed $(1,1)$ equation then gives
\[
 K_{11}=\frac{s yX^{1/2}-y^{-1}X^{-1/2}}{\epsilon_*(q-q^{-1})},
\]
which reduces to $\kappa_s(X)$ with the declared $w$ branch. **Change part 1 to “Use node 0 and the node-1 spectator equation ...”.** The existing final body link `eq-an-k-matrix` corresponds. No heading link. With that change the outline is complete.

**5.2 `exr-lec5-serre-ej` — correct.** The bulk Serre relation removes the $e_i^2$ terms. For a component $v$ with $ve_j=q^{-1}e_jv$, the $v^2$ coefficient is $q^{-2}-(q+q^{-1})q^{-1}+1=0$. The mixed Cartan terms cancel using the same commutation rules. The surviving $e_i\varphi_i$ terms are $-q^2D_ie_j+e_jD_i$; commuting $k_i^2$ past $e_j$ reduces them to $qc_ie_j$, as calculated explicitly in §6. Solvable from the stated relations. No heading link; body `prp-serre-obstruction` corresponds.

**5.3 `exr-lec5-cyclic-unitarity` — wrong transpose identity, otherwise correct.** $Z$ is a signed permutation matrix, so $Z^TZ=1$ for either $s$ and $Z^T=Z^{-1}$, not $sZ^{-1}$. Expanding the polynomial in $Z$ reproduces the matrix; evaluating it at $\zeta^3=s$ gives the stated rational eigenvalue. At reflected rapidity $X'=w^2/X$ and transposed eigenvalue $\zeta'=\zeta^{-1}$, multiply and use $\zeta^3=s$; P1 §2.9 gives the resulting scalar $U_s$ and independent checks of both signs. Remove the extra $s$ in part 3; no other change to the target is needed. No heading link; the body section link `sec-an-nondiagonal-k` corresponds. The transpose identifies the reverse reflection matrix in this fixed dual basis, not in every arbitrary gauge.

### Lecture 6

**6.1 `exr-lec6-bootstrap-a2` — wrong fusion-angle label, correct bootstrap and doubling identity.** The mass triangle gives $m^2=2m^2(1+\cos u)$, so $u=2\pi/3$ and $\bar u=\pi/3$. The constituent shifts are therefore the printed $\pm i\pi/3$. Here $S_{11}=\{1\}_B$ is correct. Set $b=B/2$. Shifted reflection factors multiply to $(4-b)(2-b)/[(3)(1)(2-b)(-b)]=(4-b)/[(3)(1)(-b)]$. Doubling the bulk blocks gives $S_{11}(2\theta)=-(1)(4)/[(b)(b+3)(1-b)(4-b)]$ (with $(0)(3)$ accounted for in the signs). Multiplication leaves $-(4)/[(3)(b+3)(1-b)]$; using $(3)=-1$, $(4)=(2)^{-1}$ and $(b+3)=(3-b)^{-1}$ gives $(3-b)/[(2)(1-b)]=K_1$. Thus the exercise's doubling identity is **correct**, unlike L1's use of it. P1 §2.14 and §3's independent evaluation agree. Swap the two introductory angle labels. No heading link; body `sec-uniform-boundary` is appropriate.

**6.2 `exr-lec6-kim-oneloop` — part 1 correct, part 2 has a false equality, part 3 is not well posed.** Logarithmic differentiation gives $\partial_x\log(x)=i\pi\sinh\theta/[h(\cosh\theta-\cos(\pi x/h))]$, proving part 1. Apply it to the *first* expression for $K_1^N$; the coefficient is the stated $-\beta^2\cos(\theta/2i)/(4\sin(3\theta/2i))$. The second expression is not equal to the first. Delete it. The purported self-dual product in part 3 is not self-dual and changes the first derivative; being made of $[x]$ blocks alone does not force a particular one-loop coefficient. **Replace part 3 by:** “Let $K^{\rm Kim}=[1/2][3/2]$. Prove its self-duality, verify that its derivative at $B=0$ agrees with $K_1^N$, and show their second derivatives differ. Explain why one-loop agreement cannot select between them.” The new independent script gives first-derivative difference $4.08\times10^{-56}$ and $B^2$ coefficient $0.37357808231122i$ in $K^{\rm Kim}-K_1^N$ at $\theta=0.71$; this is a concrete solvable question. Replace “true” by “proposed”, and “failed at two loops” by “differs at order $\beta^4$”. No heading link; body `sec-neumann-duality` is appropriate.

**6.3 `exr-lec6-rotor-spectrum` — energy identities correct, interpretation needs qualification.** For nonnegative $n,m$, summing geometric exponentials yields the displayed $F(n)$ and $F(m)$, independently for the two species. “At the same rapidities” means the same formula with its own index $k=0,\ldots,m-1$, not necessarily the same number of particles or simultaneous collisions. Negative labels are bootstrap-state labels, not negative particle counts. Explicitly, $F(-1)+F(1)=\sin\delta$, matching the separate diagonal binding pole. Also $F(-2)+F(2)=\sin^2(2\delta)/\sin\delta$ equals $F(1)$ at $\delta=\pi/12$, proving the energy degeneracy at $B=1/2$. This equality does not prove distinct simple poles at an exceptional coupling or the global rotor quantization. Supply the state ranges from P1 §2.24 and the explanation of negative labels in §7. No heading link; body `sec-boundary-spectrum-a2` and `eq-a2-energies` correspond. Parts 1–3 are then solvable with the stated distinctions.

## 12. What the lecture scripts actually check

**Overlap:** P1 risk 25, checked fully as a source/figure audit; risks 6,9,10,16,18 exposed specific gaps. P1 §2.25 remains the evidence for the six inspected SVGs and the complete assertion inventory. No author script was run in either phase. This section maps **every explicit script-verification sentence located in the seven lecture/index files** to the actual assertion and gives the necessary narrower wording. **Verdict: several true sample checks are advertised as proofs of different or stronger statements.** Independent evidence is P1's calculations and the additional derivations in §§1–11, not the author PASS lines.

| Location | What the script actually tests | Verdict and replacement scope |
|---|---|---|
| Index:69, “Every non-trivial calculation ... from first principles” | Selected matrix/scalar assertions and plot generation; it does not implement every derivation, exercise or attribution. | **Wrong.** Use P1 §2.25's replacement about selected formulas and numerical examples. |
| L1:162, braid conversion “checked numerically” | A permutation-symmetric six-vertex matrix with first-slot $K$; subtracts residuals of two equations that both already vanish at the sample. | **Weaker and convention-insensitive.** Say “The script tests two reflection equations in a symmetric six-vertex example; the general conversion is the algebraic permutation calculation above.” Correct the conversion first (P1 §2.11). |
| L1:264, three crossing forms “verified side-by-side” | The three scalar identities for the LY example. | **Correct as a sampled LY check**, not a derivation of equivalence for all species. Add “for the Lee–Yang example at the sampled rapidities”. |
| L1:266, figure “Verified and plotted” | Schematic plus analytic pole-position curves; neither determinant of the coordinate map nor all pole cancellations/boundary identifications is tested. | **Wrong caption assurance.** Correct the rotation, tick labels and $b=0$ identity-wall annotation (P1 §2.25); say “Schematic and candidate pole-position curves generated by ...”. |
| L1:404, LY sign/bootstrap | Direct bootstrap residual for both signs. | **Correct test of the final amplitude**, independent of the erroneous intervening proof. Say “The final bootstrap identity is also tested numerically”; use §2's complete derivation. |
| L2:172, both panels and iff rigidity | Explicit matrix tests, all eight $a_2$ sign choices, one off-branch sample, plotted slices. | **Weaker than a classification** and the orange slice does not vanish at zero. State the two slices, and cite §8's independent propagation proof for the classification. The transport picture is schematic. |
| L2:238, arbitrary complex $C_0,C_1$ and invertibility | Particular parameter/spectral samples. | **Weaker.** Replace by “The formula solves the node equations for arbitrary boundary parameters and is invertible at generic spectral parameter; the script tests selected samples.” P1 supplies symbolic evidence for the formula. |
| L2:354, every coweight soliton-preserving extension | One explicitly constructed $a_2$, $\dagger$ boundary configuration and its matrix compatibility. | **Weaker.** Say “The script illustrates this identity for one $a_2$ coweight; the root-phase argument proves the stated extension.” The general derivation is §8. |
| L3:127, mirror condition “to $10^{-14}$” | The correct mirror formula with the requisite vacuum shift, using correct code kinematics. | **Correct sampled final check; not the printed derivation.** Repair $k$ and derivative and say “At the sampled parameters, the residual is below the script's $10^{-14}$ tolerance.” A tolerance is not a uniform error bound in $t$ or coupling. |
| L3:175, sign pattern/energy/time-delay figure | Chooses whichever of three sign patterns minimizes each vacuum residual; tests EOM, BC and energy at a few modulus centers. | **Wrong as verification of the stated assignment.** It cannot check a mapping that it selects after looking at the answer. Replace the caption/legend by $(+,-,-)$ at $\chi=.25$; describe the tests as sampled EOM/BC/energy checks. The independent interval derivation is P1 §2.1. |
| L4:263, iff sign of $q$, Klein-twisted conjugation | Projects the crossed vector onto a one-dimensional span, at two signs of $q$ with fixed physical parameters. It checks proportionality, not the required scalar factor. It does not test all couplings or the full Klein-twisted theory. | **Overstated.** Say “At the sampled generic parameters, plain-conjugation matrix crossing is compatible with the physical sign and fails for the opposite sign. Full crossing additionally requires the scalar factor; the twist must act consistently on the scattering conventions.” P1 §2.16 and §5 give independent matrix evidence. |
| L4:312, claimed residual/energy panels | Actual figure contains a coideal schematic and pole trajectories. | **Wrong.** Delete the claims that crossing residuals and excitation energies are plotted, as in P1 §2.25. |
| L4:353, common eigenvalues | Substitution at selected $y_n$, yielding small residuals for the specified eigenvector/eigenvalues. | **Correct on the selected branch**, not a proof of exclusive reducibility. Say “This reproduces the eigenvalue check implemented at selected poles.” The algebraic solution is P1 §2.12. |
| L5:131, $U_s$ and matrix crossing | Explicit $3\times3$ products for both signs and matrix crossing proportionality at selected values. The mixed-module reflection-equation routine elsewhere fits a scalar. | **Correct for the advertised matrix checks at samples**, not full scalar boundary crossing or normalized mixed reflection equations. Keep “matrix” and state the range. P1 §2.9 independently verifies $U_s$; §5 explains the normalization gap. |
| L5:212, existence graph and “exact” UV–IR curve | Samples the one-parameter equal-$\hat\epsilon$ slice; plots the formula implied by the SG UV–IR relation. | **Overstated and one endpoint mislabeled.** A slice is not the mixed-parameter classification. Label the curve conditional on UV–IR and its zero by $\mu_0$, not 1 (P1 §2.17). |
| L6:119, general $c_+(\mu)$ | Constructs one $a_2$ singlet and the fusion operator; compares to the $h=3$ specialization with tolerance $10^{-13}$. | **Correct for $a_2$ at the tested parameters, weaker than all rank.** Replace parenthesis by “The $a_2$ specialization is tested in ...”. P1 independently repeats the construction; no all-rank scalar/Gamma cancellation is thereby checked. |
| L6:227, CT identity “to machine precision” | Symmetric finite differences at step $10^{-7}$, acceptance threshold $10^{-8}$; checks the residue/zero identity, not the propagator or loop counts. | **Wrong precision/scope claim.** Use “The script tests this residue identity with relative tolerance $10^{-8}$.” The independently controlled residue limit agrees; the topology remains unclear. |
| L6:229, duality and spectral figure | Plots the two chosen block formulas (related identically by $B\mapsto2-B$), energy expressions and allowed diagonal ranges. | **Correct formula plots; no independent boundary identification or rotor quantization.** Caption: “Duality of the proposed amplitudes and the rotor-like diagonal tower.” The code never integrates the classical inertia despite its header suggesting that calculation. |

The scripts' README adds claims beyond these sentences; its stale IDs, wrong $K_{-b}$ shift, wrong vacuum energy and reversed attraction labels were found only by P1's broader inventory (P1 §2.25). They should be corrected, but were left untouched. Neither phase's source audit asserts that the author scripts successfully execute in this environment. Their reported assertion counts are not scientific coverage metrics.

## 13. Overlap table

“Wrong/mixed” means at least one stated component is wrong, with surviving correct components and explicit unclear limits in its section. A broad coordinator item is not assigned a single universal truth value. P1 ranks are fixed and equal their entry numbers.

### Coordinator items

| Claim | Coordinator list | P1 list (rank; checked?) | Verdict |
|---|---|---|---|
| Status, scope, physical dictionary and CDD duality | 1 | 2,5,17,18,24,25; all checked in stated scope | Wrong overclaims; conditional product identities correct |
| Lee–Yang derivation, residues, Cardy data and CDD exercise | 2 | 3,4,21,22; checked | Mixed: doubling and bulk coupling wrong; shift and $g$ correct |
| Conservation, reflection types, rotation, crossing, bootstrap, elasticity | 3 | 11,20,21,25; checked | Wrong signs/indices/proof steps; repaired identities given |
| Mirror solitons, vacua, parity and Robin factor | 4 | 1,2,10,15,23,25; checked | Mixed: several false intermediate formulas; parity/time-delay results correct |
| Coideal perturbation, SG crossing, excited line and Schur argument | 5 | 12,16,17,21,25; checked with limits | Mixed; CPT normalization and full scalar crossing unclear |
| $a_2$ matrices, Serre, QSP, characters and dictionary | 6 | 8,9,17,18,25; checked with limits | Matrix/Serre identities largely correct; proof/scope claims wrong |
| Scalar fusion, continuation, amplitudes, spectra, rotor and CT | 7 | 2,6,7,13,14,17,18,24,25; checked with limits | Mixed; global rotor and diagram topology unclear |
| Sklyanin, classical rigidity, normalization, Dirichlet, images and coweights | 8 | 5,19,25; checked | Transport/matrices and $d_4$ signs correct; scope/class-count claims wrong |
| Prior art and named discoveries | 9 | 17,18,22,25; mostly incomplete historically in P1 | Mixed; conjectures overstated, genuine prior art retained |
| Remaining citations | 10 | 3,17,18,22,25; selected sources checked | Mixed; four unavailable originals explicitly unclear |
| All 18 exercises | 11 | 4,8–15,19,20,23,24; relevant parts checked | Mixed; per-exercise verdicts in §11 |
| Script-verification sentences and figures | 12 | 25 and 6,9,10,16,18; checked | Overstated coverage and several wrong captions |

### Every P1 risk entry

| Claim | Coordinator list | P1 entry/rank | Verdict |
|---|---|---|---|
| Modulus intervals and vacuum signs | 4,12 | 1 | Wrong assignments/legend; integer-index formula correct |
| Classical inertia and global rotor | 1,7,11 | 2 | Inertia correct; asserted global quantization not established |
| $K_b^{(1)}=K_{b+2}$ | 2,10 | 3 | Correct |
| LY doubling proof and $f_x$ blocks | 2,11 | 4 | Wrong |
| Dirichlet Cartan equation and scope | 1,8 | 5 | Wrong lecture proof/scope |
| CT arithmetic/residue/diagram | 7,12 | 6 | Arithmetic wrong, residue correct, topology unclear |
| Imaginary-to-real block continuation | 7 | 7 | Correct after using periodicity |
| Serre expansion and Cartan obstruction argument | 6,11 | 8 | Relation correct; omitted scalar/inverse argument wrong |
| Matrix unitarity and cyclic shift | 6,11,12 | 9 | $U_s$ correct; $Z^T=sZ^{-1}$ wrong |
| Mirror perfect square and kinematics | 4,11,12 | 10 | Wrong intermediate identities; repaired result correct |
| Braid and conjugating equation domains | 3,11,12 | 11 | Wrong printed convention/domain |
| SG common eigenline | 5,11 | 12 | Specified line/eigenvalues correct; exclusivity/determinant wrong |
| One-loop exercise | 7,9,11 | 13 | Coefficient correct; alternative equality/ansatz wrong |
| $a_2$ fusion angle | 7,11 | 14 | Wrong angle label; bootstrap correct |
| Uniform Robin reflection | 4,10 | 15 | Rational expression correct; $C=0$ block equality wrong |
| Sign of $q$ and crossing | 5,12 | 16 | Matrix obstruction correct; general/full claim overstated |
| Neumann/uniform status and UV–IR | 1,6,7,9 | 17 | Wrong promotion of conditional identifications |
| Singlet fusion and scalar factors | 1,7,12 | 18 | Checked $a_2$ identities correct; all-rank scalar derivation unclear |
| Classical matrices and mixed rigidity | 8 | 19 | Matrices correct; omitted mixed-case proof supplied |
| Continuity equation and Euclidean rotation | 3 | 20 | Wrong |
| Pair-only boundary state and pole interpretations | 3,5,10 | 21 | Overstated; zero-rapidity term and other pole mechanisms needed |
| GZ LY sign and residue normalization | 2,10 | 22 | Sign quotation correct; lecture $\Gamma^2$ wrong |
| Tau positivity and breather/self-conjugate exercises | 4,11 | 23 | Root proof wrong; positivity/zero-energy family correct with endpoint limits |
| Excited energies and physical ranges | 7,11 | 24 | Energy identities correct; physical range/interpretation incomplete |
| Script/figure assertions and README | 12; README specifics not expressly listed | 25 | Coverage/captions/README wrong in specified places |

## 14. What each list found

This comparison uses **claims**, not which phase happened to compute them first. Every one of my 25 ranked risks falls within the coordinator's unusually broad items. It would therefore be misleading to report many wholly disjoint blind discoveries. The overlap does not make the blind phase redundant: its list and evidence existed before the coordinator prompts were read.

**Wrong findings reached through both lists:** the LY doubling/proof and block exercise; LY bulk-coupling normalization; wrong conservation and rotation steps; braid/twisted-domain mismatch; Dirichlet anticommutator and scope; mirror $k$/derivative and vacuum-sign table; Robin $C=0$ block formula; positivity/root and breather-endpoint mistakes; common-eigenline determinant/exclusivity wording; Serre obstruction's missing scalar and inverse argument; cyclic transpose; false one-loop exercise equality and proposed self-dual ansatz; fusion-angle label; CT arithmetic and precision claim; promoted physical identifications; missing spectrum ranges; misleading script/figure assurances. Neither list found a counterexample to the LY excited shift, inertia coefficient, continuation, matrix unitarity or checked energy identities.

**Wrong findings reached only through the blind list's extra inventory:** the script README's stale/nonexistent target IDs, $K_{-b}$ in place of the actual shifted LY factor, incorrect stationary-vacuum energy, and reversed attractive/repulsive sign descriptions. The coordinator requested lecture verification sentences, not a README audit. Figure-axis and legend problems were first noticed blind, but fall under coordinator item 12's broad figure coverage and are therefore counted as overlap, not inflated into exclusive findings.

**Additional wrong findings exposed by the coordinator's list:** the lower-index order in the non-diagonal pair crossing equation; L4's wrong bulk crossing argument; the extra charge-normalization factor of two; the reflection/pair residue Jacobian sign; $S_{b\bar a}$ in the diagonal bootstrap; the faulty maximal-rapidity limit and equal-mass/equal-charge conflation; generic “spin-1” currents and the reversed incoming chirality; the Schur normalization argument; node-0-only attribution in the Toeplitz exercise; the unrestricted statement of Kolb's parameter condition and unqualified affine universal-$K$ citation; the coweight/inversion class count and spectral-parameter statement; the book's assertion that missing root normalization destroys node locality; the boundary-state tensor-space error in L6; the parity-only uniqueness in Exercise 4.1; the wrong book target in the mirror-square exercise; and the historical group claim that the early minimal factors all matched one-loop Neumann scattering. These are additional to P1's frozen verdicts, including cases P1 had mentioned only as unverified leads.

The blind list was effective at local algebra, explicit counterexamples, branch/index errors, plotted conventions and the unadvertised README. The coordinator's list was better at auditing provenance, status labels, comparisons between distant lecture passages and book objects, and hypotheses hidden behind broad explanatory sentences. Its leading questions also required restraint: the general CDRS attribution, all 32 $d_4$ signs, coweight nature of the parity shift, CT state name $b_{-1,2}$, inertia integral and several energy identities survive independent checks. A suspicion is not a defect.

## 15. Claims left unclear and verification limits

The requested items and all 18 exercises have been addressed, but this is **not a certificate of every formula at all ranks or all couplings**. The following claims remain unclear, rather than being assigned an unsupported correct/wrong verdict:

- The full first-order boundary-OPE normalization for every non-simply-laced/dual algebra, including every image-field coefficient, topological-charge rescaling and logarithm branch in $q^{T_j}$. §5 checks the contact-distribution sign, SG weights, coproduct and the stated simply-laced normalization, and identifies the remaining gap.
- Full **scalar** SG crossing-unitarity at generic parameters and every exceptional point; the independent tests establish the matrix obstruction/proportionality and convention repairs, not a complete scalar solution.
- An analytical all-rank proof of the non-diagonal scalar Gamma product, its full telescoping fusion calculation, and a unique physical CDD choice. The book/lecture equations were compared, the $a_2$ matrix factor independently computed, and the CDD nonuniqueness demonstrated. No unproved scalar claim was used to certify a physical dictionary.
- A complete general-rank quantum intertwiner classification proof beyond the book's sketch; the explicit $a_2$ mixed-case propagation and sampled larger-rank Toeplitz tests are evidence of their stated scope. Nor does the character-image calculation prove injectivity of every universal presentation or the applicability of every completed universal $K$ to every physical evaluation module.
- The separate internal-propagator/loop counts of the proposed boundary Coleman–Thun graph. Its printed arithmetic is wrong and its residue identity is correct; those facts do not uniquely reconstruct a diagram. DTW's zero-reduction mechanism alone does not supply that reconstruction.
- Quantum gluing and endpoint conditions for the classical modulus intervals. The inertia is computed, the leading spectrum agrees, but a global circle Hilbert space is not derived. Similarly no full classification of all regular complex breather parameters or completeness of the boundary spectrum at exceptional rational couplings is claimed.
- Exact source scope of the **four unavailable non-arXiv originals**: Cherednik (1984), Sklyanin (1988), Bernard–LeClair (1991), Letzter (1999). Abstracts, metadata and later papers are not substituted for their full text. Historical priority for the universal “only later” bulk-particle/soliton chronology, and a full gauge match to Gandenberger's general-rank 1999 matrices, also remain unclear.
- No direct two-loop field-theory calculation, all-rank Kim perturbation calculation, or independent physical proof of the Neumann/coideal dictionary was carried out. Equality of first-order coefficients and differences at second order have the narrower meaning stated in §§1,7,11.
- The author scripts were **read, not executed**; source figures were inspected in P1 without regeneration. No rendered HTML/PDF audit was performed. No Git commit identity can be independently recovered from this export, which lacks `.git`.

## 16. Other observations — unverified

These are editorial or mathematical leads, not additional wrong verdicts:

- Several headings still type chapter/statement numbers into link labels. The target audit above does not ensure those labels survive future renumbering.
- The precise relation between the antilinear classical double row and a universal quantum construction deserves an explicit account of real forms and completions; this report verifies the stated classical identities, not that entire correspondence.
- An expanded teaching treatment could distinguish matrix equations, scalar analyticity, boundary identification and state existence in separate derivations. That is a recommendation about exposition, not evidence for an extra mathematical error.

Validation: both independent phase-2 scripts completed successfully after the final additions. `python3 check-refs.py lecture-notes/half-line` passed all seven structural checks (8 files, 102 IDs, 351 bibliography keys). The frozen P1 file still has the SHA-256 recorded in this report header. These checks do not certify the unresolved physical or historical claims.

Phase 2 stops with this report. All proposed source changes remain recommendations for the editor to verify and integrate; the frozen phase-1 report is unchanged.

---

*Integrated 2026-10-25: see the `## Integration` section of `referee_lectures_halfline_combined.md`.*
