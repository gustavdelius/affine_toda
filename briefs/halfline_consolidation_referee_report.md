<!-- Referee report on the half-line consolidation round (brief: briefs/halfline_consolidation_referee_prompt.md). Written 2026-10-13 by Claude (Opus 5.5) in place of GPT6-Astra, which was out of quota; the user asked for this. Independent checks only: scripts in briefs/halfline_consolidation_referee_scripts/. No book file was edited, git was used read-only. -->

# Referee report: half-line consolidation round (Part VII, 2026-10-13)

**Integration, 2026-10-14.** This historical report has been independently
re-derived for PR #14; see [the integration record](halfline_consolidation_integration.md)
for every item's disposition. In particular, B2's extra $d_5$ pole is verified,
but its proposed spinor-channel explanation is rejected: the relevant channels
decouple. B4's finite-distance ODE values do not independently determine the
singular $\chi$-path limit. The book retains only the qualified results recorded
there.

**Branch** `half-line-consolidation`. **Referee** Claude Opus 5.5, standing in for GPT6-Astra. This breaks the house rule that Claude drafts and another model family referees. The user asked for it because Codex was out of quota. Treat the verdicts as a strong cross-check, not as model-diverse review.

**Method.**
- Every numerical claim was recomputed with my own code, at points I chose: complex rapidities, several couplings.
- The book's library modules were used only where the brief allows: `boundary_common` for representations, the coproduct and the dense intertwiner solver, and `smatrix_simply_laced` for the S-matrix blocks. I wrote my own weight-filtered R-check solver for $d_n$.
- None of the book's check scripts was used as evidence.
- I read one of the book's check scripts, `boundary_krein.py` check 7, for B5, and only after my own computation had disagreed with the text. I did not read the coordinator's scripts. The agent reports (A3, A5, A6) were read as records.
- The literature was read in the arXiv TeX sources:
  - hep-th/9501098, 9612138, 9809140, 9505138, 9506031, 0308004;
  - hep-th/9806003, 9404108, 9909174, 9503083.
- Equation numbers were counted in the source.

## Summary

| Item | Verdict | Label OK? | Action needed |
|---|---|---|---|
| A1 `prp-no-dirichlet` | CONFIRMED WITH CHANGES | yes | Qualify "not integrable at all" (lead-in and open problem 7); tighten two proof steps; "the $C_j$ cannot be sent to infinity" is false for some non-simply-laced nodes |
| A2 `prp-soliton-preserving-half` | CONFIRMED (theorem); the "Half-solitons at the wall" paragraph is **WRONG as stated** for $n\ge3$ | yes | The half-soliton picture needs a single soliton carrying a weight of the class. For $a_3$ two of the seven nontrivial classes have none, for $a_4$ five of 15, for $a_5$ 21 of 31. Fix ch. 30, ch. 31 (third soliton) and `tbl-boundary-types` |
| A3 App. D, Delius 1998b (1.4)–(1.5) | CONFIRMED (a fair statement of an error) | — | Inherit the A2 qualification in "the static solutions at the wall are half-solitons" |
| A4 `prp-dn-uniform-vacuum` | CONFIRMED | yes | None needed. Optional: one sentence on the function space; open problem 2 should say "closed forms", since existence and uniqueness hold for every $n$ |
| A5 attributions (a)–(d) | CONFIRMED | — | (c): Kim's figure 7 shows $d_2$–$d_6$, so $d_7$ is his stated rule |
| B1 `prp-dn-nondiagonal-scalar`, `eq-dn-neumann-amplitude` | CONFIRMED WITH CHANGES | yes | **Ξ is not fixed by one loop.** Three other CDD factors satisfy every axiom and reproduce Kim–Cho at $O(\beta^2)$. Ξ is singled out only by exact duality. Fix ch. 33, ch. 35 and the "same boundaries" sentence; one minor imprecision about $\rho_\varsigma$ |
| B2 `prp-dn-uniform-classical` | CONFIRMED WITH CHANGES | yes | $K^{\rm cl}_a=\prod_bS_{ab}(\theta,1)$, diagonal: confirmed for $d_4,d_5,d_6$ (own ODE); no bound state: confirmed for $d_4,d_5$ ($d_6$ not scanned). **The pole statement is false for $d_5$**: particle 3 has a second physical-strip pole, not at $p=\frac i2m_2$ |
| B3 `prp-sp-scalar`, `eq-sp-breather` | CONFIRMED | yes | none |
| B4 `prp-solitonic-linear` and the endpoint bullets | CONFIRMED (proposition and middle-soliton bullet, by a different method); stationary-pair bullet partly confirmed | yes | none; see §B4 for what I could not reproduce |
| B5 boundary Krein, two and three solitons | CONFIRMED WITH CHANGES | yes | The four thresholds and the mirror convention are right. For $V\bar V$ the breaking appears first at slow rapidity about 1.5–4, **not** for $\theta_2\to0$; fix the sentence |
| B6 $a_n$ $s=-1$ | CONFIRMED | yes | none |

Other problems (§O): a missing credit to Castro-Alvaredo and Fring; "equals Kim's amplitude" should be "agrees at one loop"; the $a_3$ CDD statement (outside this round's diff) has the same one-loop loophole as Ξ.

---

## A1. `prp-no-dirichlet` — CONFIRMED WITH CHANGES

**The proof, step by step (re-derived by hand).**
- **Dirichlet forces the Cartan condition.** With $\vec\phi(0,t)=\vec\phi_0$, $\partial_x\vec\phi(0)$ is a free dynamical variable. The compatibility condition is affine in it, so the Cartan and step parts vanish separately.
- **Transposition reduces to the antipode form.** For transposition, $\rho(X)^{\rm T}=-\rho^*(X)$ turns the condition into $K\rho(a_t(\lambda))=\rho^*(a_t(1/\lambda))K$.
- **The node coefficients drop out.** The weight shifts $\pm\vec\alpha_j$, $j=0,\dots,r$, are distinct exactly when $\theta$ is not simple, i.e. rank $\ge2$. So each node separates, and $c_je^{\beta\vec\alpha_j\cdot\vec\phi_0/2}$ drops out. This leaves $K\rho(E_{\pm j})=\lambda^{\mp2}\rho'(E_{\pm j})K$.
- **The cycle condition.** $E_{\vec\alpha_0}\propto E_{-\theta}$ is a nested commutator of $h-1$ generators $E_{-\alpha_j}$, $j\ge1$, since ${\rm ht}\,\theta=h-1$ for every simple $g$. Hence $(\lambda^{2h}-1)\rho'(E_{-\theta})K=0$.
- **Why the image must be trivial.** For generic $\lambda$, $\rho'(E_{-\theta})$ kills ${\rm Im}\,K$, which is $\rho'$-invariant. The annihilator of an invariant subspace is an ideal of the simple $g$, so it is all of $g$, and ${\rm Im}\,K$ is a trivial summand. The text's "which a representation without trivial summands excludes" skips this ideal argument.
- **Why $a_1$ escapes.** For $a_1$, $\alpha_0=-\alpha_1$, so the node coefficients do not drop out, and the GZ Dirichlet $K$ survives.
- **Antilinear $\sigma$.** With complex $\partial_x\vec\phi=u+iv$, the $u$ and $v$ parts force $KH=HK$ and $KH=-HK$, so $KH=HK=0$. The step terms then have disjoint weight support and give $\sigma(E_{\pm j})K=0$, hence $K=0$. This holds for every antilinear $\sigma$, also one composed with a Weyl twist.

**Scope.**
- **Conjugations.** Absorbing a conjugation into $K_+$ works for any matrix $M$, and $\sigma(X)=MX^{\rm T}M^{-1}$ is legitimate for any $M$ in part (a) of `prp-sklyanin`. The argument also goes through for $\sigma=-\tau$ or $\tau\circ{}^{\rm T}$ with $\tau$ any automorphism of $g$, including outer ones, with $\rho'=\rho\circ\tau$. So the theorem holds for every linear anti-automorphism.
- **Representations.** "Representation" means an evaluation representation. Tensor products of evaluation representations at different points are also excluded, but by an irreducibility argument for loop modules, not by the cycle argument as written. I suggest saying "evaluation representation".

**Numerics** (`a1_no_dirichlet.py`). These are representations the book's script does not use: $b_2$ (5), $c_3$ (6), $g_2$ (7, checked to close on 14 dimensions), ${\rm Sym}^2$ of $a_2$ and $a_3$, and the $d_4$ adjoint (28). Antipode, transposition, Hermitian and $-$conjugation, generic complex $\lambda$: nullity 0 in every case. The control $a_1$ has nullity 1 for linear $\sigma$.

**Spin ±2** (`a1_spin2_flux.py`, sympy, $a_2$).
- I solved for the unique spin-2 current myself.
- The flux of $Q_2+cQ_{-2}$ at a Dirichlet wall has a $\dot u$-term with curl $\propto(c-1)$.
- At $c=1$ the remaining flux is a nonzero cubic plus a linear term in $u=\partial_x\vec\phi(0)$.
- So no combination survives, as stated.

**Loopholes.**
- **A field-dependent $K$** is excluded by the hypothesis, which is stated.
- **The lead-in is stronger than the theorem.** "Dirichlet conditions are not integrable at all" claims more than the theorem, which is for Sklyanin with a field-independent $K$, plus spin 2 for $a_n$. Open problem 7 likewise says "Dirichlet on all components is not classically integrable … So classically the walls are integrable only for sine-Gordon".
- **"For the other theories the $C_j$ cannot be sent to infinity" is false as a general statement.** In `tbl-durham`, $b_n$ has the short end $C_n$ arbitrary, and $c_n$ has $C_0$, $C_n$ arbitrary when the short nodes vanish. What is true is that they cannot all be sent to infinity.
- **The partial-Dirichlet remark is right.** $\vec\phi(x)=M\vec\phi(-x)$ extends $C^1$ across $x=0$, so it is a weak solution. The $PM$-even densities give one conserved combination at each spin. Note that for $a_2$ no Sklyanin $K$ exists for these conditions (A3 agent); the text does not claim one.

**Prior art.**
- **Corrigan, hep-th/9612138, §4**, has the claim without proof. Quote: "setting the field to a specific value at the boundary will not be compatible with integrability in the sense described. In other words, Dirichlet boundary conditions are not generally permitted", and "the Neumann condition is possible in all cases (but not Dirichlet)". So "Corrigan remarked this without proof" is accurate.
- **BCDR hep-th/9501098 and CDRS hep-th/9404108** work with $\partial_1\phi=-\partial\mathcal B/\partial\phi$ and do not discuss Dirichlet.
- **Delius hep-th/9809140** calls (1.5) "the integrable generalization of the Dirichlet condition". Nothing contradicts the theorem.

**Replacement text.**
- ch. 30, lead-in: "For the other theories the $C_j$ cannot all be sent to infinity, and Dirichlet conditions are not integrable, at least not in Sklyanin's sense with a field-independent $K$ (@prp-no-dirichlet); for $a_n^{(1)}$ not even the spin-2 charges survive. Corrigan remarked this without proof [@corrigan1996]."
- Proof: replace "So the image of $K_+$ is invariant, and on it $\rho'$ is equivalent to $\rho\circ\gamma$." with "So $\rho'$ leaves the image of $K_+$ invariant." Replace "unless $E_{-\vec\theta}$ annihilates the image, which a representation without trivial summands excludes" with "unless $E_{-\vec\theta}$ annihilates the image; since $g$ is simple, the annihilator of an invariant subspace is an ideal, so the image would be a trivial summand".
- Open problem 7: "Dirichlet on all components is not classically integrable for $n\ge2$ in Sklyanin's sense, and no spin-2 charge survives (@prp-no-dirichlet). So classically the walls are integrable, in this sense, only for sine-Gordon; …"

## A2. `prp-soliton-preserving-half` — theorem CONFIRMED; "Half-solitons at the wall" WRONG as stated

**The proof (re-derived).**
- **Cartan part.** At imaginary coupling with $\sigma=\dagger$ and real $\lambda$, the Cartan part is $\frac{i\beta}2[K,\partial_x{\rm Re}\vec\phi\cdot H]-\frac\beta2\{K,\partial_x{\rm Im}\vec\phi\cdot H\}$, as stated.
- **Step terms.** The step terms of node $j$ give $w_jKE_j=\bar w_jE_jK$ and the same with $E_{-j}$. Hence $\kappa_j=\kappa_j^{-1}=\bar w_j/w_j=e^{-i\beta\vec\alpha_j\cdot{\rm Re}\vec\phi}$.
- **Every coweight works.** For $K=e^{i\pi\vec\nu\cdot H}$, $\kappa_j=(-1)^{\vec\alpha_j\cdot\vec\nu}$, which matches for ${\rm Re}\vec\phi\in\frac\pi\beta\vec\nu+\frac{2\pi}\beta\Lambda_W^\vee$.
- **Uniqueness for $a_n$.** In the defining representation $K$ is diagonal. If one entry is nonzero, all are, by the cycle, and $\kappa_j=\pm1$.
- **Reflection equation.** I checked the term-by-term argument that every constant diagonal $K$ solves the classical reflection equation, for the matrix units of $r$, $r^{t_1t_2}$, $r^{t_2}$, $r^{t_1}$. It is correct.
- **Is it a legitimate boundary condition?** Yes. It is mixed Dirichlet (Re)/Neumann (Im), variational with no boundary term. As for $K_+=1$, only the real part of the holomorphic energy is conserved, $\frac d{dt}{\rm Re}E=0$; this was already the case before this round.

**Numerics** (`a2_half_coweights.py`).
- Compatibility of $e^{i\pi\vec\nu\cdot\vec H}$ for all $2^n$ classes, $n=2,3,4$ (defining), in $\Lambda^2$ of $a_3$, and in the $d_4$ vector: residual $\le6\times10^{-15}$.
- Solution space of a general $K$ (defining representation, $n=2,3$): dimension 1 on the half lattice, 0 off it.
- **Three-soliton images.** Here $\xi_2=-\bar\xi_1-\ln A_{aa}(2\theta)-\ln A_{ab}(\theta)$ and ${\rm Re}\,\xi_3=-\ln A_{ab}(\theta)$.
  - $\tau_j(x)/\overline{\tau_j(-x)}=c_jF(x,t)$ with $|c_j|=1$ to $10^{-15}$, so $\vec\phi(x)+\overline{\vec\phi(-x)}$ is constant. Cases: $a_2$ ($a,b$)=(1,1),(1,2); $a_3$ (1,1),(2,1),(1,2); $a_4$ (2,1), at regular phases found by scan.
  - $\beta{\rm Re}\vec\phi(0)/\pi$ is a coweight, constant in $t$, not in $2\Lambda_W^\vee$.
  - Without the $-\ln A_{ab}$ corrections the image fails.

**The error.** The paragraph says: "In the other classes … The static solutions are halves of static solitons. … Each class contains a weight of the representation $b$ and one of $h-b$. So there are two such half-solitons of equal mass".
- **Two premises check out.**
  - A static soliton with ${\rm Re}\,\xi=0$ is self-mirror, and ${\rm Re}\,E(x\le0)=M_b/2$, by $e(x)=\overline{e(-x)}$.
  - Each class mod $2\Lambda_W^\vee$ contains exactly one complementary pair $\{\epsilon_S,\epsilon_{S^c}\}$. This follows from the parity vector $([k\in S]+[k+1\in S])\bmod2$.
- **But the inference needs that weight to be a topological charge of a single soliton, and often it is not** (the missing-charge problem):

  | algebra | nontrivial classes | classes with no single-soliton charge |
  |---|---|---|
  | $a_2$ | 3 | 0 |
  | $a_3$ | 7 | 2: those of $\epsilon_1+\epsilon_2$ and $\epsilon_1+\epsilon_4$, i.e. the alternating patterns $(\pm,\mp,\pm,\mp)$ |
  | $a_4$ | 15 | 5 |
  | $a_5$ | 31 | 21 |

  (Check 4 of `a2_half_coweights.py`. A species-$b$ soliton reaches only $h/\gcd(b,h)$ charges.)
- **Consequences.** For those classes there is no static half-soliton of this kind. A three-soliton image with a self-mirror single soliton cannot be built either.
- **What is known for them.** Static two-soliton configurations could provide halves; this is the same structure as `tbl-missing-vacua`.
- **"The static solutions are"** also asserts more than was shown: "are given by", not "are all".

**Replacement text** (ch. 30, "Half-solitons at the wall"; replace from "The static solutions are halves of static solitons." to "differ by the charge of a soliton."):
> When a single soliton carries a weight of the class, halves of static solitons satisfy the condition. A soliton of species $b$ at rest with $\operatorname{Re}\xi=0$ satisfies $\vec\phi(x)+\overline{\vec\phi(-x)}=2\operatorname{Re}\vec\phi(0)$. Here $\frac\beta\pi\operatorname{Re}\vec\phi(0)$ is its topological charge, a weight of the $b$-th fundamental representation, and the real part of its energy on $x\le0$ is half the soliton mass. Each class contains a weight of $\Lambda_b$ and its complement in $\Lambda_{h-b}$. When these are soliton charges, there are two such half-solitons of equal mass, whose vacua at $x\to-\infty$ differ by the charge of a soliton. This holds for every class of $a_2^{(1)}$. It fails for the classes whose weights are missing charges (@sec-topological-charges): for $a_3^{(1)}$ the two classes of $\epsilon_1+\epsilon_2$ and $\epsilon_1+\epsilon_4$, for $a_4^{(1)}$ five classes of fifteen.

Also qualify:
- ch. 31: "For the boundaries of @prp-soliton-preserving-half with $\vec\nu\notin2\Lambda_W^\vee$ whose class contains a soliton charge, the image solution contains a third, self-mirror soliton …".
- `tbl-boundary-types`: "off the vacua (@prp-soliton-preserving-half) a static half-soliton where a soliton carries the class".
- The "Three soliton image" sentence in ch. 30: the same qualification.

**Prior art.** None found in Delius 1998b, BCDR or Corrigan.

## A3. Appendix D, Delius 1998b — CONFIRMED

Read in hep-th/9809140:
- **(1.4)**: $\phi^*(-x)=\phi(x)+2\pi i\lambda$, "$\lambda\in$ coweight lattice", justified by "adding $2\pi i$ times a coweight to a soliton solution gives an equivalent soliton solution". That shift maps $\lambda\to\lambda-2\mu$, so it produces only $2\Lambda$ from $\lambda=0$.
- **(1.5)**: ${\rm Im}\phi/2\pi\in$ coweight lattice, $\partial_x{\rm Re}\phi=0$. "The boundary condition corresponding to (1.4) is (1.5)."
- At $x=0$, (1.4) gives ${\rm Im}\phi(0)=-\pi\lambda$, which is (1.5) only for $\lambda\in2\Lambda$.
- **§5, (5.1)–(5.2)**: "If one imposes our new boundary condition (1.5) then these equations are solved by $K_+=1$."
- For $K=e^{i\pi\lambda\cdot H}$ I checked by hand that (5.1) and (5.2) hold separately. With ${\rm Im}\,\phi(0)=-\pi\lambda$, $e_j=e^{\alpha_j\cdot\phi/2}$ is real where $\alpha_j\cdot\lambda$ is even, where $\kappa_j=+1$ and the commutator in (5.1) vanishes. It is imaginary where $\alpha_j\cdot\lambda$ is odd, where $\kappa_j=-1$ and the anticommutator in (5.2) vanishes. For complex $K$ only (5.1)$+i$(5.2) is the actual condition; it holds too.

**The entry is a fair statement of an error.** One change: "the static solutions at the wall are half-solitons" → "static solutions at the wall can be halves of solitons" (see A2).

## A4. `prp-dn-uniform-vacuum` — CONFIRMED

**Convexity.**
- $\int\frac12\vec\phi'^2$ is convex.
- $V=\sum n_j(e^{\vec\alpha_j\cdot\vec\phi}-1)$ has Hessian $\sum n_j\vec\alpha_j\vec\alpha_j^{\rm T}e^{\cdots}$, which is positive definite because the $\vec\alpha_j$ span.
- $\mathcal B$ is convex for $C_j\ge0$.
- So the energy is strictly convex on the convex set of finite-energy fields.

**Critical point ⇒ global minimum.**
- Finite energy forces $\vec\phi\to0$ and $\vec\phi\in L^2$ (since $V\simeq\frac12\vec\phi^{\rm T}M^2\vec\phi$ near 0).
- So $E(\phi+t\eta)$ is a finite convex function of $t\in[0,1]$, and its one-sided derivative is $\phi'\eta|_{-\infty}^0+\nabla\mathcal B\cdot\eta(0)+\int(-\phi''+V')\eta=0$.
- Dominated convergence is justified by $\int|\eta||V'(\phi)|\le\|\eta\|_2\|\phi\|_2\cdot{\rm const}$.
- The argument is legitimate.

**Invariance.** Affine diagram automorphisms of $d_n$ preserve the marks and act orthogonally, so they map the $C\equiv+1$ problem to itself. The invariant subspace is one-dimensional for $d_4$ ($\vec\alpha_2$) and for $d_5$ ($\vec\alpha_2+\vec\alpha_3$).

**Closed forms (by hand).**
- The exponents $-f$/$2f$ ($d_4$) and $-f$/$f$ ($d_5$) are right.
- The reduced equations are right: $f''=2(e^{2f}-e^{-f})$ and $f''=4\sinh f$.
- So is $u'=(u-1)\sqrt{2u(u+2)}$.
- The boundary conditions give $f'(0)=2e^{-f/2}-\sqrt2e^f$ and $f'(0)=2e^{-f/2}-\sqrt2e^{f/2}$. With the first integrals, $u(0)=9/8$ and $e^{f(0)}=4-2\sqrt2$, and $f'(0)>0$ in both cases.

**Energies** (`a4_dn_vacuum.py`).
- By mpmath quadrature: $-0.0695214587295478$ and $-0.0927092091312040$. These equal $16\sqrt2-8-6\sqrt6$ and $32(1+\sqrt2)\sin\frac\pi8-24-4\sqrt2$ to $10^{-25}$; the $d_5$ form I also reduced by hand.
- A full-field boundary-value solve (scipy `solve_bvp`, no symmetry assumed, three random starts) gives the same energies and boundary values to $10^{-9}$. For $d_6$ it gives $-0.103251515$, against the book's $-0.1032515$.

**Literature.**
- Not found in BCDR, CDRS, Fujii–Sasaki hep-th/9503083, Corrigan's review, or Bowcock–Perkins.
- BP §3 prove only that the energy is bounded below for $A_i\le1$ (their $A_i=-C_i$), by a Bogomolny bound, not uniqueness.
- The argument is elementary, and no citation is required. BP §3 could be cited for the stability context.

**Optional.**
- Existence also holds for every $n$, by the direct method: $E$ is convex, bounded below and coercive in $H^1_{\rm loc}$. So open problem 2's "the vacua for $n\ge7$" should read "closed forms of the vacua for $n\ge7$ (existence and uniqueness follow from convexity)".

## A5. Literature attributions — CONFIRMED

**(a) Kim–Cho, hep-th/9505138.**
- Counting the equations in the source, (29) is $K_1=1+\frac{i\beta^2}{24}\big(\frac{\sh}{\ch}-\frac{\sh}{\ch-1}+\frac{2\sh}{2\ch+\sqrt3}-\frac{2\sh}{2\ch-1}\big)$, and (40) is the one-loop $K_2$, with $B=\frac{\beta^2/2\pi}{1+\beta^2/4\pi}$.
- Kim hep-th/9506031 refers to Kim–Cho for $d_4$. His own one-loop calculation is $a_3$ only (appendix).
- The correction is right, and no stale `kim1995` $d_4$ citation remains.

**(b) Kim's $d_n$ vector conjecture.**
- In Kim's labels the vector is $K_{n-2}$, the end of the long arm.
- Figure 7: "only the middle number is two". So $K=[\frac{h-1}2]\prod_{k\ {\rm odd}\le2h-3}[\frac k2]$. Kim's $[x]$ is CAF's $\|2x\|$.
- This equals CAF's $\tilde R_1(\theta+i\pi)=\|h-1\|\prod_{k=1}^{h-1}\|2k-1\|$ (their (dn1)). CAF state it themselves: "(KimCF) $R^K(\theta,B)=\tilde R(\theta+i\pi,B)$", and "For higher ranks only a solution for $j=1$ was also reported in [Kim], which once again coincides with ours subject to (KimCF)".
- My numbers agree to $10^{-40}$, $n=4..7$ (`a5_one_loop.py`).
- Figure 7 prints the rows $d_2$–$d_6$, so "$d_7$" uses Kim's stated rule.

**(c) Appendix D sentence.**
- At complex $\theta$, $K^N_1$ and $K^N_2=K^N_1(\theta\pm\frac{i\pi}6)S_{11}(2\theta)$ ($S_{11}=\{1\}\{5\}$) reproduce Kim–Cho (29), (40) at $O(\beta^2)$ to $10^{-14}$.
- Kim–Cho's exact $K_1=[\frac12][\frac32][\frac52]^2[\frac72][\frac92]$ and Kim's vector rule agree with $K^N_1$ at $O(B)$ and differ at $O(B^2)$, $n=4..7$.
- "Which they call uniquely determined by duality" is fair: "The exact boundary reflection matrix is determined uniquely if we assume the strong-weak coupling 'duality'."

**(d) Gandenberger and CDRS.**
- **Gandenberger hep-th/9806003.** By the equation count, (5.16)–(5.18) are $K_B^{(+)},K_B^{(-)},K_B^{(d)}$, and (6.4) their continuation.
  - (5.17) $=-\langle\frac\lambda2\rangle\langle-\frac\lambda2-\frac12\rangle\langle-\frac32\lambda+\frac12\rangle$, with his $\langle a\rangle=\sin\frac\pi{3\lambda}(\mu+a)/\sin\frac\pi{3\lambda}(\mu-a)$. This is `eq-an-minus-breather` term by term.
- **CDRS hep-th/9404108.** By the equation count, (5.7) is the "minimal hypothesis" $-(1)(2+\frac B2)/(\frac B2)$ for $A_1=A_2=-2$, i.e. $C\equiv-1$.
  - (5.16) is the general $a_n$ expression $K^0_a=\prod_k\langle k-\frac12\rangle/\langle h-k+\frac12\rangle$, "following the bootstrap".
  - Their summary calls the picture "tentative".
  - So "conjectured for $C\equiv-1$, for general $a_n$" is a fair description.

## B1. `prp-dn-nondiagonal-scalar`, `eq-dn-neumann-amplitude` — CONFIRMED WITH CHANGES

**(a), (b): matrix level** (`b1_dn_matrices.py`, $n=4,5$, $\omega=0.37,0.71$). $K$ was built from the book's coideal with the library solver; normalizations, crossing vector and breather are my own.
- **Corner entry.** $K_{1,2n}/K_{11}=O(\delta)$ as $y^h\to\varsigma$, and $O(1)$ as $y^h\to-\varsigma$, for both $s$. So the data depend on $\varsigma=(-1)^ns$. ($y^h=\pm1$ exactly is degenerate, $x=1/x$, with a two-dimensional intertwiner space.)
- **Unitarity.** $K(-\theta)K(\theta)=U_\varsigma\mathbf1$ to $3\times10^{-14}$.
- **Crossing.** $k(\theta)=\Psi\check R_{\bar V\bar V}(2\theta)k(-\theta)$, with $\bar V$ the Chevalley-twisted representation and the $V\to\bar V$ intertwiner. $\Psi=\Psi_\varsigma$ to $4\times10^{-14}$. Control: $\hat\epsilon=0$ gives $\Psi=1$.
- **Singlet coefficient.** The $B_1$ singlet is mapped to $c^d_\varsigma$ times the singlet, to $6\times10^{-14}$. Control: $\hat\epsilon=0$ gives 1. The construction was calibrated on $a_2$, $a_3$, where it reproduces the established $c_\pm$ exactly.

**(c), (d): functional equations and breathers** (`b1_dn_scalar.py`).
- mpmath at 30 digits; the infinite products are summed as series of log-Gamma differences with `nsum` (Richardson for the $1/k^2$ tails, Shanks for the alternating $g_\varsigma$).
- $n=4,5,6$, $\omega=0.37,0.83$, two complex $\mu$.
- Residual $\le9\times10^{-30}$ for both unitarity equations, both crossing equations, both breather identities (with $c^d$ and Ξ as stated), and the $\hat\epsilon=0$ control (`eq-dn-scalar`, `eq-dn-breather`).

**(e): continued amplitudes**, at $\omega=-2/B$, $n=4..7$, $B=0.31,1.37,0.6+0.2i$.
- The closed forms continue to $K_1|_{B\to2-B}$, $K_1/K_1^{\rm cl}$ and $K_1$ to $10^{-30}$.
- Unitarity, crossing-unitarity with $S_{11}$ (BCDS $d_n$ blocks) and bootstrap closure (the $(1,a)$ and $(2,a-2)$ routes agree) hold to $10^{-29}$ for $K^N$, $K^{(-)}$ and the minimal (no-Ξ) amplitude.
- All tend to 1 as $B\to0$.

**(f): one loop** ($d_4$, `b1_dn_scalar.py` part B, `a5_one_loop.py`).
- With Ξ: $K(0)=1$ and $d\log K/dB=$ Kim–Cho to $10^{-14}$.
- Minimal: off by 3.45.
- $K^{(-)}_1(\theta)$ misses, $K^{(-)}_1(\theta+i\pi)$ matches.

**(g): assessment, and the finding.** One loop does not fix Ξ.
- Write $P_x=(\frac T2-x)_I(\frac T2+x)_I$, so $\Xi=P_0P_{\omega/2}P_{(T-\omega)/2}$.
- Replacing $P_0\to P_{1/2}$, or $P_{\omega/2}\to P_{\frac12+\frac\omega2}$, or both, gives amplitudes that:
  - satisfy unitarity, crossing-unitarity and bootstrap closure ($n=4,5,6$, $10^{-30}$);
  - reproduce Kim–Cho's $O(\beta^0)$ and $O(\beta^2)$ exactly ($2\times10^{-14}$);
  - differ from $K_1|_{B\to2-B}$ by 4–27% at $B=0.6$.
  (`b1_cdd_alternatives.py`.)
- So Ξ is singled out only by the requirement that $K^N_1$ be exactly the dual of $K_1$ (the $d_n$ analogue of DG's (3.20) conjecture). The A6 agent noted this ambiguity in its report, but it did not reach the book.
- There is no soliton-bootstrap constraint for the $d_n$ vector alone, because $V\otimes V$ contains no $V$. The classical limit cannot decide, since all candidates tend to 1.
- I found no derivation of Ξ.

**Replacement text.**
- ch. 33, after (v): "The CDD factor Ξ is unitary and crossing-symmetric by itself, so the functional equations do not fix it. Kim and Cho's one-loop $d_4^{(1)}$ amplitudes require a CDD factor, and Ξ reproduces them (@sec-dn-uniform), but one loop does not single it out: replacing $(\frac T2)_I^2$ by $(\frac T2-\frac12)_I(\frac T2+\frac12)_I$, or $(\frac T2\mp\frac\omega2)_I$ by $(\frac T2\mp\frac{\omega+1}2)_I$, gives amplitudes that satisfy every axiom and agree with Kim and Cho to $O(\beta^2)$. Ξ is the choice for which the continued amplitude is exactly the dual of the uniform one."
- ch. 33, "So $d_n$ has the same three solutions as $a_n$, and they describe the same boundaries." → "…, and, given @cnj-reflection-continuation and the choice of Ξ, they describe the same boundaries."
- ch. 35: "given @cnj-reflection-continuation and the CDD factor Ξ, which is fixed by one loop and not derived" → "…, which one loop does not fix uniquely and which is not derived. What singles out Ξ is that it makes $K^N_1$ exactly the dual of $K_1$ for every $n$."
- **Minor (ch. 33).** "The functions $\rho_\varsigma$ are the $a_n$ factors of @prp-an-nondiagonal-scalar without their CDD factor" is exact for $\varsigma=-1$ ($\rho_-=A^{(-)}/A$). For $\varsigma=+1$, $\rho_+=(-\frac T2)_I^2\,A^{(+)}/(A\cdot{\rm CDD})$. The extra $(-\frac T2)_I^2$ is itself a CDD factor: $(-\frac T2)_I$ is unitary with crossing ratio $-1$. Say so, or drop "are".
- **Labels.** [Numerical] is right for all three propositions.

## B2. `prp-dn-uniform-classical` — CONFIRMED WITH CHANGES (one pole statement is wrong for $d_5$)

**Method** (`b2_dn_classical.py`, `b2_quick.py`, `b2_quick6.py`).
- **Vacua.** $d_4$ and $d_5$ from the closed forms of `prp-dn-uniform-vacuum`, which I checked against the field equation and the boundary condition. $d_6$ from a full-field boundary-value solve with no symmetry assumed.
- **Integration.** The linearized equation $\varepsilon''=(M^2(x)-\Omega^2)\varepsilon$ is integrated forward from $x=-18$. It starts from the exact asymptotic modes: in- and outgoing waves in open channels, decaying exponentials in closed ones. So no channel is integrated against its growth.
- **Boundary.** $\varepsilon'(0)=-R\varepsilon(0)$, $R={\rm Hess}\,\mathcal B(\psi(0))$.
- **Comparison.** $\prod_bS_{ab}(\theta_a,1)$ from the book's block content (`smatrix_simply_laced`), particles matched by mass.

**Results.**
- **The reflection is diagonal and equals $\prod_bS_{ab}(\theta,1)$** in every open channel, at frequencies from the lowest gap to above all thresholds:

  | algebra | deviation | off-diagonal | particles |
  |---|---|---|---|
  | $d_4$ | $2.3\times10^{-11}$ | $1\times10^{-15}$ | the degenerate triplet |
  | $d_5$ | $2.7\times10^{-11}$ | $4\times10^{-14}$ | spinors included |
  | $d_6$ | $1.8\times10^{-11}$ | $1\times10^{-13}$ | spinors included |

- **$B\to0$ limit.** For the vector particle this is the $B\to0$ limit of `eq-dn-uniform-amplitude`: $-(\frac h2-1)(\frac h2+1)/((1)(h-1)(\frac h2)^2)$.
- **No bound state** (`b2_bound_fast.py`, `b2_bound_d5_part.py`). The bound-state determinant $D(\Omega^2)=\det[\varepsilon'(0)+R\varepsilon(0)]$ was taken over the solutions that decay at $-\infty$.
  - $d_4$: 45 points on $-4<\Omega^2<m_{\min}^2$, no sign change, $\min|D|=1.36$.
  - $d_5$: 22 points on $-1.4\le\Omega^2\le1.1706$ ($m_{\min}^2=1.1716$). $D$ is positive throughout and falls to 0.04 only at threshold. At $\Omega^2=1$, $D=1.84$, so the $d_5$ pole below is not a bound state.
  - $\Omega^2\le0$ is excluded by the strict convexity of `prp-dn-uniform-vacuum`, since the linearized operator with its boundary condition is the Hessian of a strictly convex functional.
  - I did not scan $d_6$.
- **Poles**, from the block content and checked by direct evaluation near each pole.
  - $d_4$: every physical-strip pole sits at $p=\frac i2m_2$, i.e. $2|p|=\sqrt6$.
  - $d_6$: at $2|p|=m_2,\ m_4,\ m_2+m_4$ (the last ones double), as stated.
  - **$d_5$: the statement is false for particle 3.** $\prod_bS_{3b}(\theta,1)$ has, besides its pole at $p=\frac i2m_2$ ($\theta=\frac{i\pi}8$), a simple pole at $\theta=\frac{3i\pi}8$ with a small residue ($|\prod S|\approx8.1\times10^{-3}/\epsilon$). There $p_3=i\,m_3\sin\frac{3\pi}8=i(1+\sqrt2)$, which is not $\frac i2m_2=i$. The frequency is $\Omega=m_3\cos\frac{3\pi}8=1$, exactly the frequency at which the spinor channels have their own pole $p_s=\frac i2m_2$ ($\Omega^2=m_s^2-1=1$). So it is a pole that particle 3 inherits from the spinor channel's tail pole through the multichannel coupling, not a pole at $p_3=\frac i2m$ for some tail rate $m$.
  - (`b2_poles_d5.py`, `b2_pole_z3.py`. The A6 report's "all physical-strip poles … for $d_4$ and $d_5$ sit at $p=\frac i2m_2$" missed this small-residue pole.)

**Replacement text** (ch. 35, after `prp-dn-uniform-classical`): replace "For $d_4^{(1)}$ and $d_5^{(1)}$ every physical-strip pole of $\prod_bS_{ab}(\theta,1)$ sits at $p=\frac i2m_2$." with
> For $d_4^{(1)}$ every physical-strip pole of $\prod_bS_{ab}(\theta,1)$ sits at $p=\frac i2m_2$. For $d_5^{(1)}$ so do all but one: particle 3 also has a pole at $\theta=\frac{3i\pi}8$, at the frequency $\Omega=m$ at which the spinor channels sit at $p=\frac i2m_2$; it is inherited from that channel through the coupling of the linearized equations.

The label [Numerical] is right.

## B3. `prp-sp-scalar`, `eq-sp-scalar-eqs`, `eq-sp-breather` — CONFIRMED

(`b3_sp_scalar.py`; $\omega=0.41$, $\eta=0.3+0.7i$.)
- **Sanity of the $K$ convention.** The de Vega–González-Ruiz $K_V$ solves the braid reflection equation ($10^{-15}$). $\check R_{V\bar V}$ is 1 on every transmission ($10^{-15}$).
- **Crossing functions**, $(n,l)=(2,1),(3,2),(3,0),(2,0)$. The matrix parts give $\Phi_1=-(-q)^{1-l}s(\mu)$ and $\Phi_2=-(-q)^{l-1}s(\mu)$ to $2\times10^{-15}$, and they are independent of $\eta$.
- **Scalar factors.** `eq-sp-scalar` solves all four equations at complex $\mu$, $(n,l)=(2,1),(3,2),(3,1)$, to $3\times10^{-25}$.
- **Breather.**
  - The singlet in $V(\theta_a)\otimes\bar V(\theta_b)$ is mapped by $(1\otimes aK_V)F(2\mu)\check R(1\otimes\tilde bK_{\bar V})$ to the singlet in $\bar V(-\theta_b)\otimes V(-\theta_a)$.
  - The coefficient is $e^{-2\pi il/h}K_0E_\eta$ to $4\times10^{-15}$, for $(n,l)=(2,1),(3,2),(2,0),(3,1)$, with the normalizations stated in ch. 35. The phase is included.
- **Continuation of $K_0$.** Checked by hand: $(x)_I\to(2x/\omega)$, $\omega\to-2/B$.
- I did not check the soliton bootstrap, the $n=1$ GZ identification or the $B\to0$ statement about $L_1$.

## B4. `prp-solitonic-linear` — CONFIRMED (different method); endpoint bullets partly confirmed

The book's script and the A5 agent use Bowcock's exact Hirota dressing. I integrated the linearized equation directly (`b4_solitonic_linear.py`, DOP853), starting at $x=-16$ from the exact asymptotic modes, with the boundary Hessian at $x=0$.

**The proposition.**
- Coverage:
  - every regular continued vacuum: $a_2$ ($a=1$), $a_3$ ($a=1$), $a_4$ ($a=1,2$), $a_5$ ($a=1,2$), three $\chi$-intervals each, with the pattern recovered from the boundary condition;
  - the self-conjugate family of $a_3$ and $a_5$ ($d=0.3, 0.8$).
- Frequencies: $\Omega=m_1\cosh0.83$, 1.9, 2.7.
- $K$ is diagonal (off-diagonal $\le4\times10^{-8}$, limited by the ODE) and $K_{cc}=K^{(b)}_c(B=0)$ to $5\times10^{-12}$, with $b=a$, and $b=h/2$ for $C\equiv-1$.

**Middle-soliton endpoint** ($a_3$, Hessian of $(+,-,+,-)$). $|K_{11}+(1)(3)|$ and the off-diagonal entries are $5.0\times10^{-3}$ and $5.0\times10^{-4}$ at $1-d=10^{-2},10^{-3}$, linear in $1-d$, as stated. At $10^{-4}$ the ODE floor ($2\times10^{-4}$) takes over.

**Stationary-pair endpoint** ($a_5$, $(+,+,-,-,-,-)$, $\Omega=m_1\cosh0.6$).
- Along the $d$-path, $|K_{15}|=0.9278,\ 0.9365,\ 0.9367$ at $e=10^{-2},10^{-3},10^{-4}$. This is consistent with the book's 0.937 (the A5 agent's 0.93744) within my ODE's floor of about $10^{-3}$ near the wall layer.
- Along the $\chi$-path my ODE loses accuracy in the wall layer: $|K|>1$ at $\delta=10^{-3}$. At $\delta=10^{-2}$ it gives 0.86 and 0.58, already far from DG's diagonal value.
- So the path dependence is confirmed: the $d$-path limit differs from both the diagonal value and the $\chi$-path values at $\delta=10^{-2}$. The book's number 0.719 for the $\chi$-path limit is **not independently reproduced** by my method.
- I did not check the exceptional $a_4$ background.

**"Pair endpoints are path-dependent" and the withdrawal of "BP Table 2 corrected" are justified.**

## B5. Boundary Krein: two and three solitons — CONFIRMED WITH CHANGES

**Construction** (`b5_boundary_krein.py`, my own braid-form double-row transfer matrices).
- Representations and R-check come from the library.
- The walls are U ($K=1$) and N (`eq-an-k-matrix`, $s=+1$), with $K_{\bar V}=CK_VC$ and each matrix part normalized to $K^\dagger CK=C$.
- Every $T_j$ preserves the Krein form ($10^{-14}$) and commutes with $Z^{\otimes N}$ ($5\times10^{-15}$). Of the two candidates $Z\otimes Z$ and $Z\otimes Z^{-1}$ on $V\otimes\bar V$, only $Z\otimes Z$ commutes, so "sector $\zeta=1$" is unambiguous.
- I read the book's check 7 only after my threshold computation disagreed with the text (below). It uses the same construction, so the disagreement is not about conventions.

**Mirror convention.**
- With the stated mirror rule, $[T_i,T_j]=0$ to $2\times10^{-14}$ (VV, $V\bar V$, $VV\bar V$; walls UN, NN, NU).
- With the right wall's matrix on $\bar V$ they fail by $O(1)$ whenever the left wall is N.
- For a U left wall the two rules coincide, because $CKC=K=1$. So "with the right wall's matrix on $\bar V$ instead they do not" is true for a Neumann left wall.

**Thresholds** (bisection, sector $\zeta=1$).
- **VV.** Fast rapidity 10 or 12, slow one $10^{-4}$–$10^{-2}$: UN 0.6982, NN 0.5362, as stated. The breaking does first appear for $\theta_2\to0$.
- **$V\bar V$.** In the same limit $\theta_2\to0$ I get 0.6977 (UN) and 0.5358 (NN), essentially the VV values. But a fine grid at $\omega=0.672$ (UN) and 0.505 (NN) already shows $|\Lambda|\ne1$ at slow rapidity 1–2 (`b5_explore3.py`). With the slow rapidity on $[10^{-4},4]$, both assignments and fast rapidity 8–12, the thresholds are 0.6667 (UN) and 0.5001 (NN), with the first breaking at (4, 12) (`b5_vvbar_threshold.py`). This agrees with the book's numbers. With both rapidities up to 20 (`b5_vvbar_extended.py`):
  - at $\omega=0.660$ (UN) and $0.495$ (NN) the levels stay on the unit circle ($3\times10^{-11}$, $10^{-12}$);
  - at $0.672$ and $0.505$ they leave it, $\max\big||\Lambda|-1\big|=0.049$ at $(20, 1.5)$ and $0.062$ at $(20, 2.0)$.
  So the thresholds hold, and the slow rapidity at which the breaking appears is about 1.5–2 for a fast rapidity of 20.
- **So the numbers are right, but the description is wrong for $V\bar V$.** "Above it $|\Lambda|\ne1$ appears first for $\theta_2\to0$ with $\theta_1$ large; the threshold is that of the limit $\theta_1\to\infty$" holds for VV only. For $V\bar V$ the breaking first appears with one rapidity large and the other of order a few units. The book's script brackets the thresholds correctly ($\theta_2\le3$), but it does not locate them.

**Uniform walls.** UU is unbroken for VV and $V\bar V$ on a $7\times7$ grid $10^{-4}\le\theta\le12$ at $\omega=0.15, 0.45, 0.9, 1.3$. UN at $\omega=0.45$ is unbroken on the same grid ($10^{-14}$).

**Three solitons** (`b5_three.py`). VVV in the sector $\zeta=1$, 45 rapidity triples, $\omega=0.1, 0.2, 0.35, 0.5, 0.9, 1.3$: broken for UU, UN and NN at every coupling, with $\max\big||\Lambda|-1\big|$ between 0.07 and 8.7. A handful of triples is not enough: with three triples NN at $\omega=0.2$ looked unbroken.

**Replacement text** (ch. 35, second bullet):
> With one or two Neumann walls they are phases only below a threshold in the coupling: $\omega\approx0.698$ (UN, VV), $0.666$ (UN, $V\bar V$), $0.536$ (NN, VV) and $0.501$ (NN, $V\bar V$). For VV the breaking appears first for $\theta_2\to0$ with $\theta_1$ large, and the threshold is that of the limit $\theta_1\to\infty$, reached by $\theta_1\approx10$. For $V\bar V$ it appears first with one rapidity large and the other of order a few units; in the limit $\theta_2\to0$ the $V\bar V$ thresholds would be close to the VV ones.

The label [Numerical] is right.

## B6. $a_n$, $s=-1$ — CONFIRMED

(`b6_an_minus.py`, mpmath, 30 digits.)
- **Scalar factor.** $A^{(-)}$ of `eq-an-nondiagonal-scalar-minus` satisfies $A(\mu)A(-\mu)U_-=1$ and $A(\frac T2-\mu)/A(\frac T2+\mu)=F_{11}(2\mu)/\Psi_-$, for $n=1..4$, $\omega=0.43,0.91$, complex $\mu$: $10^{-30}$.
- **Breather identity** `eq-an-minus-breather` with $c_-$: $n=1..4$, $10^{-30}$.
- **Continuation.** The closed form continues to $-(1)(h-\frac B2)(-1+\frac B2)=-(1)(h-1)K_1$.
- **Axioms.** $K^{(-)}_a=-(a)(h-a)K_a$ satisfies unitarity, crossing-unitarity ($S_{a,h-a}$), $K_a=K_{h-a}$ and the bootstrap at every fusing $a+b\to a+b$, for $n=2..7$ and $B=0.37,1.2+0.3i$: $10^{-30}$. It tends to 1 as $B\to0$.
- **Duals.** $-(a)(h-a)K^N_a=$ CDRS (5.16) for every $a$, $n\le7$. For $n=2$ it is CDRS (5.7).

## O. Other problems noticed

1. **Missing credit to Castro-Alvaredo and Fring** (ch. 35, `sec-dn-uniform`, unchanged text next to the new material).
   - The text says CAF constructed $\hat R$ "without a boundary condition in mind". But CAF §5 (after (Rhat)) state $\hat R_i(\theta,B=0)=\prod_jS_{ij}(\theta,B=1)$ and say that $\hat R$ "tends to a theory with fixed boundary conditions for $B\to0$". In their conclusion "fixed" means the $\pm1$ boundaries with the classical limit (l5), not Dirichlet.
   - So the new `prp-dn-uniform-classical` remark "The right-hand side is also CAF's $\hat R_a$ at $B=0$" should credit them for that formula. What is new is the identification with the linearization about the $C\equiv+1$ kink.
   - Suggested: "…one of the generic solutions they constructed. They noted that it tends to $\prod_jS_{1j}(\theta,1)$ as $B\to0$ and associated it with 'fixed' boundary conditions [@castro_alvaredo_fring2004, sec. 5]; @prp-dn-uniform-classical identifies this limit as the uniform boundary."
2. **ch. 35, "$K^{(-)}_1$ of $a_3^{(1)}$ equals Kim's amplitude at $\theta+i\pi$".** It agrees with Kim's one-loop (33) at $\theta+i\pi$ (checked). It differs from his exact self-dual amplitude by 9–16% at $B=0.6$. Say "agrees at one loop with Kim's (33) at $\theta+i\pi$".
3. **ch. 33, the $a_3$ CDD factor** (outside this round's diff, same issue as B1). "For $n=3$ it is required by Kim's one-loop amplitudes": $P_{\frac12+\frac\omega2}$ in place of $P_{\omega/2}$ also reproduces Kim's (33) at $O(\beta^2)$, and it satisfies DG's soliton-bootstrap condition $\sigma(\mu)=\prod_k\sigma(\mu+T-k\omega)$ (`b1_an_cdd_check.py`). Say "a CDD factor is required, and $(\frac T2\mp\frac\omega2)_I$ is the one that gives exact duality".
4. **The semiclassical statement** in ch. 35 ("Neither the classical limit nor the semiclassical soliton time delay distinguishes the two matrices. For $d_4$ both give half the phase of $F(2\theta)$") was not checked by me.
5. **Cross-references.** All `@`-references and citation keys added in this round resolve.

## Scripts

All are in `briefs/halfline_consolidation_referee_scripts/`. Run them from that directory; they find `textbook_code/` by relative path. Each ends with a PASS/FAIL line, except the exploration helpers.

| script | item | runtime |
|---|---|---|
| `a1_no_dirichlet.py`, `a1_spin2_flux.py` | A1 | 1 min, 1 min |
| `a2_half_coweights.py` (`a2_regular_scan.py`: phase scan) | A2 | 2 min (10 min) |
| `a4_dn_vacuum.py` | A4 | 5 s |
| `a5_one_loop.py` | A5, B1(f) | 1 s |
| `b1_dn_matrices.py` (calibration `b1_calib_an.py`) | B1 (a), (b), (iv) | 20 s |
| `b1_dn_scalar.py`, `b1_cdd_alternatives.py`, `b1_an_cdd_check.py` | B1 (c)–(g), §O.3 | 5 min, 2 s, 1 s |
| `b2_dn_classical.py` (library), `b2_quick.py`, `b2_quick6.py`, `b2_poles_d5.py`, `b2_pole_z3.py`, `b2_bound_fast.py`, `b2_bound_d5_part.py` | B2 | 10–40 min each. The full `b2_dn_classical.py` main is too slow; use the pieces |
| `b3_sp_scalar.py` (log `b3_sp_scalar.log`) | B3 | 2 min |
| `b4_solitonic_linear.py` (log `b4_solitonic_linear.log`; `p3` argument for the endpoint part) | B4 | 15 min + 22 min |
| `b5_boundary_krein.py`, `b5_vvbar_threshold.py`, `b5_vvbar_extended.py`, `b5_three.py`, `b5_explore*.py` | B5 | 20–60 min each. `b5_boundary_krein.py` ends with FAIL for the two reasons explained in §B5: its $V\bar V$ bisection runs in the limit $\theta_2\to0$, which gives 0.698/0.536 rather than the book's values, and its three-soliton test uses only three triples. `b5_vvbar_threshold.py`, `b5_vvbar_extended.py` and `b5_three.py` settle both |
| `b6_an_minus.py` | B6 | 3 min |
