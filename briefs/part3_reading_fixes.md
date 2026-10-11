# Brief: fixes and gaps from the Part III reading (Chapters 9–11)

Status: done 2026-10-10 (uncommitted, branch `reading-fixes`; see Resolution at the end). Written 2026-10-10 after a line-by-line reading of Part III by a reader with graduate-level mathematical physics and QFT, who knows the Part I and Part II material as written.

Files: `part3-classical/part.qmd`, `09-lagrangian.qmd`, `10-solitons.qmd`, `11-reality.qmd`.

Scope: fix the items below in place. Items A1–A3 and A5 need a check against the literature or a computation before any change; say in your report what you checked and what you could not confirm. Do not edit anything else without checking with the author. The Part III results I checked by hand are listed in section D.

## A. Errors and likely errors

**A1. Coincident orbits and the degree of the one-soliton $\tau$-function (`10-solitons.qmd`, Section `sec-vertex-solitons`, first bullet).** The text says "a coincident orbit $O$ has degree $\lvert O\rvert n_j^\vee$, e.g. 2 at the long end nodes of $c_n^{(1)}$". In `tbl-foldings` (`06-affine-algebras.qmd`) the long end nodes of $c_n^{(1)}$ have label 1, since they are fixed by the folding reflection, so their orbit size is 1. The orbits of size 2 are the short nodes, which carry label 2. Fix: "e.g. 2 at the short nodes of $c_n^{(1)}$", and check the degree of the species-$a$ soliton at a short node with the explicit $c_2^{(1)}$ example from `06-affine-algebras.qmd`.

**A2. Sign of the field in the real-coupling check (`11-reality.qmd`, Section `sec-theta-classical`, paragraph after the $\Theta$-invariant field discussion).** The text says "$-i\vec\phi$ solves the real-coupling equation". Let $\chi=-i\vec\phi$. The imaginary-coupling equation gives $\partial_+\partial_-\chi=(m^2/\beta)\sum_jn_j\vec\alpha_je^{-\beta\vec\alpha_j\cdot\chi}$. Setting $\psi=-\chi=i\vec\phi$ gives $\partial_+\partial_-\psi=-(m^2/\beta)\sum_jn_j\vec\alpha_je^{\beta\vec\alpha_j\cdot\psi}$, which is the real-coupling equation (`eq-toda-eom`). So the solution of the real-coupling equation is $i\vec\phi$, or $-i\vec\phi$ with the field sign reversed. Fix the sentence, and check the same point in `exr-theta-solutions`.

**A3. Time advance or time delay (`10-solitons.qmd`, Section `sec-multi-solitons`, paragraph "Time delays").** The text says that for real rapidities $0<A_{ab}<1$, "so the displacement is real and positive, a time advance", and then calls $d\delta_{ab}/d\theta=(2h/\beta^2)\ln A_{ab}$ "the time delay", although $\ln A_{ab}<0$. These two statements use opposite signs. Decide which convention is intended, state it once, and check it against the sine-Gordon classical shift (solitons are shifted forward on collision, which is a time advance in the usual sign convention). Also check the prefactor $2h/\beta^2$ in $d\delta/d\theta$ against the semiclassical formula cited, [@fring1994].

**A4. Forward references that depend on unproven hypotheses (`10-solitons.qmd`, Section `sec-vertex-solitons`, first bullet).** The bounds on products of different $F$'s are "hypothesis (D) of `sec-hirota-form-dn1-en1`, which is verified only numerically". Several results of Part III are stated as if established. Check which claims in Part III (in particular the degree claim and the multi-soliton formula for all algebras) depend on (D), and add a sentence to each such claim that says so.

**A5. The $e_8^{(1)}$ mass ratios (`09-lagrangian.qmd`, Section `sec-masses-couplings`, the table).** The entry reads "$m_1:m_2:\dots=1:2\cos\frac\pi5:2\cos\frac\pi{30}:\dots$". The ratio $m_2/m_1=2\cos\frac\pi5$ agrees with the Perron–Frobenius eigenvector in `sec-perron-frobenius` of Part II. *Partly resolved on 2026-10-10:* the third entry $m_3/m_1=2\cos\frac\pi{30}$ is confirmed by the fusing angles of `13-simply-laced.qmd` (Part IV), where the poles of $S_{11}=\{1\}\{11\}\{19\}\{29\}$ at $i\pi/15$ give $m_3=2m_1\cos\frac\pi{30}$ through the mass-triangle formula. Remaining: the table lists only three of the eight masses. Either complete it from the same calculation or mark it as partial.

**A6. Cubic coupling is for real coupling only (`09-lagrangian.qmd`, Section `sec-masses-couplings`).** The formula $C_{abc}=m^2\beta\sum_jn_j(\vec\alpha_j\cdot\vec e_a)(\vec\alpha_j\cdot\vec e_b)(\vec\alpha_j\cdot\vec e_c)$ follows from $V_{\rm real}$. For $V_{\rm imag}$ the cubic term picks up a factor $i^3$ and a sign. State "real coupling" before the formula. (I checked that the stated formula reproduces $|C_{111}|=3\beta$ for $a_2^{(1)}$ and the area rule $(4\beta/\sqrt h)\Delta$ with $h=3$, $\Delta$ the equilateral triangle of side $\sqrt3$, so the real-coupling formula is consistent.)

**A7. Boost compensation in the principal gradation (`09-lagrangian.qmd`, Section `sec-lax-pair`, first bullet).** The text says a boost is "compensated by $\lambda\mapsto e^\vartheta\lambda$". With $\partial_\pm\mapsto e^{\mp\vartheta}\partial_\pm$, the component $A_+$ must scale by $e^{-\vartheta}$, so $\lambda$ must scale by $e^{-\vartheta}$ unless the sign convention for $\vartheta$ is reversed. The statement also needs the gauge transformation by the principal grading automorphism. Fix the sign and add the gauge transformation, or say "up to a gauge transformation". Low priority.

**A8. Selection-rule logic (`11-reality.qmd`, Section `sec-reality-energy-momentum`, first bullet).** The text says "If the asymptotic data are $\Theta$-invariant ... then every conserved charge satisfies $Q=\epsilon_Q\bar Q$ ... in particular $E$ and $P$ are real." It then says a single soliton with real rapidity is not $\Theta$-invariant when $\bar a\ne a$, yet its $E$ and $P$ are real "because $M_{\bar a}=M_a$". So the reality of $E$ and $P$ has two different sources. State the proposition so that its hypothesis covers both cases, or add the second argument explicitly.

**A9. Unsupported claim about spin-2 charge (`11-reality.qmd`, Section `sec-reality-energy-momentum`).** "In $a_2^{(1)}$ the spin-2 charge, normalized to be real on single solitons, is imaginary on breathers." No computation or citation is given. Either derive it (for example from the spin-2 current of $a_2^{(1)}$ and a breather $\tau$-function) or cite the source, and mark it as an unverified claim if neither is possible.

**A10. Antilinear symmetry of the real-coupling theory (`11-reality.qmd`, `exr-theta-solutions`).** The exercise asks for "the antilinear map that is" a symmetry of the real-coupling equation, but no answer is given. Provide one (the natural candidate is complex conjugation composed with the Chevalley involution of the affine diagram, which is a diagram automorphism, and it reduces to $\vec\phi\mapsto-\vec\phi$ for $a_1^{(1)}$), and verify it on $a_2^{(1)}$. Check that the diagram involution maps $\sum_j n_j e^{\beta\vec\alpha_j\cdot\vec\phi}$ to itself.

**A11. Dorey rule: the roots and the projection (`09-lagrangian.qmd`, Section `sec-dorey-rule`).** Two sentences are not checked in the text: "each $\vec\gamma_a$ is a root" (true for the $A_2$ example, which I checked: $\vec\gamma_1=\alpha_1+\alpha_2$), and "the projections of the roots onto this plane form triangles whose sides are proportional to $m_a$, $m_b$, $m_c$". Check the second against Dorey (1991), and state the first only for simply-laced $\hat g$.

## B. Gaps for a reader following Parts I–II

**B1. Chevalley generators of the loop algebra.** `09-lagrangian.qmd`, Section `sec-lax-pair`, uses $E_{\pm\vec\alpha_j}$ for the affine simple roots, including $j=0$, with the normalization $[E_{\vec\alpha_j},E_{-\vec\alpha_j}]=\frac{2\vec\alpha_j}{\vec\alpha_j^2}\cdot\vec H$. A reader needs one sentence saying that $E_{\pm\vec\alpha_0}$ are the loop-algebra generators $E^{\mp\vec\theta}\otimes\lambda^{\pm1}$ of `sec-affine-dynkin` (Part II), and that $[E_{\vec\alpha_i},E_{-\vec\alpha_j}]=0$ for $i\ne j$ is used.

**B2. Drinfeld–Sokolov procedure.** Named in `09-lagrangian.qmd` (Section `sec-conserved-charges`) and not explained. Give one sentence: the gauge transformation that brings the Lax connection into the principal Heisenberg form, so that the spectral parameter becomes the grading variable.

**B3. Hirota's bilinear derivative and the convention for $\partial_\pm$.** Chapter 10 uses $D_x^2 f\cdot g$ and $\partial_+\partial_-\ln\tau=-\frac12(D_x^2-D_t^2)\tau\cdot\tau/\tau^2$. This is correct for $\partial_+\partial_-=\partial_t^2-\partial_x^2$, which is the convention of chapter 9 (see C1), and it is stated without that convention. Say so once.

**B4. Improvement term and virial relation.** `11-reality.qmd` says the energy-momentum tensor is a total derivative "up to an improvement term". The improvement term is not defined, and `exr-soliton-mass` uses $E=2\int V\,dx$ for static solitons without derivation. Add one sentence each: the improvement term is the term that makes the stress tensor traceless or symmetric (state which), and $E=2\int V$ follows from the first integral $\frac12\phi'^2=V$ of the static equation.

**B5. Embedded sine-Gordon kinks.** `10-solitons.qmd`, opening paragraph. The phrase is used without definition. Give the example in one line: for $a_{2n-1}^{(1)}$ at $\operatorname{Im}\xi=\pi/2$, $\tau_{2k+1}=\bar\tau_{2k}$ and the field is real.

**B6. Forward references.** Several results are used before their proof: Dorey's rule (chapter 9, cited to Dorey 1991 and not reproved in Part III), Lie duality (chapter 10, `sec-soliton-masses`, with the example `exr-lie-duality`), the breather and excited-soliton claims (chapter 10, `sec-breathers`), and the $\Theta$ selection rule's role in the quantum theory (chapter 11, `thm-krein-transfer`). State at each place whether the claim is proved in Part III, cited, or deferred.

**B7. Missing-charge problem.** `11-reality.qmd`, Section `sec-missing-charges`, gives three responses without saying which the book adopts. Add one sentence pointing to `sec-which-saddles-contribute` and `cnj-missing-charge` as the place where the choice is made.

## C. Notation clashes

**C1. $\partial_\pm$.** Part I (`02-integrable-qft.qmd`, `sec-factorization`) uses $x^\pm=x^0\pm x^1$, so $\partial_\pm=\tfrac12(\partial_0\pm\partial_1)$. Chapter 9 uses $\partial_\pm=\partial_t\pm\partial_x$, so $\partial_+\partial_-=\partial_t^2-\partial_x^2$. The conservation-law computation of Part I uses the first convention and the spin-3 current of chapter 9 uses the second. Both are internally consistent, but the same symbols mean different things. Fix: adopt one convention in Part I and III (and in the chapter 9 sign of the Lagrangian derivative), or state the convention at the first use in each part.

**C2. $\Theta$.** The antilinear operator $\Theta$ (chapter 9, `sec-real-imaginary`; chapter 11, `sec-theta-classical`) and the spin-$(s-1)$ current component $\Theta_{s-1}$ of Part I (chapter 2, `sec-factorization`). Also used in chapter 11, `prp-hirota-energy`, for the charges $\epsilon_Q$. Rename the current component as in the Part I brief (C4 there).

**C3. $\xi$.** Chapter 10 uses $\xi$ for the complex position constant of a Hirota soliton (`eq-an-soliton`). Part I (chapter 3) and Part II use $\xi$ for the renormalized breather coupling $\xi=\pi\beta^2/(8\pi-\beta^2)$ and $\xi=\pi/\lambda$. The two are used in the same chapter 10 (for breathers, "$\xi$'s"). Rename the Hirota constant, for example $\zeta$, throughout chapters 10 and 11.

**C4. $u$.** Chapter 10 says "$u$ here is $\frac\pi2$ minus the $u$ of `sec-sg-classical`". This is correct, but the same symbol is then used for two parameters in two chapters. Rename the chapter 10 parameter, for example $\vartheta$ or $u'$, and say so.

**C5. $E$.** Three uses: the Hirota exponential $E=\exp(m_a(x\cosh\theta-t\sinh\theta)+\xi)$ (chapter 10, `eq-an-soliton`), the energy $E$ (`eq-soliton-mass`, `prp-hirota-energy`), and the root generators $E_{\pm\vec\alpha_j}$ (chapter 9, `sec-lax-pair`). Rename the Hirota exponential, for example $\mathcal E$ or $e_a$, or the energy to $\mathcal E$.

**C6. $A$.** The Lax connection $A_\pm$ (chapter 9, `eq-lax-pair`) and the interaction coefficient $A_{ab}(\theta)$ (chapter 10, `eq-interaction-coefficient`). Rename the interaction coefficient to $\kappa_{ab}$, or the connection to $\mathcal A_\pm$.

**C7. $M$.** The mass matrix $M^2$ (chapter 9), the classical soliton mass $M_a$ (chapter 10), and the quantum soliton mass $M$ (Part I, chapter 3). Keep $M$ for the soliton mass, and write the matrix as $\mathsf M^2$ or $\mathbb M^2$.

**C8. $\omega$.** The $h$-th root of unity $\omega=e^{2\pi i/h}$ (chapter 10, `eq-an-soliton`) and the frequency $\omega$ of Part I (chapter 3, `eq-sg-fluctuation`). Rename the root of unity, for example $\varpi$.

**C9. $a$.** The species index $a\in\{1,\dots,n\}$ (chapters 9–11) and the algebra labels $a_n^{(1)}$, $a_{2n-1}^{(2)}$. Consider writing the species index as $\mathsf a$, or $i$, in chapter 10 and 11 only where the algebra name is nearby.

**C10. $\lambda$.** The Lax spectral parameter $\lambda$ (chapter 9, `eq-lax-pair`), the sine-Gordon parameter $\lambda=8\pi/\beta^2-1$ (Part I, chapter 3, and Part II, chapter 7), and the loop variable of Part II. Rename the Lax spectral parameter, for example $z$. This clash is the same as C3 in the Part II brief.

**C11. $k$.** The twist order in $kh$ and $X_N^{(k)}$ (chapter 9, `sec-conserved-charges`) and the level of Part II. Same as C1 in the Part II brief.

**C12. $c_j$ and $C_{abc}$.** The Lax-pair coefficients $c_j=\sqrt{n_j\vec\alpha_j^2/2}$ (chapter 9) and the cubic couplings $C_{abc}$ (chapter 9), against the Casimir $C(\vec\lambda)$ of Part II. Low priority.

## D. Items checked and found correct (do not re-check unless you change them)

- Mass matrix $M^2=m^2\sum_jn_j\vec\alpha_j\vec\alpha_j^{\mathsf T}$ from either potential (the linear term vanishes by $\sum_jn_j\vec\alpha_j=0$); $V_{\rm imag}$ has a minimum at $\vec\phi=0$ and $\operatorname{Re}V_{\rm imag}\ge0$ with equality exactly on the lattice $\frac{2\pi}{\beta}\Lambda^\vee$.
- The field equations (`eq-toda-eom`) follow from the potentials; the imaginary one is the real one with $\beta\to i\beta$.
- $V_{\rm real}$ has a unique vacuum at $\vec\phi=0$, by convexity and Jensen's inequality.
- Lax pair (`eq-lax-pair`): the $E_{\pm\vec\alpha_j}$ terms cancel, and the Cartan part gives exactly (`eq-toda-eom`), with $c_j^2\cdot2\vec\alpha_j/\vec\alpha_j^2=n_j\vec\alpha_j$.
- The spin-3 current (`eq-spin3-current`) is conserved for sinh-Gordon, checked by direct differentiation; the charge $Q_3=\int(T_4+\Theta_2)$ has the sign consistent with the convention of Part I, chapter 2.
- Spin sets: $a_1$ odd integers, $a_n^{(1)}$ all residues except $0$ mod $h$, and the dual pairs $(b_n,a_{2n-1}^{(2)})$, $(c_n,d_{n+1}^{(2)})$ all odd integers.
- Sum rule $\sum_am_a^2=2hm^2$ for simply-laced, and the $d_n^{(1)}$ masses of the table: checked for $d_4^{(1)}$ (sum $12$) and $d_5^{(1)}$ (sum $16$); the central node of $d_4^{(1)}$ is $\sqrt3$ times the legs, as in `exr-d4-masses`.
- $c_n^{(1)}$ masses $2m\sin\frac{\pi a}{2n}$, $a=1,\dots,n$, by folding $a_{2n-1}^{(1)}$; checked for $c_2^{(1)}$ from $a_3^{(1)}$.
- Cubic couplings of $a_2^{(1)}$: only $C_{111}$ and $C_{222}$ are nonzero, $|C_{111}|=3\beta$, consistent with the area rule (`eq-area-rule`) with $h=3$; the sinh-Gordon cubic coupling vanishes.
- Dorey's rule for $A_2$: the orbits of $\gamma_1=\alpha_1+\alpha_2$ are $\{\alpha_1+\alpha_2,-\alpha_1,-\alpha_2\}$ and its negative; they reproduce $C_{111}\neq0$, $C_{222}\ne0$, $C_{112}=0$.
- Hirota ansatz (`eq-hirota-ansatz`): $e^{i\beta\vec\alpha_k\cdot\vec\phi}=\prod_j\tau_j^{-a_{jk}}$; the bilinear equation (`eq-hirota-bilinear`) follows from the field equation, with the constant $-1$ allowed because $\sum_jn_j\vec\alpha_j=0$.
- Single soliton (`eq-an-soliton`): static check gives $\frac12D_x^2\tau\cdot\tau=m_a^2\omega^{ja}E$ and $\tau_j^2-\tau_{j-1}\tau_{j+1}=4\sin^2\frac{\pi a}h\,\omega^{ja}E$, so $m_a=2m\sin\frac{\pi a}h$.
- Sine-Gordon dictionary: $\vec\phi=\frac{i}{\beta}\sum_j\frac{2\vec\alpha_j}{\alpha_j^2}\ln\tau_j$ reproduces $\phi_+=\frac{2i}{\beta_{\rm SG}}\ln(\tau_1/\tau_0)$ of Part I, chapter 3, and $\tau_{0,1}=1\pm iE$ at $\operatorname{Im}\xi=\pi/2$.
- Topological charges: the species-1 soliton of $a_n^{(1)}$ realizes all $h$ weights; the species-2 soliton of $a_3^{(1)}$ realizes 2 of 6, and the species-1 soliton of $a_2^{(1)}$ realizes all 3 (counts $h/\gcd(a,h)$ against $\binom ha$).
- Soliton mass $M_a=2hm_a/\beta^2$ at sine-Gordon: $M=4m_0/\beta^2=8m_0/\beta_{\rm SG}^2$, consistent with Part I.
- Interaction coefficient $A_{ab}$: for sine-Gordon $A_{11}=\tanh^2(\theta/2)$; $0<A_{ab}<1$ for real rapidities; $\overline{A_{ab}(\theta)}=A_{h-a,h-b}(\bar\theta)$ (`exr-theta-interaction`).
- Fusion identity $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$ with $m_a=\sin\frac{\pi a}h$ (real part gives $\sin\frac{\pi(a+b)}h$, imaginary part vanishes).
- Breather energy $E=2M_a\cosh\theta_0\cos u$ is real; the breather data are $\Theta$-invariant.
- $\Theta$ maps solutions to solutions: $\overline{\omega^{ja}}=\omega^{j(h-a)}$ and $(a,\theta,\xi)\mapsto(\bar a,\bar\theta,\bar\xi)$.
- Isolated zeros of $\tau_j$ in the real $(x,t)$ plane: one complex condition, two real conditions.
- Gap ids H6 and H11 of `28-twelve-gaps.qmd` agree with the status claimed in chapter 11.
- All 59 `@` labels and all 28 citation keys used in Part III resolve, and `textbook_code/toda_classical.py` exists.

## E. Verification checklist for the agent who makes the fixes

1. Re-run the label check: every `@sec-`, `@eq-`, `@prp-`, `@exr-`, `@cnj-` and `@thm-` in `part3-classical/` resolves to an anchor elsewhere in the book (exclude `.claude/` and `_site/`).
2. For A1, recompute the degree of the species-$a$ soliton at a short node of $c_2^{(1)}$ from its $\tau$-function before editing the sentence.
3. For A2, A3 and A6, check the sign and the factor by direct computation with `textbook_code/toda_classical.py` or a short script.
4. For each notation rename (C3–C10), grep the whole book (excluding `.claude/` and `_site/`) for the old symbol in the same meaning, and change only those occurrences.
5. Report what changed, what was left, and anything you could not verify.

## Resolution (2026-10-10)

Scratch checks are in `checks.py` (A1, sympy), `e8.py` (A5, A11, numpy), `sgcheck.py` (A3, mpmath) and `labels.py`, in the agent's scratch directory.

- **A1 (claim wrong; clarified).** The $O$ in the text is an orbit of parent *solitons*, not of diagram nodes. A coincident pair of species $a$, $2n-a$ in $c_n^{(1)}$ has degree $2=\lvert O\rvert n_j^\vee$ at every node, including the long end nodes, because every dual label of $c_n^{(1)}$ is 1. Checked for $c_2^{(1)}$ from $a_3^{(1)}$: the pair (1,3) gives $\tau_0=1+2E+\frac12E^2$, $\tau_{1,3}=1+\frac12E^2$ and $\tau_2=1-2E+\frac12E^2$, and the species-2 soliton has degree 1 everywhere, including the short node. Both solve the bilinear equations exactly. The sentence was rewritten to say what $O$ is, and the explicit example was added.
- **A2 (changed).** Confirmed: $i\vec\phi$, not $-i\vec\phi$, solves the real-coupling equation. `exr-theta-solutions` contains no such sign. Its new answer uses $\vec\phi\mapsto i\vec\phi$.
- **A3 (text was consistent; conventions now stated).** $\Delta x>0$ (forward shift) and $\Delta t=-\Delta x/v<0$ (advance) agree with Fring–Johnson–Kneipp–Olive (hep-th/9405034, eqs. 2.5, 2.10, 4.10): $E_1\Delta x=-p_1\Delta t=-(2h/\beta^2)\ln X$, and the time delay is negative. The text now defines $\Delta t$ and $S=e^{i\delta}$, and derives $d\delta/d\theta=p_1\Delta t=(2h/\beta^2)\ln A$ from Wigner's relation. That paper states Eisenbud–Wigner only qualitatively, so this derivation is my own. Prefactor checked against sine-Gordon: the small-$\xi$ limit of $-i\,d\ln S_0/d\theta$ is $(2/\xi)\ln\tanh\frac\theta2=(16/\beta_{\rm SG}^2)\ln\tanh\frac\theta2$, using $\int_0^\infty\frac{dt}t\tanh\frac{\pi t}2\cos\theta t=-\ln\tanh\frac\theta2$ (mpmath). This equals $(2h/\beta^2)\ln A_{11}$.
- **A4 (changed).** The only Part III statements that depend on (D) are the degree statements: the exact degree $n_j$ for $d,e$ and the degree of coincident orbits. They now say so, and "Other algebras" no longer asserts the degree as fact. The multi-soliton formulas for other algebras come from the literature, not from (D), and the text now says that too.
- **A5 (changed).** All eight $e_8^{(1)}$ masses are now in the table. They agree with the eigenvalues of $\mathsf M^2$ to $10^{-14}$.
- **A6 (changed).** "Expanding $V_{\rm real}$", the normalization $V\supset\frac16C_{abc}\phi_a\phi_b\phi_c$, and a line saying the imaginary-coupling couplings carry an extra factor $i$, since $-(i\beta)^3=i\beta^3$.
- **A7 (partly wrong).** The sign depends on whether the boost acts on fields or on coordinates. For the active boost $\vec\phi_\vartheta(x^\pm)=\vec\phi(e^{\mp\vartheta}x^\pm)$ one finds $\mathcal A_\pm[\vec\phi_\vartheta](x;\lambda)=e^{\mp\vartheta}\mathcal A_\pm[\vec\phi](x';e^\vartheta\lambda)$, so $\lambda\mapsto e^{\vartheta}\lambda$ is correct and agrees with `sec-gradations` (multiply degree $s$ by $e^{s\vartheta}$). The brief's $e^{-\vartheta}$ is the passive version. In the principal gradation no gauge transformation is needed, because $\lambda$ multiplies every $E_{\vec\alpha_j}$; one is needed only in other gradations. The text now says all this.
- **A8 (changed).** A second bullet gives the weaker condition for $E$ and $P$: the multiset $\{(M_{a_i},\theta_i)\}$ is closed under $\theta\mapsto\bar\theta$. It covers single solitons with $\bar a\ne a$.
- **A9 (derived, with a stated assumption).** Freeman (hep-th/9408092, eqs. 3.47–3.49 and Conclusions; PDF read) computes the charges of *single* solitons only. The spin-$N$ values form the left eigenvector of the Cartan matrix of $g$ with eigenvalue $2(1-\cos\pi N/h)$, real up to an overall phase. For $A_2$, $N=2$, this is $(1,-1)$, so $q_2=-q_1$. Assuming additivity over constituents at complex rapidity, a breather has $2iq_1e^{2\theta_0}\sin2u$, which is imaginary, consistent with $\epsilon_Q=-1$. The text gives this derivation and flags the additivity assumption. Side finding, not edited: `prp-hirota-energy` cites [@freeman1995] for "the same holds for all higher charges" of multi-soliton solutions, but Freeman covers single solitons only.
- **A10 (changed; the brief's candidate was wrong).** The real-coupling equation has real coefficients, so plain complex conjugation $\vec\phi\mapsto\bar{\vec\phi}$ is the antilinear symmetry, possibly composed with a diagram symmetry. It is exactly the image of $\Theta$ under $\vec\phi\mapsto i\vec\phi$. The Chevalley involution is not a diagram automorphism, and no extra map is needed. The answer is now in the exercise.
- **A11 (partly changed).** $\vec\gamma_a$ is a root for every $g$: it is $\vec\alpha_a$ for white $a$ and $s_{\rm white}\vec\alpha_a$ for black $a$. The text now gives this one-line proof and restricts the partition statement to simply-laced $g$. The projection statement was confirmed in Dorey II (hep-th/9110058, §2, around eq. 2.14: fusing angles are relative angles of projections into the $\omega^1$ eigenspace), and [@dorey1992] was added to the citation. It was also checked numerically for $E_8$: the projected lengths equal the masses, and all 13440 root triples have angles that are multiples of $\pi/30$. I could not access Dorey 1991 itself.
- **B1 (changed).** $E_{\vec\alpha_j}=e_j$, $E_{-\vec\alpha_j}=f_j$ of `sec-affine-dynkin`, with $e_0,f_0$ carrying the loop variable. $\lambda$ is the principal-gradation spectral parameter. **C10** is handled in the same place by a note, as in Part II, not a rename: $\lambda$ is unrelated to the sine-Gordon $\lambda$.
- **B2 (changed).** One sentence on Drinfeld–Sokolov.
- **B3 (changed).** The convention is stated at the Hirota-derivative identity in ch. 10.
- **B4 (changed, differently).** The phrase "up to an improvement term" is replaced by the explicit formula $T_{\mu\nu}=(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$, $C=-\frac2{\beta^2}\sum_j\frac2{\vec\alpha_j^2}\ln\tau_j$, with a note that $C$ involves $\ln\tau_0$ and is not local in $\vec\phi$. The formula is confirmed via Corrigan's review (hep-th/9412213, §4), which attributes it to OTU; OTU itself was not accessed. Its trace part and the static energy were checked by hand from the bilinear equations. The hint of `exr-soliton-mass` now derives $E=2\int V$ from the first integral, and $V$ from the bilinear equations. (The brief located B4 in ch. 11; both places are in ch. 10.)
- **B5 (changed slightly).** The example was already present. Added the definition $\vec\phi=\vec v\,\varphi$.
- **B6 (changed).** Status sentences were added: Dorey's rule is quoted without proof, Lie duality is checked only for $c_n^{(1)}$ (otherwise quoted from kneipp1996), breathers are quoted (quantization deferred to Parts V and VI), excited solitons are not treated, and `thm-krein-transfer` is proved in Part VII and used here only as motivation.
- **B7 (changed).** The text now says the book adopts the third response, made in `sec-which-saddles-contribute` and stated in `cnj-missing-charge`.
- **C1 (note).** Ch. 9 now says that $\partial_\pm$ is twice $\partial/\partial x^\pm$ of `sec-factorization`, and that the factor drops out of conservation laws. Part I was not edited.
- **C2 (nothing to do).** Already done by the Part I fixes ($\mathcal X_2$ in ch. 9). No $\Theta_{s}$ current remains in the book.
- **C3, C5, C8 (note plus minimal rename).** The Hirota $\xi$, $E$ and $\omega$ are also used in Part V (ch. 17–19), so a rename would touch about 30 lines across two Parts. The coupling $\xi$ never occurs in Part III, so the "same chapter" clash claimed in C3 does not exist (the "$\xi$'s" of the breathers are Hirota constants). A one-paragraph notation note follows `eq-an-soliton`. The real in-chapter clash, energy $E$ against the Hirota $E$ in ch. 10, was removed by dropping the symbol $E$ for the energy in the breather display and in the `exr-soliton-mass` hint.
- **C4 (not changed).** The text already flags the shift ("$u$ here is $\frac\pi2$ minus the $u$ of …"). A rename would touch ch. 10 and ch. 11 for little gain.
- **C6 (renamed).** The Lax connection is now $\mathcal A_\pm$ (6 occurrences, ch. 9 only). $A_{ab}$ is kept. Note that Part V uses $\mathcal A$ for the fluctuation operator, which has no $\pm$ index.
- **C7 (renamed).** The mass matrix is now $\mathsf M^2$ in ch. 9 (4 lines) and `12-perturbation-theory.qmd` (1 line, same meaning). Part V's $M(x)$, $M_0$ is left as is.
- **C9, C12 (not changed).** Low priority. No local ambiguity was found.
- **C11 (nothing to do).** The level is now $\ell$ (Part II fixes).
- **Verification.** All `@` labels and citation keys in ch. 9–11, `part.qmd` and ch. 12 resolve. A scratch `quarto render --to html` finished with exit 0 and no errors or warnings in the log, and there is no `?@` anywhere in `_site`.

