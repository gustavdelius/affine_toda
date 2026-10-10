# Brief: fixes and gaps from the Part I reading (Chapters 1–4)

Status: done 2026-10-10 (uncommitted, branch `part1-reading-fixes`; see Resolution at the end). Written 2026-10-10 after a line-by-line reading of Part I by a reader with graduate-level mathematical physics and QFT, but no prior knowledge of integrable QFT, affine Toda theory, Kac–Moody algebras or quantum groups.

Files: `part1-background/part.qmd`, `01-introduction.qmd`, `02-integrable-qft.qmd`, `03-sine-gordon.qmd`, `04-non-hermitian.qmd`.

Scope: fix the items below in place. Do not edit anything else without checking with the author. Verify each fix by reading the surrounding text and by re-running any formula check you rely on. The Part I formulas that I checked by hand are listed at the end, so you need not re-derive them.

## A. Errors to fix

**A1. Proof of Proposition `prp-krein-spectrum`, item 2 (`04-non-hermitian.qmd`, the proof, item 2).** The text reads "$\lambda[x,y]=[Hx,y]=[x,Hy]=\mu[x,y]$". The form $[\cdot,\cdot]$ is antilinear in its first argument, so the first equality must be $[Hx,y]=\bar\lambda[x,y]$. The conclusion $(\bar\lambda-\mu)[x,y]=0$ is correct, so only the first displayed equality needs changing.

**A2. Same chapter, the paragraph after the proposition.** "Part 4 says that time evolution is Krein-unitary" should refer to item 4 of `prp-krein-spectrum`.

**A3. Convergence region of $S_0$ (`03-sine-gordon.qmd`, the sentence after `eq-sg-s0`).** The text says the integral converges for $|\operatorname{Im}\theta|<\min(\xi,\pi)$. For large $t$ the integrand behaves like $e^{(|\operatorname{Im}\theta|-\xi)t}/t$, so the convergence region is $|\operatorname{Im}\theta|<\xi$. The stated region is a sufficient but non-sharp condition, and it is misleading when $\xi>\pi$ (repulsive regime, $\lambda<1$). Replace $\min(\xi,\pi)$ by $\xi$, and check the estimate. The sign and normalization of the formula itself were not re-derived; check them against Zamolodchikov–Zamolodchikov (1979) if you change anything else in it.

**A4. Attribution to check (`04-non-hermitian.qmd`, the paragraph on the scaling Lee–Yang model and the Bethe–Yang discussion in `02-integrable-qft.qmd`).**
- In ch. 2 (end of "Finite volume: the Bethe–Yang equations"), complex finite-volume energies in non-unitary integrable theories "were first predicted in exactly this way [tw1999; tw2002]". The bibliography entries are Takács–Watts, "Non-unitarity in quantum affine Toda theory and perturbed conformal field theory" (1999) and "RSOS revisited" (2002). Check whether the 1999 paper makes this prediction and whether priority belongs to someone else. Soften the sentence if it does not.
- In ch. 4 (end of the Lee–Yang section), "it is a restriction of $a_2^{(2)}$ affine Toda theory too [smirnov1991]". The bibliography title of `smirnov1991` is "Exact S-matrices for $\phi_{1,2}$-perturbated minimal models". Check that this paper supports the $a_2^{(2)}$ claim, or cite the correct source.

## B. Gaps for a reader without Part II

**B1. The affine Toda theory is never defined in Part I.** Chapter 1 says real and imaginary coupling refer "to the two potentials of (-@eq-imaginary-potential)" (`01-introduction.qmd`, Conventions). Chapter 3 derives the $a_1^{(1)}$ potential from the same equation ("The imaginary-coupling potential ... becomes"). The label `eq-imaginary-potential` is in `part3-classical/09-lagrangian.qmd`, Part III. The first use of "real coupling", "imaginary coupling" and "affine Toda field theory" therefore comes before any definition of the theory. Fix: state the general potential in one or two lines in ch. 1, in the form $V=\sum_i(\text{Kac label})\,(\dots)$ or whatever Part III uses, or add an explicit pointer that the potential is taken on trust until `eq-imaginary-potential`. Use the exact form from Part III so the notation matches.

**B2. Part II terms used in Part I without definition.** Each needs at least a one-line gloss or a forward pointer. Items:
- Chapter 1: "affine Kac–Moody algebra $\hat g$"; "Dorey's fusing rule", "Perron–Frobenius masses", "Lie duality" (listed as features, not explained); "W-algebra minimal models"; "quantum-group restriction".
- Chapter 2: "exponents of the affine algebra, repeated with period $h$ for untwisted and $kh$ for twisted"; "Cartan matrix"; "Hopf algebra structure" (in the crossing paragraph); "charge conjugation matrices".
- Chapter 3: "Kac labels $n_0=n_1=1$" and the roots $\alpha_1=\sqrt2$, $\alpha_0=-\sqrt2$; "affine Dynkin diagram"; "Perron–Frobenius eigenvector"; "dual algebra" and $U_q(\hat g^\vee)$ in the table; "evaluation representation"; "IRF (RSOS) basis"; "Lie dual" in the duality remark.
- Chapter 4: "PT-type" is clear to a mathematical physicist, but the sentence "$\Theta=PK$ ... in the field representation" needs the general potential (see B1).

The dependency on Part II is acceptable, but it should be stated. Fix: add a short paragraph to ch. 1 (after "How this book is organised") saying which Part I claims depend on Part II, and that they can be taken on trust on a first reading. Alternatively, add a glossary box.

**B3. Chapter 2, spin formula.** The statement "the spins of the conserved charges are the exponents of the affine algebra, repeated with period $h$ ..." is only checked for $a_1^{(1)}$ in ch. 3 (spins $1,3,5,\dots$). Add a remark that for $a_1^{(1)}$ this gives $m=1$ and $h=2$, to make the example concrete.

## C. Notation clashes

Fix these by renaming in the affected chapter, or by adding a note at the first occurrence. Choose the fix that changes the fewest occurrences across Parts I–VII; grep before renaming.

**C1. $s$** is used for the Mandelstam variable (`eq-mandelstam-s`), the spin label in $Q_s$ (ch. 2, factorization), and the $s$-channel (ch. 2, residues and bound states). Suggestion: keep $s$ for Mandelstam, write the spin as $\mathsf s$ or $j$ in $Q_j$, and the bound-state label as $c$ only.

**C2. $h$** is used for the Coxeter number (ch. 2, "period $h$"), the conformal weight $h=-\tfrac15$ and $h_{\min}$ (ch. 4, Lee–Yang and exr-ceff), and the magnetic field in Fisher's Lagrangian $i(h-h_c)\phi$ (ch. 4). Chapter 4 uses $h$ in two senses. Suggestion: magnetic field $H$ or $\mathfrak h$; keep $h$ for the Coxeter number and use $\Delta$ or $h_{\rm w}$ for weights, after checking Part II and Part VII usage.

**C3. $g$** is used for the dimer gain/loss parameter (ch. 4, `eq-dimer`) and for the cubic coupling $ig\phi^3$ in Fisher's Lagrangian (ch. 4), both in the same chapter. It is also the symbol for the Lie algebra $\hat g$ (ch. 1, ch. 2). Suggestion: dimer gain/loss as $\gamma$.

**C4. $\Theta$** is the PT-type operator (ch. 1, ch. 4, Section `sec-antilinear`), and also the spin-$(s-1)$ component of a conserved current $\Theta_{s-1}$ (ch. 2, `sec-factorization`). Suggestion: rename the current component to $\mathcal X_{s-1}$, or write $\Theta$ for the operator only.

**C5. $\lambda$** is an eigenvalue (ch. 4, Krein and dimer sections) and the sine-Gordon parameter $\lambda=8\pi/\beta_{\rm SG}^2-1$ (ch. 3). The two uses are in different sections. Note this in the ch. 3 sentence that introduces $\lambda$.

**C6. $\beta$ and $\beta_{\rm SG}$.** Both appear in ch. 3 with the dictionary $\beta_{\rm SG}^2=2\beta^2$. This is consistent, but a reader will need the dictionary stated before first use in the $a_1^{(1)}$ section (it is currently on the line after the Lagrangian). Check the order.

## D. Items checked and found correct (do not re-check unless you change them)

- `exr-strip`: the map $s(\theta)=m_a^2+m_b^2+2m_am_b\cosh\theta$ sends the open strip bijectively onto the cut plane, and the two lines $\operatorname{Im}\theta=0,\pi$ map to the two edges of the cuts.
- Braiding unitarity from the exchange relation, `eq-braiding-unitarity`; $S(\theta+2\pi i)=S(\theta)$ for self-conjugate particles.
- $f_x$ in tanh form; its poles at $i\pi x$, $i\pi(1-x)$ and zeros at their negatives.
- $f_{2/3}(\theta)=f_{2/3}(\theta-i\pi/3)f_{2/3}(\theta+i\pi/3)$ (since $f_x$ depends only on $\sin\pi x$, $f_{2/3}=f_{1/3}$).
- The residue of $f_x$ at $i\pi x$ is $2i\tan\pi x$, so the sign rule in `exr-fx` holds, and the Lee–Yang residue is $i$ times a negative number.
- Bethe–Yang free-fermion quantization in `exr-free-fermion`.
- Sine-Gordon: the dictionary $\beta_{\rm SG}^2=2\beta^2$; the kink, $M_{\rm cl}=8m_0/\beta_{\rm SG}^2$; the breather with $\omega=m_0\cos u$, $E_b=2M\sin u$; the fluctuation operator and zero mode $\operatorname{sech}$; DHN masses $m_n=2M\sin(n\xi/2)$ with $\xi=\pi\beta^2/(8\pi-\beta^2)$.
- The one-loop check in `exr-sg-expansion`: $m_1=M\xi(1+O(\beta^4))=m_0+O(\beta^4)$.
- Sinh-Gordon continuation: $B=2\beta_{\rm SG}^2/(8\pi+\beta_{\rm SG}^2)$ and $B\to-2\xi/\pi$ under $\beta_{\rm SG}^2\to-\beta_{\rm SG}^2$; the duality $B\to2-B$.
- $\lambda=1$ gives $S_0=S_T=-1$, $S_R=0$; crossing relations for $S_0,S_T,S_R$.
- RSOS: $\xi=\pi p/(p'-p)$; for $(2,5)$, $\xi=2\pi/3$, one breather of mass $\sqrt3M$.
- Chapter 4: Krein adjoint, $J$-self-adjointness and the pseudo-Hermitian statements; the dimer eigenvalues, eigenvectors, Krein norms $\pm2s\sqrt{s^2-g^2}$, and $\Theta$-action in each regime; exr-self-orthogonal; exr-jordan-determinant (double eigenvalue 2, $\det A=4$, integral $=\pi$); exr-ceff ($c_{\rm eff}=1-6/pp'=2/5$ for $\mathcal M(2,5)$).
- All 64 `@` labels and all 57 citation keys used in Part I resolve.

## E. Verification checklist for the agent who makes the fixes

1. Re-run the label check: every `@sec-`, `@eq-`, `@prp-`, `@exr-` and `@thm-` in `part1-background/` resolves to a `{#...}` anchor elsewhere in the book (exclude `.claude/`).
2. Render or parse the Quarto files to confirm no broken syntax was introduced.
3. For each notation rename (C1–C5), grep the whole book (excluding `.claude/` and `_site/`) for the old symbol in the same meaning, and do not change other uses.
4. Report what changed, what was left, and anything you could not verify.

## Resolution (2026-10-10)

- A1, A2: fixed as stated.
- A3: **not changed — the brief's estimate was wrong.** For large $t$, $\sinh\frac{(\pi-\xi)t}2\sim\pm\frac12e^{|\pi-\xi|t/2}$, so the integrand decays like $e^{(|\operatorname{Im}\theta|-\min(\xi,\pi))t}$. For $\xi>\pi$ the rate is $\pi$, not $\xi$. The text's $\min(\xi,\pi)$ is sharp; checked numerically for $\xi=\pi/2$ and $3\pi$.
- A4: Takács–Watts 1999 (§2, eqs. 2.9–2.10) does make the Bethe–Yang argument, but credits Kausch–Takács–Watts (`tkw1997`, hep-th/9605104) for the first complex levels, which they found in $M_{3,14}+\Phi_{1,5}$ and confirmed by TCSA. The ch. 2 sentence now credits tkw1997 first and drops "first predicted". In ch. 4, $\phi_{1,3}\equiv\phi_{1,2}$ in $\mathcal M(2,5)$, so Lee–Yang is a $\phi_{1,2}$ perturbation, and those are $a_2^{(2)}$ restrictions [smirnov1991; efthimiou1993]. The sentence now says this. Smirnov's paper itself was not accessible, so whether it treats $\mathcal M(2,5)$ explicitly is unverified.
- B1: ch. 1, §why-atft-matters, now displays $V_{\rm real}$, $V_{\rm imag}$ verbatim from Part III with a gloss; Conventions, ch. 3 and ch. 4 point to it.
- B2: glossary paragraph "What Part I borrows from Part II" at the end of §reader-guide, plus local glosses and pointers in ch. 2 and 3. Ch. 3's "Lie dual" in the duality remark is replaced by "dual algebra $\hat g^\vee$", because Part II (§affine-folding) reserves "Lie duality" for $g^{(1)}\to(g^\vee)^{(1)}$.
- B3: added to ch. 2.
- C1: note at the first spin $s$ in ch. 2; no rename (spin $s$ is also used in Part III, §conserved-charges).
- C2: the magnetic field is renamed (see C3); a note in ch. 4 says $h$ is a conformal dimension there. $h$ for weights is kept, because Part VI uses it too.
- C3: the cubic coupling is renamed $\kappa$ and the field term $i\nu\phi$ (1 occurrence each); the dimer $g$ is kept (about 10 occurrences).
- C4: the current component is renamed $\Theta_{s-1}\to\mathcal X_{s-1}$ in ch. 2, and $\Theta_2\to\mathcal X_2$ in `part3-classical/09-lagrangian.qmd` (same meaning).
- C5: note added where $\lambda$ is defined in ch. 3.
- C6: the dictionary $\beta_{\rm SG}^2=2\beta^2$ is now stated in prose before eq-sg-lagrangian.
