# Brief: fixes and gaps from the Part IV reading (Chapters 12–15)

Status: done 2026-10-10 (uncommitted, branch `reading-fixes`; see Resolution at the end)

Files: `part4-real-coupling/part.qmd`, `12-perturbation-theory.qmd`, `13-simply-laced.qmd`, `14-non-simply-laced.qmd`, `15-beyond-smatrix.qmd`.

Scope: fix the items below in place. Items A2 and A3 need the literature (Delius–Grisaru–Zanon 1992, Corrigan–Dorey–Sasaki 1993) before any change, and the report must say what was checked. Do not edit anything else without checking with the author. The Part IV results I checked by hand or by script are in section D.

## A. Errors and likely errors

**A1. Exponents in Dorey's formula (`13-simply-laced.qmd`, Section `sec-dorey-formula`, after `eq-dorey-formula`).** The text says the exponents "then come out as integers" after reducing with $\{x+2h\}=\{x\}$ and $\{-x\}=\{x\}^{-1}$. The individual exponents are half-integers. For $A_2$, $S_{11}$ has exponents $\tfrac12,0,-\tfrac12$ on $\{1\},\{3\},\{5\}$. Since $\{5\}=\{-1\}=\{1\}^{-1}$ for $h=3$, the product is $\{1\}^{1/2}\{1\}^{1/2}=\{1\}$, which is correct. Likewise $S_{12}$ has exponents $\tfrac12,-\tfrac12,0$ on $\{2\},\{4\},\{6\}$, and $\{4\}=\{2\}^{-1}$ gives $\{2\}$. Fix: say that the individual exponents are half-integers, and that they combine into integer powers of the blocks after reduction. Give the $A_2$ example in one line.

**A2. Floating-Coxeter table and the parameter range (`14-non-simply-laced.qmd`, Section `sec-floating-coxeter`, `tbl-floating-h`, and Section `sec-nsl-masses`).** Three statements cannot all hold with $B\in[0,2)$ from `eq-b-of-beta`:
- the text gives $H=2n+B+O(B^2)$ (Section `sec-nsl-masses`), but the table gives $c_n^{(1)}$: $H=2n+\tfrac12B$ with range $2n\le H\le2n+2$. With $B\in[0,2)$, $2n+\tfrac12B$ lies in $[2n,2n+1)$;
- $g_2^{(1)}$: $H=6+B$ with range $6\le H\le12$. With $B<2$ this gives $[6,8)$;
- $f_4^{(1)}$: $H=12+\tfrac32B$ with range $12\le H\le18$. With $B<2$ this gives $[12,15)$.

The stated ranges correspond to $B$ up to 4 for $c_n$ and $f_4$ and up to 6 for $g_2$, and these are not consistent with each other. The short-root normalization sentence (after the table) gives $2n+B$, $6+3B$, $12+3B$, which is a different set of coefficients again. Also, the dual endpoint of $c_n^{(1)}$ is given as $d_{n+1}^{(2)}$ with masses $\sin\frac{\pi a}{2n+2}$, that is $H=2n+2$. But `tbl-foldings` (`06-affine-algebras.qmd`) gives $d_{n+1}^{(2)}$ the values $h=n+1$ and $h^\vee=2n$. Fix: decide the normalization of $B$ for each algebra (state it), recompute the endpoints from the dual algebra's $h^\vee$ or from the literature, and make the table, the one-loop sentence and the exercise `exr-floating-endpoints` agree.

**A3. The $H-12=3B_{\rm tree}$ sentence (`14-non-simply-laced.qmd`, after `tbl-floating-h`).** "The relation $H-12=3B_{\rm tree}$ confirmed in @sec-sm-f4 refers to the parameter of the CDS blocks, which is $\beta^2/4\pi$ in our normalization." Two parameters are called $B$ here ($B_{\rm tree}$ and $B$) and the statement that $B_{\rm tree}=\beta^2/4\pi$ conflicts with `eq-b-of-beta`, which gives $B=\beta^2/2\pi+O(\beta^4)$. Fix after A2: state the relation between the two parameters once, with its normalization.

**A4. Status of the non-simply-laced S-matrices (`13-simply-laced.qmd`, Section `sec-real-coupling`, versus `14-non-simply-laced.qmd`, Section `sec-nsl-status`).** The summary's title is "the real-coupling S-matrices are known for every algebra". Chapter 14 states that "the non-simply-laced S-matrices are conjectures". Make the summary say which S-matrices are conjectural. (Chapter 13 itself says the simply-laced ones are exact; keep that.)

**A5. "Every block tends to 1" (`13-simply-laced.qmd`, Section `sec-strong-weak`).** At $B\to0$ and $B\to2$ only the building blocks $\{x\}$ tend to 1. The blocks $(x)$ are $\theta$-dependent and do not. Change "every block" to "every $\{x\}$".

**A6. Ambiguous coupling label in `exr-tree-residue` (`12-perturbation-theory.qmd`).** "compute the tree amplitude for $1+1\to1+1$ from the cubic coupling $C_{11\bar1}$". With $\bar1=2$ this is $C_{112}$, which is zero (Part III, `exr-shg-cubic`, and my check there). The vertex used is $C_{111}$, with particle 2 as the exchanged bound state. Fix the label. Also the phrase "compare with the residue of the exact S-matrix $S_{11}=\{1\}$ of @sec-simply-laced expanded to first order in $B$" contains a stray brace and the $B$-expansion is not stated. I verified the residue $i\beta^2/6$ for $a_2^{(1)}$ with $|C_{111}|=3\beta$ and $m_1=\sqrt3$, $u=2\pi/3$.

**A7. Double-pole claim (`12-perturbation-theory.qmd`, Section `sec-coleman-thun`).** "In affine Toda theory, at each double pole of the exact S-matrix there is an on-shell box diagram built from the three-point couplings" uses "each", which is a strong claim. Check the statement in Braden, Corrigan, Dorey and Sasaki (1991) and replace "each" with what is shown. Also check "up to order 12 in $e_8^{(1)}$, contributing diagrams have been identified" against Corrigan–Dorey–Sasaki (1993), since the number is specific.

**A8. Ultraviolet theory of the real-coupling theory (`15-beyond-smatrix.qmd`, Section `sec-tba`).** "For the real-coupling affine Toda theories the ultraviolet limit is $r$ free bosons, and the TBA gives $c_{\rm eff}=r$." The ultraviolet theory is the Toda conformal field theory, which is not $r$ free bosons (its central charge is $r+12\rho^2$-type and is not $r$). The TBA plateau gives $c_{\rm eff}=r$. Rewrite the sentence to say that the TBA gives $c_{\rm eff}=r$ and that the ultraviolet CFT is not free.

**A9. Plateau equation in the hint of `exr-tba-uv` (`15-beyond-smatrix.qmd`).** "the plateau equations $x_a=1+x_a$ have no finite solution". This is not a plateau equation as written. Taking $\varepsilon_a$ constant in the central region and using $N_{ab}=\frac1{2\pi}\int\varphi_{ab}=\delta_{ab}$ (see A10), the TBA of `eq-tba` gives a plateau equation of the form $\varepsilon_a=-\sum_bN_{ab}\ln(1+e^{-\varepsilon_b})$ up to the sign convention of the kernel. With $x_a=e^{-\varepsilon_a}$ this becomes an equation of the form $x_a=(1+x_a)^{\pm1}$. Derive the exact plateau equation from `eq-tba` with its signs, state it, and then take the limit. The hint's "$x_a=1+x_a$" must be replaced by the correct equation.

**A10. Integral $N_{ab}=\frac1{2\pi}\int\varphi_{ab}$ (`15-beyond-smatrix.qmd`, `exr-tba-uv` hint).** The value $N_{ab}=\delta_{ab}$ is stated without proof. For $\varphi=-i\frac{d}{d\theta}\ln S$, $N_{ab}$ is the total winding of $S_{ab}$ between $\theta=-\infty$ and $+\infty$ divided by $2\pi$, which depends on the behaviour of $S_{ab}$ at $\pm\infty$. State and check it on one affine Toda S-matrix (for example $a_2^{(1)}$, $S_{11}=\{1\}$).

**A11. Summary sentence on the rigorous construction (`13-simply-laced.qmd`, Section `sec-real-coupling`, "Rigour", and `15-beyond-smatrix.qmd`, Section `sec-rigorous-real`).** The same paragraph (Albeverio and Høegh-Krohn; "extending it to several fields is expected to be routine, but has not been written down") appears twice. Keep one copy and cross-reference it. Also check the claim "the Wightman axioms and a mass gap for a single field in the super-renormalizable range" against Albeverio–Høegh-Krohn (1974), since the single-field construction there is for a different class (the potentials are polynomial-exponential in the real coupling with a smooth cut-off). State the precise hypotheses or label the claim as expected.

## B. Gaps for a reader following Parts I–III

**B1. Mass–coupling relation and bulk energy have no formula (`15-beyond-smatrix.qmd`, Section `sec-mass-coupling`).** The section says the relation "follows exactly from Fateev's results" and gives no formula. Give at least the form of the relation (the mass $m$ as a function of the coupling and the physical mass), or say that no formula is given and point to where it is.

**B2. The floating Coxeter number beyond one loop (`14-non-simply-laced.qmd`).** The text uses $H(\beta)$ and the one-loop value $H=2n+B$. A reader needs one sentence saying how $H$ is determined beyond one loop (by the bootstrap with the given masses, or by the literature), and that the form of $H(\beta)$ is not derived in this book beyond one loop.

**B3. Blocks with numerical arguments.** `13-simply-laced.qmd` uses $(x)$ both as a $\theta$-dependent function and as a number ($(0)=1$, $(h)=-1$, $(x+2h)=(x)$). Add one line: the number $(x)$ is the value of the block at $\theta=0$ or the identity used in the reduction, and state which.

**B4. Terms used without definition.** Each needs a gloss: "Watson's equations" and "cyclicity" (defined in `15-beyond-smatrix.qmd`, Section `sec-form-factors`, but the reader needs the form of Watson's equation with the S-matrix); "kinematic poles" (described, not explained); "$N_{ab}$" and "plateau" in the TBA (see A9–A10); "Wightman axioms" (name only); "ODE/IM" (named, not defined: say that the Schrödinger/Lax connection is the linear problem); "Hermitian analyticity" and "braiding unitarity" are from Part I and are fine.

**B5. Form factors of Part IV.** The form-factor axioms in `15-beyond-smatrix.qmd` are listed without the residue normalization. The reader cannot reproduce the $(n-2)$-particle residue formula. Add a pointer to Smirnov's book (already cited) as the source of the normalization.

**B6. Lie duality and the non-simply-laced dual pair.** `14-non-simply-laced.qmd` refers to "Lie duality" in the soliton masses (`sec-soliton-masses`, Part III) and to the duality $\hat g\to\hat g^\vee$ of Part II. Both are used without a sentence distinguishing them. Keep the distinction that Part III, chapter 10, makes (the solitons of $X^{(1)}$ have the particle masses of $(X^\vee)^{(1)}$, and the S-matrix duality exchanges $\hat g$ and $\hat g^\vee$), and give the reader the example $c_n^{(1)}$, $b_n^{(1)}$ in one line.

## C. Notation clashes

**C1. The sinh-Gordon/ATFT coupling is written $b$, $\beta$, and $B$.** Chapter 12 uses $b$ for the sinh-Gordon coupling with $b=\sqrt2\beta$ (Section `sec-tree-level`). Chapter 13 uses $\beta$ and $B(\beta)$ (`eq-b-of-beta`), and chapter 14 uses $B$ and $B_{\rm tree}$ (see A3). Part I (chapter 3) uses $\beta_{\rm SG}$ and $B$ for sinh-Gordon. Suggestion: use $B$ for the S-matrix parameter everywhere, and write the tree-level example with $\beta$ and $B$ only.

**C2. $R$.** The circumference $R$ in the TBA (`15-beyond-smatrix.qmd`, Section `sec-tba`) and the R-matrix $R$, $\check R$ (Part II). Rename the circumference $R$ to $L$ or $\ell$ in Part IV, or the R-matrix to $\mathsf R$.

**C3. $\varphi$.** The field $\vec\phi$ (Parts I, III, IV) and the scattering kernel $\varphi_{ab}(\theta)=-i\frac{d}{d\theta}\ln S_{ab}$ (`eq-tba`). Rename the kernel to $K_{ab}$ or $\Phi_{ab}$.

**C4. $H$.** The floating Coxeter number $H$ (chapter 14) and $H'$ (the second Coxeter number for $f_4$, chapter 14), and the Hamiltonian $H$ of Part I (chapter 4) and Part VII. Rename the floating Coxeter number to $H_{\rm f}$ or keep $H$ and rename the Hamiltonian in Part VII.

**C5. $\varepsilon$.** The TBA pseudo-energy $\varepsilon_a$ (`eq-tba`) and $\hat\varepsilon$ of Part I (chapter 2 exercises, where it is a symbol in the Bethe–Yang discussion, if any). Check before renaming.

**C6. $c$.** The central charge $c_{\rm eff}$ (chapter 15) and the algebra family $c_n^{(1)}$ (chapters 13–14), and the Part II use of $c$. Low priority; keep, but note in the text that $c_{\rm eff}$ is a central charge.

**C7. $m$ and $m_a$.** The Lagrangian parameter $m$ (Part III) and the particle mass $m$ (chapter 12 and 13). Consistent with Part III ($m_a=2m\sin\frac{\pi a}{h}$) but the sinh-Gordon example in chapter 12 uses $m$ for the particle mass and $m^2/b^2$ in the potential, which is the Lagrangian parameter. State the convention at the start of chapter 12.

## D. Items checked and found correct (do not re-check unless you change them)

- Chapter 12: $n$-point couplings $C_{a_1\dots a_n}=m^2\beta^{n-2}\sum_jn_j\prod_k(\vec\alpha_j\cdot\vec e_{a_k})$ reproduce the cubic coupling of Part III for $n=3$ and the mass matrix for $n=2$.
- Normal ordering of $e^{\beta\vec\alpha\cdot\vec\phi}$ and the power count $\int d^2k/(k^2-m^2)^n$, convergent for $n\ge2$.
- The S-matrix normalization $S=1+i\mathcal M/(4m_am_b\sinh\theta)$ (`eq-s-from-m`); the $\varphi^4$ example: $\mathcal M=-m^2b^2$, $S=1-ib^2/(4\sinh\theta)$, which agrees with the expansion of the exact sinh-Gordon $S$ with $B=\frac{b^2/4\pi}{1+b^2/8\pi}$.
- Tree-level exchange: $\mathcal M\supset-C^2/(s-m_c^2)$, the expansion $s-m_c^2\simeq2im_am_b\sin u\,(\theta-iu)$, and the residue $\frac{i C^2}{8m_a^2m_b^2\sin^2u}$. The area rule gives $2m_am_b\sin u=4\Delta$ and residue $i\beta^2/(2h)$. Checked for $a_2^{(1)}$: $|C_{111}|=3\beta$, residue $i\beta^2/6$ with $m=\sqrt3$, $u=2\pi/3$.
- Coleman–Thun power count: $P=4$, $L=1$ gives a double pole; $P=1$, $L=0$ gives a simple pole.
- Blocks: $(x)(-\theta)=(x)^{-1}$, $(h)=-1$, $(x+2h)=(x)$, $(x)(i\pi-\theta)=-(h-x)(\theta)$ (`exr-block-identities`), derived by direct computation.
- $\{x\}$: unitarity $\{x\}(\theta)\{x\}(-\theta)=1$; crossing $\{x\}(i\pi-\theta)=\{h-x\}(\theta)$; free limit $B=0$; $B\to2-B$ invariance; pole and zero positions of $\{x\}$ in the physical strip.
- $B(\beta)=\frac{1}{2\pi}\frac{\beta^2}{1+\beta^2/4\pi}$ (`eq-b-of-beta`): $B(4\pi/\beta)=2-B(\beta)$ (`exr-dual`, checked directly: both sides equal $8\pi/(4\pi+\beta^2)$); the self-dual point $\beta^2=4\pi$ gives $B=1$; range $[0,2)$.
- $a_2^{(1)}$: $S_{11}=\{1\}$ has its only pole at $2\pi i/3$; $S_{12}=\{2\}$ has its pole at $i\pi/3$; the bootstrap holds (`exr-a2-bootstrap`).
- Dorey's formula reproduces $\{1\}$ and $\{2\}$ for $A_2$ after reduction (see A1).
- $e_8^{(1)}$: $S_{11}=\{1\}\{11\}\{19\}\{29\}$ has poles at $2\pi i/3$ ($\{19\}$), $2\pi i/5$ ($\{11\}$) and $i\pi/15$ ($\{1\}$). These give bound states with masses $m_1$, $m_2=2m_1\cos\frac\pi5$ and $m_3=2m_1\cos\frac\pi{30}$. The third ratio is the one I flagged as unverified in the Part III brief (A5 there); it is now confirmed by these fusing angles, so that item in the Part III brief can be closed.
- Ch. 14 $c_n^{(1)}$: $m_a\propto\sin\frac{\pi a}{H}$ is exact by construction; the one-loop expansion of $\sin\frac{\pi a}{H}/\sin\frac{\pi n}{H}$ has coefficient $-\frac{\pi a}{4n^2}\cot\frac{\pi a}{2n}$ in $B$ (`exr-cn-one-loop`).
- TBA equations (`eq-tba`) and $E_0(R)$ have the standard form; Watson's equation and cyclicity give $F_{\min}(i\pi-\theta)=F_{\min}(i\pi+\theta)$ (`exr-watson`).
- $e_8$ and $c_n$ masses are the Perron–Frobenius values in Part II and III.
- All 35 `@` labels and all 31 citation keys used in Part IV resolve.
- `textbook_code/smatrix_simply_laced.py`: the unitarity, $S_{ab}=S_{ba}$, $B\to2-B$ and crossing checks PASS for $e_8,e_7,d_4,d_6,a_4$; $e_8$ $S_{11}=\{1\}\{11\}\{19\}\{29\}$ PASS; the $a_n$ closed form satisfies the bootstrap for $n=2..5$ PASS; the Dorey formula equals the $a_n$ closed form for all $a,b$ PASS. The bootstrap at every fusing of $e_7$ and $e_8$ passes: 112 fusings for $e_7^{(1)}$ and 224 for $e_8^{(1)}$, matching chapter 13.

## E. Verification checklist for the agent who makes the fixes

1. Re-run the label check on `part4-real-coupling/` (exclude `.claude/` and `_site/`).
2. The script `textbook_code/smatrix_simply_laced.py` was run to completion on 2026-10-10 (exit 0, 11 PASS, 0 FAIL). Its bootstrap checks confirm the fusing counts of chapter 13 (112 for $e_7^{(1)}$, 224 for $e_8^{(1)}$). Rerun it only if the S-matrix code changes.
3. For A2, recompute the endpoints from the literature before editing the table.
4. For A9–A10, check the TBA plateau equation and $N_{ab}$ with a short numerical computation, not only by hand.
5. For each notation rename (C1–C7), grep the whole book (excluding `.claude/` and `_site/`) for the old symbol in the same meaning.
6. Report what changed, what was left, and anything you could not verify.

## Resolution (2026-10-10)

- **A1.** Done (ch. 13). Exponents are half-integers that combine into integer powers; $A_2$ example added. Verified with a script (Dorey formula exponents for $A_2$: $\tfrac12,0,-\tfrac12$ for $S_{11}$).
- **A2.** Done (ch. 14). The brief's diagnosis was right: the old table mixed $B_{\rm book}=\beta^2/2\pi$ (range beyond 2) with ranges valid for the S-matrix parameter $B\in[0,2]$. Now: $B\in[0,2]$ is always the block parameter, $H=h+\tfrac12(h_\vee-h)B$ (new `eq-floating-h`), and a new column gives $B$ in terms of the book coupling, $B(\beta_*)$ with $\beta_*^2=\beta^2/k$ for $c_n^{(1)},g_2^{(1)},f_4^{(1)}$ and $\beta_*=\beta$ for $b_n^{(1)}$; one-loop column $2n+\beta^2/4\pi$, $2n-\beta^2/4\pi$, $6+\beta^2/2\pi$, $12+3\beta^2/4\pi$ (unchanged from before, and consistent with Part VI chs. 23–26). Endpoints $h_\vee=2n+2,2n-1,12,18$ explained as $kh$ (resp. $h$ for $a_{2n-1}^{(2)}$), not $h^\vee$ of `tbl-foldings`. Duality paragraph gives the twisted members' $B$: $2-B=B(\beta^\vee_*)$ with $\beta^\vee_*=\beta^\vee$ for $d_{n+1}^{(2)},d_4^{(3)},e_6^{(2)}$ and $\beta^{\vee2}/2$ for $a_{2n-1}^{(2)}$. Literature: DGZ (hep-th/9201067) read in full text: their $a_{2n-1}^{(2)}$ Lagrangian has longest root squared length 4, $b_n^{(1)}$ and $d_{n+1}^{(2)}$ ("$d_n^{(2)}$", obtained from $b_n$ by dropping particle $n$) have longest roots of squared length 2, $g_2^{(1)}$ has long roots of squared length 2 with $H=6+B'$; $c_n^{(1)}$: $H=h+B$, "too many higher order poles unless $H=h\pm B$"; §9 states $c_n\leftrightarrow d^{(2)}$ and $a_{2n-1}^{(2)}\leftrightarrow b_n$ under $\beta\to4\pi/\beta$. CDS (hep-th/9304065): $H=6+3B$, $12+3B$ with $0\le B\le2$, "not intended to imply $B$ has the form (1.4)", so CDS do not normalize $B$. All four pairs are consistent with the coroot form of strong–weak duality; ch. 24 needed no change. Exercise `exr-floating-endpoints` updated.
- **A3.** Done (ch. 14). $B_{\rm tree}$ (ch. 26) is the value of $B$ fixed by the tree amplitude, $=\beta^2/4\pi$ at leading order for $f_4^{(1)}$; the old "$\beta^2/4\pi$" was right but now stated with its normalization.
- **A4.** Done (ch. 13 summary title and paragraph; part page wording).
- **A5.** Done (ch. 13).
- **A6.** Done (ch. 12): $C_{111}$ in the basis of `exr-shg-cubic`, $B$-expansion stated. The "stray brace" is not stray: $\{1\}$ is the block. Residue $i\beta^2/6$ of $\{1\}$ at first order in $B$ rechecked by hand.
- **A7.** Done (ch. 12). "each" removed. CDS p. 3 confirms: detailed checks for second and third order poles [BCDS 1991], contributing diagrams identified for each higher pole up to order twelve in $e_8^{(1)}$. BCDS 1991 itself not accessed.
- **A8.** Done (ch. 15). TBA gives $c_{\rm eff}=r$; UV theory is Toda CFT with $c>r$.
- **A9.** The brief's claim is wrong: with the signs of `eq-tba` the plateau equation is $x_a=\prod_b(1+x_b)^{N_{ab}}$, which for $N=\delta$ is exactly $x_a=1+x_a$. The hint now derives it explicitly. Numerical $a_2^{(1)}$ TBA ($B=0.6$): $c_{\rm eff}=1.856,1.953,1.986$ at $mL=10^{-2},10^{-4},10^{-8}$, $\varepsilon(0)\to-\infty$ slowly.
- **A10.** Done (hint of `exr-tba-uv`). Winding computed numerically: $N(\{1\})=1$, $N(\{x\})=0$ for $2\le x\le h-1$ ($h=3$, $h=30$, several $B$), so $N_{ab}=\delta_{ab}$ for $a_n^{(1)}$; stated for $a_n$ only.
- **A11.** Done. Ch. 13 now cross-references `sec-rigorous-real`. Ch. 15 states AHK's hypotheses as recalled (single field, positive superposition of Wick-ordered exponentials with $|\alpha|$ below a bound, plus a free field of positive mass) and labels the Toda case as expected. The AHK paper could not be accessed; the precise $\alpha$-range is deliberately not quoted.
- **B1.** Done: form of the relation described, formula not reproduced, pointers given.
- **B2.** Done (ch. 14, "Beyond one loop" paragraph).
- **B3.** Done: the "numbers" are identities between functions of $\theta$ ($(0)\equiv1$, $(h)\equiv-1$), not values at a rapidity.
- **B4.** Done: Watson's equation and cyclicity written out, kinematic poles explained, Wightman axioms and ODE/IM glossed, plateau and pseudo-energy defined.
- **B5.** Done: residue axiom for one self-conjugate particle, normalization pointer to Smirnov.
- **B6.** Done (ch. 14, $c_n^{(1)}$ example). Also fixed "Lie dual" in ch. 13 `sec-strong-weak`, which meant $\hat g^\vee$.
- **C1.** Done: ch. 12 sinh-Gordon example rewritten with $\beta$ and $B$ only ($b$ removed).
- **C2.** Done: TBA circumference renamed $L$ (matching Part I's Bethe–Yang), Euclidean length $L'$.
- **C3.** Done: kernel renamed $K_{ab}$ (no other use in the book).
- **C4.** Note added at the first occurrence (ch. 12); no rename (floating $H$ is used throughout Part VI).
- **C5.** No change: no $\hat\varepsilon$ in Part I; the $\varepsilon$'s of Part VII are in other chapters.
- **C6.** Done: one-line note that $c_{\rm eff}$ is a central charge.
- **C7.** Done: convention sentence after `eq-n-point-couplings`; sinh-Gordon example now uses $m_1=2m$.
- Label check: 0 unresolved in `part4-real-coupling/`. Scratch render: exit 0, no `?@` in Part IV HTML.
