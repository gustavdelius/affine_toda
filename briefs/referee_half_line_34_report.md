## Referee report: chapter 34 (`part7-half-line/34-quantum-symmetric-pairs.qmd`), branch claude/affine-toda-boundary-dae69e

All 15 items were addressed; items 9, 14 and 15 only partly (details under each). No repository file was edited. Scripts are in `/tmp/claude-1000/-home-gustav-Git-affine-toda--claude-worktrees-affine-toda-boundary-dae69e/e8814860-aa66-4a52-8bab-4efe4e4077d9/scratchpad/referee_34/`, listed at the end.

**Headline errors:**
- **Physical Kolb parameter.** It is $c^{\rm Kolb}=-q_j^{-1}$, not $-q_j$. This is wrong at line 71, line 152 and table rows 257–263, which also makes the table internally inconsistent.
- **Universal K-matrix (Balagović–Kolb).** They give a $\tau\tau_0$-universal K-matrix in general; (1.1) holds as written only when $\tau\tau_0={\rm id}$.
- **Admissibility.** Its "only if" direction can be proved; text below.
- **Serre obstruction.** I redid the full computation by hand, and it is now a complete proof; text below.
- **Sine-Gordon analogy for (iii).** The submodule of `prp-one-dim-modules` (iii) is not the analogue of the sine-Gordon bound-state submodules. It is the kernel of $K$.
- **"Open: not matched" claims that can now be settled.** RV (D.2) matches `eq-an-k-matrix` exactly for both $s$. The $d_n$ anti-diagonal $K$ is $K=1$ in twisted form. RV Rem. 3.13(ii) asserts the intersection property.
- **New structural fact.** $B_{\pm\epsilon_*}$ is the image of the standard QSP $B_0$ under a character twist.

---

### 1. `prp-serre-obstruction` — **correct; now fully proved by hand**

I computed all weight components for $a_{ij}=-1$ directly, without Kolb (7.10) or Rem. 5.18. The result is exactly (-@eq-serre-obstruction), for arbitrary $c_i,c_j,\hat\epsilon$.

- The bracket vanishes iff $(1-q)(1-q^{-1})\hat\epsilon_i^2=qc_i$. For $c_i=q^{-1}$ this is $\hat\epsilon_i^2=\epsilon_*^2$.
- In Kolb's variables ($c^{\rm K}=-c$, $s=\hat\epsilon$) the bracket is $q[c^{\rm K}-(q^{-1}-1)^2s^2]$. Both forms are correct.
- Kolb's Rem. 5.18 concerns expansions for $s\in\mathcal S$, so it is a weak justification for $s$-independence. The direct computation removes the need for it.

**Replace the proof block, lines 111–118**, with:

> Write $b_i=e_i+\varphi_i+\hat\epsilon_ik_i$ with $\varphi_i=c_if_ik_i$. For $a_{ij}=-1$ we have $k_ie_j=q^{-1}e_jk_i$, $k_i\varphi_j=q\varphi_jk_i$, $\varphi_ie_j=q^{-1}e_j\varphi_i$, and $k_ie_i=q^2e_ik_i$, $k_i\varphi_i=q^{-2}\varphi_ik_i$, $\varphi_ie_i=q^2e_i\varphi_i-q^2D_i$ with $D_i=c_i(k_i^2-1)/(q-q^{-1})$. $F_{ij}(b_i,b_j)$ is linear in $b_j$; take its three terms in turn.
>
> *The term $\hat\epsilon_jk_j$.* $F_{ij}(b_i,k_j)=\big(b_i^2-[2]\,b_i\sigma(b_i)+\sigma(b_i)^2\big)k_j$ with $\sigma(x)=k_jxk_j^{-1}$, so $\sigma(e_i)=q^{-1}e_i$, $\sigma(\varphi_i)=q\varphi_i$, $\sigma(k_i)=k_i$. A product $XY$ of two terms of $b_i$ comes with $1-[2]s_Y+s_Xs_Y$, where $\sigma(Y)=s_YY$. This vanishes for $e_ie_i$ and $\varphi_i\varphi_i$. The pairs $e_ik_i,k_ie_i$ and $\varphi_ik_i,k_i\varphi_i$ cancel after reordering. As an example, the coefficient of $\hat\epsilon_je_i^2k_j$ is $1-[2]q^{-1}+q^{-2}=0$. What remains is
> $$\hat\epsilon_j\Big[(1-q^2)e_i\varphi_i+(1-q^{-2})\varphi_ie_i+(2-[2])\hat\epsilon_i^2k_i^2\Big]k_j=\hat\epsilon_j\Big[-q\,c_i(k_i^2-1)+(1-q)(1-q^{-1})\hat\epsilon_i^2k_i^2\Big]k_j ,$$
> since $(1-q^2)e_i\varphi_i+(1-q^{-2})\varphi_ie_i=c_i(1-q^2)[e_i,f_i]k_i$.
>
> *The term $e_j$.* In $F_{ij}(b_i,e_j)$ the part quadratic in $e_i$ is the quantum Serre relation. The parts quadratic in $\varphi_i$ or in $k_i$, and those linear in $\hat\epsilon_ik_i$, reorder to single monomials whose coefficients sum to zero ($q+q^{-1}-[2]$, $1-[2]q+q^2$, $2+q^2+q^{-2}-[2]^2$). The part linear in both $e_i$ and $\varphi_i$ leaves, after moving $\varphi_i$ to the right, $-q^2D_ie_j+e_jD_i=q\,c_ie_j$.
>
> *The term $\varphi_j$.* In the same way only the part linear in both $e_i$ and $\varphi_i$ survives. The part quadratic in $\varphi_i$ vanishes because $F_{ij}(f_ik_i,f_jk_j)=F_{ij}(f_i,f_j)\,k_i^2k_j=0$. The surviving part gives $D_i\varphi_j-q^2\varphi_jD_i=q\,c_i\varphi_j$.
>
> Adding the three terms gives (-@eq-serre-obstruction). At $s=0$ the weight-$\pm\alpha_j$ part is Kolb's (7.10), $F_{ij}(B_i,B_j)=-q\,c_iB_j$, transported by (-@eq-kolb-dictionary). In Kolb's conventions the result reads $F_{ij}(B_i,B_j)+q\,c_iB_j=q\,s_j\big(c_i-(q^{-1}-1)^2s_i^2\big)K_i^{-2}K_j^{-1}$. Both forms are checked on $V\otimes\bar V\otimes V$ of $a_2$, on $V\otimes\bar V$ of $a_3$ and on $V\otimes V$ of $d_4$, with random $c$ and $\hat\epsilon$, at the physical and a generic $q$ (`code/qsp_checks.py`).

**Line 101:** change **[Proposition, sketch]** to **[Proposition]**.

### 2. `prp-admissibility` — **the sketch is right and can be made a proof of "obstruction ≠ 0 ⇒ no K"**

The argument holds for the vector multiplets of $a_n$ ($n\ge2$) and $d_n$, for $q$ not a root of unity and generic $y$. Keep the [Numerical] label for the full equivalence: the "if" direction rests on `prp-an-k-existence` and `prp-dn-k`, which are themselves numerical or sketched.

**Replace line 137** with:

> Checked in `code/qsp_checks.py`. The "only if" direction can be proved for the vector multiplets of $a_n^{(1)}$, $n\ge2$, and $d_n^{(1)}$, when $q$ is not a root of unity. Suppose the obstruction of the pair $(i,j)$ does not vanish, so that $\kappa=k_i^2k_j\in B$. Since $(2\alpha_i+\alpha_j,\alpha_i)=3$, the three elements $\kappa^rb_i\kappa^{-r}=q^{3r}e_i+q^{-3r}c_if_ik_i+\hat\epsilon_ik_i$, $r=0,1,2$, put $e_i$ and $f_ik_i$ separately into $B$. Hence $[e_i,f_ik_i]=(k_i^2-1)/(q-q^{-1})$ and $k_i^2$ are in $B$. Conjugation by $k_i^2$ multiplies $e_m$, $f_mk_m$, $k_m$ of a neighbour $m$ by $q^{-2},q^2,1$. So $e_m,f_mk_m,k_m^2\in B$, and along the connected diagram this holds for every node. For $a_n$, $K:V(\theta)\to\bar V(-\theta)$ must intertwine every $k_m^2$, which acts on $\bar V$ as $k_m^{-2}$ does on $V$. So $K$ maps the weight-$\mu$ vector to a vector of weight $-\mu$, and for $n\ge2$ the vector representation has no such vector: $K=0$. For $d_n$, $K:V(\theta)\to V(-\theta)$ (@prp-dn-k) commutes with every $k_m^2$, so it is diagonal, $K={\rm diag}(d_\mu)$. Then $K\,yE_m=y^{-1}E_m\,K$ gives $d_{\mu+\alpha_m}=y^{-2}d_\mu$ along every non-zero entry of every $E_m$. The $E_m$ act irreducibly (@exr-generic-irreducibility), so some word in them of length $L\ge1$ has a non-zero diagonal entry, and $K=0$ unless $y^{2L}=1$. The converse is the existence part of @prp-an-k-existence and @prp-dn-k.

### 3. "Not a QSP" — **correct, with a nuance worth adding**

- **Kolb.** His (5.11) reads $\mathcal S=\{s\mid s_i\ne0\Rightarrow i\in I_{\rm ns},\ a_{ij}\in-2\mathbb N_0\ \forall j\in I_{\rm ns}\setminus\{i\}\}$. For $(\emptyset,{\rm id})$ on $a_n$ ($n\ge2$), $d_n$ or $e_n$ every node has a neighbour with $a_{ij}=-1$, so $\mathcal S=\{0\}$.
- **Later literature.** None of it enlarges $\mathcal S$.
  - Appel–Vlaar 2022 §6.8: $\Sigma_q$ is the same set (confirmed).
  - Balagović–Kolb use $\mathcal S$.
  - Kolb–Yakimov suppress $s$ altogether (their §1.1).
  - Regelskis–Vlaar keep Kolb's $\mathcal S$ but treat the q-Onsager constraint explicitly: Rem. 3.13(ii), the third bullet before Rem. 9.2, and App. D. They say such algebras are "inequivalent to any QP algebra" when $I_{\rm nsf}\ne I_{\rm ns}$, which is the case here ($I_{\rm nsf}=\emptyset$). So RV App. D is the later paper that treats $\hat\epsilon_j=\pm\epsilon_*$; the book already cites it.
- **New fact (my derivation; injectivity not checked).**
  - Kolb's Thm 7.1 presents $B_0$ by the relations $F_{ij}(b_i,b_j)=qc_ib_j$ (adjacent) and $[b_i,b_j]=0$ (non-adjacent).
  - So $\chi(b_j)=\beta_j$ is a character of $B_0$ exactly when $\beta\in\{0\}\cup\{\pm\epsilon_*\}^N$ (for $c=q^{-1}$).
  - Since $B_0$ is a left coideal, $({\rm id}\otimes\chi)\circ\Delta$ is an algebra map $B_0\to U_q$ sending $b_j\mapsto b_j+\beta_jk_j$. Its image is $B_{\hat\epsilon=\beta}$.
  - This is exactly the Dobson–Kolb mechanism (arXiv:1804.02912 §3.5, $\varphi_s=(\chi_s\otimes{\rm id})\circ\Delta$), used there for $s\in\mathcal S$. It is also the book's own `prp-boundary-module`.
  - So $B_{\pm\epsilon_*}$ is not a QSP coideal, but it is a homomorphic image (presumably isomorphic) of the standard QSP $B_0$, embedded by a character twist.

**Add after line 89:**

> As algebras they are nevertheless close to $B_0$. By Kolb's presentation [@kolb2014, Thm. 7.1], $b_j\mapsto\beta_j$ defines a one-dimensional representation $\chi$ of $B_0$ exactly when $\beta\in\{0\}\cup\{\pm\epsilon_*\}^{r+1}$. Then $({\rm id}\otimes\chi)\circ\Delta$ maps $b_j$ to $b_j+\beta_jk_j$, so $B_{\pm\epsilon_*}$ is the image of the standard coideal $B_0$ under this character twist. This is the construction by which Dobson and Kolb obtain the nonstandard coideals from the standard ones for $s\in\mathcal S$, and it is @prp-boundary-module read backwards.

(Citing Dobson–Kolb needs a new bib entry: Dobson, Kolb, "Factorisation of quasi K-matrices for quantum symmetric pairs", Selecta Math. 25 (2019) 63, arXiv:1804.02912.)

### 4. Baseilhac–Belliard prior art — **correct, one sentence misleading**

What I verified in their paper:
- Coproduct $\Delta(e)=e\otimes q^{-h/2}+q^{h/2}\otimes e$ and $\Delta(f)=f\otimes q^{-h/2}+q^{h/2}\otimes f$, with the same $q$ as the book.
- With $e^{\rm book}=e\,q^{h/2}$, $f^{\rm book}=q\,f\,q^{-h/2}$, $k=q^h$ one gets $A_j=b_j$ ($c_j=\bar c_j=1$, $w_j=\hat\epsilon_j$). In general $c^{\rm book}=q^{-1}\bar c/c$ and $\hat\epsilon=w/c$.
- Their Prop. 2.1 constraint $w_i\big(w_j^2+c_j\bar c_j/(q+q^{-1}-2)\big)=0$ gives $w^2=\epsilon_*^2$. Their $\rho^0_{ij}=c_i\bar c_i$ matches $F_{ij}=qc^{\rm book}_ib_j$.
- §3.1 gives $\hat\epsilon_j=\pm i/(q^{1/2}-q^{-1/2})$ for $a_n$ ($n>1$), $d_n$, $e_{6,7,8}$, and they also treat non-simply-laced and twisted algebras. Their Prop. 2.1 assumes all $w_i\ne0$; mixed zero patterns (consequence 2) are new in the book.

So `prp-an-k-existence` is a corrected re-derivation:
- of Delius–MacKay §4.2 (the intertwiner analysis, which had the $|\hat\epsilon_i|=1$ slip, already documented at ch. 32 lines 318–319);
- with the correct value first stated by BB, who themselves say they agree with [DeM, DeG].

Ch. 32 line 253 already credits BB and is fine. **Line 130**, replace the last sentence with:

> Their derivation came after the intertwiner analysis of Delius and MacKay [@delius_mackay2003, sec. 4.2], whose condition $|\hat\epsilon_i|=1$ is a slip for this one (@sec-boundary-qg-notes), and it was the first to state the value in the normalization $c_j=q_j^{-1}$. @prp-an-k-existence rederives it from the intertwining equation.

### 5. `prp-hypothesis-ii` — **correct**

- **Their hypotheses.** AV 2203.16503 Thm 6.1.1: a modified (lower) nilpotent subalgebra is generated by $\tilde F_i$ with $\tilde F_i-F_i\in\hat U^{\ge0}=\hat U^+\hat U^0$, for all $i$. Every irreducible finite-dimensional $U_q(L\mathfrak g)$-module is then generically irreducible for the principal grading shift. No parameter condition is needed.
- **The book's coideal is of their type.** $\omega_q(b_j)-f_j=c_je_jk_j^{-1}+\hat\epsilon_jk_j^{-1}\in\hat U^+\hat U^0$, for any $c_j$ (even 0) and any $\hat\epsilon$.
- **Tensor products.** The principal shift is a Hopf automorphism, so $W(z)=V_\mu(a_1z)\otimes V_\nu(a_2z)$ for fixed generic $a_1/a_2$.
- **Burnside and Schur step.** Words in $y^{-1}\pi(b_j)$ have leading parts that are words in $\pi(e_j)$. Those span ${\rm End}\,W$ because $\hat U^+$ acts irreducibly (AV Prop. 6.2.x and Cor. 6.2.5). So a minor is a non-zero Laurent polynomial, and Schur gives ${\rm End}_B=\mathbb C$ off a proper analytic subset. Uniqueness of $K$ is ${\rm Hom}_B$ between two generically irreducible modules, which is at most one-dimensional.
- **Existence.** "Uniqueness, not existence" is right. The ch. 32 sentence (line 115) is correct.

### 6. Dictionary — **$c^{\rm book}=-c^{\rm Kolb}$, $\hat\epsilon=s$ correct; physical value wrong**

$\omega_q(F-cEK^{-1}+sK^{-1})=e-cfk+sk$, so $c^{\rm book}=q^{-1}$ gives **$c^{\rm Kolb}=-q^{-1}$**.
- Check: with $c^{\rm K}=-q^{-1}$, $c^{\rm K}=(q^{-1}-1)^2s^2$ gives $s^2=\epsilon_*^2$.
- With $-q$ it gives $s^2=q^2\epsilon_*^2$. Numerically the obstruction at $s=\epsilon_*$ is $3\times10^{-17}$ for $-q^{-1}$ and $1.1$ for $-q$.
- The script used random $c$, so it never tested this value.

Invariants: the torus maps $(c,s)\to(c/t^2,s/t)$, so $s^2/c$ is invariant; the $\omega_q$ and antipode dictionaries both give $\hat\epsilon^2/c^{\rm book}=-s^2/c^{\rm K}$. All consistent.

Fixes:
- **Line 71:** "…is $c_j^{\rm Kolb}=-q_j^{-1}$."
- **Line 152:** "These are Kolb's parameters $c_j^{\rm Kolb}\ne-q_j^{-1}$".
- **Lines 257–261:** $c=-q$ → $c=-q^{-1}$ (five places).
- **Line 263:** $c_j\ne-q_j$ → $c_j\ne-q_j^{-1}$.

**Also correct line 152 on how $K$ transforms.** $K\mapsto\pi_{\bar V}(t)K\pi_V(t)^{-1}$ holds only when $\prod_jt_j^{n_j}=1$. In general the gauge must be combined with a spectral shift. Example: $t_0=-1$ maps $s=+1$ to $s=-1$. I verified $K_{s=-1}(y)=\lambda^{N/2}D_\lambda K_{s=+1}(\lambda y)D_\lambda^{-1}$ with $D_\lambda={\rm diag}(\lambda^i)$ and $\lambda^N=-1$; this is the same as $\Phi\to\Phi+i\pi/2$, sinh ↔ cosh in ch. 33. Insert:

> …and maps a reflection matrix $K$ to $\pi_{\bar V}(t)\,K\,\pi_V(t)^{-1}$ when $\prod_jt_j^{n_j}=1$. A general $t$ also shifts the spectral parameter, $y\to\lambda y$ with $\lambda^h=\prod_jt_j^{n_j}$. For example $t_0=-1$ maps the $s=+1$ coideal to the $s=-1$ one, and $K_{s=-1}(y)$ is a diagonal gauge transform of $K_{s=+1}(\lambda y)$ with $\lambda^{n+1}=-1$. Physically the two remain distinct, since a shift of rapidity is not a symmetry of the boundary.

### 7. `prp-q-dolan-grady` — **correct**

- **Kolb (7.19)** reads $\ldots=q(q+q^{-1})^2c_0(B_1B_0-B_0B_1)$. Transported by $\omega_q$ with $c^{\rm book}=-c^{\rm K}$, it gives $\rho=qc(q+q^{-1})^2$. At $c=q^{-1}$ and $q=-e^{-i\pi\lambda}$, $\rho=4\cos^2\pi\lambda$.
- **Nested-commutator form.** I checked by hand that the left side equals $[b,[b,[b,b']_q]_{q^{-1}}]$.
- **Onsager limit.** $\rho=4$, and $A=2b$ gives 16 on both sides. The `exr-dolan-grady-limit` identity holds because $B_0^2=1$.
- **BB (1.1)** uses $[X,Y]_q=XY-qYX$; the cubic is identical.
- **Baseilhac–Koizumi 2005** use $[X,Y]_q=q^{1/2}XY-q^{-1/2}YX$, so their $q$ is the book's $q^2$. Their $\rho=(q^{1/2}+q^{-1/2})^2k_+k_-$ has the same structure.
- **BK 2003 (19)** (checked on the rendered PDF): $[\hat E_+\hat E_-,\hat E_-\hat E_+]=c^2(q+q^{-1})^2(\hat E_+^2-\hat E_-^2)$.
  - Their Heisenberg realization satisfies (19) exactly, and also q-DG with $\rho_{\hat E}=c^2(q+q^{-1})^2$. This is the book's $\rho$ after $\hat E=c\hat\epsilon$, because they normalize $\hat Q=c\,b$.
  - But (19) is not implied by q-DG. It fails (residual about 10) on $V(y)$ and on $V(y_1)\otimes V(y_2)$ restricted to $B$, where q-DG holds to $10^{-15}$. So it is a condition tied to their ansatz (18).

**Line 185**, replace "We have not checked their (19) against (-@eq-q-dolan-grady)." with:

> Their (19), $[\hat E_+\hat E_-,\hat E_-\hat E_+]=c^2(q+q^{-1})^2(\hat E_+^2-\hat E_-^2)$, belongs to their ansatz (18) for $K$ and is stronger than the q-Onsager relations: $V(y)$ restricted to $B$ satisfies (-@eq-q-dolan-grady) but not (19). Their Heisenberg realization satisfies both, the q-Dolan–Grady relations with $\rho=c^2(q+q^{-1})^2$, which is (-@eq-q-dolan-grady) with $c_j=q^{-1}$ after their normalization $\hat Q=c\,b$, $\hat E=c\,\hat\epsilon$.

**Line 173:** $\texttt{baseilhac\_koizumi2005}$ (hep-th/0507053) supports "XXZ symmetry". Kolb credits the name to a different paper: Baseilhac–Koizumi, "A new (in)finite-dimensional algebra for quantum integrable models", Nucl. Phys. B 720 (2005) 325, hep-th/0503036. Add that entry, or cite it for "named it".

### 8. `prp-one-dim-modules` — **(i)–(iii) correct; the interpretation at line 150 is wrong**

- **What (iii) is.** At $y_*^{N}=s(-q)^{N/2}$ I found:
  - $K(y_*)$ is unique with rank $N-1$, and its kernel is the submodule $\mathbb C_{\hat\epsilon}$.
  - $\bar V(1/y_*)$ contains no $\mathbb C_{\hat\epsilon}$ ($a_2$, $a_3$, all roots, both $s$).
  - Symbolically ($N=3,4$), $\det K\propto(X^N-s(-q)^N)(X^N-s)^{N-1}$. So $y_*$ is a simple zero of $\det K$.
  - The partner point $X^N=s$ is where $K$ has rank one. For $s=+1$ that is $\theta=+i\pi/2$, where $K\propto$ the all-ones matrix (also true for the sine-Gordon Neumann point).
- **Physical reading.** By unitarity the scalar factor must be singular at $\theta_*$ or $-\theta_*$. The rank-one point at $i\pi/2$ is the Ghoshal–Zamolodchikov structure of a boundary-crossing pole (boundary coupling to a soliton–antisoliton pair). It is not a boundary bound state.
- **Sine-Gordon contrast.** In sine-Gordon the submodules at bound-state poles have $\beta\ne\hat\epsilon$ and lie in the physical strip. Statement (ii), every $(\beta_0,\beta_1)$, is correct.

**Replace line 150 from "The submodule in (iii)…" to the end of the paragraph** with:

> The submodule in (iii) is not the analogue of the sine-Gordon submodules at the bound-state poles, which have $\beta\ne\hat\epsilon$. Here $\beta=\hat\epsilon$, and since $\bar V(1/y_*)$ contains no copy of $\mathbb C_{\hat\epsilon}$, the matrix (-@eq-an-k-matrix) has rank $n$ at $y_*$ with kernel $\mathbb C_{\hat\epsilon}$. Indeed $\det K\propto\big(X^{N}-s(-q)^{N}\big)\big(X^{N}-s\big)^{N-1}$ (checked for $N=3,4$). By unitarity $K$ has rank one at the mirror point; for $s=+1$ this is $\theta=i\pi/2$, where $K$ is proportional to the matrix with all entries 1. This is the structure of a boundary-crossing pole, with the boundary coupling to a soliton–antisoliton pair. Whether the scalar factor has the pole has not been analysed.

**Line 139:** remove the meta reference "The referee of @sec-boundary-qg asked which $\beta$ are possible." Replace with "Which $\beta$ are possible?"

### 9. `tbl-physics-qsp` — **rows correct except the $c$ values (item 6); three "open" entries can be settled or advanced**

- **$d_n$, $\hat\epsilon=0$.** Checked numerically on $d_4$:
  - In twisted form, $K_{\rm tw}:V(y)\to\bar V(1/y)$ is exactly the identity.
  - The book's anti-diagonal $K$ is $J$, the $U_q(\hat g)$-isomorphism $\bar V(1/y)\cong V(1/y)$.
  - So the Satake data are $(\emptyset,{\rm id})$, the same as $a_n$. The anti-diagonal shape is an artefact of writing the target as $V$.

  **Line 261** Satake-data cell: "$(\emptyset,{\rm id})$ on $d_n^{(1)}$; $K=1$ as a map $V\to\bar V$, anti-diagonal after $\bar V\cong V$".
- **RV (D.2).** It equals `eq-an-k-matrix` exactly, for both $s$:
  - Gauge ${\rm diag}(X^i)$ from principal to homogeneous gradation, times $X^{N/2}$.
  - Then $u=sX^N$ and $q_{\rm RV}^{-1/2}=i(-q)^{1/2}$, i.e. $q_{\rm RV}=q^{-1}$.
  - Residual about $10^{-16}$ for $N=3,4,5$; $s=-1$ is $u\to-u$.

  **Line 253**, replace "We have not matched their (D.2) to it in their conventions." with:

  > Their bare (D.2) is exactly (-@eq-an-k-matrix): after the gauge transformation ${\rm diag}(X^i)$ to the homogeneous gradation, it is their matrix with $u=sX^{n+1}$ and $q_{\rm RV}=q^{-1}$, for both signs $s$ (checked for $n=2,3,4$).

  Table line 260: add "matches RV (D.2)".
- **Soliton-preserving row.** Advanced, not settled.
  - AV 2025 §7.7–7.8: $\tau=\eta_0$ gives the standard reflection equation and covers RV A.3a and A.3b.
  - RV (9.3)–(9.5) for A.3 with $M_2=0$ (that is, $r=\ell$) give ${\rm diag}\big(u\tfrac{\lambda\mu u-1}{\lambda\mu-u}\,[i\le\ell],\ 1,\ \tfrac{\lambda\mu u-1}{u(\lambda\mu-u)}\,[i\ge t+1-\ell]\big)$.
  - With the third block empty this is exactly de Vega–González-Ruiz's matrix, with $u=e^{2\theta}$ and $\lambda\mu=e^{2\xi}$; the ratio is $e^{2\theta}\sinh(\xi+\theta)/\sinh(\xi-\theta)$.
  - I did not check which $(\ell,r,t)$ are allowed, nor book conventions.

  **Line 264** status: "shape matches RV A.3 with $M_2=0$ ($\tau=\eta_0$); parameter ranges and book conventions not checked".
- **Line 251.** Add that RV's Thm 9.3 is "based on computational work for low $n$, thus only a conjecture for general $n$" (their Rem. 9.2(i)). At line 11, "classification" → "a (partly conjectural) classification".

### 10. `eq-chevalley-theta` — **correct**

$J E_jJ^{-1}=-E_j$ for $j\ge1$ and $(-1)^nE_0$ for node 0. The torus element that flips every $e_j$ has $\prod_jt_j=(-1)^{n+1}$. For even $n$ this cannot be a diagonal matrix and needs the derivation part (the loop rotation $y\to-y$). So it is a real sign, not a gradation artefact, and it is consistent with "an element of the extended torus".
- $e+f+ih$ is nilpotent: trace 0, determinant 0. The limit $b_j-\hat\epsilon_j\to e_j+f_j\pm i(\hat\epsilon_j/\epsilon_*)h_j$ checks out.
- There is no conflict with classical $C_j=\pm1$: the classical limit is $\beta\to0$ at oscillating $q$, together with the non-analytic $g^{1/2}$ normalization of `prp-dm-first-order`.
- **Line 285:** "$|C_j|=1$" → "$C_j=\pm1$" (ch. 30 table says $C_j=\pm1$).

### 11. Quoted theorems — **one real slip**

- **BK (1.1)** is quoted accurately. But Balagović–Kolb's $\mathcal K=\mathfrak X\xi T_{w_X}^{-1}T_{w_0}^{-1}$ is a **$\tau\tau_0$-universal** K-matrix (their Def. 4.12 and Thm 9.5): $\mathcal Kb=\tau\tau_0(b)\mathcal K$ and $\Delta(\mathcal K)=(\mathcal K\otimes1)\hat R^{\tau\tau_0}(\mathcal K\otimes1)\hat R$. For $(\emptyset,{\rm id})$ on $sl_N$, $\tau_0\ne{\rm id}$.
  - **Line 191:** after "(up to a twist by a diagram automorphism)" add "; when the twist $\varphi$ is non-trivial the first $\hat R$ in (-@eq-universal-k-coproduct) is replaced by $\hat R^\varphi=({\rm id}\otimes\varphi)$-twisted braiding [@balagovic_kolb2019, Def. 4.12]".
  - **Line 201:** "has a universal K-matrix" → "has a $\tau\tau_0$-universal K-matrix [@balagovic_kolb2019, Thm. 9.5], a universal one in the sense of (-@eq-universal-k-coproduct) when $\tau\tau_0={\rm id}$,".
- **AV 2022 and AV 2025** are quoted correctly: Thms 4.2.1, 5.1.2(a)–(c) and 6.1.1; the §6.8 parameter set; the standard/transposed statement of §7.
- **Line 210, Kolb–Yakimov:** "a universal K-matrix without parameter constraints" → "a weakly universal K-matrix without the parameter constraints needed for the bar involution, for standard coideals ($s$ suppressed) [@kolb_yakimov2020, sec. 1.1]".

### 12. Labels

| Result | Current | Verdict |
|---|---|---|
| `prp-admissibility` | [Numerical] | Keep, and add the proof of the "only if" direction (item 2). |
| `prp-serre-obstruction` | [Proposition, sketch] | Upgrade to [Proposition] with the full proof above. |
| `prp-universal-k`, `prp-trig-k` | [Established] | Appropriate for published theorems once `prp-universal-k` is corrected (item 11). They are not verified here beyond checking the statements in the papers. |

### 13. Exercises — **all eight well posed and solvable**

- **`exr-omega-coalgebra`:** fine.
- **`exr-kolb-dictionary`:** $S(e)=-k^{-1}e$ and $S(f)=-fk$ are correct, and $S(B_i)=c_i(e_i-c_i^{-1}f_ik_i+c_i^{-1}s_ik_i)$.
- **`exr-serre-zero-weight`:** coefficients $1-q^2$ and $1-q^{-2}$ confirmed. The hint (palindromic antisymmetry, Kolb Lemma 5.11 step 4) is correct.
- **`exr-classical-limit-eps-star`:** fine apart from $|C_j|=1$ → $C_j=\pm1$.
- **`exr-dolan-grady-limit`:** true, since $B_0^2=1$ in spin ½.
- **`exr-torus-dressing`:** add the condition $\prod t_j^{n_j}=1$, or the spectral shift (item 6).
- **`exr-one-dim-sg`:** the failing relations are $F_{20}(b_2,b_0)=b_0$ and $F_{21}(b_2,b_1)=b_1$.
- **`exr-generic-irreducibility`:** the Burnside determinant argument works, and $E_0,\dots,E_n$ generate all $e_{ab}$.

### 14. Bibliography

- **Correct** (DOIs checked against arXiv metadata where available): `kolb2014`, `balagovic_kolb2019`, `kolb_yakimov2020`, `appel_vlaar2022` (Represent. Theory 26, 764–824, DOI ert/623), `appel_vlaar2025` (DOI JEMS/1686; volume and pages not yet in the entry), `baseilhac_belliard2010` (LMP 93, 213–228), `baseilhac_koizumi2003`, `baseilhac_koizumi2005` (J. Stat. Mech. P10005), `dolan_grady1982`, `onsager1944`, `noumi1996` (title matches arXiv).
- **Not independently verified:** `letzter1999` (J. Algebra 220, 729, 1999) — I believe it is correct.
- **`regelskis_vlaar2016`:** no journal reference exists. arXiv lists only report number DMUS-MP-16/05, and AV 2025 still cites it as arXiv 2016. The `@misc` entry is right.
- **`mackay_short2003`:** add "Erratum ibid. 245 (2004) 425".
- **`noumi1996` attribution (lines 11 and 306):** it is a single-author paper. Write "Noumi [@noumi1996], with Sugitani and Dijkhuizen in companion papers".
- **Line 306, "Kolb … noted that … generalized q-Onsager algebras … are examples":** this is a fair paraphrase of Kolb §1, but it refers to the abstract algebras, i.e. the standard ($s=0$) coideals. Add "(as algebras; the embeddings with $\hat\epsilon_j=\pm\epsilon_*$ are the character twists of sec-letzter-kolb)".

### 15. Open items in `sec-qsp-notes` (line 312) — **partly no longer open**

- **Intersection property.** RV Rem. 3.13(ii) state, without proof ("it appears"), that the BB prescription gives coideals satisfying their (3.16), $B\cap U_q(\mathfrak h)=U_q(\mathfrak h)^{\theta_q}$. Kolb uses $s\in\mathcal S$ only through (5.21), which holds at the q-Onsager point by `prp-serre-obstruction`. Together with the character-twist map $B_0\to B_{\pm\epsilon_*}$, this makes the property very likely; I have not checked it.
- **Universal K-matrix.** Still open: RV Rem. 9.2(ii) defer it, and I found no later construction. The Dobson–Kolb twist $\varphi_s=(\chi_s\otimes{\rm id})\Delta$ is a natural route.
- **$q\to1$ limit.** Still open.
- **Later work found:** Belliard–Fomin, arXiv:1106.1317 (J. Phys. A 45 (2012) 025201), constructs dynamical K-matrices for $O_q(a_2^{(1)})$. Relevant to "Dynamical boundaries"; not read.

**Replace line 312** with:

> - **Are the $s=\pm1$ coideals quantum symmetric pairs?** Not in Kolb's sense, because $s\notin\mathcal S$, but they are images of the standard coideal $B_0$ under the character twist $({\rm id}\otimes\chi)\circ\Delta$ (@sec-letzter-kolb), and they escape the obstruction that motivates $\mathcal S$. Regelskis and Vlaar state that such coideals satisfy the intersection property $B\cap U^0=\mathbb C$ [@regelskis_vlaar2016, Rem. 3.13(ii)], without proof. Whether they admit a universal K-matrix, and what their $q\to1$ limit is, are open [@regelskis_vlaar2016, Rem. 9.2].

---

### Scripts

All are in the scratch directory above; each runs in under 10 s.

| Script | What it checks | Output |
|---|---|---|
| `bk03_eq19.py` | BK 2003 Heisenberg realization satisfies their (19) and q-DG | (19) residual 0; $\rho=c^2(q+q^{-1})^2$ |
| `bk03_eq19_generic.py` | Is (19) implied by q-DG? On $V(y)$ and $V\otimes V$ restricted to $B$ | q-DG $\sim10^{-15}$; (19) residual $\sim10$ |
| `onedim_meaning.py` | $a_2$, $a_3$, both $s$, all roots $y_*$: $\dim{\rm Hom}_B=1$, rank of $K$ is $N-1$; $\mathbb C_{\hat\epsilon}\subset V(y_*)$ but $\not\subset\bar V(1/y_*)$ | as stated |
| `det_an_k.py` | sympy factorization of $\det K$ for $N=3,4$ | $(X^N-s m^{2N})(X^N-s)^{N-1}$ |
| `dn_antidiag2.py` | $d_4$: $K_{\rm tw}=1$; $J$ unique and anti-diagonal; $K_{\rm book}\propto J$ | all True |
| `dn_antidiag.py` (end) | Kolb obstruction at $s=\epsilon_*$ | $3\times10^{-17}$ for $c=-q^{-1}$; $1.10$ for $c=-q$ |
| `rv_d2.py` | RV (D.2) against `eq-an-k-matrix`, $N=3,4,5$, $s=\pm1$ | residual $\le3\times10^{-16}$ |

Fetched for this pass: Kolb–Yakimov 1901.00490, Noumi math/9503224, MacKay–Short hep-th/0104212 and Dobson–Kolb 1804.02912 (`.pdf` and `.txt` in the same directory).

## Integration (2026-10-11)

Integrated in the main session.

**Integrated:**
- **Item 1.** The full hand proof replaces the sketch, and the label is now [Proposition] (changed to [Theorem] when chapter 1 was updated, to match the label list there). Each step was re-derived in the main session: $\varphi_ie_i=q^2e_i\varphi_i-q^2D_i$; the weight-zero bracket $-qc_i(k_i^2-1)$; $-q^2D_ie_j+e_jD_i=qc_ie_j$; $D_i\varphi_j-q^2\varphi_jD_i=qc_i\varphi_j$.
- **Item 2.** The proof of the "only if" direction is added after `prp-admissibility`, with one correction: the element in $B$ is the q-commutator $e_i(f_ik_i)-q^{-2}(f_ik_i)e_i=(k_i^2-1)/(q-q^{-1})$, not $[e_i,f_ik_i]$. "(in the principal gradation)" was added for the $d_n$ step.
- **Item 3.** Character-twist paragraph after the parameter-set bullets, with a new bib entry `dobson_kolb2019` (verified via Crossref). Injectivity is stated as unchecked.
- **Item 4.** Baseilhac–Belliard sentence replaced.
- **Item 6.**
  - $c^{\rm Kolb}=-q^{-1}$ at line 71, in the torus paragraph and in table rows 257–263.
  - The torus paragraph now has the condition $\prod t_j^{n_j}=1$ and the spectral shift.
  - `exr-torus-dressing` was adjusted.
- **Item 7.**
  - The BK 2003 (19) paragraph replaced; marked as a numerical check outside `code`.
  - "Named it" now cites `baseilhac_koizumi2005a`. The referee's arXiv id hep-th/0503036 is wrong; the paper is math-ph/0503036, NPB 720 (2005) 325, verified via the arXiv API.
- **Item 8.** The interpretation of `prp-one-dim-modules` (iii) is replaced: it is the kernel of $K$, with the structure of a boundary-crossing pole. The meta sentence at line 139 is removed.
- **Item 9.**
  - RV (D.2) matched.
  - $d_n$ Satake cell.
  - Soliton-preserving status advanced.
  - RV Thm. 9.3 marked partly conjectural.
- **Item 10.** $|C_j|=1\to C_j=\pm1$ in `exr-classical-limit-eps-star`.
- **Item 11.** $\tau\tau_0$-universal K-matrix (Balagović–Kolb Def. 4.12, Thm. 9.5); Kolb–Yakimov wording.
- **Items 14–15.**
  - Noumi attribution.
  - Kolb's "examples" qualified.
  - MacKay–Short erratum in the bib note (from the arXiv journal_ref).
  - Open-items bullet rewritten, with RV Rem. 3.13(ii) and 9.2.
- **Scripts.** `qsp_checks.py` §6 (23 PASS) adds:
  - the physical Kolb parameter, a control that $-q$ fails;
  - RV (D.2) $=$ `eq-an-k-matrix`, $N=3,4,5$, $s=\pm1$;
  - $d_4$: $K_{\rm tw}=1$ and the anti-diagonal $K$ is $J$;
  - $\det K$ factorization. As a check it needs the monomial prefactor $X^{-N^2/2}$, which the report omits; the zeros are as stated.

**Not integrated:**
- The numerical check of BK (19) (`bk03_eq19*.py`) and the one-dimensional-submodule analysis (`onedim_meaning.py`) are not ported to `code`.
- `appel_vlaar2025` volume and pages are not yet available.
- `letzter1999` was not independently verified.
- Belliard–Fomin (arXiv:1106.1317) was not added.
- The parameter ranges of the soliton-preserving row are not checked.
