# Brief: referee chapter 34 (quantum symmetric pairs and universal K-matrices)

You are a mathematical-physics referee. The repository is `/home/gustav/Git/affine_toda`; the Part VII work is on branch `claude/affine-toda-boundary-dae69e`, a Quarto book on affine Toda field theory. Do not edit any file in it. Work in your own scratch directory, and deliver your full report in your final message.

## Context

Part VII, "Affine Toda theory on the half-line" (`part7-half-line/`, chapters 28–35), was drafted on 2026-10-10 and 2026-10-11. This brief covers chapter 34 (`34-quantum-symmetric-pairs.qmd`, `sec-qsp`), which places the boundary coideals of chapter 32 in the mathematical theory of quantum symmetric pairs. Its sections:
- symmetric pairs and Satake diagrams, and the classical limit of the boundary coideal (`sec-symmetric-pairs`, `eq-chevalley-theta`);
- Letzter–Kolb coideals and the dictionary with the book's coideal (`sec-letzter-kolb`, `eq-kolb-generators`, `eq-omega-map`, `eq-kolb-dictionary`, `prp-omega-dictionary`);
- the Serre relation of the boundary coideal and its obstruction term (`prp-serre-obstruction`, `eq-serre-obstruction`), admissibility (`prp-admissibility`) and one-dimensional modules (`prp-one-dim-modules`);
- the q-Onsager algebra as the sine-Gordon boundary algebra (`sec-q-onsager`, `prp-q-dolan-grady`, `eq-q-dolan-grady`);
- the universal K-matrix (`sec-universal-k`, `prp-universal-k`, `eq-universal-k-coproduct`);
- Appel–Vlaar trigonometric K-matrices and generic irreducibility (`sec-trigonometric-k`, `prp-trig-k`, `prp-hypothesis-ii`);
- Regelskis–Vlaar generalized Satake diagrams and the table of physical boundaries (`sec-generalized-satake`, `tbl-physics-qsp`).

The chapter's main new claims:
- **The $s=\pm1$ coideals are not quantum symmetric pairs in Kolb's sense.** The coideal of $a_n^{(1)}$ ($n\ge2$) or $d_n^{(1)}$ with $\hat\epsilon_j=\pm\epsilon_*$, which belongs to the Neumann boundary for $s=+1$, lies outside Kolb's parameter set. Its Serre obstruction nevertheless vanishes there, and only there (besides $\hat\epsilon=0$).
- **Hypothesis (ii) of `thm-k-reflection` holds for all $\hat\epsilon$**, by Appel and Vlaar's generic irreducibility theorem transported by $\omega_q$, so $K$ is unique whenever it exists.

Conventions (chapter 32):
- $b_j=e_j+c_jf_jk_j+\hat\epsilon_jk_j$ with $c_j=q_j^{-1}$ in the parity-symmetric normalization. It is a left coideal for the Part VI coproduct $\Delta(e)=e\otimes1+k\otimes e$, $\Delta(f)=f\otimes k^{-1}+1\otimes f$.
- $q=-e^{-i\pi\omega}$ (`eq-q-physical`).
- $\epsilon_*^2=1/((1-q)(1-q^{-1}))$ (`eq-an-eps-condition`).
- Kolb's coideals are right coideals for the same coproduct; the map $\omega_q$ ($e\leftrightarrow f$, $k\to k^{-1}$) relates the two.

Find any object with `grep -rn '{#<id>}' part*/`.

Sources:
- Already fetched as text in `scratchpad/papers/` of the session (the path is given in your prompt): Kolb arXiv:1207.6036, Balagović–Kolb 1507.06276, Appel–Vlaar 2007.09218 and 2203.16503, Regelskis–Vlaar 1602.08471, Baseilhac–Belliard 0906.1215, Baseilhac–Koizumi hep-th/0507053 and hep-th/0208005, Delius–MacKay hep-th/0112023, Delius–MacKay–Short hep-th/0109115.
- Fetch with `curl -sL https://arxiv.org/pdf/<id>`: Kolb–Yakimov 1901.00490, Noumi math/9503224, MacKay–Short hep-th/0104212.
- Letzter's 1999 J. Algebra paper has no arXiv version (math/9812044 is a different paper).

Experience in this repository is that every adversarial referee pass finds real errors.

## Already verified (do not redo unless you find a reason)

**By script** (`textbook_code/qsp_checks.py`, 19 checks, all PASS, under a second):
- **q-Dolan–Grady.** The sine-Gordon coideal satisfies $F_{ij}(b_i,b_j)=qc_i[2]^2[b_i,b_j]$ for generic $\hat\epsilon$ and $c$, at physical and generic $q$, on spin $\frac12\otimes\frac12$, $1$ and $\frac12\otimes1$. A control with $\rho\to1.1\rho$ fails. $\rho=(q+q^{-1})^2=4\cos^2\pi\omega$ at the physical $q$. At $q=1$ the Onsager–Dolan–Grady relation holds with $A=2b_0$, $A^*=2b_1$.
- **Sine-Gordon uniqueness.** $K:V(y)\to V(1/y)$ exists and is unique for generic $(\hat\epsilon_0,\hat\epsilon_1)$.
- **Kolb's conventions.** Kolb's $B_i=F_i-c_iE_iK_i^{-1}+s_iK_i^{-1}$ is a right coideal ($a_1$, $a_2$). $\omega_q$ satisfies $\Delta\circ\omega_q=(\omega_q\otimes\omega_q)\circ\Delta^{\rm op}$, and $\omega_q(B_i)=e_i-c_if_ik_i+s_ik_i$, so $c^{\rm book}=-c^{\rm Kolb}$ and $\hat\epsilon=s$. Kolb's (7.19) (q-Onsager) and (7.10) ($a_{ij}=-1$, $s=0$) hold.
- **Kolb's Remark 5.12, made explicit:** $F_{ij}(B_i,B_j)+qc_iB_j=qs_j(c_i-(q^{-1}-1)^2s_i^2)K_i^{-2}K_j^{-1}$ ($a_2$, $a_3$, generic $s$).
- **Serre obstruction, book form.** $F_{ij}(b_i,b_j)=qc_ib_j+\hat\epsilon_j((1-q)(1-q^{-1})\hat\epsilon_i^2-qc_i)k_i^2k_j$ for $a_{ij}=-1$ ($a_2$, $a_3$, $d_4$; random $\hat\epsilon$, $c$; physical and generic $q$). $b_ib_j=b_jb_i$ for $a_{ij}=0$.
- **Admissibility.** "All Serre obstructions vanish" $\Leftrightarrow$ "$K$ exists", for $a_2$, $a_3$, $d_4$ and 8 $\hat\epsilon$ patterns each.
- **One-dimensional modules.**
  - In $V(y)$ there is none for $\hat\epsilon=0$.
  - For $s=\pm1$ there is exactly $\mathbb C_{\hat\epsilon}$, at $y^{n+1}=s(-q)^{(n+1)/2}$ ($n=2,3$).
  - The $a_2$ generalized q-Onsager algebra has one-dimensional modules only for $\beta\in\{0\}\cup\{\pm\epsilon_*\}^3$.
- **Table rows.** For $a_n$, $\hat\epsilon=0$ gives $K=1$ and $s=\pm1$ a unique non-diagonal $K$ ($n=2,3$). For $d_4$, $\hat\epsilon=0$ gives a unique anti-diagonal $K$, and $\hat\epsilon_j=\pm\epsilon_*$ a unique $K$.
- **Classical limit.** $\theta$ fixes $e_j+f_j$ ($a_2$, $a_3$, $a_4$; node 0 needs the loop rotation $y\to-y$ for even $n$).

**By hand (drafting agent):**
- The weight-zero part of `prp-serre-obstruction`.
- The limit $b_j-\hat\epsilon_j\to e_j+f_j\pm i(\hat\epsilon_j/\epsilon_*)h_j$ as $q\to1$.
- Torus dressing maps $B_{c,\hat\epsilon}$ to $B_{c/t^2,\hat\epsilon/t}$.

**Quoted, not checked:**
- `prp-universal-k` (Balagović–Kolb) and `prp-trig-k` (Appel–Vlaar).
- The DMS (2.3) branching-graph rule.
- The soliton-preserving row of `tbl-physics-qsp` (de Vega–González-Ruiz against RV family A.3).

## To be checked

Report on each item: correct / wrong (with the correction) / unclear (with what would settle it). Items are ordered roughly by risk.

1. **`prp-serre-obstruction` and `eq-serre-obstruction`.** The weight $\pm\alpha_j$ terms rely on Kolb (7.10) and his Remark 5.18 plus numerics, not on a hand computation. Redo the full computation for $a_{ij}=-1$ in the book's conventions, or find an independent derivation (Baseilhac–Belliard's relations are one). Check that the obstruction coefficient vanishes exactly at $\hat\epsilon_i^2=\epsilon_*^2$ when $c_i=q_i^{-1}$, and the claimed equivalent in Kolb's variables, $c=(q^{-1}-1)^2s^2$.
2. **`prp-admissibility`, beyond numerics.** The text sketches why a non-vanishing obstruction rules out a non-zero $K$. Is the sketch correct? Can it be made a proof (for instance: $K$ intertwines the coideal into the commutant, so the obstruction $k_i^2k_j$ term must act by zero on $V$, which it cannot)?
3. **The answer "not a QSP".** Is it right that Kolb's parameter set $\mathcal S$ is $\{0\}$ for $(X,\tau)=(\emptyset,\mathrm{id})$ on $a_n^{(1)}$ ($n\ge2$), $d_n^{(1)}$ and $e_n^{(1)}$? Check against Kolb's definition (his (5.11)) and against the later literature (Balagović–Kolb, Kolb–Yakimov, Regelskis–Vlaar Appendix D, Appel–Vlaar), which may have enlarged the parameter set. If any later paper does treat $\hat\epsilon_j=\pm\epsilon_*$, say so.
4. **Prior art: Baseilhac–Belliard.** The chapter credits them with the condition $\hat\epsilon^2=\epsilon_*^2$ (their Prop. 2.1 and sec. 3.1). Check the claim and their conventions: their coproduct, the map $A_i=b_i$, $w=\hat\epsilon$, and which algebras they cover. Is the book's `prp-an-k-existence` (ch. 32) then a rediscovery? If so, check that ch. 32 now says so (one sentence was added after the proof).
5. **`prp-hypothesis-ii`.** The claim is that Appel–Vlaar 2203.16503, Thm 6.1.1 / Cor 6.1.2 (generic irreducibility under any "modified nilpotent subalgebra") applies to every book coideal $b_j=e_j+c_jf_jk_j+\hat\epsilon_jk_j$ after $\omega_q$, with no parameter condition. Check:
   - the hypotheses of their theorem;
   - that the book's coideal is of their type after $\omega_q$;
   - the Burnside and Schur step from irreducibility of $V\otimes W$ under $B$ to uniqueness of $K$;
   - whether "uniqueness, not existence" is stated correctly.
   This also changed a sentence in ch. 32 after `thm-k-reflection`.
6. **`prp-omega-dictionary` and `eq-kolb-dictionary`.** Is $c^{\rm book}=-c^{\rm Kolb}$, $\hat\epsilon=s$, with physical value $c^{\rm Kolb}=-q$, right? The chapter says only $s^2/c$ is invariant across Kolb, Baseilhac–Belliard and Regelskis–Vlaar. Is that consistent with the torus dressing $B_{c,\hat\epsilon}\to B_{c/t^2,\hat\epsilon/t}$?
7. **`prp-q-dolan-grady`.** Check $\rho_i=qc_i(q+q^{-1})^2$ against Kolb Example 7.6 and Baseilhac–Koizumi, including their normalizations. Check the claim that at the physical $q$, $\rho=4\cos^2\pi\omega$, and the Onsager limit with $A=2b_0$. Is Baseilhac–Koizumi 2003 eq. (19) consistent? The drafting agent did not check it.
8. **`prp-one-dim-modules`.**
   - Is (iii) right that $V(y)$ contains $\mathbb C_{\hat\epsilon}$ at $y^{n+1}=s(-q)^{(n+1)/2}$?
   - What is its physical meaning? Is it a pole or zero of the scalar factor of the $s=\pm1$ $K$, or a boundary bound state?
   - Compare with ch. 32's `prp-boundary-module` and the sine-Gordon $\beta_{1,0}$ found there. In sine-Gordon $\beta$ is unconstrained; is the contrast stated correctly?
9. **`tbl-physics-qsp`.** Check every row against its stated support. In particular:
   - **The $d_n$ rows.** $\hat\epsilon=0$ gives an anti-diagonal $K$. Which Satake diagram or diagram involution does that correspond to?
   - **The soliton-preserving row**, which is marked open: can it be matched with RV family A.3?
   - **Regelskis–Vlaar (D.2)**, which they attribute to Gandenberger, against `eq-an-k-matrix`.
10. **`eq-chevalley-theta` and the classical limit.** Is the identification of $\theta$ correct? Is the need for $y\to-y$ on node 0 for even $n$ a gradation artefact or a real sign? The chapter says $e+f+ih$ is nilpotent, so the $s=\pm1$ coideals have no fixed-point classical limit. Is that consistent with the classical boundary conditions $C_j=\pm1$, which do exist (see `exr-classical-limit-eps-star`)?
11. **Quoted theorems** (`prp-universal-k`, `prp-trig-k`, `eq-universal-k-coproduct`). Are they quoted accurately, with the right hypotheses, in the book's left-coideal conventions? Check the coproduct formula of the universal K-matrix especially: a left/right or $\Delta/\Delta^{\rm op}$ slip is the most likely error.
12. **Status labels.**
    - Is `prp-admissibility` correctly [Numerical]?
    - Should `prp-serre-obstruction` be [Proposition, sketch], given that part of it rests on numerics?
    - Are `prp-universal-k` and `prp-trig-k` correctly [Established] although quoted unchecked?
13. **Exercises.** For every `exr-` (8 of them), check that it is well posed and solvable with the text, and that no hint is wrong. Pay particular attention to `exr-classical-limit-eps-star`, `exr-serre-zero-weight` and `exr-generic-irreducibility`.
14. **Bibliography.** Check the new entries in `references.bib`: `letzter1999`, `noumi1996`, `kolb2014`, `balagovic_kolb2019`, `kolb_yakimov2020`, `appel_vlaar2022`, `appel_vlaar2025`, `regelskis_vlaar2016`, `baseilhac_koizumi2005`, `baseilhac_belliard2010`, `mackay_short2003`, `dolan_grady1982`, `onsager1944`. Check authors, titles, journal data and DOIs, and that each citation supports its sentence. Note that `regelskis_vlaar2016` has no journal reference.
15. **Prior art and later work.** Are the open items in the notes (`sec-qsp-notes`) still open?
    - the intersection property of the $s=\pm1$ coideals;
    - their universal K-matrix (RV Remark 9.2);
    - their $q\to1$ limit.
    Search Kolb's and Regelskis–Vlaar's citations, and Baseilhac's later work on generalized q-Onsager algebras.

## Deliverable

One entry per item: a verdict, the justification, and, for errors, the corrected text ready to paste, with file and line. Include any scripts you wrote, either inline or as a list of what they check with their output. Do not mark anything verified on the strength of a citation alone.
