# Referee prompt: half-line open problem 12b
Tier: B
Brief: briefs/halfline_open_12_rank_one.md (part (b) only)
Report: briefs/halfline_open_12b_report.md
Scripts: briefs/halfline_open_12b_scripts/. Run them from a scratch copy.
- **Imports.** `p12b_common.py` is shared. Three scripts insert `/home/gustav/Git/affine_toda/textbook_code` into `sys.path`; that is read-only use, and acceptable.
- **Runtimes.** `p12b_an_bootstrap.py` takes about 26 minutes. The others take seconds to minutes.
- **Usage.** `python3 p12b_spectrum_poles.py <lambda> <eta/eta_N> <maxdepth>`.

Scope: part (b) of problem 12, the boundary sine-Gordon excited states, plus proposed text for `prp-bsg-excited` and `prp-bsg-spectrum`. The worker was Claude Sonnet. Part (a) and the lead are tag 12a. Tier B because the report claims analytic proofs for all $n$ (claims 1–3) and finds gaps and misprints in BPTT, BPT and Mattsson–Dorey (claims 9–11).

## Claims to referee

1. [Theorem] **The functional equation.** The bulk amplitude $a=S_{++}$ (BPTT 3.1, the book's $S_0$) obeys $a(u-\pi/\lambda)=\tan\frac u2\tan(\frac\pi{2\lambda}-\frac u2)\,a(u)$. The proof shifts $v=\lambda u/\pi$ by 1 in the Γ-product. Consequence: $a(u+x-\frac\pi\lambda)a(u-x+\frac\pi\lambda)/(a(u+x)a(u-x))=\{2\lambda x/\pi-1\}$.
2. [Theorem] **The two forms of $a_n$ agree for every $n$.** The ratio form telescopes to the product form $\prod_{l\le n}\{2(\eta/\pi-l)\}$, because $\nu_n=\nu_{n-1}-\pi/\lambda$.
3. [Theorem, conditional] **The $n$-th bootstrap follows from the $n=0$ one.** $a(u-\nu_n)Q(\eta)b(u+\nu_n)=Q(\bar\eta)a_n$ holds for all $n$ once the $n=0$ identity holds, because $\lambda\nu_n=\lambda\nu_0-n\pi$ makes the sines in $b/a$ cancel. The $n=0$ identity itself is numerical.
4. [Numerical] **The iterated scalar bootstrap** reproduces BPTT §3.8's closed forms (even and odd depth) for depth ≤4 and labels ≤4: 54 cases, 1e-19. $a_n$ is unitary and crossing-symmetric for $n\le7$.
5. [Numerical] **The full 2×2 matrix bootstrap.** It is derived from the residue of the braid-form reflection equation at a pole of $K$. The residue is rank one, and the new reflection matrix is $XK(\bar\eta)X$ on $|n\rangle$, $K(\eta)$ on $|n,m\rangle$, and so on. Coverage: 120 steps ($n\le5$, depth ≤3) at four parameter points, to 1e-14. Negative controls fail as they should.
6. [Numerical] **Unitarity and crossing-unitarity of the full $K^S$** for 14 states of depth ≤4, at 3e-15, so no extra CDD factor is needed.
7. [Numerical] **Breather reflection.** $B_1$ and $B_2$ on 10 excited states give $K_{B_m}(x_S)\prod b_{mn_i}$, generalizing BPTT (3.9). The breather decouples from the soliton multiplet.
8. [Numerical / literature] **Comparison with TCSA.**
   - BPT Table 5.2's 7 predicted values are reproduced, but this tests only $|0\rangle$.
   - At $\lambda=7$, $\eta=0.9\eta_N$, $\vartheta=0$, the rule's allowed states and C-parity sectors agree with BPT (4.6)–(4.8) and the Fig. 4.3 captions, and $|1,1\rangle$ is excluded.
   - The identified-by-TCSA set is listed in the report's table.
   - Zamolodchikov's energy gives $E(\bar\eta,\vartheta)-E(\eta,\vartheta)=M\cos\nu_0$ identically, so BPTT's conjecture $m_{|0\rangle}(\eta)=m_{|\rangle}(\bar\eta)$ is exact.
9. [Literature claim] **The "U-class" poles: a gap in the pole bookkeeping of BPTT, Mattsson–Dorey and the book.**
   - The blocks of $a_m(\bar\eta)$ at even label positions put poles at $\nu_0,\nu_{-1},\dots,\nu_{-m}$, and $\sigma(\bar\eta)$ has b-type poles at $\nu_{-j}$.
   - These lie in the physical strip iff $\eta<\eta_N-j\pi$.
   - They are not explained in the sources, and no TCSA point lies in that region.
   - The worker flags them as "unexplained in the literature, not an inconsistency".
10. [Literature correction] **Mattsson–Dorey (hep-th/0008071) Table 7** prints "$\nu_0$: order $k$". Their own §5.3 example and the worker's count give order 2 for $n\ge1$. MD Tables 6 and 7 list $\nu_{-N}$ only with order 1 for charge 1.
11. [Literature claim] **BPT's text and captions disagree.** The lists (4.7), (4.8) in BPT (hep-th/0108157) omit $|0,1\rangle$, $|0,3\rangle$ and $|1,2\rangle$, which the Fig. 4.3 captions include.
12. [Text] **The proposed replacement text.**
    - `prp-bsg-excited` becomes [Numerical] for $n=0$, with all $n$ following analytically.
    - `prp-bsg-spectrum` becomes [Conjecture] for the rule, with the TCSA support listed.
    - A new paragraph is added on the U-class poles.

## Coordinator notes (suspected weak points)

- **Claim 9 is the one that matters most.**
  - The worker notes that MD dispose of some poles with "the same mechanisms with $\xi\to\pi(\lambda+1)-\xi$". Check whether that substitution, which is the map $\eta\to\bar\eta$ in the book's notation, maps the U-class poles onto poles that are already explained, such as the $w$ tower or $\nu_j$ with $j\ge0$ for the dual parameter. If it does, claim 9 collapses to "explained by MD's duality remark".
  - Check also whether zeros elsewhere in $Q^S$ cancel any U-class poles. The net order counts; recompute net orders independently at two or three (state, position) pairs inside the region $\eta<\eta_N-\pi$.
  - Read Mattsson–Dorey's full Coleman–Thun analysis (Dirichlet case) and BPTT §5 before accepting "not in the sources".
- **Claims 1–3.** The algebra is short, so redo it by hand or in sympy. Check that $b=(\sin\lambda u/\sin\lambda(\pi-u))\,a$ is the correct relation between $S_T$ and $S_0$ in the book's conventions (`sec-bsg-parameters`, the S-matrix chapter), and that $\nu_n=\nu_0-n\pi/\lambda$ matches the book's definition of $\nu_n$.
- **Claim 5.** Check that the residue method is the right derivation of the boundary bootstrap: GZ/BPTT use the boundary-bound-state fusion. In particular, check the claim that the residue is rank one, and that the ordering of $R$ and $K$ (braid form) matches the book.
- **Claim 8, the energy identity.** Verify $E(\bar\eta,\vartheta)-E(\eta,\vartheta)=M\cos\nu_0$ symbolically from BPT (3.4).
- **Claims 10 and 11.** Read MD's Tables 6 and 7, §5.3, and BPT (4.7), (4.8) and the Fig. 4.3 captions in the arXiv sources yourself.
- **The C-parity sector rule** ("parity of the sum of labels") is taken from BPT §4.3.2 and holds only at $\vartheta=0$. Check that no claim uses it at $\vartheta\ne0$.
- **Prior art for claims 1–3.** The all-$n$ telescoping may already be in BPTT, Mattsson–Dorey, or Mattsson's thesis. Check before the book credits it as new.
