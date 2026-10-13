<!-- Opus referee follow-up, 2026-10-12. Model: Claude Opus 5.5, effort max. Prompt: briefs/halfline_open_12b_opus_prompt.md. Combined session with 06_01 (see halfline_open_06_01_referee_opus.md). Saved verbatim by the coordinator (section 12b of the hand-back). Work directory: /tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/ref_12b_0601/ (W). -->

Headline (from the hand-back): Gemini misquoted Mattsson–Dorey. MD never say that negative-index poles stay out of the physical strip; they list the ν_{−N} poles. The literature correction that does stand is smaller: MD's pole-order tables are wrong. All items were checked; only the later-literature search (item 1) is left UNCHECKED. "MD" line numbers refer to bsgdht.tex (arXiv hep-th/0008071, source dated 28 Nov 2001).

# 12b

## Item 1 — U-class poles: Gemini's "literature correction" REFUTED; a narrower correction stands (pole orders); poles "unexplained but consistent"

**What MD actually say.**
- Gemini's "md.txt l.1298" is MD bsgdht.tex l.1876–1878: "(For the charge 0 case, there are poles at w_x for x ≤ 0, but none of these are in the physical strip.)"
- That sentence is about w_x, x ≤ 0, on charge-0 states only. It is true in MD's domain 0<ξ<π(λ+1)/2 (l.372), because w_0 = π−ξ/λ+π/(2λ) > π/2 there.
- MD do list the ν_{−N} poles:
  - l.984–986: the barred P^± have "common simple poles at ν_0 and ν_{−N}";
  - Table 7 row 1 (l.1884): "order 1 at ν_{−1}, ν_{−2}, …".
- BPTT list them too: bsgrevision.tex l.437–440, "poles of σ(η̄) … of the form (spole.b) are located at ν_k, k=0,−1,…".

**Exact parameter map.**
- ξ_MD = 4πφ_0/β (l.306) equals the book's η (29-boundary-sine-gordon.qmd l.213). λ is the same in both (4π²/β² = π(λ+1)/2).
- ν_n and w_n are defined identically (MD l.859, 990; BPTT l.437).
- MD's domain coincides with BPTT's fundamental range 0<η≤η_0=η_N (BPTT l.352).
- The U-region π/2 < η < η_N − jπ lies inside this domain. For j=1 it is non-empty when λ>2.

**Bar map (Gemini's claim CONFIRMED).**
- ν_j(ξ̄) = π+π/λ−ξ/λ−(2j+1)π/(2λ) = w_j(ξ) exactly, for all integers j. So ν_{−j} ↔ w_{−j}.
- MD's charge-0 analysis gives no mechanism for w_{−j}, since it puts them outside the strip. Their transfer "same mechanisms with ξ→π(λ+1)−ξ" (l.1977–1979) therefore does not cover ν_{−j}.
- Structural fact for the book: ν_{−j} = π − w_j. The a_m(η̄) factors are crossing-symmetric, so each pole at w_j has a mirror at π − w_j. The U-poles are the crossing images of w_j poles once w_j > π/2, i.e. η < η_N − jπ.

**The real error: MD's pole-order tables.**
- I recomputed the orders independently from MD's own Γ-product (`ref12b_md_orders.py`): 2 parameter points, 13 states, P^+ and P^−.
- The order is measured as −d log|P|/d log ε between ε=1e-7 and 1e-9 at 40 digits. Every value is an integer to 4 decimals and matches my analytic count.
- General result for P^±_{|1;n_1…n_{2k−1}⟩} (MD's charge-1 states, labels n_1…n_{2k−1}):
  - order at ν_0 is 2k − δ_{n_1,0};
  - order at ν_{−j} is 1 + Σ_{i even} c_j(n_i), with c_j(n) = 2, 1, 0 for j<n, j=n, j>n.
- For charge-0 P^−: ν_0 has order 2k − δ_{n_1,0}, and ν_{−j} has order −1 + Σ c_j.
- Examples: |1;0,5,1⟩ has order 3 at ν_{−1}; |1;1,6,2⟩ has 4 at ν_0 and 3 at ν_{−1}.
- So MD Table 6 (l.1855–1870) omits the poles at ν_0 and ν_{−j}. Table 7's "ν_{−N}: 1" fails from level 3. Its "ν_0: k" is item 3.

**Other literature.**
- BPTT §3.5 (l.445–449) explains "the ν-type poles" by Fig. 3.2(c), order 2, "except ν_n", order 1. The simple poles at ν_{−j}, j≥1, are covered by neither clause, yet BPTT's Conclusions (l.876–879) claim every pole is explained.
- BPT does not treat these poles: all its TCSA points lie outside the region.
- Mattsson's thesis repeats MD (thesisht.tex l.4652–4654).
- Later papers: UNCHECKED.

**First look at a Coleman–Thun explanation (not settled).**
- A crossed soliton emission at ν_{−j} from a state S needs an intermediate state of energy e_S + M cos w_j. That is the formal state |S, j⟩, which the ordering rule excludes because w_j > π/2.
- Example: e_{|0⟩} − M cos ν_{−1} = e_{|⟩} + M cos ν_0 + M cos w_1 = e_{|⟩} + m_{B1} sin(η/λ). This is the energy of |0,1⟩, the B_1 bound state, which exists only for η > λπ/2 > η_N − π.
- The extended Fig. 3.2(c) would have order 2 against the measured 1. That is the worker's claim; I did not check the diagram's geometry.

**Proposed book text** (paragraph after `prp-bsg-spectrum`):
> "The closed forms for all states follow from the bootstrap and are unitary and crossing-unitary, but their pole content is not fully explained. The factors $a_m(\bar\eta,u)$ at even label positions are crossing symmetric, so each pole at $w_j$ comes with one at $\pi-w_j=\nu_{-j}$; with the poles of $\sigma(\bar\eta,u)$ at $\nu_0,\nu_{-1},\dots$ [@bptt2002, sec. 3.5] this puts poles at $\nu_{-j}$ ($j\ge1$) in the physical strip whenever $\frac\pi2<\eta<\eta_N-j\pi$. On $|n_1,\dots,n_k\rangle$ with $k$ odd their order is $1+\sum_{i\ \rm even}c_j(n_i)$, $c_j(n)=2,1,0$ for $j<n$, $j=n$, $j>n$. Mattsson and Dorey list them with order 1 and Bajnok et al. list those of $\sigma(\bar\eta,u)$, but neither gives a Coleman–Thun diagram; the obvious candidate, soliton emission at $\nu_{-j}$, would need the state $|\dots,j\rangle$, which the ordering rule excludes. No TCSA data exist in this region. These poles are unexplained but consistent."

For Appendix D: list MD's Tables 6 and 7 pole orders (the formulas above) as a literature correction.

## Item 2 — BPT's lists: claim 11 CORRECTED (Gemini right on content; should have been PARTIAL)
- (4.8) is at uv_ir3.tex l.1118–1121 and contains |1,2⟩. Only |0,1⟩ and |0,3⟩ are missing; (4.7) at l.1115–1118 is complete for two labels.
- Both missing states are TCSA lines in the odd caption (l.1157–1158, "∘ stand for |0,1⟩, |0,3⟩ and |1,2⟩").
- Both are allowed by the rule (ν_0 > w_1 and ν_0 > w_3 from (4.6), l.1109) and lie in the odd sector by label-sum parity.
- |0,1⟩ is the lightest odd-sector state: E/M = 0.445, below |1⟩ at 0.588 (`ref12b_bpt_lists.py`).
- Further omissions, not identified in TCSA: the allowed three-label states |0,2,2⟩ and |0,3,3⟩ (even) and |0,2,3⟩ and |1,3,3⟩ (odd).
- Every state BPT do list is allowed and in the right sector, and |1,1⟩ is excluded.
- Book text: "BPT's lists (4.7)–(4.8) omit |0,1⟩ and |0,3⟩, which their Fig. 4.3 caption identifies as TCSA lines (and four allowed three-label states)."

## Item 3 — MD Table 7 "ν_0: order k": the worker's claim is CONFIRMED; the entry is an error
- k is defined by "a state with 2k−1 indices" (l.1874–1875). The true order is 2k − δ_{n_1,0}.
- MD's own §5.3 (l.1494–1497) gives a double pole at ν_0 on |1;1⟩ (k=1). My measurements give 2 on |1;1⟩ and |1;2⟩, 4 on |1;1,6,2⟩ and 3 on |1;0,6,2⟩.
- No reading of k fits these. "k" is right only for |1;0⟩.
- Mattsson's thesis (l.4661) prints "ν_0, ν_{−1}, …: 1" instead, which is also wrong.
- Correct fix: drop the separate ν_0 row and start the "2k" row at ν_0; when n_1 = 0, ν_0 = ν_{n_1} has order 2k−1. Table 6 needs the same fix.

## Scripts (12b part)
- W/ref12b_md_orders.py and .out
- W/ref12b_bpt_lists.py and .out
- Sources: W/dl_md/x/bsgdht.tex, W/dl_thesis/x/thesisht.tex, W/dl_bptt/bptt_utf8.tex (from W/dl_bptt/x/bsgrevision.tex), W/dl_bpt/bpt_utf8.tex (from W/dl_bpt/x/uv_ir3.tex)
