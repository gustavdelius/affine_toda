# Opus referee prompt: half-line open problem 12b (tier-B follow-up)

Written 2026-10-12 by the Claude coordinator, after reading the Gemini verdict `briefs/halfline_open_12b_referee.md`.

**What the Gemini pass established.** It independently established the following, so do not redo them:
- claims 1–3 (the functional equation, the telescoping and the sine cancellation; sympy);
- the net pole orders 1, 2, 3 at $\nu_{-1}$ on $|0\rangle$, $|0,5\rangle$, $|0,5,1\rangle$ (its own code);
- the energy identity $E(\bar\eta)-E(\eta)=M\cos\nu_0$.

**Gemini's own Opus items.** Its item 1 (an analytic proof of the $n=0$ base identity) is research. Its item 2 is folded into item 1 below only as a first look.

Read:
- the brief `briefs/halfline_open_12_rank_one.md`, part (b);
- the report `briefs/halfline_open_12b_report.md`;
- the Gemini verdict;
- the scripts in `briefs/halfline_open_12b_scripts/`;
- the papers, fetched from arXiv yourself: Mattsson–Dorey hep-th/0008071, BPTT hep-th/0106070, BPT hep-th/0108157. Also Mattsson's Durham thesis, if accessible.

## Items

1. **The U-class poles: a literature correction to Mattsson–Dorey?** Gemini reports that MD state, near line 1298 of their text, that poles with negative index never enter the physical strip. It says this fails for $\eta<\eta_N-\pi$, where $\nu_{-1}\in(0,\frac\pi2)$ and the net pole orders on $|0\rangle,|0,5\rangle,|0,5,1\rangle$ are 1, 2, 3.
   - Settle whether MD's statement is false *within MD's own parameter domain*. Map their $\xi$ (Dirichlet, $\phi_0$) to the book's $\eta$ exactly, and check whether MD restrict $\xi$ to a range, for example $\xi<\frac{\pi(\lambda+1)}2$ or a "fundamental range" in which the statement holds.
   - Check Gemini's claim that $\eta\to\bar\eta$ maps $\nu_{-j}\to w_{-j}$, and that this does not land on an explained pole.
   - Check whether BPTT §5, BPT, or later papers (Bajnok and collaborators on boundary sine-Gordon, Mattsson's thesis) treat these poles.
   - Take a first look at whether a boundary Coleman–Thun diagram could produce a pole of the measured order. A full resolution is not required.
   - Conclude with exact wording for the book: "literature correction", "unexplained but consistent", or "explained by …".
2. **Claim 11, BPT's lists.** The worker says (4.7) and (4.8) omit $|0,1\rangle$, $|0,3\rangle$ and $|1,2\rangle$. Gemini says $|1,2\rangle$ *is* in (4.8), so only $|0,1\rangle$ and $|0,3\rangle$ are missing, yet it marked the claim CONFIRMED. Read BPT's TeX source and settle this. Also check whether the omission matters: are those states in the captions as TCSA lines, and are they allowed by the rule?
3. **Claim 10, MD Table 7.** Quick check only. MD's Table 7 says "order $k$" at $\nu_0$, while their §5.3 text has a double pole for $k=1$. Is the table entry a misprint, or does "$k$" mean something else in their notation?

For each item give the verdict (CONFIRMED / CORRECTED / REFUTED / UNCHECKED), the evidence with TeX line or equation references, and book text. Give script paths.
