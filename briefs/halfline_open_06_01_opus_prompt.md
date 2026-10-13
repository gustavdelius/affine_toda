# Opus referee prompt: half-line open problems 06 and 01 (tier-B follow-up)

Written 2026-10-12 by the Claude coordinator, after reading the Gemini verdict `briefs/halfline_open_06_01_referee.md`.

**What the Gemini pass established.** It independently established the following, so do not redo them:
- the energy conversion and the BP translation for $a_2$;
- the bound $E_P\ge-4\#\{C_j=+1\}$ and its saturation;
- the DG (5.7) agreement, with its own shooting solver at $5\times10^{-13}$ on six $a_3$/$a_4$ vacua;
- the $a_3$ alternating result.

**Gemini's own Opus items** (perturbation theory about singular backgrounds; an asymptotic analysis of the $a_4$ exceptional background) are research problems, not referee items.

Read:
- the briefs `briefs/halfline_open_06_missing_vacua.md` and `briefs/halfline_open_01_negative_product_vacua.md`;
- the report `briefs/halfline_open_06_01_report.md`;
- the Gemini verdict;
- the scripts in `briefs/halfline_open_06_01_scripts/`;
- Bowcock–Perkins hep-th/9909174. Fetch the TeX source, not a PDF text dump.

## Items

1. **The $a_5$ diagonal count: 12/9 (worker) or 6/15 (Gemini)?**
   - The worker's table marks the orbit of $C=(++----)$ (BP $A=--++++$, singular $2\bar2$) as diagonal, "= DG $b=2$". Gemini reads BP's Table 2 as grouping four rows (`--++++`, `--++-+`, `--+-++`, `--+--+`) under one "not diagonal" label. The worker said it did not recompute BP's $K$ entries.
   - Settle this by computation. Construct the singular $2\bar2$ background for $(++----)$ as the limit of regular configurations (see `exc4*.py`, `t7.py`, `tauvac.py` for the method). Compute the linearized $K$ as $\tau_j(0)\to0$, with `linrefl2.py` and with an independent check.
   - Then read BP's Table 2 in the TeX source and say which reading of the layout is right.
   - Report the corrected counts for `tbl-missing-vacua`.
2. **Do the proposed book texts state the singular nature honestly?** Read what BP say about these boundary-singular backgrounds: do they call them vacua, do they discuss regularity, and do they compute $K$ about them as the limit of regular configurations? Then judge the worker's proposed text for `sec-boundary-missing-charges`, `tbl-boundary-types` and `sec-real-coupling-vacua`.
   - The book must not say that DG's amplitudes are "confirmed" by a configuration whose field diverges at the wall without saying so.
   - Propose exact corrected wording and labels, for example "[Numerical] energy infimum, not attained".
3. **Claim 8: "$\prod C=-1$ vacua are not limits of regular solutions of the same pattern".** The worker's argument is that the regular neighbours of the endpoint have a fixed product of the signs at $j\pm1$, and flipping one gives $\prod C=-1$. Check the argument. Gemini confirmed this claim without examining it.

For each item give the verdict (CONFIRMED / CORRECTED / REFUTED / UNCHECKED), the evidence, and book text. Give script paths.
