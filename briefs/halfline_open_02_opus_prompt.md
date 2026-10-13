# Opus referee prompt: half-line open problem 02 (tier-B follow-up)

Written 2026-10-12 by the Claude coordinator, after reading the Gemini verdict `briefs/halfline_open_02_referee.md`.

**What the Gemini pass established.** It independently established the following, so do not redo them:
- from CAF's TeX source, that their (cuR) has $S_{ii}(2\theta)$ while their (dn2) satisfies the relation only with $S_{ss'}$ for odd $n$ (30 digits);
- the sign chain $F(0)=-1$, with no $(-T/2)_I$ prefactor in $A$;
- the $d_4$/$d_5$ kink reductions and the boundary condition.

**The coordinator's check of the crossing form.** Starting from GZ's $K_a(\frac{i\pi}2-\theta)=S_{a\bar a}(2\theta)K_{\bar a}(\frac{i\pi}2+\theta)$, bulk crossing $S_{a\bar a}(i\pi-x)=S_{aa}(x)$ and unitarity give $K_a(\theta)K_{\bar a}(\theta+i\pi)=S_{aa}(2\theta)$. So the worker's and Gemini's form appears correct.

Read:
- the brief `briefs/halfline_open_02_beyond_an.md`;
- the report `briefs/halfline_open_02_report.md`;
- the Gemini verdict;
- the scripts in `briefs/halfline_open_02_scripts/`.

## Items

1. **Claim 7, before it enters the book's literature-corrections appendix.**
   - (i) Check that CAF's bulk $S_{ij}$ uses the same crossing convention as the book ($S_{ab}(i\pi-\theta)=S_{a\bar b}(\theta)$), and the same labelling of $s$ and $s'$ for odd $n$. A different convention in CAF's own S-matrix could reconcile (dn2) with (cuR).
   - (ii) Search for an erratum, a later version (check arXiv versions v2 and later), or a later paper by Castro-Alvaredo, Fring or others that corrects or comments on (dn2).
   - (iii) Draft the Appendix D entry.
2. **Claim 3, the bulk test that nobody ran.** Check that the new $d_n$ vector–vector scalar factor $F$, combined with the $d_n$ vector R-matrix, reproduces the bulk $S_{11}$ of $d_n$ for the lowest breathers. Compare $F$ with the published $d_n$ vector soliton S-matrix (Johnson 1997, Gandenberger–MacKay–Watts, or the book's Part VI) after converting conventions.
3. **The Fateev–Onofri reference.** The brief cites Fateev–Onofri hep-th/0207152 (IJMPA 18 (2003) 879); Gemini's verdict cites "Fateev & Onofri (2000, hep-th/0008149)". Identify which paper contains (r1c) and the $B_r$, $C_r$ families the worker used (see `briefs/halfline_open_02_scripts/fo_cmp.py` and the worker's `dl_fo` notes), and give the correct citation.

For each item give the verdict (CONFIRMED / CORRECTED / REFUTED / UNCHECKED), the evidence, and book text. Give script paths.
