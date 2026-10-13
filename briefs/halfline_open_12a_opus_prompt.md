# Opus referee prompt: half-line open problem 12a (tier-B follow-up)

Written 2026-10-11 by the Claude coordinator, after reading the Gemini verdict. The coordinator replaced Gemini's "Items for an Opus pass" with the items below. Gemini's items 2 and 3 were tangential (item 3 is a new research question, recorded as a follow-up). The real weak points were not on its list.

You are an independent adversarial referee. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it, and do not spawn subagents. Work in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/ref12a/`. Read:

- the brief `briefs/halfline_open_12_rank_one.md`;
- the worker report `briefs/halfline_open_12a_report.md`;
- the Gemini verdict `briefs/halfline_open_12a_referee.md`;
- the scripts in `briefs/halfline_open_12a_scripts/`.

The worker is Claude Opus, the same model as you, so distrust steps that merely look plausible. Fetch the arXiv TeX sources yourself and read the equations directly.

## Items

1. **Claim 9, the AKR $\Delta c_{\rm boun}$ misprint.** The worker and Gemini describe it differently.
   - Worker: AKR have $\sin(\pi(1-B)/2h)$ where Fateev has $\cos(\pi(1-x)/2h)$, and a factor $1/\pi$ is missing. At $h=2$ the error is a factor $\pi\cot(\pi(1+B)/4)$.
   - Gemini: in addition, the numerator is $\cos(\pi/2h)$ versus Fateev's $\sin(\pi/2h)$, and there is a "scale $B$ vs $x$" issue with $x=B/2$.
   - At $h=2$ a numerator $\cos\leftrightarrow\sin$ slip is invisible, so the worker's $h=2$ check could not have caught it.

   Settle exactly what AKR print, for general $h$ (AKR hep-th/0110218; Fateev hep-th/0103014, eq. "eQ"), and how AKR's $B$ relates to Fateev's $x$ (or $p$). Write the correct formula. Then decide whether the misprint is purely typographical or propagates into AKR's numerical tables and figures: compare their quoted numbers for $h=3,4,5$ with both formulas. Produce exact wording suitable for the book's literature-corrections appendix.

2. **Claim 4, "for $n\ge2$ the only integrable value".** The proposed book text says that for $n\ge2$ there is no free boundary parameter and that the uniform coupling is fixed at $\mu_B^2=\frac\mu2\cot\frac{\pi b^2}2$. Gemini's support is that Fateev derives this value from $W$-symmetry and null vectors in *conformal* Toda theory. Determine precisely:
   - what Fateev and AKR state;
   - whether the step from "preserves $W$-symmetry" to "is the integrable affine Toda boundary" is argued anywhere, or only assumed;
   - how it squares with the classical fact that integrable $a_n^{(1)}$ boundaries with $n\ge2$ have $|C_j|=1$, with no free parameter.

   Propose a corrected sentence and status label for the book text.

3. **Claim 10, the minor misprints.**
   - BRZ (arXiv:0710.4789) eq. "Bfugacity": $\cos(\eta_a p)$ versus $\cosh(p\eta)$. Confirm that only the cosh reading reproduces their Table 1, by your own numerics at one table entry.
   - Takács (arXiv:0801.0962) eqs. "zEFmatch" and "zetathetamatch": Gemini's explanation is vague. State exactly where a factor $i$ is missing and show that the two slips cancel in $E=2\xi\eta/\pi$.

## Not in scope

Claims 1–3 and 5–8 are confirmed by the worker and, independently, by Gemini; in particular Gemini re-derived the sign of the continued $\mu$ from AKR's mass–coupling relation. Do not redo them unless item 1 or 2 shows they depend on something wrong.

## Report

For each item give the verdict (CONFIRMED / CORRECTED / REFUTED / UNCHECKED), the evidence with equation numbers and line references in the TeX sources, and proposed text. Give the full paths of your scripts. Deliver everything in your final message.
