# Referee prompt: half-line open problem 12a
Tier: B
Brief: briefs/halfline_open_12_rank_one.md
Report: briefs/halfline_open_12a_report.md
Scripts: briefs/halfline_open_12a_scripts/ (`rank_one_lead_check.py`, `shg_btba_uniform.py`, which takes about 3 minutes, and their outputs `lead_out.txt`, `btba_out.txt`; they import only numpy and mpmath)

Scope of the report: the lead ("A lead found while writing this brief") and part (a) of problem 12. Part (b) is a separate report (tag 12b). Tier B because the report corrects published papers (claims 9, 10) and the book (claim 11).

## Claims to referee

1. [Derivation] **The convention chain.** With $b^2=\beta^2/4\pi=B/(2-B)$, AKR's $a_1$ action maps to FZZ/BRZ form by $\phi=\sqrt2\phi'$, with $M^\pm=\mu_B$. With $M=M_0\cosh\pi bs$ and $M_0=\sqrt{\mu/\sin\pi b^2}$, AKR (10) is $s=ib/2$.
2. [Derivation] **The sign correction to the brief.** Under the continuation $b^2\to-1/(\lambda+1)$ the continued $\mu$ is negative, so $\mu_{SG}=-2\mu$ and $\tilde\mu=-2\mu_B$. The signs cancel, giving $(\tilde\mu/\tilde\mu_{\rm crit})^2=+\cos^2\frac\pi{2(\lambda+1)}$, and BPT (3.2) at $\phi_0=0$ gives $\eta=\pi/2$, $\vartheta=0$, i.e. $\hat\epsilon=0$.
3. [Derivation] **The sinh-Gordon relation.** $\cos\frac{\pi E}{2-B}=M/M_0$ (equivalently Chenaghlou–Corrigan's $E=2a(1-B/2)$ with $C=\cos\pi a=M/M_0$), so $s=ib/2$ gives $E=B/2$ and $s=0$ gives $E=0$. Also: classically $C_0=C_1=1$ corresponds to $M/M_0\to1$.
4. [Literature] **The origin of AKR (10).** It is Fateev's (hep-th/0103014, eq. "MND") quantum-modified boundary coupling, from null-vector equations and W-symmetry, and for $n\ge2$ it is the only integrable value.
5. [Numerical] **AKR's TBA.** AKR's $\delta_1$ table is reproduced at $s=ib/2$ (12 entries, to 6 digits). At $s=0$, $\delta_1$ shifts by 0.46–8.9.
6. [Numerical] **The boundary energies.** Fateev's (+) boundary energy at $h=2$ equals the BPT/BRZ sinh-Gordon boundary energy at $E=p=B/2$ (15 digits at five values of $p$). Its dual matches $E=1-p$. The energies at $E=B/2$ and $E=0$ differ by $-m\pi p/16+O(p^2)$.
7. [Numerical] **The worker's sinh-Gordon boundary TBA.**
   - It is validated against BRZ Table 1.
   - With $E=B/2$ input it matches the boundary Liouville phase at $s=ib/2$ to about $10^{-11}$ and misses $s=0$ by 0.025–0.20 (at $b^2=0.8086, 0.4286$).
   - The $E=0$ runs mirror this.
8. [Status] **Proposed status for `cnj-bsg-uv-ir`.** "[Numerical] in sinh-Gordon form ($0<b^2<1$); [Conjecture] for sine-Gordon". In addition, the BPT analytic route needs a factor ½ in $\log g_S$, and arXiv:0801.0962 is by Takács alone, with only an $O(g^2)$ check.
9. [Literature correction] **AKR's printed $\Delta c_{\rm boun}$ is wrong.** It has $\sin(\pi(1-B)/2h)$ where Fateev has $\cos(\pi(1-x)/2h)$, and it lacks a factor $1/\pi$; at $h=2$ it is off by $\pi\cot(\pi(1+B)/4)$.
10. [Literature correction] **Two printed slips.**
    - BRZ eq. (Bfugacity) prints $\cos(\eta_a p)$ where (R1) has $\cosh(p\eta)$; only cosh reproduces their Table 1.
    - Takács's (zEFmatch) and (zetathetamatch) each lack a factor $i$, and the two cancel.
11. [Book correction] **Self-duality.** Ch. 29's "Self-duality" paragraph calls $B=2$ the self-dual point; by `eq-b-of-beta` it is $B=1$.
12. [Text] **The proposed replacement text** for the "Rank one" paragraph, for open problem 12(a), and for the `cnj-bsg-uv-ir` status paragraph. Check correctness, labels and citations, and the proposed bib entries (`brz2008`: NPB 796 (2008), page "622?" unverified; `fzz2000`; `fateev2001`).

## Coordinator notes

- **The crux is claim 2.** The sign of the continued $\mu$ decides whether the AKR point is $\hat\epsilon=0$. Re-derive it from AKR's mass–$\mu$ relation yourself, and check that it really puts the sine-Gordon vacuum at $\phi=0$. Also check the identifications $\mu_{SG}$ and $\tilde\mu$, including the factors of 2 from $\cosh$ being two exponentials.
- **Claim 3, "classically $C=1$ ↔ $M/M_0\to1$".** The classical boundary coupling is $\mu_B^2=\mu/\pi b^2$, and $\sin\pi b^2\approx\pi b^2$. Check that $M_0$ reduces to this with AKR's normalization ($\mu_B^2$ versus $M^2$, and the $\frac\mu2$ in (10)).
- **Claim 7.** The TBA code and its Liouville-phase extraction were written by the worker. Write your own check of at least one table entry, for example the BRZ Table 1 validation point, rather than relying only on a rerun.
- **Claim 4, "for $n\ge2$ the only integrable value".** It is unclear whether this is Fateev's statement or the worker's inference. Find out which.
- **Claim 12, unlisted provenance.** The proposed text cites `textbook_code/shg_boundary_tba.py`, which does not exist yet; it would be the worker's `shg_btba_uniform.py`. The worker mentions a factor ½ "confirmed privately by Zamolodchikov". Check where BPT say this.
- **Prior art.** BRZ (arXiv:0710.4789) was not in the brief, so check whether BRZ or a later paper already states the conclusion for the uniform point $s=ib/2$ explicitly. Search, for example, for Bajnok–Rim–Zamolodchikov follow-ups and for Ahn and collaborators on boundary affine Toda.
