# Brief: exchange statistics of affine Toda solitons in the Bethe–Yang equations

You are a mathematical-physics research assistant. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files.

## The question

The paper Part I (`foundations/*.qmd`) classifies sectors of imaginary-coupling affine Toda field theory as $\Theta$-unbroken or $\Theta$-broken. A sector is unbroken when every multi-particle Bethe–Yang transfer matrix has eigenvalues of modulus 1 for real rapidities (Section `sec-krein-unitarity`).

For several self-conjugate theories the verdict depends on how soliton labels are exchanged in the Bethe–Yang equations:

- with the ordinary flip $P$ (bosonic exchange), or
- with a $\mathbb Z_2$-graded flip (some soliton states fermionic).

A change of gradation is a similarity and changes no eigenvalue. A $\mathbb Z_2$ grading of the flip does change eigenvalue moduli. The current findings (Section `sec-sector-status` of the paper, paragraphs "Fundamental kinks of $a_2^{(2)}$", "Spinor solitons" and "Caveat: exchange statistics"):

| theory, sector | bosonic exchange | graded exchange |
|---|---|---|
| $a_2^{(2)}$, $N$ fundamental kinks | unbroken iff $\lvert\psi\rvert\le\pi/N$ | unbroken iff $\lvert\psi\rvert\le2\pi/N$ (charged kinks fermionic) |
| $c_2^{(1)}$ spinor | two-soliton amplitude broken | unbroken up to three solitons |
| $c_3^{(1)}$ spinor | broken | broken for every one of the 512 $\mathbb Z_2$ bicharacter twists |
| $a_n^{(1)}$ fundamental | three-colour sectors broken | the $\mathbb Z_n$ symmetry allows only a uniform grading, so the result does not depend on this |

Here $q^2=-e^{i\psi}$ and $\psi=2\pi^2/3\xi$ reduced to $(-\pi,\pi]$, in the conventions of Takács–Watts (hep-th/9810006), with $q=i\,e^{i\pi^2/3\xi}$ and $x=e^{2\pi\theta/\xi}$.

**Goal.** Decide which convention the field theory selects, first for $a_2^{(2)}$ and then for the $c_2^{(1)}$ and $c_3^{(1)}$ spinor solitons (and, if time permits, the $a_{2n-1}^{(2)}$ spinor). Then state the resulting classification.

## Read first

- Part I (`foundations/*.qmd`): Section `sec-krein-unitarity` (Propositions `prp-krein-unitarity` and `prp-krein-fundamental`, the ladder of conditions), Section `sec-sector-status` (the table, Proposition `prp-scalar-breathers` and the paragraphs named above), and open problem 6 (the list item "The sector classification for self-conjugate theories" in `foundations/12-open-problems.qmd`).
- Part II (`smatrices/*.qmd`): Section `sec-crossing-sign` in `smatrices/02-construction.qmd` ("The crossing sign": the numerical crossing test of the $c_n^{(1)}$ spinor returns $(-1)^n$) and open problem 4 there (the list item "Crossing convention" in `smatrices/09-discussion.qmd`). The grading question may be related.
- Code: `foundations_code/sectors/` (`a22*.py` for the Izergin–Korepin amplitude and Takács–Watts' cubic, `graded*.py`, `ikgraded.py`, `spin*.py`, `lib.py`), `foundations_code/krein_bethe.py`, and `soliton_rmatrix_code/smatrix/` for the Part II amplitudes (see its `README.md`). Copy code to your scratch directory before running it.

## Suggested approach

1. **Cheap analytic test.** In an interacting integrable theory the two-particle S-matrix at zero relative rapidity is minus the exchange operator, $S(0)=-\mathcal P$, in the statistics the Bethe–Yang equations must use; for sine-Gordon solitons $S(0)=-P$. Compute $S(0)$, including the scalar factor, for the $a_2^{(2)}$ kink amplitude (Smirnov; Takács–Watts eqs. 5.1–5.5) and for the $c_2^{(1)}$ and $c_3^{(1)}$ spinor amplitudes of Part II. Is it $-P$, $-P_{\rm graded}$ or neither? Check whether the answer is consistent with the crossing sign $(-1)^n$.
2. **Literature.** Search for derivations of Bethe–Yang equations with an explicit statistics phase from a lattice or from non-linear integral equations: for the Izergin–Korepin / Bullough–Dodd / $a_2^{(2)}$ case (for example work by Dunning, by Ahn, Bajnok, Palla and Ravanini, by Takács and Watts, or by Vernier, Jacobsen and Saleur on $a_2^{(2)}$ chains), and for twisted or spinor cases. Record exactly what they establish.
3. **Decisive numerical test, if steps 1–2 do not settle it.** Compare the low-lying finite-volume levels of a regularization whose scaling limit is imaginary-coupling $a_2^{(2)}$ Toda with Bethe–Yang predictions in both conventions, at a coupling where they differ, that is, $\cos\psi<0$ for two kinks. Options are the staggered or light-cone Izergin–Korepin vertex model in the right regime, a non-linear integral equation from the literature, or published TCSA data for an RSOS restriction where the kink content survives. Two-kink states at moderate volume suffice.
4. Report which convention holds, for which theories, and with what status, and how the Section `sec-sector-status` table changes.

## Known traps

- The Bethe–Yang operator is the R-form $R=P\check R$, with $R(1)\propto P$, not the braided $\check R$, whose eigenvalues are trivially phases.
- Gradations are similarities. Only a grading of the flip, or a change of the physical normalization at $\theta=0$, can change eigenvalue moduli.
- Twisted algebras: the literature normalizes the longest root to squared length $2k$, so its coupling is $\beta^2/k$ in the paper's normalization (longest roots of squared length 2).
- Prior art: Takács and Watts (1999) found much of the $a_2^{(1)}$ and $a_2^{(2)}$ picture. Their Bethe ansatz (their eq. 2.9) has only a scalar statistics phase. Check what they and Smirnov assume before claiming novelty.
- In either convention, eigenvalues must come in pairs $(s,1/\bar s)$ (Proposition `prp-krein-unitarity`). Use this as a consistency check.
- Usage is limited: run computations incrementally, save outputs, and keep the scope realistic. A clean answer for $a_2^{(2)}$ alone is valuable.

## Deliverable

In your final message:

- a summary;
- the method, with conventions stated explicitly;
- results with status labels in the paper's convention: [Theorem] only with a complete proof, [Proposition, sketch], [Established] for literature results, [Numerical], [Conjecture];
- the numerical precision and the scripts used;
- proposed replacement text for the Section `sec-sector-status` caveat, the affected table rows and the open problem 6 item, in the paper's style (short declarative sentences, LaTeX in `$...$`, no hype);
- every doubt, stated plainly. If the question cannot be settled, say what would settle it.
