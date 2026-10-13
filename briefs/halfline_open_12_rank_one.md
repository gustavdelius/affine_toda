# Brief: rank one on the half-line — the uniform sinh-Gordon boundary and the higher boundary sine-Gordon states (half-line open problem 12)

You are a mathematical physicist with expertise in boundary sine-Gordon and sinh-Gordon theory, boundary CFT and the boundary TBA. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it, including the lead in "A lead found while writing this brief".

## The problem

Open problem 12 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`) has three parts:

> (a) For sinh-Gordon, the uniform boundary sits at $E=\frac B2$ by `cnj-reflection-continuation` and at $E=0$ by Chenaghlou and Corrigan's conjecture extended to $C_0=C_1=1$. The two are reconciled if the continued Zamolodchikov relation (`cnj-bsg-uv-ir`) holds: then $\hat\epsilon=0$ is $C_0=C_1=\cos\frac{\pi B}{2(2-B)}$, a boundary parameter renormalized at $O(B^2)$ (`sec-uniform-boundary`). The two differ at $O(\beta^4)$, so a two-loop calculation would test this.
> (b) For sine-Gordon, the boundary CDD factors $a_n$ of the excited states (`eq-bsg-excited-k`) have been checked here only for $n\le2$ (`prp-bsg-excited`), and the higher states of `prp-bsg-spectrum` are extrapolated.
> (c) The relations between reflection parameters and Lagrangian parameters rest on `cnj-bsg-uv-ir`.

## What the book already has

Read `sec-bsg-parameters`, `sec-bsg-breathers`, `sec-bsg-spectrum`, `sec-shg-boundary-breathers` (ch. 29), `sec-uniform-boundary` (ch. 35, paragraph "Rank one") and `exr-shg-uniform`. Scripts: `textbook_code/boundary_sine_gordon.py` (GZ matrix, excited-state bootstrap for $n\le2$, breather factors), `shg_boundary_breathers.py` (Corrigan–Delius WKB, the $E(\epsilon,\beta)$ fit), `particle_reflection.py` (`K_shg`), `boundary_coideal.py` (the map $\hat\epsilon(\eta,\vartheta)$).

## A lead found while writing this brief (verify it first)

Ahn, Kim and Rim (`akr2002`, hep-th/0110218, eq. (10)) use the action $\int\frac1{8\pi}(\partial\phi)^2+\mu\sum_ie^{be_i\cdot\phi}+\mu_B\sum_i\eta_ie^{be_i\cdot\phi/2}$ and fix the boundary coupling of their $(+)$ boundary ($\eta_i=1$) by

$$\mu_B^2=\frac\mu2\cot\frac{\pi b^2}2 ,$$

citing Corrigan–Dorey–Rietdijk–Sasaki and the FZZ boundary Liouville paper (hep-th/0001012). They then check DG's uniform amplitude (their (12)) by the boundary TBA for $a_{n-1}^{(1)}$, $n\ge3$; sinh-Gordon is excluded "because of its two continuous parameters". A rough continuation, done while writing this brief, compares this with BPT's $\tilde\mu_{\rm crit}=\sqrt{2\mu/\sin(\beta_{\rm SG}^2/8)}$ for $n=1$. It uses $b^2\leftrightarrow-\beta_{\rm SG}^2/8\pi$ (the scaling dimension of $e^{ibe\cdot\phi}$ is $2|b|^2$, that of $\cos\beta_{\rm SG}\phi$ is $\beta_{\rm SG}^2/4\pi$), $\mu_{\rm SG}=2\mu$ and $\tilde\mu=2\mu_B$ (from $\cosh$ = two exponentials). The result is

$$\Big(\frac{\tilde\mu}{\tilde\mu_{\rm crit}}\Big)^2=-\cos^2\frac{\pi b^2}2\ \longleftrightarrow\ \cos^2\frac\pi{2(\lambda+1)} \quad(\text{up to the sign from }\mu\to-\mu),$$

which is exactly the $\hat\epsilon=0$ point $\tilde\mu/\tilde\mu_{\rm crit}=\cos\frac\pi{2(\lambda+1)}$ of `cnj-bsg-uv-ir` at $\phi_0=0$. If this survives a careful check, the "uniform" boundary that AKR test non-perturbatively for $n\ge3$ is already the renormalized point. At $n=1$ that point would carry $E=\frac B2$, and the classical $C_0=C_1=1$ ($\tilde\mu=\tilde\mu_{\rm crit}$) would carry $E=0$. This is the reconciliation stated in the book, now with independent support.

## Literature

- Bajnok–Palla–Takács (`bpt2002`) and Bajnok–Palla–Takács–Tóth (`bptt2002`): the UV–IR relation, the boundary energy, TCSA, and the full excited-state spectrum. Also BPT, "Finite size effects in boundary sine-Gordon theory", hep-th/0108157, and the review "Boundary sine-Gordon model", hep-th/0211132 (not in `references.bib`).
- Fateev–Zamolodchikov–Zamolodchikov, boundary Liouville I, hep-th/0001012 (not in `references.bib`): the $\cosh(\pi bs)$ parametrization of the boundary cosmological constant, from which both AKR's (10) and the UV–IR relation descend.
- Chenaghlou–Corrigan, hep-th/0002065 (`chenaghlou_corrigan2000`): one loop, first order in $C_0-C_1$. Chenaghlou, hep-th/0108077 (IJMPA 16 (2001) 4613; not in `references.bib`): second order in $C_0-C_1$, still **one loop** in the bulk coupling. Also hep-th/0004121, "On the quantum reflection factor for the sinh-Gordon model with general boundary conditions" (not in `references.bib`; check its authors and content). No genuine two-loop boundary calculation was found in the search.
- Corrigan–Taormina (`corrigan_taormina2000`): WKB for $C_0\ne C_1$, $N_0=\frac12$ from imaginary stability angles.
- Bajnok–Böhm–Takács, hep-th/0309119 (NPB 682 (2004) 585) and the boundary reduction formula, hep-th/0207079 (not in `references.bib`): perturbation theory about Neumann, Landau equations and Cutkosky rules with a boundary. This is the framework a two-loop calculation would use.
- Kormos–Takács, boundary form factors of exponential operators in sinh-Gordon, arXiv:0801.0962 (not in `references.bib`): the exact relation between bootstrap and Lagrangian parameters "in exact form in the perturbed boundary CFT framework", as their abstract states.

## Suggested plan of attack

1. **Check the lead.** Redo the continuation of AKR's (10) to $n=1$ with all factors of 2, $i$ and signs tracked, and compare with FZZ's parametrization and BPT's (3.2)–(3.3). Decide whether the AKR $(+)$ point is $\hat\epsilon=0$. Also check whether AKR's TBA would distinguish $\mu_B$ at the AKR value from $\mu_B$ at the classical critical value: the boundary energy depends on $\mu_B$.
2. **Part (a) without two loops.** The boundary ground-state energy is known exactly as a function of the boundary parameters (BPT; for $n\ge3$ AKR (55)). Expand the exact energy at the uniform point, with either identification, to $O(b^4)$. Compare with a perturbative calculation of the boundary energy, which is a vacuum diagram and much easier than a two-loop reflection factor. Alternatively, run a boundary TBA (sinh-Gordon with $F=0$, $E=\frac B2$ versus $E=0$) and compare the ground-state energy with a Hamiltonian-truncation or TCSA computation of the sinh-Gordon model on an interval at $C_0=C_1=1$.
3. **Part (a) with two loops (optional).** If step 2 is inconclusive, set up the two-loop reflection factor around the uniform boundary in the Bajnok–Böhm–Takács framework. Isolate only the $O(\beta^4)$ coefficient that distinguishes $E=\frac B2$ from $E=0$: by `exr-shg-uniform` the two factors differ by a known function at $O(B^2)$.
4. **Part (b).** Extend `prp-bsg-excited` to $n=3,4,5$ and to the two-label states $|n_1,n_2\rangle$ of `prp-bsg-spectrum`. Check the boundary bootstrap with the full matrix reflection factors, the ratio and product forms of $a_n$, and the pole bookkeeping (each pole a state or a Coleman–Thun diagram). Do this at several generic $(\lambda,\eta,\vartheta)$. Compare low levels with BPT's published TCSA data where they quote numbers.
5. **Part (c).** Collect the evidence for `cnj-bsg-uv-ir` (FZZ, BPT sec. 3.5, Kormos–Takács). State precisely which part is derived (for example, on the line $\phi_0=0$ from the boundary energy and the FZZ one-point function) and which is conjectured. Propose a status label.

## Known traps

- The node crossing $\hat\epsilon_1\leftrightarrow C_0$ and $\hat\epsilon_0\leftrightarrow C_1$ in sine-Gordon (`sec-bsg-parameters`, "Node labels"), and the sign of $F$ that goes with it.
- GZ's boundary $\xi$ clashes with the book's $\xi=\pi/\lambda$; use BPT's $(\eta,\vartheta)$.
- (`eq-shg-reflection`) depends on $E^2$ only, so $E=\pm\frac B2$ cannot be distinguished; one loop cannot distinguish $E=0$ from $E=\frac B2$ either.
- `cnj-shg-e` needs $\frac12<a<1$ ($-1<\epsilon<0$); $C=1$ is outside the WKB regime, so it is reached only by continuation.

## Deliverable

- A verdict on the lead, with the full convention chain.
- For (a): which identification the evidence supports, with the computation (energy expansion, TBA or two-loop) and its precision.
- For (b): a table of checked excited states and the remaining unexplained poles, if any.
- For (c): a proposed status for `cnj-bsg-uv-ir` and replacement text for the "Rank one" paragraph of `sec-uniform-boundary` and for `prp-bsg-excited` and `prp-bsg-spectrum`.
- Scripts.

## Feasibility (assessment when the brief was written, 2026-10-11)

**High for (b)** and for the convention check of the lead. These are numerics and bookkeeping with existing scripts. **Medium for (a)** through the boundary energy or a TBA, **low** through a genuine two-loop reflection factor. **Low for (c)** as a proof; a sharper status statement is feasible.
