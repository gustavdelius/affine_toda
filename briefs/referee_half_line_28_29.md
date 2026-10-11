# Brief: referee chapters 28–29 (scattering with a boundary; boundary sine-Gordon)

You are a mathematical-physics referee. The repository is `/home/gustav/Git/affine_toda` (branch `claude/affine-toda-boundary-dae69e`, PR #6), a Quarto book on affine Toda field theory. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files.

## Context

Part VII, "Affine Toda theory on the half-line" (`part7-half-line/`, chapters 28–35), was drafted on 2026-10-10 by agents working from the boundary literature. This brief covers:

- **Chapter 28** (`28-boundary-scattering.qmd`, `sec-boundary-scattering`): the scattering axioms with a boundary. It covers the reflection equation in S/K, braid and soliton-conjugating forms; boundary unitarity, the boundary state and crossing-unitarity; boundary bound states and the boundary bootstrap; boundary Coleman–Thun; boundary Bethe–Yang; a map of g-function, boundary TBA and TCSA; and boundary Lee–Yang as the worked example.
- **Chapter 29** (`29-boundary-sine-gordon.qmd`, `sec-boundary-sine-gordon`): the rank-one template. It covers the boundary sine-Gordon action, the soliton reflection matrix at the physical $q$, the parameter relations, breather reflection factors and their continuation to sinh-Gordon, the boundary bound-state spectrum, and the sinh-Gordon boundary breathers with their WKB quantization.

Find any object with `grep -rn '{#<id>}' part*/`. Book conventions are in `sec-conventions` (ch. 1), `sec-boundary-potential` (ch. 30) and `sec-qg-coproduct` (ch. 20). Specifically:
- Longest roots have length² 2, and $\beta_{\rm SG}^2=2\beta^2$.
- Part VI uses the coproduct $\Delta(e)=e\otimes1+k\otimes e$ and the physical $q=-e^{-i\pi\omega}=e^{-4\pi^2i/\beta^2}$.
- Boundary parameters are written $C_j$ (classical) and $\hat\epsilon_j$ (quantum coideal).
- Ghoshal–Zamolodchikov's boundary $\xi$ is replaced by Bajnok–Palla–Takács–Tóth's $(\eta,\vartheta)$, because $\xi=\pi/\lambda$ is taken.

The sources are mostly arXiv papers; fetch them with `curl -sL https://arxiv.org/pdf/<id>`. The main ones:
- Ghoshal–Zamolodchikov hep-th/9306002; Ghoshal hep-th/9310188; Fring–Köberle hep-th/9304141.
- Dorey–Pocklington–Tateo–Watts hep-th/9712197; Dorey–Tateo–Watts hep-th/9810098.
- Mattsson–Dorey hep-th/0008071; Bajnok–Palla–Takács–Tóth hep-th/0106070; Bajnok–Palla–Takács hep-th/0108157.
- Saleur–Skorik–Warner hep-th/9408004; MacIntyre hep-th/9410026.
- Corrigan–Delius hep-th/9909145; Chenaghlou–Corrigan hep-th/0002065; Baseilhac–Delius–George nlin/0201007.

Experience in this repository is that every adversarial referee pass finds real errors. The drafting itself already corrected several published formulas, so look hard at the corrections too.

## Already verified (do not redo unless you find a reason)

**By script** (`textbook_code/`, see its README; every script prints PASS lines):

- `boundary_axioms.py`:
  - the braid and S/K forms of the reflection equation agree on the sine-Gordon $K$, with a random-$K$ control;
  - the book's diagonal crossing-unitarity form (`eq-boundary-crossing-diagonal`) is equivalent to Corrigan–Delius (2.7) and Delius–Gandenberger (3.17);
  - `prp-ly-reflection`: the boundary Lee–Yang factors satisfy unitarity, crossing-unitarity with $f_{2/3}$ and the bootstrap. The sign of GZ (3.51) fails the bootstrap; residues and $g$ are as stated; the excited-boundary pole and gap are right;
  - the free-boson boundary Bethe–Yang check.
- `boundary_sine_gordon.py`:
  - the coideal $K$ equals the GZ/BPT matrix (`eq-bsg-k-matrix`, `eq-bsg-eta-theta`);
  - the reflection equation holds at $q$, and also at $-q$, so it cannot detect the sign of $q$;
  - unitarity and crossing-unitarity hold with the book's $S_0$. Crossing-unitarity with $C=\sigma_x$ holds only at the physical $q$;
  - the Neumann point is $\hat\epsilon_0=\hat\epsilon_1=(2-q-q^{-1})^{-1/2}$;
  - GZ (5.23)'s normalizer diverges, so BPT's form is used;
  - the breather factors from the soliton $K$ equal Ghoshal's (`prp-bsg-breather-k`), and their continuation equals Corrigan–Delius (2.5) with $E=B\eta/\pi$, $F=iB\vartheta/\pi$ (`eq-shg-reflection`);
  - the first excited boundary state (`prp-bsg-excited`, state $|0\rangle$ only);
  - the limits of the UV–IR relation (`cnj-bsg-uv-ir`).
- `shg_boundary_breathers.py`:
  - the Corrigan–Delius τ-functions solve the field equation and the boundary condition with $C_0=C_1=\varepsilon$;
  - the regular region holds;
  - the energy (3.7) and action (3.9) agree;
  - the WKB level spacings equal the bootstrap ones iff $E=2a(1-B/2)+(1-2N_0)B/2$;
  - $N_0=\frac12$ is needed for the Neumann limit;
  - the excited factors follow by bootstrap.
- `boundary_coideal.py` (used for `eq-reflection-twisted`): the twisted reflection equation holds for the $a_n$ solutions, and the crossing-unitarity index order is right.

**By hand** (worth a second look):
- Chapter 28:
  - the derivation of the reflection equation from factorization;
  - the boundary-state form of crossing-unitarity and its index order $k(\theta)=\sum K_{b\bar a}(i\pi/2-\theta)v_a\otimes v_b$;
  - the one-point residue normalization.
- Chapter 29:
  - the dictionary $C_{0,1}=(\beta_{\rm SG}^2M_0/4m_0)e^{\pm i\beta_{\rm SG}\phi_0/2}$ (also checked in `boundary_classical.py`);
  - the status of the Dirichlet relation, which GZ themselves call a conjecture.

## To be checked

Report on each item: correct / wrong (with the correction) / unclear (with what would settle it). Items are ordered roughly by risk.

1. **Corrections to the sources.** Each of the following is stated in the text; confirm or refute it independently:
   - the overall sign of GZ (3.51) ($R_{(1)}$ of boundary Lee–Yang);
   - the divergence of GZ (5.23);
   - the GZ (5.9) typo;
   - the claim that GZ present $\xi=4\pi\phi_0/\beta$ as a conjecture;
   - that Corrigan–Delius's $N_0=\frac12$ is fitted.
2. **Status labels.** Is anything overclaimed? In particular:
   - `cnj-bsg-uv-ir`: Al. Zamolodchikov's relation. What exactly is its status today, and has later work (e.g. TBA/NLIE, TCSA, the boundary sine-Gordon literature after 2002) established it?
   - `cnj-shg-e`, and its relation to Corrigan–Taormina's sinh-Gordon relation, which was not compared;
   - `prp-bsg-spectrum` [Established]: the BPTT/Mattsson–Dorey spectrum, quoted rather than derived;
   - `prp-boundary-elastic`.
3. **The BPTT factor $a_n$ for $n\ge1$.** The drafting agent could not reproduce BPTT's printed $a_n(\eta,u)=\prod_{l=1}^n\{2\eta/\pi-l\}$ from $a(u+\nu_n)a(u-\nu_n)/(a(u+\nu_0)a(u-\nu_0))$, and the chapter therefore does not print it. Settle it: is it a reading or block-convention issue, or an error in BPTT?
4. **Breather ordering.** In Part VI conventions the breather singlet sits in $V(\theta-iu_n/2)\otimes V(\theta+iu_n/2)$, the smaller imaginary part on the left, and this ordering reproduces Ghoshal exactly. Is this consistent with the ZF ordering of ch. 2 and the opposite coproduct of ch. 20, or does it signal a convention slip somewhere?
5. **Statements given only qualitatively** (ch. 28):
   - the residue normalization for boundary bound states;
   - "$(Q_s-Q_{-s})|B\rangle=0$ up to $i^s$": which combinations of charges survive on the half-line, and with which sign convention for $s$;
   - the boundary Coleman–Thun rule $P-2L$, quoted from Dorey–Tateo–Watts;
   - the paraphrased Mattsson–Dorey lemmas.
6. **The spectrum section** (`sec-bsg-spectrum`): the $\nu_n$ towers, which states exist for which $(\eta,\vartheta)$, and the excited-boundary reflection factors beyond $|0\rangle$. The breather factors on excited states are unchecked.
7. **The Neumann limit claim**: a tower of poles at $\theta=i(\pi/2-n\pi/\lambda)$ with energies $M\sin(n\pi/\lambda)$. Not checked by script.
8. **Claims about Saleur–Skorik–Warner and MacIntyre** (images still work for the general boundary; semiclassical agreement with GZ). The drafting agent took these from the abstracts only. Check them against the papers.
9. **Optional:** Delius–Gandenberger (3.16) at $n=1$ against Ghoshal/Corrigan–Delius at the uniform boundary.
10. **The finite-size map** (`sec-boundary-finite-size`): are the attributions for the g-function, boundary TBA and TCSA right? The cited works are Affleck–Ludwig, LeClair–Mussardo–Saleur–Skorik, Dorey–Fioravanti–Rim–Tateo and Dorey–Pocklington–Tateo–Watts. Is anything important missing?
11. **Exercises.** For every `exr-` in the two chapters: is it well posed and solvable with the text, and is any hint wrong?
12. **Bibliography.** These were added for these chapters: `gz1994`, `ghoshal1994`, `cherednik1984`, `sklyanin1988`, `fring_koberle1994a`, `sasaki1993`, `ssw1995`, `macintyre1995`, `mattsson_dorey2000`, `bptt2002`, `bpt2002`, `dptw1998`, `dtw1999`, `lmss1995`, `dfrt2004`, `affleck_ludwig1991`, `chenaghlou_corrigan2000`, `corrigan_delius1999`, `bdg2002`, `baseilhac_koizumi2003`. Check authors, titles, journal data and DOIs (`references.bib`), and whether each citation supports the sentence it is attached to.
13. **Prior art and later work.** Is anything presented as new here (the Lee–Yang sign, the divergence of GZ (5.23), the crossing-unitarity test of the sign of $q$) already in the literature?

## Deliverable

One entry per item, each with:
- a verdict;
- the justification;
- for errors, the corrected text ready to paste, with the file and line.

Include any scripts you wrote, inline or as a list of what they check with their output. Do not mark anything verified on the strength of a citation alone.
