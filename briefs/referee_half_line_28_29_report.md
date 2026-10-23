# Referee report: chapters 28–29 (scattering with a boundary; boundary sine-Gordon)

This is the report of the 2026-10-11 referee pass on `referee_half_line_28_29.md`. The referee read both chapters and the cited source texts (GZ, Ghoshal, DPTW, DTW, MD, BPTT, BPT, SSW, MacIntyre, CD), and also fetched Corrigan–Taormina (CT), hep-th/0008237. Line numbers refer to the files at commit 64ed4af: ch28 is `28-boundary-scattering.qmd` and ch29 is `29-boundary-sine-gordon.qmd`. The scripts are listed at the end. The integration status is recorded at the end.

**Most serious findings**
1. The boundary CDD factor is $f_x(\theta)$, not $f_x(2\theta)$. This affects ch28 l.128 and `exr-boundary-cdd`.
2. Ch29 l.98 puts the poles of $R_0$ in the wrong place. They are at $u=n\pi/2\lambda$, not at $\pi/2-n\pi/2\lambda$.
3. Zamolodchikov's UV–IR relation was first printed by BPTT (sec. 4), not by BPT.
4. The $N_0=\frac12$ status is outdated: Corrigan–Taormina obtain it automatically.
5. The comparison with Corrigan–Taormina is now done. Their relation is exactly the analytic continuation of Zamolodchikov's relation.
6. Several attributions overclaim:
   - the $P-2L$ rule, attributed to DTW;
   - the residue sign;
   - the Mattsson–Dorey lemmas;
   - the sinh-Gordon Coleman–Thun example;
   - MacIntyre's "boundary breathers".
7. Item 3: BPTT is correct as printed. It has $\{2(\eta/\pi-l)\}$, and pdftotext dropped the parentheses.

## Item 1. Corrections to the sources
**1a. Sign of GZ (3.51): correct.**
- **Bootstrap.** In GZ's own conventions the $+$ sign gives unitarity, crossing-unitarity and the bootstrap to about 1e-25. GZ's $-$ sign fails the bootstrap with residual 3.3.
- **Residue.** At $i\pi/2$ it is $ig^2/2$, as GZ (3.55) requires.
- **Not new.** DPTW (2.8) print $+$ without remarking on it.

**1b. Divergence of GZ (5.23): correct (also analytically).**
- **Why it diverges.** $\log\Pi(x,\pi/2)\sim\frac\lambda4\log l_{\max}$, which gives an infinite constant independent of $u$.
- **New: GZ's normalizations are inconsistent with each other.**
  - With $\sigma$ normalized as (5.24) prints it, (5.22) violates (5.20) by the factor $k^2/\cos^2\xi$. The working prefactor is $\pm1/k$.
  - Equivalently, (5.22) and (5.26) are right if (5.24) reads $\sigma(x,u)\sigma(x,-u)=\cos^2x/(\cos(x+\lambda u)\cos(x-\lambda u))$, which is BPTT's normalization.
- **Replacement text for ch29 l.297, after its first sentence:** "The convergent form is that of [@bptt2002] used above, normalized by $\sigma(x,u)\sigma(x,-u)=\cos^2x/(\cos(x+\lambda u)\cos(x-\lambda u))$. This is also the normalization that GZ's own prefactors in (5.22) and (5.26) require: with the normalization printed in [@gz1994, eq. (5.24)], (5.22) violates (5.20) by the factor $k^2/\cos^2\xi$."

**1c. GZ (5.9) typo (stated in ch. 32 l.305): incomplete.**
- All four operators on the right-hand sides of both lines of (5.9) should carry $-\theta$.
- **Replacement text:** "In their (5.9) all four operators on the right-hand sides should carry rapidity $-\theta$, as in their (3.22)."

**1d. GZ present $\xi=4\pi\phi_0/\beta$ as a conjecture: correct.** GZ add that it is conjectured only for $M=\infty$.

**1e. "$N_0=\frac12$ is fitted": correct for CD, but outdated.** Corrigan–Taormina (hep-th/0008237, footnote near their (5.28)) get it automatically by including the imaginary solutions of the stability-angle equation. They also extend the WKB match to $C_0\ne C_1$. See item 2b.

## Item 2. Status labels
**2a. `cnj-bsg-uv-ir`: the label is right, but the attribution is wrong and the status incomplete.**
- **Replacement text for ch29 l.112:** "Bajnok, Palla, Takács and Tóth proposed the relation for general $(M_0,\phi_0)$ from symmetry arguments and the Dirichlet limit [@bptt2002, sec. 4], and Al. B. Zamolodchikov presented the same relation, with the normalization $\tilde\mu_{\rm crit}$, in unpublished seminars; it was printed in that form by Bajnok, Palla and Takács [@bpt2002, eqs. (3.2), (3.3)]."
- **Add after "…to numerical accuracy [@bpt2002]." (l.123):** "On the line $\phi_0=0$ they also showed that the relation follows from the boundary energy and the Fateev–Zamolodchikov–Zamolodchikov expectation value of $e^{i\beta_{\rm SG}\phi(0)/2}$ [@bpt2002, sec. 3.5]. The TCSA comparison used the relation as input, so it tests it only indirectly."
- **Later literature: unclear.** Settle with an INSPIRE search of papers citing hep-th/0108157.

**2b. `cnj-shg-e`: the label is right. The CT comparison is now done.**
- **CT's relation.** Their (2.10) reads $E=(a_0+a_1)(1-\frac B2)$, $F=(a_0-a_1)(1-\frac B2)$, with $\epsilon_j=\cos\pi a_j$. CT derive $E$ by WKB and obtain $F$ from the symmetry $a_j\to-a_j$.
- **Continuation.** Continuing `cnj-bsg-uv-ir` with $\lambda\to-2/B$, $\eta=\pi E/B$, $\vartheta=-i\pi F/B$ gives exactly CT's relation, with $C_{0,1}=(\tilde\mu/\tilde\mu_{\rm crit})e^{\pm i\beta_{\rm SG}\phi_0/2}$ (residual 1e-16).
- **Replacement for the last sentence of l.252:** "Corrigan and Delius fitted $N_0=\frac12$; Corrigan and Taormina later found that it arises automatically when the imaginary solutions of the stability-angle equation are kept [@corrigan_taormina2000]."
- **Replacement for the l.255 bullet:** "Corrigan and Taormina extended the boundary breathers and the WKB matching to $C_0\ne C_1$ and found $E=(a_0+a_1)(1-\frac B2)$, supporting the generalization above; $F$ follows from the symmetry $a_j\to-a_j$ [@corrigan_taormina2000]. The relation is exactly the continuation $\lambda\to-2/B$ of @cnj-bsg-uv-ir, with $C_{0,1}=(\tilde\mu/\tilde\mu_{\rm crit})e^{\pm i\beta_{\rm SG}\phi_0/2}$, as Bajnok, Palla, Takács and Tóth observed [@bptt2002, sec. 4]."
- **Ch29 l.298:** replace "as that paper says" with "as that paper says; Corrigan and Taormina later derived it [@corrigan_taormina2000]".
- **New bib entry:** `corrigan_taormina2000` (Corrigan, E. and Taormina, A., J. Phys. A 33 (2000) 8739, hep-th/0008237; to be verified).

**2c. `prp-bsg-spectrum` [Established]: slightly overclaimed.**
- **Replacement for the first sentence of l.201:** "The low-lying part of this spectrum was confirmed against truncated conformal space data, mainly for $\phi_0=0$ [@bpt2002]. Bajnok et al. account for every pole by a bound state or a boundary Coleman–Thun diagram, following the Dirichlet analysis of [@mattsson_dorey2000]; the reflection factors on the higher states are extrapolated from the low cases [@bptt2002]."

**2d. `prp-boundary-elastic`: acceptable.**

## Item 3. BPTT $a_n$: a reading problem; BPTT is correct
- **Formula.** The PDF prints $a_n=\prod_{l=1}^n\{2(\frac\eta\pi-l)\}$, with parentheses that pdftotext dropped.
- **Checks.** The ratio $a(u+\nu_n)a(u-\nu_n)/(a(u+\nu_0)a(u-\nu_0))$ equals $\prod\{2\eta/\pi-2l\}$ for $n=1,2,3$. The full bootstrap holds for $n=1,2$. BPTT's $b^n_k$ product (3.9), with upper limit $\min(n,k)$, equals $a_k(u+\frac\pi2-u_n)a_k(u-\frac\pi2+u_n)$.
- **Insert after `eq-bsg-excited-k` (l.175), replacing "with a boundary CDD factor $a_n$, $a_0=1$":** "with the boundary CDD factor $a_n(\eta,u)=\frac{a(u+\nu_n)a(u-\nu_n)}{a(u+\nu_0)a(u-\nu_0)}=\prod_{l=1}^n\big\{2\big(\tfrac\eta\pi-l\big)\big\}$, where $a=S_0$ and $\{y\}=\frac{(\frac{y+1}{2\lambda})(\frac{y-1}{2\lambda})}{(\frac{y+1}{2\lambda}-1)(\frac{y-1}{2\lambda}+1)}$ with $(x)=\sin(\frac u2+\frac{\pi x}2)/\sin(\frac u2-\frac{\pi x}2)$".
- **Also:** update `prp-bsg-excited` to say "$a_1,a_2$ checked".

## Item 4. Breather ordering: consistent
- **Check.** $\check R(\theta_{12})$ with the Part VI coproduct has a rank-1 residue at $\theta_{12}=+iu_n$. Its image is the singlet in $V(\theta-\frac{iu_n}2)\otimes V(\theta+\frac{iu_n}2)$ (overlap 1.000000000000).
- **Suggested text for l.141, after "...no such vector.":** "It is the image of the residue of $\check R(\theta_1-\theta_2)$ at the physical pole $\theta_1-\theta_2=iu_n$."

## Item 5. Statements given only qualitatively (ch. 28)
**5a. Residue normalization (l.140): wrong.** The correct factor is $\frac i2gg$ (GZ (3.57), (3.62)), and the positivity claim is in neither GZ nor DTW. Replacement text: "Near the pole $K_{a\alpha}^{b\beta}(\theta)\simeq\frac i2\frac{g_{a\alpha}^{\gamma}g_{\gamma}^{b\beta}}{\theta-iv_{a\alpha}^\gamma}$ [@gz1994, eqs. (3.57), (3.62)]. In a unitary theory with real couplings the diagonal residues are $\frac i2$ times a positive number. This is a heuristic test; in non-unitary models such as boundary Lee–Yang the couplings are imaginary and their signs are part of the data [@dtw1999]."

**5b. Line 98: imprecise.** Replacement text: "...become the conditions $(Q_s-Q_{-s})|B\rangle=0$ [@gz1994, eq. (2.21)]: the rotation multiplies $Q_{\pm s}$ by $i^{\mp s}$, a relative sign $(-1)^s=-1$ for the odd spins."

**5c. $P-2L$ attributed to DTW (l.167): overclaimed.** Replace "Dorey, Tateo and Watts showed that the power counting ... [@dtw1999]" with: "Dorey, Tateo and Watts applied the bulk counting rule of @sec-coleman-thun to such diagrams, treating reflections as vertices, and found that it explains every physical-strip pole of the boundary Lee–Yang reflection factors once the order is reduced by one for each internal reflection at a zero of a reflection factor [@dtw1999]. As in the bulk, the rule is a heuristic." Also add "[@mattsson_dorey2000; @bptt2002]" to the first bullet: MD use cancellations between diagrams.

**5d. Mattsson–Dorey lemmas (l.172): the paraphrase is inaccurate.**
- **Replacement text for l.172:** "Mattsson and Dorey proved two lemmas, using only momentum conservation, that force certain poles to be bound states [@mattsson_dorey2000]. Suppose particle $a$ has a simple pole at $\theta=iu$, and let $U_a$ be the smallest angle $\pi-U_{ab}^c$ over its three-point couplings. Lemma 1, for any boundary state: if $u<U_a$, no Coleman–Thun diagram exists and $a$ binds to the boundary. Lemma 2, for the ground state, where rescattering is confined to a smaller region: $a$ must bind unless $U_a<u<\frac\pi2-U_a$, and it always binds if all its decay angles exceed $\frac\pi2$. A pole can also persist beyond the coupling at which its state leaves the spectrum; @sec-shg-boundary-breathers shows a sinh-Gordon example whose explanation is open."
- **Replacement text for ch29 l.257:** "Corrigan and Delius note that no Coleman–Thun diagram is available in sinh-Gordon theory with this boundary, since there is no three-point coupling, and suggest, following Dorey, that the explanation lies in a mechanism involving solitons, presumably via the continuation from sine-Gordon theory. Unlike the bulk sine-Gordon breathers, no double pole signals the disappearance. The puzzle is not resolved."
- **Not checked:** CD's critical couplings (6.2) and (6.4).

## Item 6. The spectrum section: correct, and more is now checked
- **Agreement with BPTT.** $\nu_n$, $w_k$, the fundamental range, the poles of $\sigma(i\vartheta,u)$, `eq-bsg-excited-k`, the sequence condition and the energies all match BPTT (3.3)–(3.8) and sec. 3.8.
- **Breather factors on excited states.** By the full matrix bootstrap, $B_1,B_2$ on $|0\rangle,|1\rangle$ equal Ghoshal's $R^{(n)}(\bar\eta)\,b^n_k$ to about 1e-15. So BPTT (3.9) is verified for $k=0,1$ and $n=1,2$.
- **Not checked:** which poles on the higher states are bound states and which are Coleman–Thun poles.

## Item 7. The Neumann tower is correct, but line 98 has an error
- **The tower.** At Neumann, $P$ and $Q$ have simple poles at $u=\frac\pi2-\frac{n\pi}\lambda$ and at $u=\frac\pi2$, with energies $M\sin\frac{n\pi}\lambda$.
- **The error.** The poles of $R_0$ are at $u=n\pi/2\lambda$, not at $\frac\pi2-\frac{n\pi}{2\lambda}$ (GZ after (5.21); BPTT sec. 3.4). The value $\frac\pi2-\frac{n\pi}{2\lambda}$ belongs to $R_0^{(n)}$, which is correctly placed at l.144.
- **Replacement for the first sentence of l.98:** "**Poles.** $R_0$ has simple poles at $u=\frac{n\pi}{2\lambda}$, $n=1,2,\dots<\lambda$, where the soliton and its reflected image fuse into $B_n$, absorbed at zero momentum (the one-point pole of (-@eq-one-point-residue)) [@gz1994; @bptt2002]. The factor $\sin\lambda u$ in $Q_0$ removes them from $Q$, and for $\vartheta=0$ the factor $\cos\lambda u$ removes the odd ones from $P^\pm$, as $C$-invariance requires. $\sigma(\eta,u)$ has the boundary-dependent poles discussed in @sec-bsg-spectrum."

## Item 8. Saleur–Skorik–Warner and MacIntyre
- **MacIntyre, l.32: correct.**
- **Notes l.293: wrong.** "Breather" does not occur in MacIntyre. Replacement text: "Saleur, Skorik and Warner [-@ssw1995] studied the classical and semiclassical theory, including boundary breathers, and MacIntyre [-@macintyre1995] derived the general integrable boundary condition from Sklyanin's formalism."
- **SSW, l.40: imprecise.** Replace "They found that the method of images ... reflection matrix of @sec-bsg-k-matrix." with: "They found that the method of images still works classically despite the non-linear boundary condition: the reflected soliton is part of a three-soliton solution on the full line, with an image soliton and a stationary soliton whose position and charge are fixed by $M_0$ and $\phi_0$. In the Dirichlet case their semiclassical phase shifts agree with the $\beta\to0$ limit of the reflection matrix of @sec-bsg-k-matrix and give $\xi_B\propto\phi_0$; they note that the images do not extend to the quantum amplitudes."
- **Line 110: unclear.** The clause about "the consistency of Mattsson and Dorey's Dirichlet spectrum" does not test the dependence on $\phi_0$. Delete it.

## Item 9. DG (3.16) at $n=1$ (optional): unchecked.

## Item 10. Finite-size map
- **The g-theorem attribution (l.188, l.273) is wrong.** Affleck and Ludwig conjectured it; Friedan and Konechny proved it (Phys. Rev. Lett. 93 (2004) 030402, hep-th/0312197).
  - Replacement for l.188: "...with a boundary entropy $\log g$ that decreases along boundary renormalization-group flows, as conjectured by Affleck and Ludwig [@affleck_ludwig1991] and proved by Friedan and Konechny."
  - Replacement for l.273: "the g-theorem was conjectured by Affleck and Ludwig [-@affleck_ludwig1991]".
  - Add `friedan_konechny2004` to the bibliography.
- **Missing (to check):** LMSS already proposed a g-function, which DFRT corrected. The text should say "obtained, correcting [@lmss1995], by DFRT".
- **Line 233, text to append:** "; the two differ in the sign of the one-point coupling, $g_{\Phi(h(0))}=-g_{\mathbb 1}$, which the infrared data must include [@dtw1999]."

## Item 11. Exercises
- **`exr-boundary-cdd` (ch28 l.252) is wrong, and so is l.128.**
  - $f_x(2\theta)$ violates $\Phi_B(\frac{i\pi}2-\theta)=\Phi_B(\frac{i\pi}2+\theta)$ (residual 2.12); $f_x(\theta)$ satisfies it (GZ (3.36)).
  - Replacement for l.128 from "exactly": "exactly the CDD ambiguity of the bulk, in the same variable $\theta$ [@gz1994, eq. (3.36)]. A boundary CDD factor is therefore a product of bulk factors $f_x(\theta)$ of (-@eq-fx-block) (@exr-boundary-cdd); for example $K_b/K_{(\mathbb 1)}=f_{(1+b)/6}(\theta)/f_{(1-b)/6}(\theta)$ in @sec-boundary-lee-yang."
  - Corrected exercise: "Show that $\Phi_B(\theta)=f_x(\theta)$ satisfies $\Phi_B(\theta)\Phi_B(-\theta)=1$ and $\Phi_B(\frac{i\pi}2-\theta)=\Phi_B(\frac{i\pi}2+\theta)$, and that $f_x(2\theta)$ does not. Where are the poles of $f_x(\theta)$ in the strip $0<\operatorname{Im}\theta<\pi/2$?" (Answer: at $i\pi x$ for $0<x<\frac12$, and at $i\pi(1-x)$ for $\frac12<x<1$.)
- **`exr-bsg-sign-q`: muddled.** Suggested wording: "Show that $\check R(-q)=(\Pi\otimes1)\check R(q)(\Pi\otimes1)$ implies that every solution $K$ at $q$ also solves the reflection equation at $-q$, and that $\Pi K(\hat\epsilon)\Pi\propto K(-\hat\epsilon)$..."
- **All the other exercises** are well posed, with correct hints.
- **Minor, `exr-ly-excited`.** For $-1<b<0$, $K_{b+2}$ also has a pole at $i(3+b)\pi/6$ in the strip, which l.233 leaves unexplained.

## Item 12. Bibliography (helper agent)
- **`chenaghlou_corrigan2000`: wrong DOI.** Use `doi = {10.1142/S0217751X0000183X}`.
- **`bdg2002`:** this is a NATO proceedings chapter in *Statistical Field Theories* (eds. Cappelli and Mussardo, Kluwer, Dordrecht 2002), pp. 129–137, doi:10.1007/978-94-010-0514-2_12.
- **`ssw1995`:** the published title is "classical and semi-classical analysis".
- **Optional:** the GZ erratum is IJMPA 9 (1994) 4353, doi:10.1142/S0217751X94002430.
- **To add:** `corrigan_taormina2000`, `friedan_konechny2004`.
- **Citations that do not support their sentences:** ch28:140, ch28:167, ch29:112 and ch29:293. See above.

## Item 13. Prior art
- **The Lee–Yang sign is not new.** DPTW print $+$, so say so.
- **The divergence of GZ (5.23)** is not remarked on in the sources read.
- **Crossing-unitarity fixing the sign of $q$** was not found in the sources read; a wider search is still needed.

## Referee scripts (summary)
| Script | What it checks | Result |
|---|---|---|
| `ref_ly_an.py` | Lee–Yang sign; GZ (3.55), (3.41); the $a_n$ ratio | as above |
| `ref_neumann_bnk.py`, `ref_r0poles.py` | the Neumann poles; $R_0$ poles at $n\pi/2\lambda$; the $b^n_k$ product | as above |
| `ref_an_bootstrap.py` | BPTT (3.4)–(3.5), $n=1,2$ | ratio 1.000000000000 |
| `ref_excited_breathers.py` | $B_{1,2}$ on $\lvert0\rangle,\lvert1\rangle$ | 1e-15 |
| `ref_breather_order.py` | residue image is the singlet | overlap 1 |
| `ref_cdd.py` | $f_x(\theta)$ is a boundary CDD factor; $f_x(2\theta)$ is not | as above |
| `ref_gz522.py` | GZ (5.22)/(5.24) normalization factor | $k^2/\cos^2\xi$ |
| `ref_uvir_shg.py` | continued UV–IR relation equals the CT/CC relation | 1e-16 |

## Integration (2026-10-11)
Integrated into `28-boundary-scattering.qmd` and `29-boundary-sine-gordon.qmd`; new checks appended to `code/boundary_axioms.py` (check 6), `boundary_sine_gordon.py` (checks 9–13) and `shg_boundary_breathers.py` (check 8). All lines PASS (16, 31 and 9; runtimes about 1 s, 55 s and 8 s). The new keys `corrigan_taormina2000` and `friedan_konechny2004` are cited but not yet in `references.bib`.

**Integrated**
- Item 11 / finding 1: boundary CDD factor is $f_x(\theta)$; ch28 text after the CDD display and `exr-boundary-cdd` corrected; check 6 of `boundary_axioms.py` (both conditions, crossing-unitarity of $K_{(\mathbb 1)}f_x$, $K_b/K_{(\mathbb 1)}$, pole positions).
- Item 7 / finding 2: poles of $R_0$ at $u=n\pi/2\lambda$; $\sin\lambda u$ removes them from $Q$, $\cos\lambda u$ the odd ones from $P^\pm$ at $\vartheta=0$; Neumann tower kept. Check 9 of `boundary_sine_gordon.py`.
- Item 2a / finding 3: UV–IR attribution (BPTT sec. 4 first, Zamolodchikov's normalization printed in BPT (3.2), (3.3)); status paragraph extended (FZZ derivation on $\phi_0=0$, TCSA tests only indirectly). Conjecture title changed to "The UV–IR relation"; notes paragraph adjusted.
- Items 1e, 2b / findings 4, 5: $N_0$ derived by CT; CT bullet replaces "we have not compared"; ch29 corrections note updated. Continuation check is check 8 of `shg_boundary_breathers.py`. The text states that the exact match with $C_{0,1}=(\tilde\mu/\tilde\mu_{\rm crit})e^{\pm i\beta_{\rm SG}\phi_0/2}$ holds for $F=(a_1-a_0)(1-\frac B2)$ when $\vartheta=-i\pi F/B$; the other sign of $F$ exchanges $C_0,C_1$.
- Item 2c: `prp-bsg-spectrum` follow-up paragraph reworded (low-lying TCSA, mainly $\phi_0=0$; higher states extrapolated). Label kept.
- Item 3: $a_n$ printed (display after `eq-bsg-excited-k`, blocks identified as `eq-blocks` with $h=1$); `prp-bsg-excited` now says $a_1,a_2$ checked; $b_{mk}$ formula of BPTT (3.9) printed and checked for $m=1,2$, $k=0,1$ (item 6). Checks 12, 13.
- Item 4: residue sentence added to `prp-bsg-breather-k`; check 11.
- Item 5a: residue normalization $\frac i2 gg$ (GZ (3.57), (3.62)), positivity as heuristic (indices relabelled so that $\beta$ stays the bound state).
- Item 5b: $(Q_s-Q_{-s})|B\rangle=0$ with GZ (2.21) and the $i^{\mp s}$ factors.
- Item 5c: DTW counting-rule paragraph replaced; the first bullet already cited MD and BPTT, and now mentions cancellations between diagrams.
- Item 5d: Mattsson–Dorey lemmas paragraph replaced (written with the book's $\pi-u_{ab}^c$); sinh-Gordon puzzle text replaced.
- Item 8: SSW sentence at l.40 replaced; MD clause at l.110 deleted; notes sentence on SSW/MacIntyre replaced.
- Item 10: g-theorem conjectured by AL, proved by Friedan–Konechny (text and notes); "obtained, correcting [@lmss1995], by DFRT"; one-point-coupling remark appended to the Lee–Yang IR paragraph.
- Item 11: `exr-bsg-sign-q` reworded; $\Pi K(\hat\epsilon)\Pi\propto K(-\hat\epsilon)$ added to `boundary_sine_gordon.py`.
- Item 1b: GZ (5.22)/(5.24) normalization note; check 10.
- Item 13: Lee–Yang sign "not new: DPTW print $+$" in `prp-ly-reflection` and the notes.
- `exr-ly-excited` (minor): one sentence noting that the pole at $i(b+3)\pi/6$ is in the strip for $-1<b<0$ and is not explained.

**Not integrated**
- 2a: later literature on the UV–IR relation (INSPIRE search of papers citing hep-th/0108157).
- 5d: CD's critical couplings (6.2), (6.4) (kept as quoted).
- Item 6: which poles on the higher excited states are bound states and which Coleman–Thun poles.
- Item 9: DG (3.16) at $n=1$.
- `exr-ly-excited`: an explanation of the extra pole.
- 1c: GZ (5.9) typo statement lives in ch. 32 (not owned here).
- 1d: no change needed. Item 12 (bibliography) and the new bib entries: for `references.bib`. Item 13: crossing-unitarity fixing the sign of $q$ needs a wider prior-art search; no novelty is claimed in the text.
- Ch. 35 still calls the relation "Zamolodchikov's relation" (not owned here).
