<!-- Worker report, 2026-10-11. Model: Claude Opus 5.5, effort xhigh. Scope: full brief halfline_open_03_neumann_soliton_route.md, including the s=-1 bonus. Interrupted once by the usage limit and resumed. Saved verbatim by the coordinator. The worker's own scripts are in briefs/halfline_open_03_scripts/; they import unchanged copies of textbook_code/*.py, which must be copied next to them. -->

## Report: problem 3, Neumann amplitudes from the s=+1 soliton reflection matrix

**Verdict.** The soliton route works. For $n=1,\dots,5$ I derived the $s=+1$ breather amplitude from (`eq-an-k-matrix`) with an explicit scalar factor. Continuing $\omega\to-2/B$ gives exactly $K^N_1$ of (`eq-neumann-amplitude`), and $K^N_a$ then follows by the bulk bootstrap. So `cnj-neumann-coideal` and `cnj-neumann-amplitude` are confirmed for $n=3,4,5$, with one proviso for odd $n$ (section 4). Nothing in the repo was edited.

The coordinator's "rounding issue" was the weight-label rounding in my weight-restricted R-check solver. It was fixed before the pause and agrees with the book's `rcheck` to 1e-15 for $n=2,3$.

### 1. Two corrections to the brief's premises

- **Boundary unitarity is not $K(\theta)K(-\theta)$.** That product is not scalar for $n\ge2$, nor for $n=1$ with $s=-1$.
  - The antisoliton is reflected by $K_{\bar V}=JKJ=K^{\rm T}$, where $J$ reverses the basis.
  - $J$ is charge conjugation. It comes from the diagram automorphism $j\leftrightarrow n{+}1{-}j$ that fixes node 0, and this automorphism preserves $B_{\hat\epsilon}$ because $\hat\epsilon_j=\hat\epsilon_{n+1-j}$.
  - I checked numerically that $K^{\rm T}$ is the intertwiner $\bar V(\theta)\to V(-\theta)$, and that the soliton bootstrap maps the fused antisoliton onto $\iota_{\rm out}K^{\rm T}$ (residual ≤5e-12).
  - So the physical amplitudes are $A\,K$ for the soliton and $A\,K^{\rm T}$ for the antisoliton, and unitarity is $A(\theta)A(-\theta)\,K(-\theta)^{\rm T}K(\theta)=1$.
- **Branch.** The book's formula needs $X^p=e^{p\omega(\theta-i\pi/2)}$ and $\epsilon_*=1/(2\cos\frac{\pi\omega}2)$, which is analytic in $\omega$.
  - `k_closed` in `boundary_coideal.py` uses `np.sqrt(-q)` and a positive $\epsilon_*$. For odd $N$ and $\cos\frac{\pi\omega}2<0$ that is the other sign class.
  - The book only checks at $\omega=0.37$, so no claim in it is wrong. Its scripts should not be reused at $\omega>1$. I used my own `k_book`.

### 2. Closed forms

Notation: $T=h\omega/2$, $\mu=-iT\theta/\pi$, DG blocks $(a)_I$.

- **Crossing function** $\Phi_n$, defined by $k(\theta)=\Phi\,\check R_{\bar V\bar V}(2\theta)\,k(-\theta)$:
  - $s=+1$: $\Phi_+=\dfrac{\sinh(\frac{i\pi\omega}2-T\theta)}{\sinh(\frac{i\pi\omega}2+T\theta)}=\dfrac{\sin\pi(\frac\omega2-\mu)}{\sin\pi(\frac\omega2+\mu)}$
  - $s=-1$: $\Phi_-=\dfrac{\cos\pi(\mu-\frac\omega2)}{\cos\pi(\mu+\frac\omega2)}$
  - It depends on $n$ only through $T$. Precision 1.3e-14, for $n=1..5$, $s=\pm1$ and $\omega\in\{0.37,1.83,2.6\}$ (also 4.3).
- **Unitarity function** $U_n$:
  - $K(-\theta)^{\rm T}K(\theta)=U\,\mathbf 1$ with $U=\dfrac{s\cosh2T\theta-\cos\pi T}{2\sin^2\frac{\pi\omega}2}$.
  - Equivalently $U_+=\sin\pi(\frac T2+\mu)\sin\pi(\frac T2-\mu)/\sin^2\frac{\pi\omega}2$ and $U_-=-\cos\pi(\frac T2+\mu)\cos\pi(\frac T2-\mu)/\sin^2\frac{\pi\omega}2$. Precision 4.4e-15.
  - **Proved.** $K=K_{ii}\mathbf 1+\sum_{m=1}^{N-1}X^{m-N/2}Z^m$, with $Z$ the twisted shift ($Z^N=s$).
  - On the $\zeta$-eigenvector ($\zeta^N=s$) the eigenvalue is $\kappa(\zeta)=\dfrac{(1-sX^N)(X\zeta+q)}{wX^{N/2}(w-w^{-1})(X\zeta-1)}$, with $w=e^{-i\pi\omega/2}$.
  - $K^{\rm T}(-\theta)$ has eigenvalue $\kappa$ at $X\to w^2/X$, $\zeta\to1/\zeta$. The product does not depend on $\zeta$.
  - The same $\kappa$ reproduces the $\det K$ formula of ch. 34.
- **These are Gandenberger's equations.** For $n=2$ the functional equations below are his (4.4) and (4.6) exactly.
- **Singlet projection.** $(1\otimes K)\check R(1\otimes K^{\rm T})$ maps the $B_1$ singlet to $c_s(\mu)$ times itself, in both orderings. Residual 1.7e-14 for $n=2..5$.
  - $c_+=-\dfrac{\cos\pi\mu\,\cos\pi(\mu+\frac\omega2)\cos\pi(\mu-T)}{\sin^2\frac{\pi\omega}2\cos\pi(\mu-\frac\omega2)}$
  - $c_-=\dfrac{\sin\pi\mu\,\sin\pi(\mu+\frac\omega2)\sin\pi(\mu-T)}{\sin^2\frac{\pi\omega}2\sin\pi(\mu-\frac\omega2)}$

### 3. Scalar factors

**Functional equations.** $A(\mu)A(-\mu)=1/U$ and $A(\frac T2-\mu)/A(\frac T2+\mu)=F_{11}(2\mu)/\Phi(\mu)$.

**Result** (precision 3.5e-8 for $n=1..5$):
- $A^{(+)}=a_0^{(+)}\,a_1\,(\tfrac T2-\tfrac\omega2)_I(\tfrac T2+\tfrac\omega2)_I$
- $A^{(-)}=-i\,(-\tfrac T2)_I\,a_0^{(-)}\,a_1$

The ingredients, with $G_j(\mu,a)=\Gamma(\mu+jh\omega+a)/\Gamma(-\mu+jh\omega+a)$:
- $a_1$ is DG's (3.9) product without its prefactor, so $A_d=(-\frac T2)_I\,a_1$.
- $a_0^{(+)}=\dfrac{\sin\frac{\pi\omega}2}{\sin\pi(\mu-\frac T2)}\prod_{j\ge1}\dfrac{G_j(\mu,\frac\omega2-\frac{3T}2)\,G_j(\mu,1-\frac\omega2-\frac T2)}{G_j(\mu,\frac\omega2-\frac T2)\,G_j(\mu,1-\frac\omega2-\frac{3T}2)}$
- $a_0^{(-)}$ is the same with $\sin\to\cos$ in the prefactor and with $+\frac12$ added to the first and third arguments and the "1" replaced by $\frac12$ in the other two.

**Checks of the scalar factors:**
- $n=1$: $A^{(+)}$ equals the GZ/BPTT factor at the Neumann point $(\eta,\vartheta)=(\frac\pi2(\omega+1),0)$. $A^{(-)}$ equals GZ at $(\frac\pi2,\frac{i\pi\omega}2)$. Precision 7e-8.
- $n=2$: $A^{(+)}$ is Gandenberger's (4.14) times his $\sigma^+=(\frac\lambda4)(\frac{5\lambda}4)$ of (5.3), to 2.5e-9. So $(\frac T2\mp\frac\omega2)_I$ is his $\sigma^+$ for general $n$.
- **Soliton bootstrap** $V^{\otimes n}\to\bar V$, measured relative to $K=1$ so that bulk factors and embedding normalizations cancel. It holds for $A^{(+)}$ at $n=2..5$ (5.8e-9 at $n=5$) and for $A^{(-)}$ at $n=2..4$.
- **The $s=-1$ sign.** $A^{(-)}$ needs no CDD factor, as Gandenberger found for $n=2$. The bootstrap fixes its overall sign for even $n$. My reconstruction of his (4.15) differs by a constant $-i$. The text extracted from his PDF is garbled there, so I did not pursue this. The breathers agree.

### 4. The CDD factor, and the proviso for odd $n$

Unitarity and crossing determine $A^{(+)}=A_d\,\rho_+\,\sigma$, where $\rho_+$ is my Gamma-product solution of the two functional equations and $\sigma$ is a CDD factor. Three things pin $\sigma$ down:

1. **The Neumann pole at $\theta=i\pi/2$.** $K$ is rank one there (all entries 1), which is the structure of a soliton absorbed at zero momentum. The Neumann boundary has degenerate vacua, and GZ's $n=1$ factor has this pole. $A_d$ has a simple zero there, so a factor $(\frac T2)_I^2$ is needed. It is present in $a_0^{(+)}$.
2. **The soliton bootstrap forces $(\frac T2-\frac\omega2)_I(\frac T2+\frac\omega2)_I$ for even $n$.** $(\frac T2)_I^2$ alone fails the bootstrap for even $n$, and this factor repairs it.
3. **For odd $n$ that factor is allowed but not forced.** It satisfies the bootstrap by itself: it is DG's $\sigma_c$ with $c=\frac{n-1}2$, up to sign. Leaving it out changes $K_1$ by DG's $k_{h/2-1}$ of (3.25).
   - That factor is $1+O(B)$ with a nonzero $O(B)$ term, so Kim's $a_3$ one-loop result requires the factor for $n=3$.
   - For $n=5,7$ only uniformity in $n$ selects it.
   - This factor also puts a pole at $\theta=i\pi(\frac12-\frac1h)$ in the physical strip for $n\ge2$. I have not explained that pole (it is in Gandenberger's $n=2$ result too).

### 5. Breathers and continuation

- **$s=+1$.** $A(\mu_a)F_{11}(2\mu)A(\mu_b)\,c_+(\mu)=(-\frac{n\omega}2)_I(\frac{n\omega}2-\frac12)_I(\frac12)_I$, with all infinite products evaluated. Precision 9.5e-8 for $n=1..5$ at $\omega=0.37$ and $0.83$.
  - Under $\omega\to-2/B$ this becomes $(-n)(n+\frac B2)(-\frac B2)=K^N_1$, exact for $n=1..7$ as a block identity.
  - $n=1$ is Ghoshal's $K_{B_1}$ at the Neumann point. $n=2$ is Gandenberger's (5.16)/(6.4).
  - The $K^N_a$ satisfy all axioms of `prp-uniform-amplitude` (i)–(iii) and tend to 1 as $B\to0$ ($n\le7$). Kim's one-loop agreement carries over, since the derived amplitude is the conjectured one.
- **$s=-1$ (bonus).** The breather amplitude is DG(3.11)$\,\times\,\sigma_1(\mu)$ for $n=1..5$.
  - Continued, this gives $K^{(-)}_a=-(a)(h-a)\,K_a^{\rm unif}$, i.e. $K^{(-)}_1=-(1)(h-\frac B2)(-1+\frac B2)$.
  - It satisfies all axioms ($n\le7$), tends to 1 as $B\to0$, equals Gandenberger's (6.4) for $n=2$, and its dual is $-(a)(h-a)K^N_a$.
- **Semiclassics.** $\frac{\beta^2}{4h}\frac{d}{d\theta}(-i\ln A^{(+)})\to\frac12\ln\big(1-\sin^2\frac\pi h/\cosh^2\theta\big)$, which is the Neumann time advance. At $\theta=0.6$ and $\omega=40,80,160$ the errors are 0.021, 0.011, 0.0055 ($n=2$) and 0.017, 0.0085, 0.0043 ($n=3$), falling like $1/\omega$.

### 6. Not proved analytically

$\Phi_n$, $c_\pm$, the identity $\rho(\mu_a)\rho(\mu_b)c_+=(-\frac12)(\frac12-T)(\frac12-T-\frac\omega2)(\frac\omega2-\frac12)$, and the soliton bootstrap are numerical, for $n\le5$ (bootstrap $n\le4$ for $s=-1$). Only $U_n$ is proved.

### 7. Proposed text

**Replacement for "Which soliton reflection matrix?" (`sec-neumann-duality`):**

> **Which soliton reflection matrix?** Delius and Gandenberger obtained the uniform amplitude from the diagonal soliton reflection matrix and the Neumann amplitude only as its dual. The soliton route needs the reflection matrix of the Neumann boundary. Following Delius and MacKay [@delius_mackay2003] one would take the diagonal one. It is not: by @cnj-neumann-coideal the Neumann boundary has $\hat\epsilon_j=\epsilon_*$ and the matrix (-@eq-an-k-matrix) with $s=+1$, with the scalar factor of @prp-an-nondiagonal-scalar. Bootstrapping it to the lowest breather gives
> $$K^{(+)}_{B_1}(\mu)=\Big(-\frac n2\omega\Big)_I\Big(\frac n2\omega-\frac12\Big)_I\Big(\frac12\Big)_I ,$$
> and $\omega\to-2/B$ turns this into $K^N_1$ of (-@eq-neumann-amplitude). The heavier particles follow by the bootstrap, as for the uniform boundary.
>
> **[Numerical]** For $n\le5$: the projection is a genuine matrix computation, since $K$ is not diagonal. The singlet is mapped to $c_+(\mu)$ times itself, and $A^{(+)}(\mu_a)F_{11}(2\mu)A^{(+)}(\mu_b)c_+(\mu)$ equals the expression above to $10^{-7}$. For $n=1$ it is Ghoshal's factor at the Neumann point, and for $n=2$ it is Gandenberger's (5.16) and (6.4).
>
> So for $n\le5$ the conjecture is derived from the soliton side, with one proviso. The scalar factor is fixed by unitarity, crossing-unitarity, the zero-momentum pole at $\theta=i\pi/2$ and the soliton bootstrap, except for odd $n$. There the factor $(\frac T2-\frac\omega2)_I(\frac T2+\frac\omega2)_I$ is allowed but not forced, and without it $K_1$ changes by DG's $k_{h/2-1}$ (3.25). For $n=3$ Kim's one-loop amplitude requires the factor; for $n\ge5$ only uniformity in $n$ selects it. The same scalar factor reproduces the classical Neumann time advance (@sec-boundary-semiclassics). This confirms @cnj-neumann-coideal and @cnj-neumann-amplitude, subject to the proviso.

**For `sec-an-nondiagonal-k`, replacing the sentence "The scalar factor ... has no closed form yet":**

> **Unitarity.** The antisoliton is reflected by $JKJ=K^{\rm T}$, with $J$ the basis reversal. $J$ is charge conjugation: it comes from the diagram automorphism $j\leftrightarrow n{+}1{-}j$, which fixes $\hat\epsilon$. Then
> $$K(-\theta)^{\rm T}K(\theta)=\frac{s\cosh2T\theta-\cos\pi T}{2\sin^2\frac{\pi\omega}2}\,\mathbf 1 .$$
> $K(\theta)K(-\theta)$ is not scalar for $n\ge2$. *Proof:* $K$ is a polynomial in the twisted shift $Z$ ($Z^N=s$), with eigenvalue $\kappa(\zeta)$ on the $\zeta$-eigenvector of $Z$. The product of the two eigenvalues does not depend on $\zeta$.
>
> **Crossing-unitarity.** $k(\theta)=\Phi(\theta)\check R_{\bar V\bar V}(2\theta)k(-\theta)$ with $\Phi=\sinh(\frac{i\pi\omega}2-T\theta)/\sinh(\frac{i\pi\omega}2+T\theta)$ for $s=+1$ and $\cosh(T\theta-\frac{i\pi\omega}2)/\cosh(T\theta+\frac{i\pi\omega}2)$ for $s=-1$. With $A\,K$ and $A\,K^{\rm T}$ as the physical amplitudes, $A(\mu)A(-\mu)=1/U$ and $A(\frac T2-\mu)/A(\frac T2+\mu)=F_{11}(2\mu)/\Phi(\mu)$. For $n=2$ these are Gandenberger's (4.4) and (4.6).
>
> ::: {#prp-an-nondiagonal-scalar}
> **[Numerical]** For $n\le5$, $A^{(\pm)}$ [the formulas of section 3 above] solve these equations. For $n=1$ they are the Ghoshal–Zamolodchikov/BPTT factors at the Neumann point and at $(\frac\pi2,\frac{i\pi\omega}2)$. For $n=2$, $A^{(+)}$ is Gandenberger's $A^+\sigma^+$. The soliton bootstrap $V^{\otimes n}\to\bar V$ holds ($s=+1$: $n\le5$; $s=-1$: $n\le4$).
> :::
>
> The CDD discussion follows (section 4 above). Also change "(checked for $n=2,3$)" in `prp-an-nondiagonal` (iii) to $n\le5$.

**Other places to update:**
- **"The third matrix" (`sec-neumann-duality`):** for general $n$, $K^{(-)}_a=-(a)(h-a)K_a$, the imaginary-coupling breather amplitude is DG(3.11)$\,\sigma_1$, all axioms hold for $n\le7$, it tends to 1, and its dual is $-(a)(h-a)K^N_a$.
- **`sec-breather-reflection`, "Other soliton reflection matrices":** both non-diagonal scalar factors are now known.
- **Last sentence of `sec-boundary-semiclassics`:** the $s=+1$ test is done (section 5 above).
- **Open problem 3:** mark it [Numerical] for $n\le5$. What remains is an analytic proof for all $n$, the odd-$n$ CDD factor, and the pole at $\theta=i\pi(\frac12-\frac1h)$.

### 8. Scripts

All in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/work/` (temp directory, so copy soon). The book scripts were copied there unchanged.
- `nk_common.py` — fast weight-restricted R-check solver, and `k_book` (the K-matrix on the analytic branch).
- `proj.py` — $B_1$ singlet projection coefficient, both orderings.
- `rhonp.py` — $\rho_+$, blocks, $c_+$, breather ratio.
- `rhom.py` — $\rho_-$ and $c_-$.
- `gz.py` — BPTT $R_0$ and $\sigma$, DG $A_1$ (numpy, Richardson extrapolation).
- `sol_boot.py` — soliton bootstrap $V^{\otimes n}\to\bar V$: fusion embeddings and $n$-soliton reflection.
- `check_all.py` — checks 1–9 (about 4 minutes; argument `full` adds the $n=5$ bootstrap).
- `closed_forms.py` — the Gandenberger-type forms of $A^{(\pm)}$, $n=1..5$.
- `particle_checks.py` and `pr_core.py` (the header of `particle_reflection.py`) — continuation, axioms of $K^{(-)}_a$ and $K^N_a$, the $k_c$ identification and its $O(B)$ term.
- `t19.py` and `boot5.log` — the $n=5$ bootstrap run.
- `t6.py` — $\Phi$, $U$ and intertwiner checks up to $\omega=4.3$.
- `fitR.py`, `fitc2.py`, `fitphi.py` — the rational fits that found the closed forms.
- `t26.py`, `t27.py` — exploration of the closed forms.
- `../g98/g98.txt` — text extracted from Gandenberger's hep-th/9806003 PDF (sec. 4, (5.2)–(5.19), app. B).

## Integration (branch half-line-open-problems-integration)

Basis: the worker report, the Gemini verdict, the Opus follow-up `halfline_open_03_referee_opus.md`, and the coordinator's own checks. The new `textbook_code/an_nondiagonal_scalar.py` (7 checks, all PASS, about 1 min) is built from the Opus follow-up scripts and the worker's closed forms, and was run by the coordinator.

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 1 (unitarity is $K(-\theta)^{\rm T}K(\theta)$; antisoliton $=JKJ=K^{\rm T}$) | **Verified, integrated**; refined as in the Opus follow-up ($JKJ$ needs $\tau$-invariant $\hat\epsilon$; $K^{\rm T}$ holds in the parity-symmetric gauge) | Check 1: $J$ is charge conjugation (exact); $K^{\rm T}$ intertwines $\bar V(y)\to V(1/y)$; $K(\theta)K(-\theta)$ not scalar | ch. 28 after `eq-boundary-unitarity`; ch. 33 procedure step 4; ch. 33 new "Unitarity" paragraph |
| 2 (`k_closed` branch) | **Verified** (Gemini; the book's statements are unaffected since the book uses $\omega=0.37$) | The new script uses its own `kmat` on the analytic branch | No book change; not mentioned |
| 3 ($U_s$; $\kappa(\zeta)$ proof) | **Verified** | $U_s$ in check 1; $\kappa(\zeta)$ against the spectrum of $K$ (coordinator, $n=2,3,4$, two $\omega$, $s=\pm1$, $\le2\times10^{-15}$) | ch. 33 "Unitarity" with the proof |
| 4 ($\Phi$, renamed $\Psi_s$ because the book's $\Phi$ is taken) | **Verified** | Check 6: closed form against the crossing matrix, $n=2,3,4$, $\omega=0.37,1.29$, $2\times10^{-14}$ | ch. 33 "Crossing-unitarity", new `eq-an-nondiagonal-scalar-eqs`; `prp-an-nondiagonal` (iii) no longer says "checked for $n=2,3$" |
| 5 ($c_+$) | **Verified** | Check 2, $n=2,3,4$ | in `prp-an-nondiagonal-scalar` |
| 6 ($A^{(+)}$), 9 (breather identity), continuation | **Verified** for $n=2,3,4$ | Check 3 (Stirling/Hurwitz tails, 30 digits: $4\times10^{-27}$); Opus integral method $10^{-36}$–$10^{-40}$; check 4 (continuation, $n=1..7$) | New `prp-an-nondiagonal-scalar` with `eq-an-nondiagonal-scalar`, `eq-an-neumann-breather`; ch. 35 "Which soliton reflection matrix?" rewritten |
| 6 ($A^{(-)}$, $s=-1$), 10 (third matrix, $K^{(-)}_a=-(a)(h-a)K_a$) | **Not integrated** | Only reruns of the worker's scripts (Gemini); no independent check | Book says these are not yet established for $n\ge3$ (open problem 3; ch. 35 line 45) |
| 7 (soliton bootstrap) | **Not integrated** as a claim | Not rechecked by the coordinator | Not mentioned |
| 8 (CDD: Kim selects the factor at $n=3$; equals $-\sigma_{(n-1)/2}$ for odd $n$; unexplained pole at $i\pi(\frac12-\frac1h)$) | **Verified** (Kim: check 5, $w/o$ factor deviates by 0.66; the $-\sigma_c$ identity by algebra) | | ch. 33 "CDD factors" paragraph. The even-$n$ bootstrap statement is left out (not rechecked). |
| 11 (semiclassics) | **Verified** | Check 7 reproduces the worker's errors exactly | ch. 33 semiclassics paragraph |
| 12 (text) | Integrated as modified above. Also: ch. 33 "Which boundary"; ch. 32 support sentence after `cnj-neumann-coideal`; open problem 3; Appendix D Delius–MacKay "Neumann and $\hat\epsilon=0$" entry (still *Disputed*). `cnj-neumann-amplitude` keeps its [Conjecture] label: the derivation assumes `cnj-neumann-coideal` and `cnj-reflection-continuation`. | | |
| Prior art | None found (Opus follow-up item 4) | | |
