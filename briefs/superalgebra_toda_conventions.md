# Conventions dictionary: affine Toda theories based on Lie superalgebras

Chapter `part4-real-coupling/16-superalgebra-toda.qmd` (`sec-super-toda`), drafted 2026-10-23. For every convention imported into the chapter: the source and its equation, the book's equation, and the transformation between them. "EH91" Evans–Hollowood NPB 352 (1991); "EH92" Evans–Hollowood PLB 293 (1992), hep-th/9207073; "EM96/97/98" Evans–Madsen hep-th/9608190, hep-th/9703065 (v2), hep-th/9712227; "DGPZ-L" Delius–Grisaru–Penati–Zanon PLB 256 (1991) 164; "DGPZ" the same authors, NPB 359 (1991) 125; "PZ" Penati–Zanon hep-th/9205045; "G98" Gandenberger hep-th/9703158; "IZ" Ito–Zhu 2206.08024. Script names refer to `code/super_toda/`.

## 1. Spacetime and light cone

| Item | Source | Book | Transformation / check |
|---|---|---|---|
| Light-cone coordinates | EH92 after (1.1): $x^\pm=\frac12(t\pm x)$; EM97, EM98 do not restate it and use the same superspace | `sec-super-superspace`: $x^\pm=\frac12(t\pm x)$, $\partial_\pm=\partial/\partial x^\pm=\partial_t\pm\partial_x$ | identical. The book's $\partial_\pm$ of `sec-lagrangians` is $\partial_t\pm\partial_x$, so $\frac12\partial_+\phi\partial_-\phi=\frac12(\dot\phi^2-\phi'^2)$ is the standard kinetic term. Check: EM97 (5.1) $L=\partial_+\phi\partial_-\phi-m^2\sinh^2\phi$ gives mass $m$, as they state |
| Overall normalization of EM Lagrangians | EM97 (4.2), (5.1); EM98 (1.1)–(1.2): $L_{\rm sinh}(\phi;m)=\partial_+\phi\partial_-\phi-\frac{m^2}2\cosh2\phi$, $L_{\rm free}(\psi)=i\psi_+\partial_-\psi_++i\psi_-\partial_+\psi_-$ | book kinetic terms $\frac12\partial_+\phi\partial_-\phi$, $\frac i2\psi_+\partial_-\psi_+$ | EM's are twice the book's; their "normalized fields" are the book's divided by $\sqrt2$ (used for the D(2,1;α) Yukawa coupling, `components.py`) |
| PZ light cone | PZ (1): $z=x^+=(x^0+x^1)/\sqrt2$, $D^2=i\partial$ | not used for formulas | only PZ's masses and roots are quoted, both convention-independent up to overall scale |
| DGPZ | DGPZ (1.5), App. A: Minkowski $\eta={\rm diag}(1,-1)$, $\gamma^0=\sigma_1$, $\gamma^1=i\sigma_2$-type, $\beta^2\mathcal L=-\frac12\phi\Box\phi+\frac12i\bar\psi\gamma\cdot\partial\psi-\dots$ | only masses, amplitudes and roots are imported | masses quoted in DGPZ units ($m_n^2=2$); couplings via $\Omega$ (§5) |

## 2. Superspace

| Item | Source | Book | Transformation / check |
|---|---|---|---|
| Grassmann coordinates, derivatives | EH92 (1.3)–(1.4), EM97 §3, EM98 §3.1: $D_\pm=\partial/\partial\theta^\pm-i\theta^\pm\partial_\pm$, $D_\pm^2=-i\partial_\pm$, $Q_\pm=\partial/\partial\theta^\pm+i\theta^\pm\partial_\pm$, $Q_\pm^2=i\partial_\pm$ | `eq-super-derivatives` | identical; $\partial/\partial\theta$ acts from the left. `superspace.py` checks the algebra |
| Other superspace conventions | EH91 (3.1): $D=\partial_\theta+\theta\partial$, $D^2=\partial$ (complexified); DGPZ (1.1): $D_\pm^2=i\partial_\pm$; IZ (2.1): $D=\partial_\theta+\theta\partial$ | — | not used for imported formulas; EH91 and DGPZ zero-curvature structure (1.1)–(1.3) is convention-independent and is re-derived in the book's conventions (`lax_currents.py`) |
| Superfield | EM97, EM98, EH92 (1.2): $\Phi=\phi+i\theta^+\psi_++i\theta^-\psi_-+i\theta^+\theta^-\sigma$ (EH92 writes $F$ for $\sigma$) | `eq-super-field` | identical |
| Berezin integral | not stated in the sources | $\int d^2\theta\,\theta^+\theta^-=i$ (`superspace.berezin`) | fixed so that EM97 (3.2) gives EM97 (3.4) with kinetic term $+\frac12\partial_+\phi\cdot\partial_-\phi$; then the fermion kinetic and Yukawa signs and $+\frac12\sigma^2$ agree with EM97 (3.4) (`components.py`, check 1) |
| Supersymmetry count | EH91, EH92, EM: "N=1", "N=2"; EM96: "(1,0)" | text: $(1,1)$ called $N=1$, $(2,2)$ called $N=2$ | — |

## 3. Action, field equation, components

| Item | Source | Book | Transformation / check |
|---|---|---|---|
| Superspace action | EM97 (3.2): $L=\frac i2D_+\Phi\cdot D_-\Phi+\sum_\tau\mu_je^{\alpha_j\cdot\Phi}+i\theta^+\theta^-\sum_\varepsilon\mu_je^{\alpha_j\cdot\Phi}$ | `eq-super-action` | coupling restored: $\Phi_{\rm EM}=\beta\Phi$, $L_{\rm book}=L_{\rm EM}/\beta^2$ |
| EH92 action | EH92 (1.9): $L=\frac i2D_+\Phi\cdot D_-\Phi-\beta^{-2}\sum\mu_je^{\beta\alpha_j\cdot\Phi}$ | — | $\mu_j^{\rm book}=-\mu_j^{\rm EH}$. Used in `n2_solitons.py` with this sign |
| Field equation | EM97 (3.1) | `eq-super-eom` | $\Phi\to\beta\Phi$; `lax_currents.py` checks that its components are the Euler–Lagrange equations of the component Lagrangian (super sinh-Gordon) |
| Component Lagrangian | EM97 (3.4) | `eq-super-components` | $\alpha\to\beta\alpha$, $\mu\to\mu/\beta^2$; checked symbolically for generic data (`components.py` 1–2) |
| Auxiliary field | EH92 (5.6): $F=-\frac\mu\beta\sum\alpha_je^{\beta\alpha_j\cdot\phi}$ | $\vec\sigma=\frac1\beta\sum_\tau\mu_j\vec\alpha_je^{\beta\vec\alpha_j\cdot\vec\phi}$ | same with $\mu^{\rm book}=-\mu^{\rm EH}$ |
| Mass matrices | — | $M_F=\sum\mu_j\alpha_j\alpha_j^{\mathsf T}G$, $M_B^2=G^{-1}{\rm Hess}V=M_F^2$ | derived (`exr-super-masses`), checked in `components.py` |

## 4. Roots, parity, invariant form, algebra names

| Item | Source | Book | Transformation / check |
|---|---|---|---|
| Node colours | EM97 §2: white even, black odd non-null, grey odd null | same words | — |
| Cartan matrix | EM97 (2.1): $a_{ij}=\alpha_i\cdot\alpha_j$ (symmetric); EH91 (3.10): $[e_i^+,e_j^-]=\delta_{ij}\alpha_i\cdot H$ | `eq-super-chevalley`, $a_{ij}=\vec\alpha_i\cdot\vec\alpha_j$ | identical; DGPZ (1.3) same |
| Sign of the form | EM97 B(1,1): $a=\begin{pmatrix}1&-1&-1\\-1&0&2\\-1&2&0\end{pmatrix}$; PZ: $\alpha_1=(1,-i)$, $\alpha_2=(0,i)$, $\alpha_0=(-1,-i)$, bilinear form, giving $-a_{\rm EM}$ | text notes the overall sign is a convention | checked in `components.py` (B(1,1) realization) |
| Super sinh-Gordon | EM97 Table 4 (2): $A(1,0)^{(2)}=C(2)^{(2)}$, $a=\begin{pmatrix}2&-2\\-2&2\end{pmatrix}$, both odd; IZ §3: $C(2)^{(2)}=\mathfrak{osp}(2|2)^{(2)}$; EH91 §7 attributes it to a deformation of $\mathfrak{osp}(1|2)$ super-Liouville | `sec-super-shg`: $\vec\alpha_1=-\vec\alpha_0$, $\vec\alpha^2=\frac12$ | $\vec\alpha_{\rm book}=\vec\alpha_{\rm EM}/2$; chosen so that $(2\vec\alpha)^2=2$ (bosonic sector in book convention). The untwisted $\mathfrak{osp}(1|2)^{(1)}=B(0,1)^{(1)}$ is excluded by the check $\{F_-,F_-\}\ne0$ (`lax_currents.py` 2) and by EM97 Table 2 (it has an even simple root) |
| $C^{(2)}(n+1)$, $A^{(2)}(0,2n-1)$ | DGPZ notation; EM97 Table 2(b): $C(m+1)^{(2)}$ (B*-B*), $A(2m+1,0)^{(2)}$ (B*-D), $m=n$ resp. $m=n-1$ | DGPZ notation in text, EM's in table caption | root data reconstructed in the orthonormal basis (`components.py` 7), DGPZ normalization: even simple roots squared length 2, odd $=1$ |
| $\mathfrak{sl}(m,m)^{(1)}$ | EH92 (2.8), (2.13), (5.1) | `sec-super-n2` | Kac $A(m-1,m-1)=\mathfrak{psl}(m|m)$; EH91 §3.2 on the centre; EM97 write $A(m-1,m-1)^{(1)}$; the projected roots are realized with the bilinear Euclidean form on traceless $(\Phi^{(1)},\Phi^{(2)})$ (`n2_solitons.py` 1, reproduces (2.8)) |
| $D(2,1;\alpha)^{(1)}$ | EM98 (2.7)–(2.12), $\epsilon_1^2=-\frac{1+\alpha}2$, $\epsilon_2^2=\frac12$, $\epsilon_3^2=\frac\alpha2$ | `sec-super-n2` | identical (`components.py` 5) |

## 5. Masses and couplings

| Item | Source | Book | Transformation |
|---|---|---|---|
| Super sinh-Gordon mass | EM98 (1.5): $\frac i2D_+\Phi D_-\Phi-m\cosh\Phi$ | `eq-ssg-action`: $\frac i2D_+\Phi D_-\Phi+\frac{4m}{\beta^2}\cosh\frac{\beta\Phi}{\sqrt2}$ | $\hat\Phi=\beta\Phi/\sqrt2$, $\mathcal L_{\rm book}=\frac2{\beta^2}[\frac i2D_+\hat\Phi D_-\hat\Phi-M\cosh\hat\Phi]$ with $M=-2m$. Bosonic truncation is the book's $a_1^{(1)}$ with mass parameter $m$ and particle mass $2m$ (`components.py` 2) |
| DGPZ coupling | DGPZ (1.5): $\beta^2\mathcal L$ with even simple roots of squared length 2; (2.38) $\gamma^2=\sin^2\pi\Omega$; text after (2.43): $\Omega=\beta^2/4\pi+O(\beta^4)$ | `eq-dgpz-amplitudes`, block parameter $B=2\Omega$ | book normalization (longest root squared length 2) divides all roots by $\sqrt2$, so $\beta_{\rm book}^2=2\beta_{\rm DGPZ}^2$ and $\Omega=\beta_{\rm book}^2/8\pi+O(\beta^4)$. DGPZ's $\{x\}$ (2.48) is the book's `eq-blocks` with $B=2\Omega$ |
| Sign of $\gamma$ | DGPZ fix only $\gamma^2$ | $\gamma=-\sin\pi\Omega$ | fixed by the tree amplitude computed from the potential (`dgpz_smatrix.py` 8) and consistent with DGPZ (2.5) |
| $N=2$ soliton masses | EH92 (5.4) | `eq-n2-masses` = `eq-soliton-mass` | EH's $\mu$ is the book's $m$; check 2–4 of `n2_solitons.py` |
| $D(2,1;\alpha)$ masses | EM98 (3.23): $m=2\mu$ | text | identical |
| Gandenberger | G98 (3.5)–(3.6): $\mu=\theta/2\pi i$, $\lambda$; special point $\lambda=1/h$, $\beta^2=4\pi h/(h+1)$ in his normalization | `sec-known-smatrices-notes` | only identities at $\lambda=1/h$ are quoted |

## 6. Scattering conventions

| Item | Source | Book | Transformation / check |
|---|---|---|---|
| Amplitudes | DGPZ (2.4), App. A (A.7)–(A.9): out-states with momentum order reversed; $D=$ out$\langle\psi\psi|\phi\phi\rangle_{\rm in}=-$out$\langle\phi\phi|\psi\psi\rangle_{\rm in}$ | `eq-dgpz-operator` | $\check S$ maps in-basis $a\otimes b$ to out-states in spatial order (braid form, as `eq-r-check` of ch. 7). Free limit: graded permutation |
| Graded permutation / tensor order | — | $\check S_{12}=\check S\otimes1$, $\check S_{23}=1\otimes\check S$ with the graded tensor product; for the even $\check S$ these are plain Kronecker products | `dgpz_smatrix.py` 3 checks the braid equation and its failure for the ungraded free limit |
| Yang–Baxter | DGPZ (2.7) (component form, with $b\leftrightarrow f$ partners) | `eq-dgpz-ybe` (braid form, ordering of `eq-ybe-braid`) | both checked for the solution |
| Unitarity | DGPZ (2.6) with real analyticity (2.30), $D^*(\theta)=-D(-\theta^*)$ | $\check S(\theta)\check S(-\theta)=1$; $\check S\check S^\dagger=1$ for real $\theta$ | checked |
| Crossing / charge conjugation | DGPZ (2.31): $T(i\pi-\theta)=T(\theta)$, $R(i\pi-\theta)/T(i\pi-\theta)=-iD(\theta)/T(\theta)$; DGPZ-L: "$D(i\pi-\theta)=iR(\theta)$" | text | the factor $i$ is the Majorana crossing phase $u(\theta+i\pi)=iv(\theta)$ with DGPZ (A.2) spinors; checked. Both particles are self-conjugate (real boson, Majorana fermion) |
| Free-fermion condition | DGPZ (B.1) printed $BF+T^2-D^2-R^2=0$ | $BF+T^2+D^2-R^2=0$ ($d_1d_2=D\cdot(-D)$) | misprint, Appendix D |
| $C^{(2)}(n+1)$ fermions | DGPZ (3.7), (3.11): $\chi=(\psi_1+i\psi_2)/\sqrt2$ ($n$ even), $\chi_\pm=(\psi_1\pm\psi_2)/\sqrt2$ ($n$ odd) | `eq-dgpz-cn` | identical |

## 7. Supercharges, central charge, currents

| Item | Source | Book | Transformation / check |
|---|---|---|---|
| Second supersymmetry | EH92 (2.1)–(2.5): $Q'_\pm=JD_\pm$, $J^2=-1$, $J$ antisymmetric, $J\alpha_j=\lambda_j\alpha_j$ | `eq-n2-condition` | identical |
| First-order soliton equation | EH92 (5.7)–(5.8), $\lambda=\sqrt{(1-v)/(1+v)}$, $\sigma=\tan(\frac{\pi a}{2m}+\frac\pi4)$ | text (not displayed) | holds componentwise with $J\alpha_j=i(-1)^j\alpha_j$; $\lambda=-\sqrt{(1-v)/(1+v)}$ for $\rho>0$ and $+\sqrt{(1-v)/(1+v)}$ for $\rho<0$ (`n2_solitons.py`). EH92 (5.3) fixes only $\rho^2$, so the sign is an orientation choice, not a correction (referee integration 2026-10-24). The chapter writes (5.9) with $c_a=i\,e^{i\pi a/m}=(1+i\sigma)/(1-i\sigma)$, regular at $a=m/2$ |
| Central charge / mass | EH92 (5.9) | `eq-n2-central` | identical; checked on the asymptotic vacua |
| Superspace conservation law | EM98 (5.44) printed $D_+J=D_-\bar J$; (5.47)–(5.51) components | `eq-ssg-currents`: $D_-J=D_+\bar J$ | EM's components (5.48)–(5.51) are those of $D_-J=D_+\bar J$ (derived in the report); partners $-iMD_+e^{-\hat\Phi}$ and $iM(\cdots)e^{-\hat\Phi}$ for $L=\frac i2D_+\hat\Phi D_-\hat\Phi-M\cosh\hat\Phi$; misprint, Appendix D |
