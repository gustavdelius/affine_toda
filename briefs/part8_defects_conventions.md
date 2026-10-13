# Part VIII (defects): conventions, names and id registry

This file is the single source of truth for the conventions of Part VIII, "Affine Toda theory with defects" (`part8-defects/`, chapters 36–43). It was written on 2026-10-11, before drafting. Every agent drafting or refereeing a chapter of Part VIII reads it first. Where a convention is marked *verified*, the named check in `textbook_code/` derives it; where it is marked *to verify*, the first agent that needs it confirms it by script and reports.

The plan of the part is recorded in the memory note `defects-part` and in the plan file `~/.claude/plans/system-reminder-you-are-operating-indexed-stroustrup.md`.

## Placement and numbering

- Part VIII sits between Part VII (the half-line, ch. 28–35) and Foundations. On 2026-10-12 Foundations became Part IX (`part9-foundations/`, ch. 44–49) and the outlook chapter 50.
- Chapters mirror Part VII: 36 ↔ 28 (axioms), 37 ↔ 29 (rank one), 38 ↔ 30 (classical conditions), 39 ↔ 31 (solitons), 40 ↔ 32 (quantum groups), 41 ↔ 33 (soliton matrices), 42 ↔ 34 (quantum symmetric pairs), 43 ↔ 35 (particles and open problems).

## Names and symbols (clash table)

| Use in Part VIII | Symbol / word | Why not the obvious alternative |
|---|---|---|
| Transmission matrix, left to right | $\mathbb T_\mu(\theta-\theta_d)$ | $T$ is already the topological charges $T_j$, the gradation $T$ and transfer matrices $T(\lambda)$; $\mathsf T$ is the transpose; $\mathcal T$ is the double-row matrix (`eq-double-row`) |
| Transmission matrix, right to left | $\bar{\mathbb T}_\mu$ | |
| Scalar particle transmission factor | $\mathbb T_a(\theta)$ | $X_{ab}$ is the fluctuation transmission factor of Part V (`eq-an-transmission`) |
| Defect state space | $\mathcal W$ | $\mathcal D$ and $D$ are used for differential operators and Hirota derivatives |
| Fields left / right of the defect | $\vec\phi$ on $x<0$, $\vec\psi$ on $x>0$ | $u,v$: $u=i\theta$ throughout Part VII; $\phi_L$, $\phi_R$: $\phi_L$ is the chiral field of ch. 20 |
| Defect potential | $\mathcal B_d(\phi,\psi)$ | it becomes a boundary potential under the doubling trick (cf. $\mathcal B$ in `eq-boundary-action`) |
| Defect matrix in the defect conditions | $\mathsf D$ (sans serif), $E=1-\mathsf D$ | plain $D$ is overloaded; sans serif follows $\mathsf M^2$ |
| Defect parameter | $\sigma=e^{-\theta_d}$, $\theta_d$ the defect rapidity | $\eta$, $\vartheta$ are the BPTT boundary parameters ($\vartheta$ is also the boost) |
| Classical defect matrix | $L_d(\lambda)$ | $K$ is Part VII's reflection matrix |
| Defect types | type I$_-$ (BCZ's original; $U_q(\mathfrak b^-)$; $\tilde L=1+\sigma\lambda^{-1}\sum_jE_{-\alpha_j}$), type I$_+$ (its time reverse; $U_q(\mathfrak b^+)$), type II (Corrigan–Zambon's, with an auxiliary field) | Both I$_\pm$ keep modified $Q_{+2}$ *and* $Q_{-2}$ (verified; contradicts BCZ 2004 §5), so the types are named by their Borel subalgebra, not by a surviving charge |
| Defect → boundary of the doubled theory | "the doubling trick" (`sec-doubling-trick`); footnote that the literature says "folding trick" | "folding" means Dynkin-diagram folding throughout the book (`sec-affine-folding`, `sec-images-cn`) and Robertson's folding of defects (`sec-defect-folding`) |
| Delfino–Mussardo–Simonetti | bib key `delfino_mussardo_simonetti1994` | `dms2001` is Delius–MacKay–Short |

## Conventions

Verified on 2026-10-11 by the conventions agent of phase 0. Status key: [V:script#n] = derived or verified by check n of that script; [hand] = derived by hand only. Scripts: SG = `defect_classical_sg.py`, AN = `defect_classical_an.py`, Q = `defect_sg_transmission.py`; helper `defect_common.py`, whose module docstring repeats these conventions.

### 1. Geometry, action, defect conditions (`eq-defect-action`, `eq-defect-conditions`)

- $\vec\phi$ on $x<0$, $\vec\psi$ on $x>0$, defect at $x=0$. Bulk Lagrangian of `eq-imaginary-potential`; longest roots have $\alpha^2=2$.
- $S=\int dt\,\{\int_{-\infty}^0\mathcal L(\phi)\,dx+\int_0^\infty\mathcal L(\psi)\,dx+[\tfrac12\phi\cdot E\dot\phi+\phi\cdot\mathsf D\dot\psi+\tfrac12\psi\cdot E\dot\psi-\mathcal B_d(\phi,\psi)]_{x=0}\}$.
- Variation gives, at $x=0$ [V:AN#2 generic $\mathsf D$ with $E$ antisymmetric; V:SG#1 for $a_1$]:
  $\partial_x\phi=E\partial_t\phi+\mathsf D\partial_t\psi-\nabla_\phi\mathcal B_d$, $\quad\partial_x\psi=\mathsf D^{\mathsf T}\partial_t\phi-E\partial_t\psi+\nabla_\psi\mathcal B_d$.
  (The bulk boundary terms are $-\phi_x\delta\phi$ from $x<0$ and $+\psi_x\delta\psi$ from $x>0$.)
- **Type I$_-$** (BCZ's original): $\mathsf D+\mathsf D^{\mathsf T}=2$, $E=1-\mathsf D$ antisymmetric. For $a_1$: $\mathsf D=1$, $E=0$.
- **Type I$_+$** is the time reverse of I$_-$: $(\mathsf D,E)\to(-\mathsf D,-E)$. The parity image $(\mathsf D,E,\mathcal B(\phi,\psi))\to(-\mathsf D^{\mathsf T},E,\mathcal B(\psi,\phi))$ is the same up to the diagram reflection $\pi\to\pi^{-1}$ [hand]. For $a_1$ the two coincide.
- The quadratic $\dot q$ terms in $dP/dt$ cancel for any antisymmetric $E$ with $\mathsf D=1-E$ [hand; implicit in V:AN#3].

### 2. The matrix $\mathsf D$ for $a_r$ (`eq-d-matrix`) [V:AN#1, $r=1,\dots,5$]

- Roots $\alpha_j=e_j-e_{j+1}$ in $\mathbb R^h$ ($j=0,\dots,r$, indices mod $h=r+1$), so $\alpha_0=-\sum_{j\ge1}\alpha_j$.
- Fundamental weights $\lambda_i$ ($i=1,\dots,r$) with $\alpha_i\cdot\lambda_j=\delta_{ij}$, and $\lambda_0=\lambda_{r+1}=0$.
- $\mathsf D=2\sum_{i=1}^r(\lambda_i-\lambda_{i+1})\lambda_i^{\mathsf T}$; then $\mathsf D\alpha_j=2(\lambda_j-\lambda_{j+1})$ and $\mathsf D\alpha_0=-2\lambda_1$.
- $\mathsf D=2(1-\pi)^{-1}$, with $\pi$ the cyclic rotation $\pi\alpha_j=\alpha_{j-1}$ ($\pi e_k=e_{k-1}$ on $\mathbb R^h$; $\pi$ is orthogonal and has no eigenvalue 1 on the root space).
- $\mathsf D^{\mathsf T}\alpha_j=-\mathsf D\alpha_{j-1}$; $\mathsf D$ is normal; $\alpha_i\cdot\mathsf D\alpha_j=2\delta_{ij}-2\delta_{i,j+1}$.
- $E=\sum_i(\lambda_{i+1}\lambda_i^{\mathsf T}-\lambda_i\lambda_{i+1}^{\mathsf T})$ (BCZ04 (3.30)).
- On the Fourier mode $k$ of $\pi$, $\mathsf D=1+i\cot(\pi k/h)$.

### 3. Defect potential (`eq-defect-potential`, `eq-sg-defect-potential`)

- **Notation.** $f_j=\exp(b\,\alpha_j\cdot(\mathsf D^{\mathsf T}\phi+\mathsf D\psi)/2)$ and $g_j=\exp(b\,\alpha_j\cdot\mathsf D(\phi-\psi)/2)$, with $b=\beta$ (real coupling) or $i\beta$ (imaginary). Here $\mathsf D$ is always the type-I$_-$ matrix.
- **Potential.** $\mathcal B_d^{\rm I_-}(\sigma)=(m/b^2)\sum_{j=0}^r[\sigma(f_j-1)+\sigma^{-1}(g_j-1)]$ and $\mathcal B_d^{\rm I_+}(\sigma)=\mathcal B_d^{\rm I_-}(1/\sigma)$. So $\mathcal B_d^{\rm imag}=\mathcal B_d^{\rm real}(\beta\to i\beta)$, as for $V$ and $\mathcal B$ in Part VII, and $\mathcal B_d(0,0)=0$.
- **Sine-Gordon** ($\beta_{SG}=\sqrt2\beta$, $m_0=2m$), both types: $\mathcal B_d=-(2m_0/\beta_{SG}^2)[\sigma(\cos(\beta_{SG}(\phi+\psi)/2)-1)+\sigma^{-1}(\cos(\beta_{SG}(\phi-\psi)/2)-1)]$ for I$_-$; for I$_+$ swap $\sigma\leftrightarrow\sigma^{-1}$ [V:SG#2, #3].
- **Integrability.** $\tfrac12|\nabla_\phi\mathcal B_d|^2-\tfrac12|\nabla_\psi\mathcal B_d|^2=V(\phi)-V(\psi)$, with prefactor exactly $m/\beta^2$ and no dependence on $\sigma$ [V:SG#2 sympy; V:AN#3 for $r=1,\dots,4$, both couplings, both types].
  - Proof sketch [hand, agrees with the numerics]: the cross terms give $(m^2/2\beta^2)\cdot2\sum_if_i(g_i-g_{i+1})$, and $f_ig_i=e^{\beta\alpha_i\cdot\phi}$, $f_ig_{i+1}=e^{\beta\alpha_i\cdot\psi}$ by $\mathsf D^{\mathsf T}\alpha_{i+1}=-\mathsf D\alpha_i$.
- **Momentum.** $P=-\int\phi_t\cdot\phi_x\,dx$ (so `eq-an-soliton` has momentum $M\sinh\theta$).
  - $P_d^{\rm I_-}(\sigma)=-(m/b^2)\sum_j[\sigma(f_j-1)-\sigma^{-1}(g_j-1)]$ and $P_d^{\rm I_+}(\sigma)=-P_d^{\rm I_-}(1/\sigma)$.
  - $E+\mathcal B_d$ and $P+P_d$ are conserved [V:AN#3 algebraically; V:SG#5 along soliton solutions, relative spread $5\times10^{-16}$].
- **Light cone.** $E+P$ (book spin $+1$) gets $2(m/b^2)\sum\sigma^{-1}(g_j-1)$; $E-P$ gets $2(m/b^2)\sum\sigma(f_j-1)$.

### 4. Defect rapidity, boosts, parity

- **Defect rapidity.** $\sigma=e^{-\theta_d}$. A boost $\theta\to\theta+\vartheta$ is compensated by $\theta_d\to\theta_d+\vartheta$: $E+P$ scales as $e^{\vartheta}$, and its defect term is $\propto\sigma^{-1}$ [hand].
- **Stored soliton.** A sine-Gordon jump $\psi-\phi=2\pi/\beta_{SG}$ stores $\mathcal B_d=M\cosh\theta_d$ and $P_d=M\sinh\theta_d$, with $M=8m_0/\beta_{SG}^2$: a soliton of rapidity $\theta_d$ [V:SG#5].
- **Soliton through the defect.** The book soliton is $\phi=(2i/\beta_{SG})\ln((1-E)/(1+E))$ with $E=\exp(m_0(x\cosh\theta-t\sinh\theta)+\xi)$; $\psi$ is the same with $E\to zE$.
  - $z_{\rm I_-}(\theta)=(e^{-\theta}+\sigma)/(e^{-\theta}-\sigma)=-\coth((\theta-\theta_d)/2)$: a soliton becomes an antisoliton for $\theta>\theta_d$.
  - $z_{\rm I_+}(\theta)=+\coth((\theta-\theta_d)/2)$: conversion for $\theta<\theta_d$.
  - Verified by root-finding and at 4 times × 3 rapidities with complex $\xi$, residual $6\times10^{-15}$ [V:SG#4].
  - At $\theta=\theta_d$ both $z$ diverge and the soliton is absorbed, $\psi\equiv2\pi/\beta_{SG}$ [V:SG#4].
- **Linear transmission** (`eq-shg-linear-transmission`). Plane waves are $e^{-i(Et-px)}$ with $E=m_0\cosh\theta$, $p=m_0\sinh\theta$, and $\mathbb T^{\rm cl}$ is the transmitted over the incident amplitude.
  - $R=0$ and $\mathbb T^{\rm cl}(\theta)=-i\sinh(\tfrac12(\theta-\theta_d)-i\pi/4)/\sinh(\tfrac12(\theta-\theta_d)+i\pi/4)$, for sinh-Gordon and the even sine-Gordon vacua [V:SG#6 sympy].
  - $|\mathbb T^{\rm cl}|=1$, $\mathbb T^{\rm cl}(\theta_d)=i$, and $\mathbb T^{\rm cl}\to1$ as $\theta_d\to+\infty$.
  - About $(\phi,\psi)=(0,2\pi/\beta_{SG})$: $\mathbb T_{\rm odd}=\overline{\mathbb T_{\rm even}}$ (BCZ05 (2.9)).

### 5. Classical defect matrix and the types (`eq-defect-lax-matrix`, `prp-defect-lax`, `prp-spin-two`)

**The defect matrix.**
- $a_t=(\mathcal A_++\mathcal A_-)/2$ as in `eq-boundary-lax-at`, and $\partial_tL_d=L_d\,a_t(\psi)-a_t(\phi)\,L_d$ (from $\Psi_\phi=L_d\Psi_\psi$ with $(\partial_t+a_t)\Psi=0$).
- $L_d=e^{-bH\cdot(E\phi+\mathsf D\psi)/2}\,\tilde L\,e^{bH\cdot(\mathsf D^{\mathsf T}\phi-E\psi)/2}$.
- In the defining representation ($E_{\alpha_j}=e_{j,j+1}$, cyclic), the field-independent $\tilde L$ is unique up to scale for $a_1,a_2,a_3$ and both couplings [V:AN#4]:
  - **I$_-$:** $\tilde L=1+\sigma\lambda^{-1}\sum_jE_{-\alpha_j}$, so $c=1$ exactly. This is BCZ04 (3.39) with $d=\sigma$ and $\lambda_{\rm BCZ}=\lambda_{\rm book}$; BCZ (3.3) is the book's $a_t$ at $\beta=m=1$, and $\hat\phi=\beta\phi$, $\hat x=mx$ maps one to the other.
  - **I$_+$:** $\tilde L^{-1}\propto1-\sigma^{-1}\lambda\sum_jE_{\alpha_j}$, i.e. $\tilde L=\sum_{k<h}(\lambda/\sigma)^kP^k$. It is more natural written as $L_d^{-1}$, the matrix taking $\phi$-side data to the $\psi$ side.
  - A non-integrable $\mathsf D$ admits no $\tilde L$.
- **Soliton ↔ Lax parameter.** The vacuum Lax solutions give the soliton exponent at $\lambda=-e^{-\theta}$ [V:AN#4]. So the off-diagonal part $\sigma/\lambda$ of $\tilde L_{\rm I_-}$ grows like $e^{\theta-\theta_d}$, matching $z_{\rm I_-}<0$ for $\theta>\theta_d$.

**Spin conventions.**
- In the book (`sec-qg-gradation`), spin $s$ goes with $e^{s\theta}$: $E+P$ has spin $+1$, and spin-$+s$ densities are built from $\partial_-=\partial_t-\partial_x$.
- BCZ04 uses $x^\pm=(t\pm x)/\sqrt2$ and $P_{\rm BCZ}=-P_{\rm book}$, so their "spin $+2$" is book spin $-2$.

**Higher charges** [V:AN#5]. **Both** the spin $+2$ and the spin $-2$ charges of $a_2$ acquire local defect corrections, for **both** types. This contradicts BCZ04 §5 and Penati–Polvara App. B; see `prp-spin-two`.
- *Tests:* the Euler operator in $t$ (residuals about $10^{-15}$) and closed-loop integrals (about $10^{-16}$).
- *Controls:* the bulk current is conserved on shell; the density is not a total $x$-derivative; spin $\pm1$ pass; a non-integrable $\mathsf D$ fails (0.15–0.48).
- *$a_3$:* same result, from a separate run of about 6 min per case (script in the session scratchpad, not in the repo).
- *Why BCZ's obstruction vanishes:* BCZ's cubic obstruction (5.11), $A(\mathsf Dv)^3+A(\mathsf D^{\mathsf T}v)^3$, vanishes identically for all $a_r$. On Fourier modes, $\mathsf D=1+i\cot(\pi k/h)$, and $\cot a\cot b+\cot b\cot c+\cot c\cot a=1$ whenever $a+b+c\equiv0$ mod $\pi$, which is what the cubic invariant selects.
- *Structural reason:* $\sum_{j=0}^rE_{-\alpha_j}$ is invertible in the defining representation, so $L_d$ is regular at both $\lambda\to\infty$ and $\lambda\to0$. Expanding the conserved generating function $\mathrm{tr}(T_+L_dT_-)$ at both ends then gives charges of both chiralities [hand, main session].
- *Consequence:* I$_-$ and I$_+$ are distinguished by the $\lambda$-structure of $\tilde L$, by the Borel subalgebra and by the direction of the soliton-to-antisoliton change, **not** by which spin-$\pm2$ charge survives.
- *BCZ's tensor (added 2026-10-12).* The symmetric tensor of BCZ04 (5.5) is $A_{abc}=2\sum_k\mu_{k,a}\mu_{k,b}\mu_{k,c}$, i.e. $A=2\sum_k(\vec\mu_k\cdot\vec u)^3$, and their density (5.2) carries $\frac13A$; the earlier "$A=\frac23\sum$" was a factor-3 slip [V:AN#11].

**Added from chapter 38 (2026-10-12).**
- **Monodromy** (`eq-defect-monodromy`): $\mathcal T_d(\lambda)=T_+(\lambda)L_d(\lambda)^{-1}T_-(\lambda)$ with $T_+=T^\psi(x_+,0)$, $T_-=T^\phi(0,x_-)$; $\tau_d=\operatorname{tr}\mathcal T_d$ is conserved. BCZ write the inverse matrix (transport from right to left).
- **The two ends of the $\lambda$-plane.** Coefficients of $\lambda^{-n}$ at $\lambda\to\infty$ have spin $-n$, those of $\lambda^n$ at $\lambda\to0$ spin $+n$. The first defect contribution at $\lambda\to\infty$ is $\propto\sigma$, at $\lambda\to0$ $\propto\sigma^{-1}$. So the defect corrections $\mathcal B_{+s}$ carry the $\sigma^{-1}g_j$ terms and $\mathcal B_{-s}$ the $\sigma f_j$ terms, as $E+P$ and $E-P$ do (§3).
- **Brackets.** The canonical local bracket is $\{a_x\overset\otimes,a_x\}=+[r,a_{x1}+a_{x2}]\delta$; with $\partial_xT=-a_xT$ the transport matrix has $\{T\overset\otimes,T\}=-[r,T\otimes T]$ (ch. 30 `prp-sklyanin`, corrected 2026-10-12; `defect_lax_rmatrix.py` check 3). The defect bracket $\{L_d\overset\otimes,L_d\}=+[r,L_d\otimes L_d]$, coefficient exactly 1, is for the defect values treated as independent coordinates that Poisson-commute with the bulk fields, with the Dirac bracket $K^{-1}$ of (`eq-dc-defect-bracket`) [`prp-defect-involution`, Numerical for $a_1$–$a_4$]; the reduction to the equal-time field theory is open. `boundary_classical.tau_bracket` uses $+[r,TT]$, harmlessly (the reflection equations and involution are homogeneous in $r$).
- **The $a_r$ delay** (`eq-dc-ar-delay`): $z_a(\theta)=\frac{e^{-(\theta-\theta_d)}-ie^{i\pi a/h}}{e^{-(\theta-\theta_d)}-ie^{-i\pi a/h}}=\frac{\mathsf a(\theta)+\varepsilon_a\mathsf a(\theta_d)}{\mathsf a(\theta)-\varepsilon_a^{-1}\mathsf a(\theta_d)}$, $\mathsf a(\theta)=-e^{-\theta}$, $\varepsilon_a=e^{i\pi(2a-h)/2h}$ (type I$_-$, no jump). It is cz2009a (4.10), BCZ04 (6.2) after $\theta\to-\theta$, and cz2007 (2.20) after $a\to h-a$ (their $\mathsf D$ is the book's $\mathsf D^{\mathsf T}$). Identities: $z_az_{\bar a}=\frac{\cosh(\theta-\theta_d)+\sin(\pi a/h)}{\cosh(\theta-\theta_d)-\sin(\pi a/h)}$, $z_{\bar a}(\bar\theta;\bar\theta_d)=\overline{z_a(\theta;\theta_d)}$; for $a=h/2$ it is the sine-Gordon $z_{\rm I_-}$ [V:AN#8].
- **Type II** (`sec-type-two-atft`, CZ 2009b units). Sine-Gordon: $\hat p=\beta_{\rm SG}p$, $\hat q=\beta_{\rm SG}q$, $\hat\lambda=\beta_{\rm SG}\Lambda$, $\hat x=m_0x/\sqrt2$, $\mathcal D=\frac{m_0}{\sqrt2\beta_{\rm SG}^2}(f+g)$, $\hat\sigma=e^{\theta_d}$ (their $\sigma=e^\eta$). The vacuum $e^{2\hat\lambda}=1/(2\cosh^2\tau)$ has two sectors: $e^{\hat\lambda}\to+1/(\sqrt2\cosh\tau)$ gives $z=z_{\rm I_-}(\theta;\theta_d+\tau)z_{\rm I_-}(\theta;\theta_d-\tau)$; $e^{\hat\lambda}\to-1/(\sqrt2\cosh\tau)$ (CZ's $z_1$) gives its inverse. Tzitzéica: $u=\beta_{\rm lit}\phi$, $\beta_{\rm lit}^2=\beta^2/2$, $\hat x=mx/\sqrt2$, $\mathcal D=\frac{\sqrt2m}{\beta^2}(f+g)$, $e^{\theta_d}=\sqrt2\hat\sigma$; three sectors $e^{\hat\lambda}\to1,e^{\mp2\pi i/3}$ with delays $z_1$, $z_{2,3}$, $z_1z_2z_3=1$ [V:AN#10].

### 6. Quantum (`eq-weyl-module`, `eq-borel-intertwiner`, `thm-t-triangle`, `prp-two-borels`, `eq-sg-transmission-matrix`) [V:Q]

- **Conventions.** Part VI throughout: $\Delta e=e\otimes1+k\otimes e$, $\Delta f=f\otimes k^{-1}+1\otimes f$; $q=q_{\rm phys}=-e^{-i\pi\lambda}$; $\lambda=8\pi/\beta_{SG}^2-1=\omega$ ($T=\lambda$ for sine-Gordon); $y=e^{\lambda\theta}$; $v_+$ is the soliton.
- **Defect space.** $\mathcal W=\mathrm{span}\{|n\rangle,n\in\mathbb Z\}$ with $k_1=q^{2n+s}$ and $k_0=k_1^{-1}$. Level zero is forced: if $k_0k_1\ne1$ there is no intertwiner [V:Q#3].
- **Weyl modules.** With $y_d=e^{\lambda\theta_d}$ and $\mu=i/(q-q^{-1})=1/(2\sin\pi\lambda)$:
  - $U_q(\mathfrak b^+)$: $e_1|n\rangle=\mu y_d|n+1\rangle$, $e_0|n\rangle=\mu y_d|n-1\rangle$.
  - $U_q(\mathfrak b^-)$: $f_1|n\rangle=(\mu/y_d)|n-1\rangle$, $f_0|n\rangle=(\mu/y_d)|n+1\rangle$.
  - $e_0e_1$ (resp. $f_0f_1$) is central, so the q-Serre relations hold [V:Q#2].
  - Only $\mu^2$ is basis-independent. The value $\mu^2=-1/(q-q^{-1})^2$ makes $\theta_d$ real and equal to BCZ's $\eta$.
- **Transmission** $\mathbb T(\theta-\theta_d):V(y)\otimes\mathcal W\to\mathcal W\otimes V(y)$ (left to right):
  - $\mathbb T(v_+\otimes|n\rangle)=A_n|n\rangle\otimes v_++C_n|n+1\rangle\otimes v_-$,
  - $\mathbb T(v_-\otimes|n\rangle)=B_n|n-1\rangle\otimes v_++D_n|n\rangle\otimes v_-$.
- **The intertwiners.** Each is unique up to a scalar on a window [V:Q#3]:
  - $\mathfrak b^-$: $A=q^{-n}$, $D=q^{n+s}$, $C=-i(y/y_d)q^{-n}$, $B=-i(y/y_d)q^{n+s}$;
  - $\mathfrak b^+$: $A=q^{n}$, $D=q^{-n-s}$, $C=i\,q^{-1-s}(y_d/y)q^{-n}$, $B=i\,q^{-1}(y_d/y)q^{n}$.
- **Triangle relation** (braid form), on $V(y_1)\otimes V(y_2)\otimes\mathcal W\to\mathcal W\otimes V(y_2)\otimes V(y_1)$: $(1\otimes\check R)(\mathbb T(y_1)\otimes1)(1\otimes\mathbb T(y_2))=(\mathbb T(y_2)\otimes1)(1\otimes\mathbb T(y_1))(\check R\otimes1)$. It holds for both Borels with `boundary_common.rcheck` (residual $1.5\times10^{-13}$) [V:Q#4]. It is the same equation as `eq-triangle-braid` of ch. 36.
- **Identification** [V:Q#5]:
  - I$_-$ ↔ $U_q(\mathfrak b^-)$ ↔ BCZ05 (3.51)/Konik–LeClair: $|C/A|$ grows through $\theta=\theta_d$, as $z_{\rm I_-}<0$ for $\theta>\theta_d$.
  - I$_+$ ↔ $U_q(\mathfrak b^+)$.
  - In $y$, the I$_-$ matrix is "$1+(y/y_d)\cdot$off-diagonal", matching $\tilde L_{\rm I_-}=1+(\sigma/\lambda)E_-$ at $\lambda\propto e^{-\theta}$.
- **Sign of $q$.** Borel intertwiners built at $q'\in\{-q_{\rm phys},1/q_{\rm phys},-1/q_{\rm phys}\}$ fail the triangle relation with the physical $\check R$ (residual $\ge0.43$), but pass with $\check R(q')$ [V:Q#6]. So the triangle relations are self-consistent at $-q$ and do not fix the sign of $q$. Unitarity, crossing or the scalar factor must, and this is still to be done (ch. 37).
- **The parameter $s$.** BCZ's $\nu=q_{\rm phys}^s$ is the $k$-eigenvalue offset. The even sector $s=0$ gives $\nu=1$; the triangle relations allow continuous $s$.
- **Fock modules (done in ch. 40).** For a shift module with $e_1$ the shift up and $e_0|n\rangle=b_n|n-1\rangle$, q-Serre forces $u_n=e_1e_0|_n$ to satisfy $u_n-[3]u_{n+1}+[3]u_{n+2}-u_{n+3}=0$, so $u_n=A+Bq^{2n}+Cq^{-2n}$ [hand; checked in `defect_universal_r.py`]. The Weyl module is $B=C=0$, and Fock modules truncate ($u_0=0$): `eq-fock-module`, `defect_common.fock_module_b`; $\mathcal F^+$ is the kind of module of BLZ's Q-operators.

**Added from chapter 40 (2026-10-12).**
- **Scope.** $q$ is generic (not a root of unity) and $q^s$ fixed, unless stated. At $q=e^{i\pi/3}$ the shift $S^3$ commutes with the $\mathcal W^-$ generators, so $(S^3\otimes1)\mathbb T$ is a second intertwiner and `prp-sg-t-unique` fails [V:`defect_universal_r.py` check 12].
- **Universal-R direction.** $\mathcal R\in U_q(\mathfrak b^-)\hat\otimes U_q(\mathfrak b^+)$ (Part VI), and $\mathbb T_{{\rm I}_+}\propto P(\pi_\theta\otimes\rho^+)(\mathcal R)$, a series in $y_d/y$, while $\mathbb T_{{\rm I}_-}\propto P(\pi_\theta\otimes\rho^-)(\mathcal R_{21}^{-1})$, a series in $y/y_d$ (`prp-two-borels`). This is the opposite of the plan's guess.
- **Universal-R normalization.** In it the matrices are $q^{\pm s/2}$ times those of `prp-sg-t-unique`; with $q^{-s/2}\Phi$ as scalar, unitarity (real $s$) and crossing hold for $s\ne0$. Unitarity of `prp-two-borels` (4) needs real $s$, since $|q^{s/2}|\ne1$ for complex $s$ at $|q|=1$.
- **L-operators** (`eq-l-operator`), $x=y_W/y$ for $\mathbb L^+$ and $x=y/y_W$ for $\mathbb L^-$: $\mathbb L^+=\begin{pmatrix}q^{\mathsf h/2}&i\kappa x\,q^{\mathsf h/2}\mathsf a_-\\ i\kappa x\,q^{-\mathsf h/2}\mathsf a_+&q^{-\mathsf h/2}(1-\kappa^2x^2Z)\end{pmatrix}$, $\mathbb L^-=\begin{pmatrix}(1-\kappa^2x^2q^2Z)q^{-\mathsf h/2}&-i\kappa x\,\mathsf a_-q^{\mathsf h/2}\\ -i\kappa x\,\mathsf a_+q^{-\mathsf h/2}&q^{\mathsf h/2}\end{pmatrix}$, for $(\pi_\theta\otimes\rho)(\mathcal R)$ and $(\pi_\theta\otimes\rho)(\mathcal R_{21}^{-1})$ up to scalars. On a Weyl module, $\mathbb L^-$ is the quantization of $L_d$ (`eq-dq-weyl-l`).
- **Defect operator.** $\mathcal D_\omega(\theta_d)=\mathrm{Tr}_{\mathcal W}[\omega^{h_{\mathcal W}}\mathbb L_{N\mathcal W}^{-1}\cdots\mathbb L_{1\mathcal W}^{-1}]$: the product is of the **inverse** L-operators. Fock traces converge at $|q|=1$, $|\omega|<1$. Weyl traces converge for no twist; they factorize by charge sector, $\langle n|M|n\rangle=q^{mn}\langle0|M|0\rangle$ for $\mathcal W^-$ and $q^{-mn}$ for $\mathcal W^+$, so the formal support is $\omega^2=q^{-m}$ for $\mathcal W^-$ and $\omega^2=q^m$ for $\mathcal W^+$. The sector twist is $\omega_m=\pm q^{-m/2}$ for $\mathcal W^-$ (written $\omega_m$ in ch. 43 to avoid the coupling $\omega$).
- **Scalar limit.** The Gamma-function ratio of `eq-sgd-scalar` tends to $e^{-i\pi(1+\lambda)/4}$ as $\theta\to-\infty$, cancelling the explicit prefactor $e^{i\pi(1+\lambda)/4}$; $|\Phi|=(1+Y^2)^{-1/2}$, so $\Phi\to1$ (`prp-defect-twist-limit`).
- **Two labels.** $\psi\to-\psi$ maps $\mathbb T_{{\rm I}_-}$ with offset $s$ to $\mathbb T_{{\rm I}_+}$ with $s_+=-1-s$ (BCZ (3.57)–(3.59)); $\mathbb T_{{\rm I}_+}\ne(1\otimes C)\mathbb T_{{\rm I}_-}$ in the same labels.
- **Translation rules.** Boos–Göhmann–Klümper–Nirov–Razumov use $\Delta(e_i)=e_i\otimes1+q^{-h_i}\otimes e_i$, Part VI with $q\to q^{-1}$; their $\mathcal R$ is the book's $\mathcal R_{21}^{-1}$. BLZ III's $B_-$ is the book's $U_q(\mathfrak b^+)$ with $q\to q^{-1}$. Weston's $U_q(\mathfrak b^+)$ is the book's $\mathfrak b^-$ through his transposition; his type-I module $r=(r_0,0,0)$ is the Weyl module with $j=-n$, $r_0=q^s$.
- **Massless defect** (`eq-dq-massless-defect`): at $m_0=0$ the type-I$_-$ conditions are $\psi_L=\phi_L$, $\psi_R=-\phi_R+$const: left-movers pass, right-movers change sign. For a compact boson this is naturally a T-duality interface $R\leftrightarrow1/R$.
- **Helper trap.** `defect_common.sg_tmatrix` returns the Weyl module with $\mu=1$, not the module ($\mu=i/(q-q^{-1})$) that its matrix intertwines (rechecked 2026-10-12: residual $10^{-15}$ with $\mu=i/(q-q^{-1})$, 3.7 with $\mu=1$); use that $\mathcal W$ only for labels and the $k$-action. The docstring now says so.

### 7. Scattering conventions fixed in ch. 36

- **Defect ZF algebra** `eq-defect-zf`. $D_\alpha$ sits between the particles on its left and those on its right. Index rule: lower indices are incoming, upper indices outgoing, each pair in spatial order, so $\mathbb T(v_a\otimes w_\alpha)=\sum\mathbb T_{a\alpha}^{\beta b}w_\beta\otimes v_b$.
- **Right to left.** $\bar{\mathbb T}$ is parametrized by $\theta>0$ for a particle of rapidity $-\theta$ arriving from the right (as in Bajnok–Simon), so it depends on $\theta+\theta_d$.
- **Crossing.** Defect crossing is `eq-defect-crossing` (Bajnok–George (15)).
- **RTT form.** $\mathbb L=P\mathbb T$ and $R=P\check R$.
- **Doubled theory.** The reflected copy scatters with $P\check SP$.

### 8. Source dictionary (`tbl-defect-dictionary`)

| Source | Their → book | Offsets / signs |
|---|---|---|
| BCZ 2003 hep-th/0305022 (`bcz2004`) | $(\phi_1,\phi_2,m,\beta,\sigma)=(\phi,\psi,m,\beta,\sigma)$, sinh-Gordon. Linearization: $\sigma=e^p$, $p=-\theta_d$, same $\theta$ ($k=2m\sinh\theta$, $e^{-i\omega t+ikx}$) | $B=\mathcal B_d+2m(\sigma+\sigma^{-1})/\beta^2$; $V=V_{\rm book}+2m^2/\beta^2$; $P_{\rm BCZ}=\int\phi_t\phi_x=-P_{\rm book}$. The soliton section (3.5)–(3.6) uses units $m_0=\beta_{SG}=1$ and $E\propto e^{x\cosh\theta+t\sinh\theta}$, so $\theta_{03}=-\theta_{\rm book}$, with $e^{i\phi/2}=(1-iE)/(1+iE)$; (3.6) becomes $z_{\rm I_-}$ after $\theta\to-\theta$. (2.13) takes the "+" sign |
| BCZ 2004 hep-th/0401020 (`bcz2004b`) | $\beta=m=1$ ($\hat\phi=\beta\phi$, $\hat x=mx$); their $D$, $E=F$ are the book's $\mathsf D$, $E$ (I$_-$); $d=\sigma$; $\tilde K=\tilde L_{\rm I_-}$, same $\lambda$; their $\Pi(i)=i+1$ ↔ $\pi\alpha_j=\alpha_{j-1}$; their (5.16) setting is I$_+$ | $B=\mathcal B_d+(m/\beta^2)h(\sigma+\sigma^{-1})$ up to scaling; their $V=\sum n_ie^{\alpha_i\cdot\phi}$ (no $-h$, no $m^2/\beta^2$); $P_{\rm BCZ}=-P_{\rm book}$; their spin $s$ is book spin $-s$. The §5 claim on spin $\pm2$ is false (both survive) |
| BCZ 2005 hep-th/0506169 (`bcz2005`) | $(u,v,m,\beta,\sigma)=(\phi,\psi,m_0,\beta_{SG},\sigma)$; $\sigma=e^{-\eta}$, $\eta=\theta_d$; same $\theta$ ($b=-\sinh\theta$); soliton $E_{\rm book}=iE_{\rm BCZ}$ ($\xi=c+i\pi/2$), same $z$ | $B=\mathcal B_d-2m_0(\sigma+\sigma^{-1})/\beta_{SG}^2$; $U=-P_d+$const; $P_{\rm BCZ}=-P_{\rm book}$; their (2.2) must read $+\delta(x)$ (the printed $-$ is a misprint) |
| BCZ 2005 quantum, Konik–LeClair | $\gamma=\lambda$ ($=\omega$); $q_{\rm BCZ}=e^{i\pi\gamma}=-1/q_{\rm phys}$; $x=y=e^{\lambda\theta}$; $\Theta=\theta_1-\theta_2$; $S_{kl}^{mn}/a=\check R[(m,n),(k,l)]$ entry by entry; their (3.9) is the book's braid triangle relation; defect label $\alpha=2n$ | (3.51) is the $U_q(\mathfrak b^-)$ intertwiner with $y_d=e^{\lambda\eta}$, $\mu^2=-1/(q-q^{-1})^2$, $\nu=q_{\rm phys}^s$, $Q^2=-q_{\rm BCZ}=1/q_{\rm phys}$. Their $B$ is the book's $C_n$ (soliton → antisoliton, $n\to n+1$) and their $C$ the book's $B_n$. Gauge invariants: $A_{n+1}/A_n=1/q$, $D_n/A_n=\nu q^{2n}$, $C_nB_{n+1}/(A_nD_{n+1})=-(y/y_d)^2$ |

**Traps.**
- **BCZ 2003 rapidity sign.** Its soliton section has $\theta>0$ moving left and the opposite kink convention, but its linearization uses the book's sign. Convert with $\theta_{03}=-\theta_{\rm book}$ only in §3 of that paper.
- **Sign of $\sigma$.** $\sigma=e^{-\eta}$ in `bcz2005`, `cz2007` and `cz2009a`, but $\sigma=e^{+\eta}$ in `cz2009b`, and $\sigma=\sqrt2e^{-\eta}$ in `cz2011` (Tzitzéica).
- **Transposed defect matrix.** cz2007 uses the book's $\mathsf D^{\mathsf T}$: their species $a$ is the book's $h-a$ (delays, unstable windows, breather factors).
- **Misprints in the $a_r$ sources (ch. 38; Appendix D).** BCZ04 (5.43) is the transpose of their (3.28); cz2007 (2.21) has $-\coth$, not $+\coth$; cz2009b (3.7): $1/\sigma$ multiplies both terms of $f$; cz2009b (4.4) has the wrong overall sign and $1\pm\rho E_0$, and the winding statement after (4.5) is reversed; cz2009b (4.11) is their $z_3$.
- **Coupling and $q$.**
  - The $a_r$ papers use $\gamma=4\pi/\beta^2-1$, $q=-e^{-i\pi\gamma}$, $x=e^{h\gamma\theta/2}$.
  - BLZ use $q=e^{i\pi\beta^2}$.
  - Penati–Polvara use $\gamma=\beta^2/2\pi$ in Euclidean signature.


### 9. Added from chapter 37 (2026-10-11)

- **Sign of $q$, settled.** Unitarity of $\mathbb T$ holds at both $\pm q$ and cannot see the sign. Crossing with the book's $C=\sigma_x$ holds at $q_{\rm phys}$ and fails at $-q_{\rm phys}$; at $-q$ it is restored only with $C'=[[0,1],[-1,0]]$ (`defect_sg_transmission.py`). This is the same mechanism as in `sec-crossing-sign`.
- **The two types for $a_1$.** Type I$_+$ is type I$_-$ followed by $\psi\to-\psi$, so $z_{\rm I_+}=-z_{\rm I_-}$ and $\mathbb T^{\rm cl}_{\rm I_+}=-\mathbb T^{\rm cl}$.
- **Bäcklund transformation.** $\mathfrak B_a$ with $a(\theta)=-e^{-\theta}$. The type-I$_-$ conditions are $\mathfrak B_a$ at $x=0$ with $a=-\sigma=a(\theta_d)$, and the delay is the Bianchi coefficient $z=(a(\theta)+a(\theta_d))/(a(\theta)-a(\theta_d))$. The interaction coefficient $A_{11}$ is its inverse square.
- **Odd sector.** Classically $\theta_d\to\theta_d+i\pi$. In the quantum theory the bootstrap gives the $s=1$ intertwiner at $\theta_d-i\pi$ (the image prescription of the bound state).
- **Poles.**
  - The scalar factor $\Phi$ (`eq-sgd-scalar`) is the Konik–LeClair minimal solution in the BCZ (3.43)–(3.44) form, with $Y=e^{\lambda(\theta-\theta_d)}$.
  - There is a resonance at $\theta=\theta_d-i\pi/(2\lambda)$ on the soliton mass shell.
  - The pole of $\bar{\mathbb T}$ is at $-\theta_d-i\pi/(2\lambda)$.
- **Renormalization (`cnj-defect-parameter`).** For sinh-Gordon $\theta_d=(1-B/2)\theta_d^{\rm cl}$; continued to sine-Gordon, $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$.
- **Names.** Λ is the type-II auxiliary field, since λ is the coupling in ch. 37.
- **Dictionary additions.**
  - KL's ν corresponds to $-\lambda\theta_d$.
  - Bajnok–Simon's $B_{\rm BS}=B/2$.
  - cz2009b has $\sigma=e^{+\eta}$.
  - Penati–Polvara's β is the book's β; their $\gamma=\beta^2/2\pi$ is defined for the sinh-Gordon coupling, so $\gamma/2=b_{\rm SG}^2/8\pi$, which continues to $-\beta_{\rm SG}^2/8\pi$ (corrected 2026-10-12; the earlier "$+\beta_{SG}^2/8\pi$" had the wrong sign).
- **Fusion onto a Dirichlet boundary.**
  - Classically it gives Ghoshal–Zamolodchikov with $C_{0,1}=\cosh(\theta_d\pm i\beta_{SG}\phi_D/2)$.
  - In the quantum theory $\bar{\mathbb T}K_D\mathbb T=K_{\rm shG}(E,F=2i\theta_d/\pi)$, i.e. $\vartheta=2\theta_d/B$.

### 10. Added from chapter 39 (2026-10-12)

- **Notation.** $\omega=e^{2\pi i/h}$ in ch. 39 only, as in Part III; in ch. 41 and 43 $\omega$ is the coupling $4\pi/\beta^2-1$ of Part VI (both meanings are book-wide; ch. 39 says so in its conventions paragraph). $\gamma_a=\pi a/h$ is not a coupling; $\hat M=2hm/\beta^2$, $M_a=2\hat M\sin\gamma_a$; a frequency is written $\nu$.
- **Strip convention.** Complex soliton rapidities are taken in $|\operatorname{Im}\theta|<\frac\pi2$; a species-$a$ soliton at $\theta+i\pi$ is the $\bar a$ soliton at $\theta$ with its profile reversed.
- **Jump class.** $p\in\Lambda_W/\Lambda_R\cong\mathbb Z_h$, defined by $\frac\beta{2\pi}(\vec\psi(0)-\vec\phi(0))=\vec\lambda_p$ mod the root lattice. A class-$p$ jump multiplies every $f_j$ by $\omega^{-p}$ (`eq-ds-class-shift`), which is $\theta_d\to\theta_d\pm2\pi ip/h$ for type I$_\mp$, i.e. $\sigma\to\sigma\omega^{\mp p}$. For $h$ even, $p=h/2$ is $\sigma\to-\sigma$, the $a_r$ odd sector [V:`defect_solitons_an.py` check 1].
- **General delay** (`eq-an-defect-delay`): $z^{\rm I_-}_{a,p}(\theta)=z_a(\theta-\theta_d-2\pi ip/h)$ and $z^{\rm I_+}_{a,p}(\theta)=z_a(\theta_d-2\pi ip/h-\theta)$, with $z_a$ of `eq-dc-ar-delay` as a function of $\theta-\theta_d$; the type-I$_+$ delays follow by time reversal, $\mathcal B_d^{\rm I_-}(\sigma)\to\mathcal B_d^{\rm I_+}(1/\sigma)$.
- **Linear factors** (`eq-an-linear-transmission`): $R_k=0$, $\mathbb T^{\rm cl}_k(\theta)=z_k(\theta-\theta_d+\frac{i\pi}2)$ for type I$_-$ and $z_k(\theta_d-\theta-\frac{i\pi}2)$ for type I$_+$, on the mode $\pi\vec\rho_k=\omega^k\vec\rho_k$ of particle $k$, where $\mathsf D=1+i\cot\gamma_k$. So $\mathbb T^{\rm cl}=z$ continued, as $X=A$ continued in Part V. cz2007's $T_k$ is the book's $\mathbb T_{h-k}$.
- **Class energies** (`eq-ds-class-energy`): energy and momentum both real iff $2p\equiv0$ mod $h$; at $\theta_d=0$ the energy alone is real in every class (corrected 2026-10-13, Codex referee; do not write \"the energy is real iff $2p\equiv0$\"). The mass shell is algebraic: complex classes are not real vectors. Corrigan and Zambon's normalization subtracts $\hat M(\cosh\theta_d,\sinh\theta_d)$. As $\theta_d\to-\infty$ (type I$_-$) the class-$p$ energy diverges like $\frac{\hat M}2e^{-\theta_d}(1-e^{-2\pi ip/h})$; only the class-0 branch tends to the rotation line $\vec\psi=\pi\vec\phi$ without subtraction.
- **Absorption** (`eq-ds-absorption`): $\theta_{\rm abs}=\theta_d\mp i(\frac\pi2-\frac{\pi(a+2p)}h)$ for type I$_\mp$, real $\theta_d$; real iff $a+2p\equiv h/2$ mod $2h$ (then exactly one class per species); the absorbed or emitted envelope moves with velocity $\tanh\theta_d$ only if $a+2p\not\equiv0$ mod $h$; for $a+2p\equiv0$ mod $h$, $\cos\chi=0$ and $|E|$ is constant (no absorption); at $\theta_d=0$ the envelope is stationary. At complex $\theta_{\rm abs}$ the configurations are meromorphic, with zeros of $\tau_j$ on real spacetime (2026-10-13). Stored energy of class $0\to a$: real for $\theta_d\ne0$ only if $a=h/2$; at $\theta_d=0$ real for every $a$, momentum only for $a=h/2$.

- **Selection rule direction** (2026-10-13). `prp-stored-charge` is for $\theta>0$, incidence from the left: the jump changes by $\frac{2\pi}\beta\operatorname{sign}(\theta)(\vec\lambda_L-\vec\lambda_R)=\frac{2\pi}\beta(\vec\lambda_{\rm in}-\vec\lambda_{\rm out})$; for incidence from the right the outgoing phase shift is $-\arg z$. The bound is $\min(a,h-a)/\gcd(a,h)$ steps; both bounds are attained over all $\theta-\theta_d$, but at fixed $\theta_d$ with $\theta>0$ fewer values may occur. For $a=h/2$: $k=1$ iff $\theta>\theta_d$ (I$_-$).
- **Reality** (2026-10-13). $\Theta$ is a symmetry iff $\sigma$ is real, $\theta_d\in\mathbb R+i\pi\mathbb Z$ (not "iff $\theta_d$ is real"). Only classes $2p\equiv0$ are fixed by $\Theta$; the breather statement needs class 0 or $h/2$. $|z|=1$ when $\cos\frac{2\pi p}h=0$.
- **Linear factors** are a theorem for every $r$ (exact elimination of the two scalar equations of `exr-ds-linear`), since 2026-10-13.

### 11. Added from chapter 41 (2026-10-12; referee integration 2026-10-13)

- **Coupling.** $\omega=4\pi/\beta^2-1$ (CZ's $\gamma$), $q=-e^{-i\pi\omega}$, $y=e^{\omega\theta}$, $X=e^{\omega(\theta-\theta_d-i\pi/2)}$, which is CZ's $\hat x$ with $\Delta=\theta_d+\frac{i\pi}2$. Ch. 39's root of unity is written $e^{2\pi i/h}$ there.
- **The $a_r$ Weyl module** (`eq-dt-weyl`): labels $N+s$ with $N$ in the root lattice and offset $s=\vec\lambda_p$ for class $p$; $f_j|N\rangle=(\mu/y_d)\,q^{(N+s)_{j+1}}|N-\alpha_j\rangle$, $\mu=-e^{-i\pi\omega/2}/(q-q^{-1})$. Orientation $+1$ is the book's type I$_-$ (rotation $\pi$); orientation $-1$ (the powers $q^{-(N+s)_{j+1}}$) is the charge-conjugate defect, cz2007's $\mathsf D^{\mathsf T}$, with the roles of $V$ and $\bar V$ exchanged. Its two-entry $\bar V$ matrix is `eq-dt-reversed`; in the basis $q^{-|N+s|^2/2}|N\rangle$ it is cz2007 (5.10) with $\hat x=qX$, so keeping $\hat x=X$ ($\Delta=\theta_d+\frac{i\pi}2$) needs $\mu/q$ for the reversed module (2026-10-13).
- **Uniqueness** (`prp-an-transmission` (3), 2026-10-13): unique up to a scalar for every $r\ge2$, both types, among charge-conserving maps, $\theta$ off the determinant zeros. Sector determinants (projected labels, independent of sector and $s$): $q^r(1+(-q^{-1})^rX^h)$ (`eq-dt-tv`), $(q-X^h)^r$ (`eq-dt-tvbar`), $q^{-r}(1+(-q)^rX_+^h)$ (`eq-dt-tplus`), $q^r(1+(-1)^rqX^h)$ (`eq-dt-reversed`). $q$ is a root of unity iff $\omega$ is rational (also $\omega=1.37$: $q^{200}=1$); then the $k_j$-equations do not force charge conservation (shift by $5\alpha_1$ at $q^5=1$), so state charge conservation explicitly. The Chevalley involution is written $\varpi$ ($\omega$ is the coupling); it gives the type-I$_+$ $V$ matrix as the inverse of `eq-dt-reversed` on the sector $-Q$, offset $-s$, $X\to-X_+$.
- **Universal-R normalization.** Diagonal entries $q^{-{\rm wt}\cdot(N+s)}$; with it the crossing scalars do not depend on $s$.
- **Universal R for $a_r$** (2026-10-13, `prp-an-transmission` (4)). $\mathcal R=q^t\big(1+(q-q^{-1})\sum_jf_j\otimes e_j+\dots\big)\in U_q(\mathfrak b^-)\hat\otimes U_q(\mathfrak b^+)$, factor $q^{c\otimes d+d\otimes c}$ dropped at level zero. Type I$_+$ ↔ $\mathcal R$ (series in $y_d/y$), type I$_-$ ↔ $\mathcal R_{21}^{-1}=\big(1-(q-q^{-1})\sum_je_j\otimes f_j+\dots\big)q^{-t}$ (series in $y/y_d$), as for $sl_2$. With $X=-(q-q^{-1})\mu y/y_d$ and $X_+=(q-q^{-1})\mu y_d/y$ the terms of order 0 and 1 of $P(\pi\otimes\rho)(\cdot)$ equal `eq-dt-tv`, `eq-dt-tvbar`, `eq-dt-tplus` and the type-I$_+$ $V$ matrix exactly, so the evaluated R-matrix is $\phi M$ with $\phi=1+O(x^2)$ ($\phi$ not computed). The triangle relation for every $r$ follows from the YBE over $\mathbb C((u))((v))$, without hypothesis (ii).
- **Crossing.** $C=1$ ($v_i^*\mapsto\bar v_i$), dual $V(y)=\bar V(ye^{i\pi\omega})$; $c_{\rm I}=1+(-q)^rX^h$, $c_{\rm II}=1-qX^h$ (`eq-dt-crossing`). **Sign of $q$:** at $-q$ crossing fails with $C=1$; for $h$ even it holds with $C'=\mathrm{diag}((-1)^i)$, for $h$ odd with no constant conjugation. Statements about the sign are relative to the book's fixed exchange and crossing conventions.
- **Kernel of the residue, class shift** (type I$_-$; corrected 2026-10-13): on the one-dimensional kernel of the residue the transmission is the class-1 matrix (offset $s+\vec\lambda_1$, $X^h\to q^2X^h$), i.e. $\theta_d\to\theta_d+2\pi i/h+2\pi ik/(h\omega)$ for any integer $k$; every $k$ has the same classical limit. Do not write "the classical limit selects $2\pi i/h$". The excited defect is the $(h-1)$-dimensional image of the residue (ch. 36), not the kernel (for sine-Gordon they give $\theta_d\mp i\pi$). The class-$p$ scalar $g(\theta-2\pi ip/h)$ ($k=0$, integer representative with $0<a+2p\le h$) is a [Conjecture], and so are the charged-class windows and the resonances at $\theta_d-i\pi/(h\omega)$.
- **Scalar factor** `eq-dt-g`: Gamma functions paired at equal $k$. cz2009a (6.10) pairs across $k$: $P_K(z)/P_K(0)=\prod_{k<K}F_k(z)/F_k(0)\cdot\Gamma(a_K-z)/\Gamma(a_K)$ exactly, a drift like $(h\omega K)^{-z}$. cz2011 (3.12) has the same structure (drift $K^{-2z}$) and its outside Gamma poles are at $\nu-i(\xi/2+\pi/3)$, $\nu-i(\xi/2+2\pi/3)$, not at the printed $\xi$.
- **Labels of the unstable states.** Windows are open, $\max(0,a-h/2)<1/\omega<a$, for the pole in $-\pi/2<\operatorname{Im}\theta_a<0$; at $\vartheta_a=\pi/2$ the pole is real (zero width). In book labels species 1 of $a_2$ is unstable for $\omega>1$ (with the classical absorption as the limit $\omega\to\infty$; do not write "includes the classical limit") and species 2 for $\frac12<\omega<2$, cz2007's (strict) window for their species 1.
- **Penati–Polvara** (2026-10-13): the anomaly is in the conservation law of the spin-3 *current* of $a_r$ (charges of spin $\pm2$; they test $Q_{+2}$), explicit for $r\le5$, with renormalizations of the potential, the derivative coupling and the improvement terms; for sinh-Gordon the spin-4 current, i.e. the spin-3 charge. Book convention: a current of leading spin $s+1$ gives $Q_s$ (ch. 36).
- **Helper trap.** `defect_common.defect_crossing_residual` uses $C^{\mathsf T}$; this matters only for a non-symmetric $C$, such as the antidiagonal conjugation of $a_2^{(2)}$.

### 12. Added from chapter 42 (2026-10-12)

- **Doubled theory.** The reflected copy carries $V'(u)=V(-u)$, acted on through $\Delta^{\rm op}$. The doubled bulk algebra is $A=U\otimes U^{\rm cop}$ with $\Delta_A(g\otimes h)=\sum(g_{(1)}\otimes h_{(2)})\otimes(g_{(2)}\otimes h_{(1)})$ (`eq-dbd-coproduct`). The doubled K is ordered $(V\otimes\mathcal W)\oplus(V'\otimes\mathcal W)$, with blocks $P\mathbb T=\mathbb L(\theta)$ and $\bar{\mathbb T}P=\mathbb L(-\theta)^{-1}$. Symmetries of the doubled boundary are left coideals of $A$, acting by $(\sigma\otimes\chi)\Delta_A$.
- $\omega_q\otimes{\rm id}$ maps $\Delta(U)$ onto Bao–Wang's right coideal $(\omega\otimes1)\Delta(U)$; $V(1/y)\circ\omega_q=\bar V(y)$, so the reflected copy is the conjugate multiplet.
- **Invertible lines** ↔ twisted diagonals $(1\otimes\psi)\Delta(U)$: the sine-Gordon twist line has $t_1=1/t_0=e^{2i\alpha}$; charge conjugation $C$ ↔ the node swap $\tau$.
- **Jump defect.** $\Delta(U_q(\mathfrak b^\mp))$ on the Weyl module $\rho$ and $(1\otimes\tau)\Delta(U_q(\mathfrak b^\pm))$ on a second module structure $\tilde\rho$ (`eq-dbd-second-module`, sine-Gordon normalization $A_n=q^{-n}$), for type I$_\mp$; the two cannot be combined into one coideal (for I$_+$ the algebra contains $1\otimes f_0$). `cnj-twisted-diagonal` is refuted for every Hopf twist (2026-10-13): with $\mathbb L=P\mathbb T=\mathcal D+y\mathcal O$ (I$_-$) or $\mathcal D+y^{-1}\mathcal O$ (I$_+$), one power of $y$ in a single generator equation has a block free of unknowns, for arbitrary $\theta$-independent operators on $\mathcal W$; Hopf automorphisms of $U_q(\hat g')$ are torus $\rtimes$ diagram automorphisms [Kolb 2014, Thm. 3.2, after Twietmeyer]. Do not write "refuted by the test" or "refuted as stated" alone.
- **Limits $\theta_d\to\pm\infty$** (2026-10-13). Matrix limits only: type I$_-$ is diagonal (on each charge sector the twist line $t_1=1/t_0=q^{-2n-s}$) as $\theta_d\to+\infty$ and a conjugation with a charge shift as $\theta_d\to-\infty$; type I$_+$ reversed. No twisted diagonal acts on all of $\mathcal W$ in either limit ($\chi(\Delta k_1)$ and $\chi(\Delta k_1^{-1})$ are forced to be the same operator). Never write that the full (twisted) diagonal is "restored".
- **RTT and the doubled RE algebra** (2026-10-13): vanishing diagonal blocks plus boundary unitarity is a quotient; the presentations are equivalent, no embedding. $\Delta(U_q(\mathfrak b))$ is a left coideal of $U\otimes U^{\rm cop}$ isomorphic to the Borel as an algebra, not a Hopf subalgebra under $\Delta_A$.
- **Coideal search** (`prp-diagonal-coideal` (4)) is a largest coideal *subspace* of the 21-dimensional $X$, not a maximality statement for coideal subalgebras ($c^{\pm1}\otimes1$, $c=k_0k_1$, enlarges $\Delta(U)$).
- **Fusion.** $K_{\rm f}=\bar{\mathbb T}(1_{\mathcal W}\otimes K)\mathbb T$ (`eq-dbd-fused-k`), $\bar{\mathbb T}(\theta)=\mathbb T(-\theta)^{-1}$. Type I$_-$ fused onto Dirichlet($\eta_D$) is Ghoshal–Zamolodchikov with $(\eta,\vartheta)=(\eta_D,\lambda\theta_d)$, scalar factor included; after cancelling $R_0\sigma(\eta_D,u)$ the scalar comparison is $\sigma(i\lambda\theta_d,-i\theta)=\frac{1+y_d^{-2}}{1+y^{-2}y_d^{-2}}\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$ (numerical). $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$ from fusion is conditional on the quantum UV–IR matching at finite coupling (`prp-defect-fusion-k` (4)); a classical limit fixes only the leading term.
- **Corrigan–Zambon 2012 conventions** (rows = incoming soliton $a$ and charge $\alpha$, $a+\alpha=b+\beta$; (2.9) is the row-incoming product $T(\theta)R(\theta)T(-\theta)^{-1}$; $Q^\alpha=q^{-\alpha/2}$, $Q^2=1/q$). Their (2.17) is right; (2.19) needs $h=-(|b_-|^2+|b_+|^2)Q^{-2}$ and the opposite sign of $s_-$ (Appendix D).
- **Appel–Vlaar 2024.** Thm. 5.6.1 acts on $M\otimes V$ for any $M$ in $\mathcal W_\theta$ (modules of $U_q(\mathfrak k)$, need not extend to $U_q(\hat g)$); the restriction hypothesis is Thm. 6.5.2, and the trigonometric tensor K-matrix is identified only up to a scalar series.
- **Classical dressing** (`eq-dbd-dressed-k`): $K_\phi(\lambda)=\sigma(L_d(1/\lambda)^{-1})K_\psi(\lambda)L_d(\lambda)^{-1}$, with the dynamical compatibility $\partial_tK=Ka_t(\lambda)+\sigma(a_t(1/\lambda))K$. `prp-fused-bdg` proves conservation, not involution: call the fused $a_n$ systems "dynamical boundaries with a conserved Lax generating function", not "integrable" (2026-10-13).
- **BDG dictionary.** The sine-Gordon type-I$_-$ defect potential is BDG's $H_{\rm osc}$ at $(\hat p,\hat q)=(\beta_{\rm SG}\psi/2-i\theta_d,\beta_{\rm SG}\psi/2+i\theta_d)$.
- **Trap.** `boundary_coideal`, `qsp_checks`, `boundary_classical`, `boundary_sine_gordon` and `defect_axioms` run their checks at import and are not import-safe.

### 13. Added from chapter 43 (2026-10-12)

- **Notation.** $\omega=4\pi/\beta^2-1$ is always the coupling; ch. 39's root of unity is written $e^{2\pi i/h}$, and the ch. 40 twist $\omega_m=q^{-m/2}$.
- **Breathers.** The lowest breather $B^{(a)}_1$ is the $U_q$ singlet with species $h-a$ at $\theta-i\bar u$ and species $a$ at $\theta+i\bar u$, $\bar u=\frac\pi2-\frac\pi{h\omega}$. It exists only in this order (checked for $a=1,r$). $\bar V\otimes V$ gives particle 1 and $V\otimes\bar V$ gives particle $r$. CZ's labels are swapped (transposed $\mathsf D$).
- **Particle factors.** $\mathbb T_a=\frac{e^{-D}+e^{i\pi a/h}}{e^{-D}+e^{-i\pi a/h}}=-e^{i\pi a/h}/(h-a)$ and $\mathbb T^+_a=e^{i\pi a/h}/(a)$, in the blocks `eq-blocks`, with $D=\theta-\theta_d$. $\bar{\mathbb T}_a(\theta)=\mathbb T_{\bar a}(i\pi-\theta)$. Hermitian analyticity holds for real $\theta_d$.
- **Transmission bootstrap (no S-factor).** For $a+b<h$: $c=a+b$, $\bar u=(\pi b/h,\pi a/h)$. For $a+b>h$: $c=a+b-h$, $\bar u=(\pi(h-b)/h,\pi(h-a)/h)$.
- **Doubled-boundary equations (from $\Delta_A$).** $e$-type: $\mathbb T[\pi(e_j)\otimes1+\pi(k_j)\otimes Y]=[Y\otimes1+\kappa\otimes\pi(e_{\psi(j)})]\mathbb T$. $f$-type: $\mathbb T[\pi(f_j)\otimes\kappa+1\otimes Y]=[Y\otimes\pi(k_{\psi(j)})^{-1}+1\otimes\pi(f_{\psi(j)})]\mathbb T$. For $a_r$ type I$_-$ the twist is $\psi=\pi$, $e_j\mapsto e_{j-1}$ (`prp-dpt-doubled`).
- **Real-coupling vacuum.** $\mathcal B_d\ge0$ for $\sigma>0$, with equality only at $\phi=\psi=0$; unbounded below for $\sigma<0$ (`prp-dpt-vacuum`).

## Status-label rules specific to Part VIII

- Every quantum $a_r$ transmission matrix is labelled as a solution of the scattering axioms ([Numerical] or [Established]); its identification with the transmission matrix of the Lagrangian defect is [Conjecture], because of the higher-spin anomalies found by Penati and Polvara (arXiv:1902.10690).
- `cnj-twisted-diagonal` is pre-registered: the expectation, stated in the text before the check, is that a jump defect preserves a diagonal *Borel* subalgebra with a dynamical boundary module $\mathcal W$, not the full twisted diagonal coideal. The script `defect_doubling.py` decides; the text reports the outcome either way. *Outcome (ch. 42):* refuted for every Hopf twist by an exact argument (2026-10-13; before that, by the least-squares test only); the diagonal Borel holds, and the opposite Borel survives twisted by charge conjugation ($\pi$ for $a_r$, ch. 43).
- The sign of $q$ is tested by unitarity, crossing and the scalar factor, not by the triangle relations alone (verified: Borel intertwiners at $-q_{\rm phys}$ satisfy the triangle relation with $\check R(-q_{\rm phys})$).
- `prp-spin-two` corrects BCZ 2004 §5 and Penati–Polvara App. B. Since the referee pass of 2026-10-12 it is labelled [Theorem] for (1) and for $a_2$ in (3) (exact, $\sigma$ symbolic), [Numerical] for $a_3$, and [Proposition, sketch] for (2) and all $a_r$; no priority is claimed.
- Unless stated otherwise, quantum statements assume $q$ generic (not a root of unity) and a fixed $q^s$.

## Id registry

All ids are already present as placeholders in the skeleton chapters, so chapters drafted in parallel can cross-reference each other. Add new ids only inside your own chapter, prefixed so that they cannot clash, and report them.

| Ch | File | Ids |
|---|---|---|
| 36 | `part8-defects/36-defect-scattering.qmd` | `sec-defect-scattering`, `sec-defect-kinematics`, `eq-defect-zf`, `prp-defect-elastic`, `sec-rt-algebra`, `prp-rt-nogo`, `sec-triangle-relations`, `eq-triangle-st`, `eq-triangle-braid`, `eq-triangle-rtt`, `prp-triangle-forms`, `sec-defect-unitarity-crossing`, `eq-defect-unitarity`, `eq-defect-crossing`, `sec-doubling-trick`, `eq-doubled-k`, `prp-doubling`, `sec-defect-bound-states`, `eq-defect-bootstrap`, `sec-defect-circle`, `prp-twist-defect`, `sec-topological-lines`, `sec-defect-finite-size`, `sec-defect-lee-yang`, `eq-ly-transmission`, `prp-ly-transmission`, `sec-defect-scattering-exercises`, `sec-defect-scattering-notes` |
| 37 | `part8-defects/37-defect-sine-gordon.qmd` | `sec-defect-sine-gordon`, `sec-backlund`, `eq-sg-backlund`, `prp-sg-backlund`, `prp-bianchi`, `sec-sg-defect-lagrangian`, `eq-sg-defect-lagrangian`, `eq-sg-defect-potential`, `eq-sg-defect-conditions`, `prp-sg-defect-classification`, `prp-sg-defect-momentum`, `sec-sg-defect-solitons`, `eq-sg-delay`, `prp-sg-delay`, `prp-half-soliton`, `sec-sg-defect-linear`, `eq-shg-linear-transmission`, `sec-type-two`, `sec-sg-transmission-matrix`, `eq-sg-transmission-matrix`, `prp-sg-transmission`, `cnj-defect-parameter`, `sec-sg-defect-spectrum`, `sec-sg-breather-transmission`, `prp-sg-breather-transmission`, `sec-defect-sine-gordon-exercises`, `sec-defect-sine-gordon-notes` |
| 38 | `part8-defects/38-defect-conditions.qmd` | `sec-defect-conditions`, `sec-defect-lagrangian`, `eq-defect-action`, `eq-defect-conditions`, `eq-defect-potential`, `eq-d-matrix`, `sec-defect-momentum`, `prp-defect-momentum`, `prp-ar-only`, `cnj-momentum-integrability`, `sec-defect-lax`, `eq-defect-lax-matrix`, `prp-defect-lax`, `prp-spin-two`, `sec-defect-r-matrix`, `eq-defect-monodromy`, `prp-defect-involution`, `sec-defect-backlund-an`, `sec-type-two-atft`, `prp-type-ii-momentum`, `sec-defect-folding`, `sec-defect-theta`, `prp-defect-theta`, `sec-defect-conditions-exercises`, `sec-defect-conditions-notes` |
| 39 | `part8-defects/39-defect-solitons.qmd` | `sec-defect-solitons`, `sec-an-soliton-defect`, `eq-an-defect-delay`, `prp-an-defect-delay`, `sec-defect-topological-charge`, `prp-stored-charge`, `sec-defect-absorption`, `sec-defect-multi`, `sec-defect-reality`, `sec-classical-transmission-factors`, `eq-an-linear-transmission`, `prp-an-linear-transmission`, `prp-defect-frozen-soliton`, `sec-defect-solitons-exercises`, `sec-defect-solitons-notes` |
| 40 | `part8-defects/40-defect-quantum-groups.qmd` | `sec-defect-qg`, `sec-cft-topological-defects`, `sec-perturbed-defects`, `sec-borel`, `eq-q-oscillator`, `eq-weyl-module`, `eq-fock-module`, `sec-t-intertwiner`, `eq-borel-intertwiner`, `thm-t-triangle`, `prp-sg-t-unique`, `prp-two-borels`, `sec-universal-r-affine`, `eq-l-operator`, `prp-l-operator`, `sec-q-operators`, `prp-defect-transfer`, `cnj-defect-tq`, `prp-defect-twist-limit`, `sec-defect-boundary-dictionary`, `sec-defect-qg-exercises`, `sec-defect-qg-notes` |
| 41 | `part8-defects/41-soliton-transmission.qmd` | `sec-soliton-transmission`, `sec-t-procedure`, `sec-an-soliton-transmission`, `prp-an-transmission`, `sec-transmission-unitarity-crossing`, `prp-transmission-crossing`, `sec-a22-transmission`, `prp-a22-transmission`, `sec-an-defect-bound-states`, `prp-a2-unstable`, `sec-transmission-semiclassics`, `prp-transmission-semiclassical`, `sec-soliton-transmission-exercises`, `sec-soliton-transmission-notes` |
| 42 | `part8-defects/42-defects-and-boundaries.qmd` | `sec-defect-boundary`, `sec-doubled-theory`, `prp-reflected-copy`, `sec-diagonal-qsp`, `prp-diagonal-coideal`, `sec-jump-coideal`, `cnj-twisted-diagonal`, `prp-diagonal-borel`, `sec-defect-boundary-fusion`, `prp-defect-fusion-k`, `prp-fused-bdg`, `sec-tensor-k`, `sec-noninvertible-boundaries`, `sec-defect-boundary-exercises`, `sec-defect-boundary-notes` |
| 43 | `part8-defects/43-particle-transmission.qmd` | `sec-particle-transmission`, `sec-breather-transmission`, `prp-lightest-breather-universal`, `sec-an-particle-transmission`, `prp-an-particle-transmission`, `sec-transmission-continuation`, `cnj-transmission-continuation`, `sec-transmission-limits`, `sec-defect-spectrum`, `sec-defect-finite-volume`, `sec-defect-open-problems`, `sec-particle-transmission-exercises`, `sec-particle-transmission-notes` |

Tables added outside the skeleton: `tbl-defect-dictionary` (ch. 38, `sec-defect-lagrangian`), `tbl-boundary-defect-algebra` (ch. 40, `sec-defect-boundary-dictionary`), `tbl-dbd-symmetries` (ch. 42) and `tbl-dpt-spectrum` (ch. 43). Equation and exercise ids added during drafting carry a chapter prefix: `sgd-` (ch. 37), `dc-` (38), `ds-` (39), `dq-` (40), `dt-` (41), `dbd-` (42), `dpt-` (43).
