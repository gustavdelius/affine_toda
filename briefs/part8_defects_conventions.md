# Part VIII (defects): conventions, names and id registry

This file is the single source of truth for the conventions of Part VIII, "Affine Toda theory with defects" (`part8-defects/`, chapters 36–43). It was written on 2026-10-11, before drafting. Every agent drafting or refereeing a chapter of Part VIII reads it first. Where a convention is marked *verified*, the named check in `textbook_code/` derives it; where it is marked *to verify*, the first agent that needs it confirms it by script and reports.

The plan of the part is recorded in the memory note `defects-part` and in the plan file `~/.claude/plans/system-reminder-you-are-operating-indexed-stroustrup.md`.

## Placement and numbering

- Part VIII sits between Part VII (the half-line, ch. 28–35) and Foundations. When the part is complete, Foundations becomes Part IX (`part9-foundations/`, ch. 44–49) and the outlook becomes chapter 50; until then the old numbers remain.
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
- **Fock modules (open).** For a shift module with $e_1$ the shift up and $e_0|n\rangle=b_n|n-1\rangle$, q-Serre forces $u_n=e_1e_0|_n$ to satisfy $u_n-[3]u_{n+1}+[3]u_{n+2}-u_{n+3}=0$, so $u_n=A+Bq^{2n}+Cq^{-2n}$ [hand]. The Weyl module is $B=C=0$, and Fock modules truncate ($u_0=0$). `fock_module_b` and its intertwiners are left to the ch. 40 agent.

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
- **Coupling and $q$.**
  - The $a_r$ papers use $\gamma=4\pi/\beta^2-1$, $q=-e^{-i\pi\gamma}$, $x=e^{h\gamma\theta/2}$.
  - BLZ use $q=e^{i\pi\beta^2}$.
  - Penati–Polvara use $\gamma=\beta^2/2\pi$ in Euclidean signature.

## Status-label rules specific to Part VIII

- Every quantum $a_r$ transmission matrix is labelled as a solution of the scattering axioms ([Numerical] or [Established]); its identification with the transmission matrix of the Lagrangian defect is [Conjecture], because of the higher-spin anomalies found by Penati and Polvara (arXiv:1902.10690).
- `cnj-twisted-diagonal` is pre-registered: the expectation, stated in the text before the check, is that a jump defect preserves a diagonal *Borel* subalgebra with a dynamical boundary module $\mathcal W$, not the full twisted diagonal coideal. The script `defect_doubling.py` decides; the text reports the outcome either way.
- The sign of $q$ is tested by unitarity, crossing and the scalar factor, not by the triangle relations alone (verified: Borel intertwiners at $-q_{\rm phys}$ satisfy the triangle relation with $\check R(-q_{\rm phys})$).
- `prp-spin-two` corrects BCZ 2004 §5 and Penati–Polvara App. B. It is labelled [Numerical] for $a_2,a_3$, with the cubic cancellation and the Lax argument by hand, and it is the first item of the ch. 38 referee brief.

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

Tables, not yet in the skeleton: `tbl-defect-dictionary` (ch. 38, `sec-defect-lagrangian`) and `tbl-boundary-defect-algebra` (ch. 40, `sec-defect-boundary-dictionary`).
