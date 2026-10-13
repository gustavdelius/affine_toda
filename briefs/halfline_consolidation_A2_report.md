<!-- Consolidation agent report A2, 2026-10-13. Model: Claude Opus 5.5, effort xhigh. Prompt: briefs/halfline_consolidation_common.md + briefs/halfline_consolidation_A2_prompt.md. Interrupted three times by usage limits and resumed. Saved verbatim by the coordinator. Scripts and logs copied to briefs/halfline_consolidation_A2_scripts/. -->

# A2 report: problem 11, scalar factors, CDD factor and breather amplitudes of the soliton-preserving (SP) matrices

All five tasks are done.
- **What holds:** the worker's crossing function, scalar factors, $n=1$ Dirichlet match and breather θ-dependence hold, to $10^{-15}$ at matrix level and $10^{-26}$ for the scalar factors.
- **What I corrected:**
  - the constant $k$ (closed form found; the worker's numerics land on the wrong branch for some $l,\omega$);
  - the sign of the $n=1$ Dirichlet relation;
  - the overall constant of the breather amplitude, which is an $l$-dependent root of unity and a normalization convention for $n\ge2$.
- **What is new:** the SP soliton bootstrap $V^{\otimes n}\to\bar V$ holds with the minimal solution. It fixes the exponential ambiguity and the overall sign, but not the CDD factor. It rules out $\sigma_A=(-\frac T2)_I^2$ for even $n$; $\sigma_B$ is allowed.
- **Continued amplitudes:** they satisfy all the particle axioms, with crossing via $S_{aa}$.

Everything below was computed in my scratch directory. No repository file was touched, and no git command was run.

## 1. Verdicts

### Task 1: crossing function and the soliton–antisoliton factor — VERIFIED
Code: book `boundary_common.rcheck`, `an_vector`, `q_phys`; $K_V$, $K_{\bar V}$ copied from `soliton_preserving.py`, identical to it at ω = 0.37 (difference 0).

**Bulk normalization.** I checked crossing of the bulk S-matrix at matrix level, with the identity pairing $v_i\leftrightarrow\bar v_i$, in both index forms:
- $\check S_{ab}^{dc}(\theta)=\check S_{\bar d a}^{c\bar b}(i\pi-\theta)$;
- $\check S_{ab}^{dc}(\theta)=\check S_{b\bar c}^{\bar a d}(i\pi-\theta)$.

Both hold for $n=1,2,3$, ω = 0.37 and 1.29, two θ, to $3\cdot10^{-15}$, in the form $\check R_{VV}(\mu)=t(\mu)\,\check R_{\bar VV}(T-\mu)$, with $t=\sin\pi\mu/\sin\pi(\omega-\mu)$. All transmission entries ($i\ne j$) of $\check R_{\bar VV}$ and $\check R_{V\bar V}$ equal 1 in rcheck normalization.

Hence the formula I use: $\check S_{\bar VV}(\mu)=S^T(T-\mu)\,\check R_{\bar VV}(\mu)$, and the same for $V\bar V$, with $S^T=t\,F_{11}$ as in `prp-an-diagonal-scalar`. In the boundary crossing this gives $F(2\mu)=S^T(T-2\mu)$, as the worker said.

**SP crossing functions.** I used (eq-boundary-crossing-unitarity) with
- $k_A(\theta)=\sum_i K_{\bar V}(\frac{i\pi}2-\theta)_{ii}\,v_i\otimes\bar v_i$,
- $k_B(\theta)=\sum_i K_V(\frac{i\pi}2-\theta)_{ii}\,\bar v_i\otimes v_i$.

Then $k_B(\theta)=\Phi_1\check R_{V\bar V}(2\theta)k_A(-\theta)$ and $k_A(\theta)=\Phi_2\check R_{\bar VV}(2\theta)k_B(-\theta)$, with $\Phi_1=-(-q)^{1-l}s$, $\Phi_2=-(-q)^{l-1}s$, $s=\sin\pi(2\mu-T)/\sin\pi(2\mu+T)$. This holds for $n=1,2,3$, all $l$, ω = 0.37 and 1.29, two η, two θ, to $3.3\cdot10^{-15}$, and is independent of η.

### Task 2: scalar factors — VERIFIED, with a closed form for $k$ and a sign correction
Code: my own Stirling–Hurwitz implementation at 30 digits (`sp_mp.py`). Its $F_{11}$ agrees with the book's `an_nondiagonal_scalar.scalar_data`.

**Equations, derived independently.** The physical amplitudes are $a\,K_V$ and $\tilde b\,K_{\bar V}$ (the worker's $b=\tilde b/\lambda$, with $\lambda=(-q)^{1-l}$).
- $a(\mu)a(-\mu)U=1$, with $U=\sinh(\eta+i\pi\mu)\sinh(\eta-i\pi\mu)$.
- $\tilde b(\mu)\tilde b(-\mu)\bar U=1$, with $\bar U$ the same function of $\bar\eta$.
- $\Phi_1(\mu)\,a(\frac T2-\mu)=F(2\mu)\,\tilde b(\frac T2+\mu)$.
- $\Phi_2(\mu)\,\tilde b(\frac T2-\mu)=F(2\mu)\,a(\frac T2+\mu)$.

These are the worker's equations.

**Structure, proved by hand.**
1. Write $a=\rho\alpha$ and $\tilde b=\rho\tilde\beta$, where $\rho\rho(-)=1$ and $\rho(\frac T2-\mu)/\rho(\frac T2+\mu)=-F/s$.
2. The crossing equations then reduce to $\tilde\beta(x)=\lambda\,\alpha(T-x)$.
3. With unitarity this gives $\alpha(x+2T)/\alpha(x)=U(x)/(\lambda^2\bar U(x+T))$.

The solution is:
- $\alpha=\lambda^{-\mu/T}\alpha_0$, with $\alpha_0=i\Lambda_{2T}(\mu;\varepsilon+2T)/[\Lambda_{2T}(\mu;\bar\varepsilon+T)\sin\pi(\varepsilon+\mu)]$;
- $\tilde\beta=\lambda^{\mu/T}\beta_0$, with $\beta_0(\mu)=\alpha_0(T-\mu)=i\Lambda(\mu;\bar\varepsilon+2T)/[\Lambda(\mu;\varepsilon+T)\sin\pi(\bar\varepsilon+\mu)]$;
- $\rho=R(\mu)(\frac\omega2)_I$, the worker's $R$.

**The constant $k$, in closed form.** $\lambda^{-\mu/T}=e^{(1-l)\omega\theta}$. In the worker's parametrization $\gamma=\lambda^{-2\mu/T}e^{k\mu}\Lambda/\Lambda$, this means $k=\frac1T\log\lambda=-2\pi i(1-l)/h$.
- The worker's principal-branch numerics equal this only when $|(1-l)\omega/2|<\frac12$.
- At $(n,l)=(3,4)$, ω = 0.37, their $k$ is off by $-2\pi i/T$. That is an extra $e^{-2\theta}$ on $a$ and $e^{2\theta}$ on $\tilde b$. (I ran their `scal.gamma_shift` only to identify the constant; it is not used as a check.)

**The solution is ambiguous, but the bootstrap fixes it.**
- The equations are exactly invariant under $a\to a\,e^{m\theta}$, $\tilde b\to(-1)^m\tilde b\,e^{-m\theta}$, $m\in\mathbb Z$ (proved).
- The soliton bootstrap (task 4) selects $m=0$. The $m=\pm1,\pm2$ members fail, at ω = 0.37 and 1.29.
- With $m=0$, the exponential centres the exponents of $K_V$ in each block: $e^{(1-l)\omega\theta}K_V$ has exponents $-(l-1),\dots,l-1$ on the block and $-(n-l),\dots,n-l$ off it.

**Numerical checks at 30 digits.** All four equations hold to $\le1.9\cdot10^{-26}$ for 12 pairs $(n,l)$ including the requested six, ω = 0.37 and 0.83, two η, two μ (`t2_scalar.log`). The book-style $G_j$ form below equals the Λ form to $6\cdot10^{-28}$.

**$n=1$ — VERIFIED up to a sign (CORRECTED).**
- With the book's $R_0$ and σ, the Dirichlet limit of (eq-bsg-k-full) is $R_0\sigma(\xi,u)/\cos\xi\cdot\mathrm{diag}(\cos(\xi\pm\lambda u))$; I checked that $\sigma(i\vartheta,u)\to1$.
- This equals $-\,a_{\min}\,(-\frac T2)_I^2\,K_V$ at $\xi_B=\pi\varepsilon-\frac\pi2$, for every η: $l=1$, ω = 0.37 and 0.83, three η, to $10^{-30}$.
- The worker's "$a/[\mathrm{GZ}]=-i(T/2)_I^2$" therefore has the wrong sign in book normalization; it is $+i(T/2)_I^2$. With the sign fixed as in (eq-sp-scalar) below, GZ $=a\,(-\frac T2)_I^2K_V$ exactly.

### Task 3: breather amplitude — VERIFIED (θ-dependence); CORRECTED (overall constant)
Code: book `rcheck`, `an_vector` and coproduct singlets for the matrix projection (as in `an_nondiagonal_scalar` check 2). For the 30-digit version I wrote an mpmath projection. It uses a closed form of the zero-weight block of $\check R_{V\bar V}$, obtained from crossing and checked against rcheck to $1.5\cdot10^{-14}$ for $n\le4$; it agrees with the book-code projection to $10^{-15}$.

- **Projection.** $B_1\to\bar B_1$ is confirmed: residual $\le2\cdot10^{-15}$.
- **Amplitude, 30 digits.**
  - With singlets normalized to 1 on (highest weight of the left factor)⊗(lowest weight of the right factor): $K_{1\to n}=e^{-2\pi il/h}\,K_0\,E_\eta$ and $K_{n\to1}=e^{2\pi il/h}\,K_0\,\bar E_\eta$.
  - This holds for $n=1,2,3,4$, all $l$, two η, three μ: the ratio is constant to $3\cdot10^{-29}$.
  - Unitarity $K_{1\to n}(\theta)K_{n\to1}(-\theta)=1$ holds to $10^{-27}$.
  - So the worker's six-block **fit $K_0$ is exact**, and their $E_\eta$ is exact.
- **The constant is a convention for $n\ge2$.** It depends on the relative normalization of the $B_1$ and $\bar B_1$ singlets (unitarity, crossing and the bootstrap are blind to it). In the worker's index-0 normalization it is a different root of unity at each $l$; for example $+K_0E$ at $l=1$. So "$-K_0E_\eta$ for all $n$" is not a normalization-free statement.
- **$n=1$.** The two singlets coincide under $\bar V\cong V$, so the constant is canonical. There $K_{B_1}=-K_0E_\eta$, and with $\sigma_A$ at both constituents it equals Ghoshal's $K_{B_1}$ (eq-bsg-b1-reflection with $\vartheta\to\infty$, $\xi_B=\pi\varepsilon-\frac\pi2$) to $10^{-30}$, for three η and two ω.
- **η drops out for $l=0,h$.** $E_\eta=1$ there. Analytically: $K_V$ then depends on η only through the scalar $\sinh(\eta\mp T\theta)$, which cancels the factor $1/\sin\pi(\varepsilon\mp\mu)$ in $a$.

### Task 4: CDD factor — neither candidate is forced; $\sigma_A$ is refuted for even $n$; the minimal solution and $\sigma_B$ are allowed
Code: book `intertwiners`, `rcheck`, `an_vector`, `cop`; my 30-digit scalars.

- **Embedding.** $\bar V(\theta)$ embeds in $V(\theta_1)\otimes\dots\otimes V(\theta_n)$ only for $\mu_k=\mu-T+k\omega$ increasing from left to right. The reversed order has no embedding for $n\ge2$.
- **Out-state.** The out-embedding has the same rapidity pattern; it equals the in-embedding to $4\cdot10^{-14}$.
- **Matrix level.** Reflecting the constituents (rcheck $VV$, $K_V$) maps $\iota$ to $P\,\iota K_{\bar V}(\theta)$, with residual $\le4.4\cdot10^{-14}$. This is an independent check of $K_{\bar V}$ and $\bar\eta$.
- **Scalar.** $\tilde b_{\min}=(-1)^{n+1}P\prod_k a_{\min}(\mu_k)\prod_{j<k}F_{11}(\mu_j+\mu_k)$ for:
  - $n=1,2,3$, all $l$, two η, three μ (ω = 0.37);
  - $n=2,3$ at ω = 1.29;
  - $n=4$ at $l=0,2,5$.

  The agreement is $\le10^{-14}$. With the sign $-i$ of (eq-sp-scalar) the bootstrap holds exactly for every $n$. The $n=1$ control is $Q=1$.
- **CDD freedom, proved.**
  - Allowed SP CDD pairs: $a\to\sigma(\mu)a$, $\tilde b\to\sigma(T-\mu)\tilde b$, with σ of period $2T$ and $\sigma\sigma(-)=1$. Single blocks are allowed, not only crossing-symmetric products.
  - A block $(x)_I$ changes the bootstrap by $\sin\frac{\pi(\mu-T-x)}\omega/\sin\frac{\pi(\mu-T+x)}\omega$. This is $+1$ for $x\in\omega\mathbb Z$ and $-1$ for $x\in\omega(\mathbb Z+\frac12)$.
  - For even $n$ a ratio $-1$ is absorbed by flipping the overall sign.
- **Numerically** (`t4_cdd`, and checks 6 and 8 of the script):
  - $\sigma_A$ passes for odd $n$ and **fails for $n=2,4$**; the ratio is θ-dependent.
  - $\sigma_B$ passes for $n=1..5$ (for even $n$ with the opposite sign).
  - The pair $((\omega)_I,-(T-\omega)_I)$ passes.
  - The pair $((-\frac\omega2)_I^2,(T+\frac\omega2)_I^2)$, which equals $\sigma_A$ at $n=1$, passes for all $n$.
- **Conclusion.** The bootstrap does not decide the CDD factor. For $n=1$ the GZ comparison requires $\sigma_A$; for even $n\ge2$ the bootstrap excludes $\sigma_A$.

### Task 5: continued amplitudes — done
Code: `particle_reflection.py` blocks and S-matrix (exec of its definitions); my continued $K_0$.

**Axioms for reflection into the conjugate, derived from ch. 28.** Write $L_a=K_{a\to h-a}$.
- Unitarity: $L_a(\theta)L_{h-a}(-\theta)=1$.
- Crossing: $L_a(\frac{i\pi}2-\theta)=S_{aa}(2\theta)L_a(\frac{i\pi}2+\theta)$. The pair in the boundary state is $aa$, so $S_{aa}$ replaces the book's $S_{a\bar a}$.
- Bootstrap: $L_c=L_a(\theta+iu_a)\,S_{a,h-b}(2\theta+iu_a-iu_b)\,L_b(\theta-iu_b)$.

**Results.**
- $K_0\to(\frac h2)(\frac h2+\frac B2)(\frac h2-1)(\frac h2-1-\frac B2)(\frac h2+1+\frac B2)/(\frac h2-1+\frac B2)$ for $n=1..7$, to $10^{-30}$.
- The worker's continued $\sigma_A$ and $\sigma_B$ block forms are confirmed to $10^{-30}$.
- For the minimal solution, $\sigma_B$ and both allowed pairs ($n=1..5$, $l\in\{0,1,n\}$, two $B$): unitarity, crossing and the bootstrap at every fusing, both $a+b<h$ and $a+b>h$, hold to $\le1.6\cdot10^{-29}$. The bootstrap-iterated $L_n$ equals the soliton-side $K_{n\to1}$.
- $\sigma_A$ breaks the particle bootstrap for $n=2,4$; this matches task 4.
- **$B\to0$:** $L_1\to(\text{phase})\,g_{h/2}\,g_{h/2-1}\,E_\eta$, with $g_x=(\sinh\theta+i\sin\frac{\pi x}h)/(\sinh\theta-i\sin\frac{\pi x}h)$, for all tested CDD factors ($10^{-12}$ at $B=10^{-12}$, $n=1..4$). This is θ-dependent for every $l$ and η when $h\ge3$; zero–pole counting shows no η can cancel it.
- So the worker's claim 10 is **VERIFIED** for these CDD factors. The CDD-independent argument remains the structural one, $a\to h-a$.

### Other worker claims — UNCHECKED (outside A2)
- Coideal identification and the Regelskis–Vlaar 2018 match.
- Claim 4: coweights outside $2\Lambda_W^\vee$ with $K_+=e^{i\pi\mu\cdot H}$. This bears on Appendix D's Delius 1998b entry (1.4).
- Claim 5: real-coupling Dirichlet is not integrable for $n\ge2$.
- The auxiliary identities $F_{11}(x-2T)=\dots$ and $c=U(\mu_a)c_0$; the second is superseded by the full amplitude check.

## 2. Check script for `textbook_code/`
File: `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a2/sp_scalar.py` (copied to `briefs/halfline_consolidation_A2_scripts/`)

It imports `boundary_common` and execs the definitions part of `particle_reflection.py`, located through `boundary_common.__file__`. Run it with `PYTHONPATH=…/textbook_code`. Its docstring names the book ids and conventions. Runtime 4.5 minutes. Output (`sp_scalar.log`):

```
PASS bulk crossing with v_i <-> vbar_i: Rc_VV(mu) = t(mu) Rc_(VbarV, VVbar)(T - mu) [4e-15], transmission entries = 1 [3e-15]; so S_(VbarV) = S^T(T - mu) Rc (n = 2, 3; omega = 0.37, 1.29)
PASS SP crossing-unitarity: kB = Phi1 Rc_VVbar(2th) kA(-th), kA = Phi2 Rc_VbarV(2th) kB(-th), Phi1,2 = -(-q)^(+-(1-l)) sin pi(2mu-T)/sin pi(2mu+T), independent of eta (n = 2, 3, all l, omega = 0.37, 1.29) [3e-15]
PASS scalar factors a, bt solve unitarity and crossing-unitarity, (n, l) = (1,1), (2,0), (2,1), (3,0), (3,1), (3,2) [3.6e-27]; G_j form = Lam form [5.6e-28]
PASS n = 1: Ghoshal-Zamolodchikov Dirichlet (book R_0, sigma; xi = pi eps - pi/2) = -a (-T/2)_I^2 K_V for all eta [6.5e-28]
PASS soliton bootstrap V^(x n) -> Vbar (n = 1, 2, 3): M iota = P iota K_Vbar [7e-15], Vbar not in the reversed order [0, 0]; bt = (-1)^(n+1) P prod a prod F11 [5.7e-16]; the e^(+-2 theta) solutions fail [1.0]
PASS CDD and bootstrap: sigma_B and the pair ((omega)_I, -(T-omega)_I) satisfy sigma_b/prod sigma_a = +-1 for n = 1..5, sigma_A = (-T/2)_I^2 for odd n; sigma_A fails for n = 2, 4 [ratio varies by 0.66]
PASS breather: zero-weight block of Rc_VVbar in closed form [2e-15]; B_1 -> Bbar_1 = e^(-2 pi i l/h) K0 E, Bbar_1 -> B_1 = e^(2 pi i l/h) K0 Ebar ((n, l) = (1,1), (2,0), (2,1), (3,2)) [3.5e-27; singlet and projection 3.8e-30]; n = 1 with sigma_A: Ghoshal K_B1 at Dirichlet [1.2e-27]
PASS continuation: K0 -> (h/2)(h/2+B/2)(h/2-1)(h/2-1-B/2)(h/2+1+B/2)/(h/2-1+B/2), n = 1..7 [1.0e-30]; particle amplitudes (minimal and sigma_B; n = 2, 3; l = 0, 1): L_n, unitarity, crossing with S_aa, bootstrap [1.0e-29]; sigma_A breaks the bootstrap for n = 2 [0.99]; B -> 0: L_1 -> g_(h/2) g_(h/2-1) E [4.1e-15]
```

Check 3's $G_j$-form test was added after this run, and a final rerun showed all 8 PASS. The docstring was then edited only.

Script's `a`/`bt` use the sign $s_-=+1$, i.e. prefactor $+i$; check 5 then shows the factor $(-1)^{n+1}$. The book formula should carry $-i$, under which the bootstrap is exact.

Suggested README row: `sp_scalar.py` | `sec-preserving-conjugating`, `eq-sp-scalar-eqs`, `prp-sp-scalar`, `eq-sp-breather`, open problem 11 | normalization $S^T(T-\mu)$ of the soliton–antisoliton S-matrix by crossing; SP crossing functions; scalar factors at 30 digits; $n=1$ = Ghoshal–Zamolodchikov Dirichlet; soliton bootstrap $V^{\otimes n}\to\bar V$; the CDD factors it allows; $B_1\to\bar B_1$ amplitude; continuation and particle axioms | about 4.5 minutes.

## 3. Proposed book text

**(a) ch. 33, `sec-preserving-conjugating`.** Replace "For $n=1$, where $\bar V\cong V$, they give Ghoshal and Zamolodchikov's Dirichlet factor." with "For $n=1$, where $\bar V\cong V$, they give Ghoshal and Zamolodchikov's Dirichlet matrix; its scalar factor is that of @prp-sp-scalar times $(-\frac T2)_I^2$." Then insert, at the end of the section (before `## Semiclassical check`):

```markdown
**Crossing-unitarity.** In the boundary state (-@eq-boundary-state) a soliton is paired with an antisoliton, so for a soliton-preserving boundary the crossed channel is soliton–antisoliton scattering. With the pairing $v_i\leftrightarrow\bar v_i$, crossing of the bulk S-matrix gives $\check S_{V\bar V}(\theta)=S^T(T-\mu)\,\check R_{V\bar V}(\theta)$, and the same for $\bar VV$, with $\check R$ equal to 1 on every transmission $v_i\otimes\bar v_j\mapsto\bar v_j\otimes v_i$, $i\ne j$, and $S^T$ the transmission amplitude of @prp-an-diagonal-scalar. Write the physical amplitudes as $a(\mu)K_V$ and $\tilde b(\mu)K_{\bar V}$, and let $k_A(\theta)=\sum_iK_{\bar V}(\frac{i\pi}2-\theta)_{ii}\,v_i\otimes\bar v_i$, $k_B(\theta)=\sum_iK_V(\frac{i\pi}2-\theta)_{ii}\,\bar v_i\otimes v_i$. Crossing-unitarity (-@eq-boundary-crossing-unitarity) consists of two equations, whose matrix parts are
$$
k_B(\theta)=-(-q)^{1-l}s(\mu)\,\check R_{V\bar V}(2\theta)\,k_A(-\theta),\qquad
k_A(\theta)=-(-q)^{l-1}s(\mu)\,\check R_{\bar VV}(2\theta)\,k_B(-\theta),\qquad s(\mu)=\frac{\sin\pi(2\mu-T)}{\sin\pi(2\mu+T)},
$$
independent of $\eta$. With $U(\mu)=\sinh(\eta+i\pi\mu)\sinh(\eta-i\pi\mu)$, $\bar U$ the same with $\bar\eta$, and $F(x)=S^T(T-x)$, the scalar factors must satisfy
$$
a(\mu)a(-\mu)U(\mu)=1,\quad \tilde b(\mu)\tilde b(-\mu)\bar U(\mu)=1,\quad
-(-q)^{1-l}s(\mu)\,a\big(\tfrac T2-\mu\big)=F(2\mu)\,\tilde b\big(\tfrac T2+\mu\big),\quad
-(-q)^{l-1}s(\mu)\,\tilde b\big(\tfrac T2-\mu\big)=F(2\mu)\,a\big(\tfrac T2+\mu\big).
$$ {#eq-sp-scalar-eqs}

::: {#prp-sp-scalar}
## Scalar factors of the soliton-preserving matrices
**[Numerical]** Let $\varepsilon=\eta/i\pi$, $\bar\varepsilon=\bar\eta/i\pi=\varepsilon+\omega l-T$, let $G_j$ be as in (-@eq-an-diagonal-scalar), and
$$
\mathcal G_{c,c'}(x)=\frac{G_0(x,c)}{G_0(x,c')}\prod_{j\ge1}\frac{G_j(x,c)\,G_j(x,1-c)}{G_j(x,c')\,G_j(x,1-c')},
$$
with $\mathcal G^{(2)}_{c,c'}$ the same with $G_{2j}$ in place of $G_j$. Then
$$
a(\mu)=-i\rho(\mu)\,e^{(1-l)\omega\theta}\,\frac{\mathcal G_{\varepsilon+2T,\,\bar\varepsilon+T}(\mu)}{\sin\pi(\varepsilon+\mu)},\qquad
\tilde b(\mu)=-i\rho(\mu)\,e^{-(1-l)\omega\theta}\,\frac{\mathcal G_{\bar\varepsilon+2T,\,\varepsilon+T}(\mu)}{\sin\pi(\bar\varepsilon+\mu)},\qquad
\rho(\mu)=\Big(\frac\omega2\Big)_I\,\mathcal G^{(2)}_{2T+\omega,\,0}(2\mu)
$$ {#eq-sp-scalar}
solve (-@eq-sp-scalar-eqs) for $n=1,2,3$ and every $l$. They satisfy the soliton bootstrap $V^{\otimes n}\to\bar V$ for $n\le3$, and for $n=4$ at $l=0,2,5$: $\bar V(\theta)$ sits in $V(\theta_1)\otimes\dots\otimes V(\theta_n)$ with $\mu_k=\mu-T+k\omega$, and reflecting the $n$ constituents with $aK_V$ and the bulk S-matrix gives $\tilde b(\mu)K_{\bar V}(\theta)$, matrix and scalar. For $n=1$, $l=1$, $a(\mu)(-\frac T2)_I^2K_V$ is Ghoshal and Zamolodchikov's Dirichlet matrix at $\xi_B=\pi\varepsilon-\frac\pi2$, the limit $\vartheta\to\infty$ of (-@eq-bsg-k-full).
:::

Checked in `textbook_code/sp_scalar.py`, the functional equations to $10^{-26}$ and the matrix parts to $10^{-14}$. The equations (-@eq-sp-scalar-eqs) are invariant under $a\to a\,e^{m\theta}$, $\tilde b\to(-1)^m\tilde b\,e^{-m\theta}$, $m\in\mathbb Z$; the soliton bootstrap selects $m=0$, and it fixes the overall sign. The factor $e^{(1-l)\omega\theta}$ centres the exponents of $K_V$: $e^{(1-l)\omega\theta}K_V$ has the exponents $-(l-1),\dots,l-1$ on the block and $-(n-l),\dots,n-l$ off it.

**CDD factors.** The remaining freedom is $a\to\sigma(\mu)a$, $\tilde b\to\sigma(T-\mu)\tilde b$ with $\sigma(\mu)\sigma(-\mu)=1$ and $\sigma$ of period $2T$ in $\mu$. Unlike the soliton-conjugating case (-@eq-an-diagonal-scalar-eqs), this includes single blocks $(x)_I$. A block $(x)_I$ in $\sigma$ multiplies the soliton bootstrap by $\sin\frac{\pi(\mu-T-x)}\omega\big/\sin\frac{\pi(\mu-T+x)}\omega$, which is $1$ for $x\in\omega\mathbb Z$ and $-1$ for $x\in\omega(\mathbb Z+\frac12)$; for even $n$ a factor $-1$ can be absorbed in the overall sign. So the bootstrap does not fix the CDD factor. It excludes, for even $n$, the factor $(-\frac T2)_I^2$ that reproduces Ghoshal and Zamolodchikov at $n=1$, while other factors that reduce to it at $n=1$ pass, such as $-[(\frac T2)_I^2(\frac T2-\frac\omega2)_I(\frac T2+\frac\omega2)_I]^{-1}$ and $\sigma=(-\frac\omega2)_I^2$.
```

**(b) ch. 35, `sec-breather-reflection`.** Insert a paragraph after "**Other soliton reflection matrices.** …":

```markdown
**Soliton-preserving matrices.** For the matrices of @sec-preserving-conjugating the lowest breather comes back as its conjugate. Reflecting the singlet in $V(\theta_a)\otimes\bar V(\theta_b)$ with $\tilde bK_{\bar V}$ and $aK_V$ of @prp-sp-scalar, with the bulk factor $F(2\mu)$ in between, gives the singlet in $\bar V(-\theta_b)\otimes V(-\theta_a)$ times
$$
K_{1\to n}(\mu)=e^{-2\pi il/h}K_0(\mu)E_\eta(\theta),\qquad
K_0=\Big(\frac T2\Big)_I\Big(\frac T2-\frac12\Big)_I(\delta)_I\Big(\frac12-\delta\Big)_I\Big(\frac12+\delta\Big)_I\Big(T-\frac12-\delta\Big)_I,\quad\delta=\frac{T-\omega}2,
$$
$$
E_\eta(\theta)=\frac{\sinh\theta-i\sin\frac{\pi A}T}{\sinh\theta+i\sin\frac{\pi\bar A}T},\qquad A=\varepsilon+\frac{T-1}2,\quad\bar A=A+\omega l-T,
$$ {#eq-sp-breather}
times $\sigma(\mu_a)\sigma(T-\mu_b)$ for a CDD factor $\sigma$. The reverse transition $\bar B_1\to B_1$ has $A\leftrightarrow\bar A$ and the phase $e^{2\pi il/h}$. The phase depends on the relative normalization of the two singlets. Here each has coefficient 1 on the highest weight of its left factor times the lowest weight of its right factor, which for $n=1$, where the two singlets coincide, is not a choice. For $l=0$ and $l=h$, $E_\eta=1$, since $\eta$ enters $K_V$ only through a scalar. For $n=1$, $l=1$ and $\sigma=(-\frac T2)_I^2$, (-@eq-sp-breather) is Ghoshal's Dirichlet $K_{B_1}$, (-@eq-bsg-b1-reflection) with $\vartheta\to\infty$. Continued by $\omega\to-2/B$,
$$
K_0\to\frac{(\frac h2)(\frac h2+\frac B2)(\frac h2-1)(\frac h2-1-\frac B2)(\frac h2+1+\frac B2)}{(\frac h2-1+\frac B2)} .
$$
The particle amplitudes $L_a=K_{a\to h-a}$ given by the bulk boundary bootstrap, $L_c(\theta)=L_a(\theta+iu_a)\,S_{a,h-b}(2\theta+iu_a-iu_b)\,L_b(\theta-iu_b)$, satisfy unitarity $L_a(\theta)L_{h-a}(-\theta)=1$ and crossing-unitarity $L_a(\frac{i\pi}2-\theta)=S_{aa}(2\theta)L_a(\frac{i\pi}2+\theta)$; for reflection into the conjugate the pair in the boundary state is $aa$, so $S_{aa}$ replaces $S_{a\bar a}$ of (-@eq-boundary-crossing-diagonal). The bootstrap closes at every fusing, and $L_n$ is the soliton-side $K_{n\to1}$, for the minimal factor and for the CDD factors allowed by the soliton bootstrap that we tested; $(-\frac T2)_I^2$ breaks it for even $n$. As $B\to0$, $L_1\to e^{-2\pi il/h}g_{h/2}\,g_{h/2-1}\,E_\eta$ with $g_x=(\sinh\theta+i\sin\frac{\pi x}h)/(\sinh\theta-i\sin\frac{\pi x}h)$, for all these factors, which depends on $\theta$ when $h\ge3$ (`textbook_code/sp_scalar.py`).
```

**(c) ch. 35, open problem 11.** Replace it with:

```markdown
11. **Soliton-preserving boundaries.** The soliton-preserving matrices now have scalar factors (@prp-sp-scalar) and a lowest-breather amplitude (@sec-breather-reflection), and they satisfy the soliton bootstrap. For $n\ge2$ they reflect particle $a$ into $h-a$, so they describe neither Dirichlet nor Neumann boundaries (@sec-preserving-conjugating). Three things are open. The first is the CDD factor, which the soliton and particle bootstraps leave free: for $n=1$, Ghoshal and Zamolodchikov's Dirichlet amplitude needs $(-\frac T2)_I^2$, which is excluded for even $n$. The second is the classical boundary that these matrices describe for $n\ge2$. The third is their coideal in the classification of Regelskis and Vlaar.
```

## 4. Contradictions with the book, Appendix D or the worker's report
- **ch. 28, below (eq-boundary-crossing-unitarity), line 142.** "with $\check S(2\theta)$ acting on $V(-\theta)\otimes V(\theta)$" is a slip: $k(-\theta)\in V(\theta)\otimes V(-\theta)$. It should read "mapping $V(\theta)\otimes V(-\theta)$ to $V(-\theta)\otimes V(\theta)$".
- **ch. 33, "they give Ghoshal and Zamolodchikov's Dirichlet factor".** Correct for the matrix; the scalar factor needs $(-\frac T2)_I^2$ (fix in 3a).
- **Appendix D.** Nothing in it is contradicted by my results. The worker's claim 4 (Delius 1998b (1.4)) was not checked here.
- **Worker's report, corrected items:**
  - $k$ is not "fixed numerically on the minimal branch". Closed form: total exponential $e^{(1-l)\omega\theta}$. Their principal-branch numerics give a different solution, excluded by the bootstrap, when $|(1-l)\omega/2|>\frac12$, e.g. $(n,l)=(3,4)$ at ω = 0.37.
  - The $n=1$ ratio has the opposite sign.
  - The breather constant is $l$-dependent and convention-dependent for $n\ge2$.
  - $\sigma_A$ is not viable for even $n$.
  - $\sigma_B$ is viable.
- **Side note for problem 3's consolidation.** `briefs/halfline_open_03_scripts/sol_boot.py` defaults to `order=+1` (largest-imaginary-part constituent on the left). For $n\ge2$ this ordering has no embedding $\bar V\to V^{\otimes n}$ in book conventions; only the reversed order has one.

## 5. Files
All in `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a2/`, copied to `briefs/halfline_consolidation_A2_scripts/`:

| File | Purpose |
|---|---|
| `sp_scalar.py` | the check script for `textbook_code/` |
| `sp_scalar.log` | its output |
| `sp_lib.py` | book-code library: SP matrices, `rcheck` cache with gesvd fallback, `tensor_rep` |
| `sp_mp.py` | 30-digit Ê, Λ, $F_{11}$, scalar factors |
| `t1_crossing.py`, `t1_crossing.log` | task 1 |
| `t2_scalar.py`, `t2_scalar.log` | task 2 equations |
| `t2_gz.py`, `t2_gz.log` | $n=1$ GZ comparison |
| `t2_gform.py` | $G_j$ form of the scalar factor |
| `t3_breather.py`, `t3_breather.log` | projection, double and 30-digit |
| `t3_amp.py`, `t3_amp_n1.log`, `t3_amp_n2.log`, `t3_amp_n3.log`, `t3_amp_n4.log` | breather amplitude, $n=1..4$ |
| `t3_ghoshal.py`, `t3_ghoshal.log` | $n=1$ breather vs Ghoshal |
| `t4_probe.py` | embedding order |
| `t4_bootstrap.py`, `t4_boot_n1.log`, `t4_boot_n2.log`, `t4_boot_n3.log` | soliton bootstrap, $n=1,2,3$ |
| `t4_boot_n4.py`, `t4_boot_n4.log` | soliton bootstrap, $n=4$ |
| `t4_boot_om.py`, `t4_boot_om.log` | soliton bootstrap at ω = 1.29 |
| `t4_cdd.py` | CDD candidates vs bootstrap |
| `t5_particles.py`, `t5_particles.log` | continued particle axioms |
| `t5_limits.py`, `t5_limits.log` | $B\to0$ limits and continued CDD blocks |
