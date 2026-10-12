# Deep referee report on chapters 36–37

Date: 2026-10-12. Model: GPT-6 (Codex). Starting commit: `2b6d6a3`; completion check: `4051990` (`git log -1 --format=%h`). Concurrent agents advanced HEAD during the review; the four audited chapter files are unchanged between these commits.

This follows `briefs/referee_defects_36_37_deep.md`. No pre-existing file was edited, no index-changing Git command was used, and Quarto was not run. All computations below are new: neither scratch script imports or runs the book's verification code. Line references refer to the chapter versions inspected during the review.

Correction targets throughout this report:

* Chapter 36: `part8-defects/36-defect-scattering.qmd`.
* Chapter 37: `part8-defects/37-defect-sine-gordon.qmd`.
* The chapter-27 comparison target is `part6-soliton-smatrices/27-rsos.qmd`; boundary crossing is in `part7-half-line/28-boundary-scattering.qmd`.

Reproduce the numerical evidence with:

```bash
python3 referee_scratch/defects_36_37_deep/independent_checks.py
python3 referee_scratch/defects_36_37_deep/supplemental_checks.py
```

The scripts use NumPy, SciPy, SymPy, and only each other's new code. Mathematical derivations below provide independent evidence for each of the eight items, including the limitations of the numerical tests. The source digests and the previous referee report were treated as claims to check, not as evidence. Source equation numbers were checked in the supplied `~/referee-copy/papers/*.txt` files.

The main results are: doubled crossing passes the requested non-diagonal check; the fixed-convention $q$-sign test passes but admits an explicit Klein-factor interpretation; the diagonal no-go proof needs an additional hypothesis; the UV–IR relation survives; the scalar needs a common-cutoff definition and the resonance discussion needs a coupling restriction; the claimed chapter-27 RSOS comparison cannot be made from what that chapter actually prints; the half-soliton identities pass; and the Lee–Yang finite-size paragraph confuses real-rapidity and crossed-channel phase properties.

## 1. Doubled crossing-unitarity with sine-Gordon transmission

**Verdict: correct.** The indices in chapters 28 and 36 fit together, including the unchanged defect indices. This assertion does not require diagonal particle scattering.

### Independent derivation

Put $L=P\mathbb T$, and write its entries as $L_{b\beta,a\alpha}=\mathbb T_{a\alpha}^{\beta b}$. In the doubled space label a particle by $A=(\ell,a)$, where $\ell=0,1$ specifies the copy, and conjugate only $a$. For each fixed input/output defect pair define the operator-valued pair coefficient

$$
k^{AB}_{\beta\alpha}(u)=K_{B\beta,\bar A\alpha}(i\pi/2-u).
$$

Because $K$ is off-diagonal in the copy label, this coefficient vanishes unless $A$ and $B$ belong to different copies. Consequently the *entire* doubled crossing equation is

$$
K_{B\beta,\bar A\alpha}(i\pi/2-u)
=K_{A\beta,\bar B\alpha}(i\pi/2+u).
$$

The same-copy blocks $\check S$ and $P\check SP$ multiply zero. With $A=(0,\bar a)$, $B=(1,b)$, and $v=i\pi/2-u$, this becomes

$$
\mathbb T_{a\alpha}^{\beta b}(v)
=\bar{\mathbb T}_{\alpha\bar b}^{\bar a\beta}(i\pi-v),
$$

exactly `eq-defect-crossing`. There is **no transpose of $\mathcal W$**. Treating all of $V\otimes\mathcal W$ as the crossed particle and taking a full transpose would give a different, incorrect equation.

### Independent numerical construction

`independent_checks.py`, function `item_1_2`, constructs:

* the chapter-37 transmission matrix on $-5\le n\le5$, including its scalar;
* $\bar L(v)=L(-v)^{-1}$, with the defect-rapidity shift applied before inversion;
* the full $44\times44$ doubled $K$;
* the $16\times16$ doubled bulk braid matrix, with the physical sine-Gordon $S_0$ evaluated independently from its integral, $\check S$ and $P\check SP$ on the same-copy blocks, and the flip on mixed blocks.

It forms the chapter-28 pair coefficient entry by entry and tests crossing for all particle indices and $-3\le\alpha,\beta\le3$. Excluding the window endpoints avoids fictitious truncation effects.

```text
lambda=1.37: doubled crossing=2.28e-15, negative control=0.0311
lambda=2.6:  doubled crossing=4.20e-15, negative control=0.0363
lambda=0.73: doubled crossing=4.50e-16, negative control=0.0416
```

The negative control multiplies one copy block of the pair coefficient by 1.07. Thus this is not a test that succeeds merely because all the tested components vanish. The proof also explains why testing the non-diagonal defect is useful even though the non-diagonal *same-copy bulk* amplitudes drop out of this particular identity.

No correction to `prp-doubling` item 3 is needed. At chapter 36, line 177, the following evidence update is ready to append:

```markdown
Item 3 also passes an independent test with the non-diagonal sine-Gordon transmission matrix, including its scalar, on $|n|\le5$ with defect indices restricted to $|n|\le3$: residuals are below $5\times10^{-15}$. The boundary pair coefficient is operator-valued on $\mathcal W$; crossing exchanges the particle indices and leaves both defect indices unchanged (`referee_scratch/defects_36_37_deep/independent_checks.py`).
```

## 2. Sign of q, ordinary gauges, and Klein factors

**Verdict: correct in the book's fixed exchange and conjugation conventions.** The assertion that only $C'$ works needs to be read as an assertion about the *displayed matrix*, not as a prohibition of changes to the exchange algebra. An explicit Klein-factor change relates the two signs.

### Independent block computation

In the fixed-total-charge block with basis $\{v_+|n\rangle,v_-|n+1\rangle\}$, the scalar-free matrix is

$$
M_n(Y)=
\begin{pmatrix}1&-iY\\-iY&1\end{pmatrix}
\begin{pmatrix}q^{-n}&0\\0&q^{n+1+s}\end{pmatrix}.
$$

For real $Y$, real $s$, and $|q|=1$,

$$
M_n^\dagger M_n=(1+Y^2)I,
\quad
M_n^{-1}=\frac1{1+Y^2}
\begin{pmatrix}q^n&0\\0&q^{-n-1-s}\end{pmatrix}
\begin{pmatrix}1&iY\\iY&1\end{pmatrix}.
$$

Thus the matrix part is unitary **after normalization**, at either sign of $q$. Its unnormalized entries are not themselves a unitary matrix.

Cross the particle indices of the inverse, leaving the defect indices in place. Write $Y'=e^{i\pi\lambda}Y$. Comparing a preserving entry and a converting entry gives respectively

$$
\Phi(\theta)\Phi(\theta+i\pi)(1+Y^2)=q^{-s},
\qquad Y'=-Y/q.
$$

The second relation holds at $q=-e^{-i\pi\lambda}$. At $q'=+e^{-i\pi\lambda}$ it has the opposite sign. No scalar can repair that relative sign. Replacing $C=\sigma_x$ by $C'=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)$ changes the sign of precisely the converting crossed entries and repairs it.

Independent output, with the scalar included:

```text
lambda=1.37, qphys:  unitarity 4.44e-16; C crossing 2.96e-13; C' crossing 1.52
lambda=1.37, -qphys: unitarity 4.44e-16; C crossing 1.52; C' crossing 2.96e-13
lambda=2.60, qphys:  unitarity 2.22e-16; C crossing 8.97e-13; C' crossing 1.99
lambda=2.60, -qphys: unitarity 2.22e-16; C crossing 1.99; C' crossing 8.97e-13
```

The conclusion about the relative sign is independent of the scalar. The scalar is needed to turn the correctly matched matrix equations into exact crossing and physical unitarity, and it passes those independent checks.

### What can and cannot be gauged away

A charge-preserving change of the one-particle and defect bases does not change the diagonal coefficient $q^{-n}$. It cannot change it to $(-1)^nq^{-n}$ for every $n$. More generally, an ordinary constant basis change transports the crossing bilinear form by congruence, which preserves whether the form is symmetric or antisymmetric. It cannot turn $C'$ into $C$. This is the same distinction made in chapter 21, `sec-crossing-sign`: a change acting on ordered multiparticle states is not a common one-particle basis change.

There is, however, a concrete Klein-factor equivalence. Let

$$
\kappa_+^2=\kappa_-^2=1,\qquad
\kappa_+\kappa_-=-\kappa_-\kappa_+,
\quad J=\kappa_+\kappa_-.
$$

Redefine the ZF generators in an enlarged Klein algebra by

$$
A'_\pm=\kappa_\pm A_\pm,\qquad D'_n=J^nD_n.
$$

For all four allowed even-sector transmissions one obtains

$$
\kappa_aJ^n=(-1)^nJ^\beta\kappa_b,
\quad (a,n)\longrightarrow(\beta,b).
$$

Hence $\mathbb T'_{an}^{\beta b}=(-1)^n\mathbb T_{an}^{\beta b}$, exactly $\mathbb T(-q)$. The script verifies this identity with explicit $2\times2$ Pauli matrices, for all four transitions and $-5\le n\le5$, with residual zero. Equivalently, for the displayed numerical matrices $L(-q)=L(q)(1\otimes\Pi_{\mathcal W})$, $\Pi_{\mathcal W}|n\rangle=(-1)^n|n\rangle$.

The same Klein transformation reverses the unlike-charge bulk *transmission* coefficient: the ordered products $\kappa_+\kappa_-$ and $\kappa_-\kappa_+$ differ by a minus sign. It leaves the bulk reflection coefficient unchanged. It also changes the annihilation/crossing pairing from symmetric to antisymmetric. Undoing the Klein transformation restores ordinary $C$, but simultaneously restores the original exchange amplitudes. It does **not** make the unchanged numerical matrix at (-q) satisfy the unchanged $C=\sigma_x$ equation.

Thus there is a genuine selection **after fixing the ZF algebra, locality/statistics convention, and crossing form**, as this book does. There is no stronger convention-independent selection excluding Klein-equivalent descriptions. Their finite-volume sectors must be transported too, as chapter 21 explains.

Suggested precision edit at chapter 37, line 204, replacing its last sentence:

```markdown
For the displayed matrix at $-q$, the off-diagonal crossing equations require the opposite sign from the diagonal ones, so no scalar factor restores crossing with $C=\sigma_x$; the antisymmetric $C'=\big(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\big)$ does. This selects $q_{\rm phys}$ in the fixed exchange conventions of this book. A Klein-factor redefinition can relate the two signs, but changes the bulk exchange amplitudes and the crossing pairing together, as in @sec-crossing-sign; it is not a common change of the one-particle basis.
```

Also replace “the matrix part is unitary at both signs” on that line by “the matrix part has $M^\dagger M=$1+Y^2$I$ at both signs.”

## 3. Reflection–transmission no-go theorem

**Verdict: wrong as a fully specified theorem; correct after adding species-preserving scalar defect amplitudes and the usual independence of asymptotic states.** Its rapidity argument and both displayed scalar reductions are correct. The hypotheses actually printed do not explicitly impose scalar/species-preserving $R$ and $\mathbb T$; absence of defect states does not impose that condition on a degenerate particle multiplet.

### Both reductions, independently

Use $A_a(x)A_b(y)=S_{ab}(x-y)A_b(y)A_a(x)$, and suppose

$$
A_a(t)D=R_a(t)A_a(-t)D+T_a(t)DA_a(t).
$$

For $t_1>t_2>0$, start from $A_a(t_1)A_b(t_2)D$. The two coefficients are:

| Out-state | Let $b$ reach $D$ first | Exchange (a,b) first |
|---|---|---|
| $A_b(-t_2)DA_a(t_1)$ | $R_b(t_2)S_{ab}(t_1+t_2)T_a(t_1)$ | $S_{ab}(t_1-t_2)T_a(t_1)R_b(t_2)$ |
| $A_a(-t_1)DA_b(t_2)$ | $T_b(t_2)R_a(t_1)$ | $S_{ab}(t_1-t_2)R_a(t_1)S_{ba}(t_1+t_2)T_b(t_2)$ |

In the bottom-right entry, after reflection of $a$, the remaining exchange is $A_b(t_2)A_a(-t_1)$, whose rapidity difference is $t_2-(-t_1)=t_1+t_2$. This explicitly establishes the (ba) order and the plus sign.

If $R_b,T_a$ are nonzero meromorphic functions, choose an open real rectangle avoiding their isolated zeros and poles and lying in $t_1>t_2>0$. On it,

$$
S_{ab}(t_1+t_2)=S_{ab}(t_1-t_2).
$$

The invertible coordinate change $u=t_1-t_2, v=t_1+t_2$ has Jacobian 2. Thus this equality holds on an open two-dimensional set in which $u$ and $v$ vary independently. Differentiating at fixed $v$ gives $S'_{ab}(u)=0$ on an interval, and the identity theorem makes $S_{ab}$ constant on its connected meromorphic domain. No asymptotic bound, positivity, parity, or second mixed equation is needed for **constancy**. The second equation gives $S_{ab}(u)S_{ba}(v)=1$. For $a=b$, the two give $S_{aa}^2=1$.

Associativity here means equality of coefficients in a nondegenerate ZF scattering algebra. Associativity of an arbitrary quotient algebra, in which the out-states can vanish or be linearly dependent, would not suffice. Relativistic difference dependence is also used explicitly; it is already implicit in the notation $S_{ab}(\theta)$.

### Trying the defect-state loophole

Two natural noncommuting constructions fail without establishing a general no-go:

1. Keep the Weyl transmission of chapter 37 and try a species-preserving reflection $R=r(\theta)\,1_V\otimes H(N)$, with nonzero diagonal $H$. Its commutator with the converting part of $T$ is generally nonzero. Nevertheless, the two-incoming-soliton channel in which both particle charges are preserved contains $S_0(t_1+t_2)-S_0(t_1-t_2)$, multiplied by a nonzero diagonal transmission/reflection coefficient. No other intermediate path has the same final charges in this channel. The attempted reflection therefore vanishes for interacting sine-Gordon.
2. Try $R=r\,1_V\otimes\sigma_x, T=t\,1_V\otimes\sigma_z$, genuinely anticommuting in the impurity space. The same-charge mixed channel now demands $S_0(t_1+t_2)=-S_0(t_1-t_2)$. Analytic continuation to $t_2=0$ forces $S_0=0$, incompatible with bulk unitarity. Allowing (r(0)) or (t(0)) to vanish does not help: the meromorphic identity is first derived at generic rapidities.

`supplemental_checks.py` evaluates the respective offending coefficients using an independent integral for $S_0$, at $(\lambda,t_1,t_2)=(1.37,0.9,0.3)$. Their absolute residuals are `0.1500080375` and `1.9943664630`, respectively. Neither vanishes. This is a concrete failed construction in a non-diagonal bulk theory, not a claim to settle all operator-valued reflection/transmission possibilities.

There is also a qualification missing from the broad CAFG summary at chapter 36, line 85. Cancelling a *matrix* transmission factor requires appropriate nondegeneracy. A direct sum of a purely reflecting defect sector and a purely transmitting sector has globally nonzero $R$ and $T$, commuting defect indices, and $RT=TR=0$. It works with an interacting bulk: in each superselection sector the corresponding mixed equations are zero. It is not simultaneous reflection and transmission within a communicating sector and not the requested genuinely mixed counterexample, but it disproves an unqualified assertion that global $R\ne0\ne T$ plus the displayed commutativity condition always forces a constant bulk matrix. Neither meromorphy nor the commutativity condition licenses cancelling a singular matrix.

I constructed no interacting, irreducible example with genuine mixing. “No interacting example is known” is not established by these calculations; the report should not promote it to a theorem.

### Ready-to-paste corrections

Chapter 36, line 70, insert immediately after “no internal states”:

```markdown
, that reflection and transmission preserve each particle label, $R_a^{\ b}=\delta_a^bR_a$ and $\mathbb T_a^{\ b}=\delta_a^b\mathbb T_a$, and that the asymptotic out-states are linearly independent
```

Chapter 36, line 85, replace the sentence beginning “Then a mixed equation” with:

```markdown
When the transmission matrix can be cancelled in the mixed relation, this commutativity condition leads to the rapidity-independent result of [@cafg2002, eqs. (22)–(24)]. Cancellation requires a nondegeneracy assumption: a direct sum of a reflecting sector and a transmitting sector has $RT=TR=0$ and is not excluded by this argument. The genuinely mixed problem concerns reflection and transmission acting within the same communicating defect sector.
```

## 4. Defect rapidity renormalization versus potential renormalization

**Verdict: correct as a conjectural UV–IR relation.** The two renormalizations concern different coordinates on the defect coupling space. Their rapidity statements do not contradict each other. Their scale statements differ if one identifies schemes and bare mass parameters without an additional matching calculation. The chapter is right to describe a tension rather than a resolved contradiction.

### Independent rederivation from Bajnok–Simon

Write $b_{\rm sh}$ for the real sinh-Gordon coupling and $g=b_{\rm sh}^2/(8\pi)$. In [Bajnok–Simon](https://arxiv.org/abs/0712.4292),

$$
B_{\rm BS}=\frac{g}{1+g},\qquad
M e^{\pm g\vartheta}=M_{\rm cl}e^{\pm\mu}\quad\text{(17)}.
$$

Dividing the two signs gives $g\vartheta=\mu$ on the real branch; multiplying them gives $M=M_{\rm cl}$ in the normalization used for the quantum potential. Equation (16) has $B_{\rm BS}\vartheta=-\theta_d$, whereas $\mu=-\theta_d^{\rm cl}$. Therefore

$$
\frac{\theta_d}{\theta_d^{\rm cl}}
=\frac{B_{\rm BS}}g=\frac1{1+g}
=1-B_{\rm BS}=1-\frac B2.
$$

On continuation $b_{\rm sh}^2\mapsto-\beta_{\rm SG}^2$, put $x=\beta_{\rm SG}^2/(8\pi)=1/(\lambda+1)$. Then

$$
\theta_d=\frac{\theta_d^{\rm cl}}{1-x},\qquad
\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}.
$$

This is precisely the chapter's formula. It uses the assumed exact boundary UV–IR relation; it is not a derivation of quantum integrability from the defect Lagrangian.

### Independent comparison of the two coupling coordinates

Write the two potential coefficients as $c_+=M_de^{\mu}$ and $c_-=M_de^{-\mu}$. They determine

$$
M_d=\sqrt{c_+c_-},\qquad \mu=\tfrac12\log(c_+/c_-).
$$

Multiplying the full potential by $\sqrt{1+g}$, as in [Penati–Polvara](https://arxiv.org/abs/1902.10690), (4.12), changes $M_d$ but leaves $\mu$, and hence the bare $\sigma$, unchanged. Thus it is not a competing formula for the rapidity renormalization.

The scale comparison is more informative than saying merely “different schemes.” From Bajnok–Simon (14), with their classical $M_{\rm cl}^{(0)}=4m_{\rm cl}/b_{\rm sh}^2$,

$$
\frac{M}{M_{\rm cl}^{(0)}}
=\sqrt{\frac{b_{\rm sh}^2}{8\sin(b_{\rm sh}^2/8)}}
=\sqrt{\frac{\pi g}{\sin\pi g}}
=1+\frac{\pi^2g^2}{12}+O(g^4).
$$

Penati–Polvara instead require

$$
\sqrt{1+g}=1+\frac g2-\frac{g^2}8+O(g^3).
$$

These series are printed independently by the scratch script. At an identically fixed bulk mass/potential normalization they disagree already at (O(g)). Bajnok–Simon explicitly allow scheme dependence of the scale after (17); their scale and the normal-ordered massless-current calculation have not been matched here. The higher-spin anomaly is a further issue, which a one-loop two-point/transmission calculation does not settle. This supports the cautious status in the chapter.

One normalization qualification is necessary. Penati–Polvara's real-coupling identity is $\gamma/2=\beta^2/(4\pi)=b_{\rm sh}^2/(8\pi)$. Upon continuation to *positive real sine-Gordon* $\beta_{\rm SG}^2$, it becomes **minus** $\beta_{\rm SG}^2/(8\pi)$. Using the symbol $\beta_{\rm SG}$ for the positive real sinh-Gordon coupling obscures this sign.

Chapter 37, line 229, replace the first two sentences with:

```markdown
Penati and Polvara computed the quantum conservation laws of the Lagrangian defect in massless perturbation theory [@penati_polvara2019]. For real sinh-Gordon coupling $b$, momentum conservation requires the whole defect potential to be multiplied by $\sqrt{1+\gamma/2}$, where $\gamma/2=b^2/(8\pi)$ [@penati_polvara2019, eq. (4.12)]. Under continuation to sine-Gordon this becomes $\sqrt{1-\beta_{\rm SG}^2/(8\pi)}$. This changes the common scale of the two exponential terms but leaves their ratio, and hence the bare $\sigma$, unchanged.
```

The remainder of the existing paragraph can retain its statement about unmatched schemes and higher-spin anomalies.

## 5. Minimal scalar, convergence, poles, and the resonance

**Verdict: wrong as literally defined, but the intended common-cutoff scalar is correct.** The auxiliary product $r$ itself diverges. Its paired quotient has a well-defined limit and satisfies the required equations. The absence of physical-strip poles is correct. The unqualified decay interpretation at all positive $\lambda$ is false.

### Deriving a minimal scalar independently

Set $z=i\lambda\theta/(2\pi)$, $Y=e^{\lambda\theta}$, and define

$$
F_k(z)=\frac{\Gamma(k\lambda+\tfrac14-z)
\Gamma((k+1)\lambda+\tfrac34-z)}
{\Gamma((k+\tfrac12)\lambda+\tfrac14-z)
\Gamma((k+\tfrac12)\lambda+\tfrac34-z)},
\quad
P(z)=\lim_{N\to\infty}\prod_{k=0}^{N-1}\frac{F_k(z)}{F_k(-z)}.
$$

The charge-block computation in item 2 fixes

$$
|\Phi(\theta)|^2=\frac1{1+Y^2},\qquad
\Phi(\theta)\Phi(\theta+i\pi)=\frac1{1+Y^2}.
$$

Factor off the simple real-axis modulus by writing $\Phi=A P(z)/(1+iY)$. Shifting $\theta\mapsto\theta+i\pi$ shifts $z\mapsto z-\lambda/2$. Pair adjacent half-step Gamma factors. At finite cutoff the products telescope, and the remaining endpoint ratio tends to one, giving

$$
P(z)P(z-\lambda/2)
=\frac{\Gamma(\tfrac14-z)\Gamma(\tfrac34+z)}
{\Gamma(\tfrac34+\lambda/2-z)\Gamma(\tfrac14-\lambda/2+z)}
=\frac{\sin\pi(\tfrac14-\lambda/2+z)}{\sin\pi(\tfrac14-z)}.
$$

The last expression equals

$$
e^{-i\pi(1+\lambda)/2}\frac{1+i e^{i\pi\lambda}Y}{1-iY}.
$$

Thus $A=e^{i\pi(1+\lambda)/4}$ gives exactly the required crossing equation. For real $\theta$, $z$ is imaginary and $P(z)^*=P(z)^{-1}$, proving physical unitarity. The overall opposite sign is also an axiom solution. More generally one may multiply by a defect CDD factor $h(\theta)h(\theta+i\pi)=1$, of unit modulus on the real axis. “Minimal” selects the indicated Gamma divisor/asymptotic prescription with no such extra factor; unitarity and crossing alone do not prove uniqueness.

As an independent check on the Gamma evaluation, its logarithm also follows directly from the integral representation of $\log\Gamma$:

$$
\log P(z)=\int_0^\infty\frac{dt}{t}\,2\sinh(zt)
\frac{e^{-t/4}(1-e^{-(\lambda+1)t/2})}
{(1-e^{-t})(1+e^{-\lambda t/2})}.
$$

This integral initially converges for $|\Re z|<1/4$; elsewhere use meromorphic continuation. The independently evaluated integral and paired Gamma product agree to $1.42\times10^{-15}$ at $(\lambda,\theta)=(1.37,0.61)$.

```text
lambda=0.73: scalar crossing 7.90e-15; real-axis unitarity 2.22e-16
lambda=1.37: scalar crossing 1.43e-13; real-axis unitarity 1.11e-16
lambda=2.60: scalar crossing 3.62e-13; real-axis unitarity 0
```

The numerical implementation derives the asymptotic coefficients from Bernoulli polynomials and sums the remaining tail with Hurwitz zeta. It does not borrow the book's evaluator.

### Comparison with Konik–LeClair and BCZ

The paired formula is BCZ (3.43)–(3.44). The numerator in Konik–LeClair (5.26) uses $3/4+k\lambda$ and $1/4+(k+1)\lambda$ instead. Its paired product $P_{\rm KL}$ obeys

$$
\frac{P_{\rm KL}(z)}{P(z)}
=\frac{\Gamma(\tfrac34-z)\Gamma(\tfrac14+z)}
{\Gamma(\tfrac14-z)\Gamma(\tfrac34+z)}
=\frac{\sin\pi(\tfrac14-z)}{\sin\pi(\tfrac14+z)}.
$$

The accompanying Konik–LeClair denominator $2\cosh((\lambda\theta-i\pi/2)/2)$, together with $e^{-\lambda\theta/2}$, is $e^{-i\pi/4}/(1-iY)$. Combining these factors reproduces the BCZ scalar exactly. The different-looking pairings are not different CDD solutions. `supplemental_checks.py` checks the finite products approaching this common limit.

For convergence, use

$$
\log\Gamma(x+a)\sim (x+a-\tfrac12)\log x-x+\tfrac12\log2\pi
+\frac{B_2(a)}{2x}+\cdots.
$$

The four shifts in $F_k$ give

$$
\log F_k(z)=\frac{\lambda+1}{4k}+O(k^{-2}).
$$

The coefficient is independent of $z$. Thus

$$
r_N(z)\sim c(z)N^{(\lambda+1)/4},
$$

whereas $\log F_k(z)-\log F_k(-z)=O(k^{-2})$ and the paired quotient converges normally away from its divisor. Numerically, at $\lambda=1.37$,

```text
r_(2N)/r_N, N=1000:  1.50785738 + 0.0000433679 i
predicted magnitude limit: 2^((lambda+1)/4) = 1.5078574026
```

The compact product in BCZ (3.47) instead has a $k$-th logarithm $-2z/k+O(k^{-2})$. Hence its phase really does drift by $-\lambda\theta\log N/\pi$, as the chapter says. Both (3.42) and (3.47) need regularization. The same computation gives, at $\lambda=1.37,\theta=0.61$,

```text
compact log-phase(2N)-log-phase(N), N=10000: -0.184381846
predicted limit:                                 -0.184385171
```

This is an omitted endpoint normalization, not an inconsistency of the desired scalar. If $C_N(z)$ denotes the right-hand side of BCZ (3.47) at cutoff $N$, the exact finite-cutoff identity is

$$
\frac{r_N(z)}{r_N(-z)}=C_N(z)\,2^{-4z}
\frac{\Gamma(N\lambda+\tfrac34-z)}{\Gamma(N\lambda+\tfrac34+z)}
\frac{\Gamma(2N\lambda+\tfrac12+2z)}{\Gamma(2N\lambda+\tfrac12-2z)}.
$$

Its endpoint factor is asymptotic to $(N\lambda)^{2z}$. This identity is also checked in the supplemental script. The chapter is right about the unregularized compact products, but wrong if its sentence means that the standalone product (3.44) converges.

### Divisor and physical strip

For $\theta=x+iu$, $z=i\lambda x/(2\pi)-\lambda u/(2\pi)$. Poles of Gamma factors can occur only at $x=0$, relative to real $\theta_d$.

For $u>0$, the numerator Gamma functions of (r(z)) have positive arguments. The possible upper-half-plane poles introduced by the reciprocal of the *denominator* of (r(-z)) begin at

$$
u=\pi+\frac{\pi}{2\lambda},
$$

already above both the physical and crossed physical strips. The elementary denominator has zeros at $u=(\pi/2+2\pi j)/\lambda$, $j\ge0$. Every such apparent pole is cancelled by $\Gamma(1/4+z)$ in (r(-z)), whose argument there is (-j). Therefore there are no poles in $0<u<\pi/2$, indeed none in $0<u<\pi$, for real defect rapidity and $\lambda>0$.

At $\theta=-i\pi/(2\lambda)$, $z=1/4$: the factor $\Gamma(1/4-z)$ gives a genuine simple pole and all possible cancelling factors are finite. Here $Y=-i$, so the first matrix in $M_n$ becomes $\left(\begin{smallmatrix}1&-1\\-1&1\end{smallmatrix}\right)$, of rank one. At the upper apparent pole $Y=i$, the scalar cancellation is finite. Direct limits at $\lambda=1.37$ give

```text
epsilon * Phi(-i*pi/(2*lambda)+epsilon) -> -0.8185369 - 0.2447944 i
Phi(+i*pi/(2*lambda)+epsilon)           ->  0.4092684 + 0.1223972 i
```

The algebraic pole does not always represent a physical decaying resonance. Its energy is

$$
\Delta E=M\cosh\theta_d\cos\alpha-iM\sinh\theta_d\sin\alpha,
\qquad\alpha=\pi/(2\lambda).
$$

The usual pole below the physical region with positive real energy and decay for $\theta_d>0$ is obtained in the attractive regime $\lambda>1$, where $0<\alpha<\pi/2$. The chapter does not impose this restriction. For example, at $\theta_d=0.7$, $M=1$:

```text
lambda=0.73: Delta E = -0.68889067 - 0.63411945 i  (negative real part)
lambda=0.30: Delta E =  0.62758450 + 0.65695276 i  (positive imaginary part)
```

Thus “for $\theta_d>0$, $\Im\Delta E<0$” is false without a coupling restriction. Poles farther down the rapidity plane require sheet/energy interpretation, not simply the same resonance assertion.

### Ready-to-paste corrections

Chapter 37, lines 188–192: keep the finite-factor expression but replace its definition/interpretation by:

```markdown
Here $r_N$ is the displayed product over $0\le k<N$, and $r/\bar r$ means the common-cutoff limit $\lim_{N\to\infty}r_N(\theta)/\bar r_N(\theta)$. The individual products diverge as $N^{(\lambda+1)/4}$; their paired quotient converges because the logarithm of each paired factor is $O(k^{-2})$. The bar replaces $i$ by $-i$ in the Gamma arguments and does not conjugate the complex rapidity.
```

Chapter 37, line 210 and the second source-correction bullet near line 320: replace the convergence claim with:

```markdown
The compact products in [@bcz2005, eqs. (3.42), (3.47)] do not converge without an endpoint normalization: their phase contains $-\lambda(\theta-\theta_d)\log N/\pi$. In (-@eq-sgd-scalar) it is the paired ratio $r_N/\bar r_N$, not either separate product, that converges. Keeping the endpoint factors in the Gamma duplication formula relates the two representations.
```

Chapter 37, line 233: replace “the defect without stored energy has no excited states” by “there are no poles corresponding to stable excited defects in the physical strip.” Chapter 37, line 243: replace the opening decay assertion with:

```markdown
In the attractive regime $\lambda>1$, the pole has positive real energy and, for $\theta_d>0$, $\operatorname{Im}\Delta E<0$, so it has the interpretation of a decaying resonance. The algebraic pole persists outside this regime, but that physical interpretation requires a separate analysis of its sheet and energy.
```

## 6. RSOS crossing and the proposed kink rescaling

**Verdict: unclear for the requested comparison with amplitudes “actually printed” in chapter 27; correct for the explicit CCK minimal-model formulas.** The named $\mathcal M_5$ or $\mathcal M(4,5)$ sine-Gordon kink amplitudes are not printed in chapter 27. That chapter gives a general restriction discussion and a different $a_2^{(1)}$, $W_3$ example. Chapter 3, `prp-sg-rsos`, identifies the models but does not supply their amplitudes either. Consequently the chapter-36 claim that its cited chapter-27 amplitudes were checked cannot be verified as stated. I did not silently substitute a different formula and call this a check of the book's formula.

### Independent check of the formulas the cited source does print

For [Copetti–Córdova–Komatsu](https://arxiv.org/abs/2403.04835), (18), (21), let

$$
\rho_{abcd}=\frac{d_ad_c}{d_bd_d},\qquad
W_{dc}^{ab}(\theta)=Z(\theta)
\left[\sqrt\rho\,\sinh(\theta/n)\delta_{bd}
+\sinh((i\pi-\theta)/n)\delta_{ac}\right],
$$

with $Z(\theta)=Z(i\pi-\theta)$. The literature/ordinary-crossing expression is

$$
\widehat W_{dc}^{ab}(\theta)=\rho^{i\theta/(2\pi)}W_{dc}^{ab}(\theta).
$$

Rotation $(a,b,c,d)\mapsto(b,c,d,a)$ sends $\rho\mapsto\rho^{-1}$. Direct substitution gives

$$
W_{dc}^{ab}(\theta)=\sqrt\rho\,W_{ad}^{bc}(i\pi-\theta),\qquad
\widehat W_{dc}^{ab}(\theta)=\widehat W_{ad}^{bc}(i\pi-\theta).
$$

Thus the *direction* of the relating factor is specified: ordinary equals $\rho^{i\theta/(2\pi)}$ times modified. For unequal dimensions, the modified formula does not also satisfy ordinary crossing.

The script enumerates all 14 allowed $A_4$ faces of $\mathcal M_5$, with heights (0,1/2,1,3/2), and tests both formulas at $\theta=0.73+0.12i$. It separately constructs all three-kink path matrices and tests the braid Yang–Baxter equation and the matrix part of unitarity. The common scalar cancels from crossing and YBE; its unitarity product is retained analytically.

```text
ordinary-crossing residual:                      1.14e-16
modified-crossing residual:                      1.14e-16
ordinary crossing applied to modified amplitude: 0.365
YBE residual, modified / ordinary:                5.72e-17 / 6.21e-17
unitarity residual after scalar-product removal: 4.34e-19 / 1.11e-16
```

The matrix normalization checked is

$$
W_0(\theta)W_0(-\theta)
=\sinh((\theta-i\pi)/n)\sinh((-\theta-i\pi)/n)I,
$$

whose reciprocal is $Z(\theta)Z(-\theta)$ in the source. No book check routine is used.

### Is the factor a legitimate state change? What about bootstrap?

For positive dimensions, let $f_{xy}(\theta)=(d_x/d_y)^{i\theta/(2\pi)}$. A process with in-state $K_{xy}(\theta_1)K_{yz}(\theta_2)$ and out-state $K_{xw}(\theta_2)K_{wz}(\theta_1)$ acquires the coefficient

$$
\frac{f_{xy}(\theta_1)f_{yz}(\theta_2)}
{f_{xw}(\theta_2)f_{wz}(\theta_1)}
=\left(\frac{d_xd_z}{d_yd_w}\right)^{i(\theta_1-\theta_2)/(2\pi)}.
$$

Whether this is $\rho^{+i\theta/(2\pi)}$ or its inverse depends on which opposite corners are designated outer vacua. This explicit in/out convention removes the ambiguity. For real rapidities the change is unitary. In unitarity and YBE, intermediate state phases telescope; the scripts verify both facts independently.

For a fusion $K_{xy}(\theta+iu)K_{yz}(\theta-iv)\to K_{xz}(\theta)$, the fusion coefficient must be transformed too. Its factor is

$$
\frac{f_{xy}(\theta+iu)f_{yz}(\theta-iv)}{f_{xz}(\theta)}
=\left(\frac{d_x}{d_y}\right)^{-u/(2\pi)}
\left(\frac{d_y}{d_z}\right)^{v/(2\pi)}.
$$

This is independent of the fused rapidity, so a covariantly transformed bootstrap remains valid, with changed height-dependent three-point couplings. The supplemental script verifies the cancellation at three rapidities. Keeping the old fusion couplings while changing only the two-body amplitudes is not a valid invariance claim. The particular unitary $\phi_{1,3}$ examples are in the repulsive regime, so they do not furnish a nontrivial breather-bootstrap test; the formula above establishes the general fusion covariance rather than pretending to test absent bound states.

This basis equivalence alone does not equate the amplitudes with a *fixed* categorical symmetry action. Transporting the states also transports the line operators. CCK's distinction holds the specified line action fixed; a rapidity-dependent rephasing generally makes that action rapidity-dependent. The chapter's warning that the symmetry action distinguishes the descriptions is therefore essential.

### Ready-to-paste correction

Chapter 36, line 259: replace the sentences beginning “The restricted S-matrices of @sec-rsos” through the numerical-check sentence by:

```markdown
For the $\phi_{1,3}$ kink amplitudes explicitly written in [@cck2024, eqs. (18), (21)], the ordinary-crossing expression is $(d_ad_c/d_bd_d)^{i\theta/(2\pi)}$ times the modified-crossing expression. @sec-rsos explains the restriction but does not print these sine-Gordon kink amplitudes. A rapidity-dependent kink rephasing produces this factor or its inverse, according to the in/out corner convention, and preserves unitarity and the Yang–Baxter equation. It also preserves the bootstrap provided the fusion coefficients are rephased consistently. The action of the symmetry lines must likewise be transported, so this change does not preserve a prescribed rapidity-independent line action. An independent calculation for all allowed $\mathcal M_5$ faces verifies both crossing rules and the relating factor (`referee_scratch/defects_36_37_deep/independent_checks.py`).
```

## 7. Half-soliton identities and the reading of BCZ's remark

**Verdict: correct for the identities and the distinction between spatial and temporal factors.** It should be described as a directional clarification of BCZ's time-delay comparison, not as proof that their temporal comparison is false.

### Independent classical computation

With $X=e^{-\theta}$, $\sigma=e^{-\theta_d}$,

$$
z_{\rm even}=\frac{X+\sigma}{X-\sigma},\qquad
z_{\rm odd}=\frac{X-\sigma}{X+\sigma}
=-\tanh\frac{\theta-\theta_d}{2}.
$$

Therefore

$$
z_{\rm odd}^2=\tanh^2\frac{\theta-\theta_d}{2}=A_{11},\qquad
z_{\rm even}^2=A_{11}^{-1}=A_{11}(\theta-\theta_d-i\pi).
$$

To check what the coefficient actually means, use the two-soliton tau ratio

$$
e^{i\phi/2}=\frac{1+iE_1+iE_2-A_{11}E_1E_2}
{1-iE_1-iE_2-A_{11}E_1E_2}.
$$

As $E_2\to0$, this is the soliton with exponential $E_1$. As $E_2\to\infty$, it is minus the soliton ratio with exponential $A_{11}E_1$; the minus sign is the background vacuum shift. At fixed time, crossing soliton 2 from left to right takes $E_2$ from 0 to infinity, so the **spatial** exponential factor is $A_{11}$.

On the worldline of soliton 1, however, $E_2\propto e^{m\cosh\theta_2(v_1-v_2)t}$. If $v_1<v_2$, soliton 1 is overtaken: $E_2$ goes from infinity to zero as time increases, and the **temporal before/after** factor is $A_{11}^{-1}$. If $v_1>v_2$, the temporal factor is $A_{11}$. This calculation directly explains both claims without appealing to the drafting script. The supplemental script evaluates the two tau limits numerically.

### Independent quantum computation

Put $X=\theta-\vartheta$, $a=\xi/2$, and

$$
t(x)=-i\frac{\sinh(x/2-i\pi/4)}{\sinh(x/2+i\pi/4)}.
$$

Using $t(x+i\pi)=1/t(x)$ and $\sinh(u+v)\sinh(u-v)=\sinh^2u-\sinh^2v$,

$$
t(X-i\pi-ia)t(X+i\pi+ia)
=\frac1{t(X-ia)t(X+ia)}
=\frac{\sinh X+i\cos a}{\sinh X-i\cos a}.
$$

This is the claimed $S_{sB_1}(X)$, for $\lambda>1$ where $B_1$ exists. Setting $a=0$ gives $t_{\rm odd}(X)^2=f_{1/2}(X)$, the classical statement. Independent numerical residuals at $\lambda=1.37,2.6,5.1$ are $2.22\times10^{-16},2.99\times10^{-16},4.04\times10^{-16}$.

### Source reading and suggested wording

BCZ Sec. 2 discusses a time delay produced by two identical defects and compares it to an overtaking collision. The slower soliton's temporal delay is indeed $\coth^2((\theta-\eta)/2)$. The chapter is right that this is the inverse of the left-to-right spatial factor across a soliton. Those are different comparisons, not contradictory results.

The current paragraph at chapter 37, line 154, mostly makes this distinction correctly. The source-correction bullet near line 321 would be more accurate as:

```markdown
Bowcock, Corrigan and Zambon's two-defect remark [@bcz2005, Sec. 2] compares temporal delays. A soliton overtaken by a faster soliton acquires the before-to-after factor $1/A_{11}$, agreeing with two even defects. A spatial left-to-right passage across a soliton instead has factor $A_{11}$, so the spatial factors are inverse. This is a directional clarification, not a failure of the overtaken-soliton time-delay comparison; @prp-half-soliton states the corresponding odd/even defect identities.
```

## 8. Defect Lee–Yang: range, analyticity, bootstrap, and finite size

**Verdict: wrong in the crossed-channel phase statement; correct for the principal range, bound-state structure, TBA driving term, and Lüscher sign with the qualifications below.** The unrestricted “iff” requires the periodic version already used correctly in the exercise.

### Independent parameter and analyticity calculation

Use

$$
S(z)=\frac{\sinh z+i\sqrt3/2}{\sinh z-i\sqrt3/2},\quad
t_b=\frac{i\pi(3-b)}6,\quad
T_b(\theta)=S(\theta-t_b),\quad\bar T_b(\theta)=S(\theta+t_b).
$$

The period $2\pi i$ gives $b\mapsto b+12$; $b\mapsto6-b$ sends $t_b\mapsto-t_b$, interchanging the two directions. For real $b$, choosing an orientation representative first gives ([-3,3]). The excitation gap $m\cos[\pi(b+1)/6]$ changes sign at $b=2$, after which the chosen reference state ceases to be the ground state. With ground/excited-state relabelling this gives the usual fundamental ground-state range ([-3,2]). This is a convention for real $b$, modulo direction reversal and state relabelling, not a statement about a general complex parameter or a fixed labelled oriented defect.

For Hermitian analyticity, $S(z)^*=S(-z^*)$, so

$$
T_b(-\theta^*)^*=S(\theta+t_b^*).
$$

Write $t_b=r+i\alpha$. Since the displayed rational function is injective as a function of $\sinh z$, equality to $\bar T_b(\theta)$ for all $\theta$ demands

$$
\sinh(\theta+r+i\alpha)=\sinh(\theta+r-i\alpha),
$$

or $2i\cosh(\theta+r)\sin\alpha=0$. Hence $\alpha\in\pi\mathbb Z$ and

$$
\boxed{\Re b\equiv3\pmod6}.
$$

On a period-12 representative these are the two lines $\Re b=\pm3$. Globally there are also $\Re b=9,15,\ldots$. In particular, “for real $b$ the factors are not phases” has the endpoint exceptions $b=\pm3\pmod{12}$.

### Independent excited-state and Coleman–Thun checks

Write $[x]=-\tanh(\theta/2+i\pi x/12)$. At $\theta=i\pi h/6$, this block has a pole for $h=6-x\pmod{12}$ and a zero for $h=-x\pmod{12}$. Consequently every pole count below follows directly from the displayed blocks, without copying a pole diagram.

Let $u_1=\pi(b+1)/6$, $u_2=\pi(b-1)/6$, and $v=\pi(3-b)/6$. Bootstrap gives

$$
T^{(1)}=T_bS(\theta+iu_1)=[b+1][b+3]=T_{b+2},
\quad
\bar T^{(1)}=\bar T_bS(\theta-iu_1)=[5-b][3-b].
$$

The energy identity is

$$
m\sin\frac{b\pi}6+m\cos u_1=m\sin\frac{(b+2)\pi}6.
$$

For $1<b<2$, the second ground-state pole gives

$$
T^{(2)}=T_bS(\theta+iu_2)=[b-1][b+1]^2[b+3],
\quad
\bar T^{(2)}=[3-b][5-b]^2[7-b].
$$

The apparently alternative excitation from state $1^+$ is the same state:

$$
T^{(1)}(\theta)S(\theta-iv)=T^{(2)}(\theta),\qquad
\cos u_1+\cos v=\cos u_2.
$$

These follow from the elementary Lee–Yang self-fusion identity. Independent residuals for the two routes at $b=1.6$ are $1.57\times10^{-16}$ for transmission and $1.11\times10^{-16}$ for energy. The state $2^+$ is not just another consecutive $b\mapsto b+2$ ground-state relabelling: its transmission has a double block.

For generic $b$ in the indicated ranges, the complete physical-strip pole list is:

| State and direction | Pole height $h=6\Im\theta/\pi$ | Range | Explanation |
|---|---|---|---|
| ground, right | $b+1$ | $-1<b<2$ | $1^+$ |
| ground, right | $b-1$ | $1<b<2$ | $2^+$ |
| $1^+$, right | $b+1$ | $-1<b<2$ | reverse/emission–reabsorption process through ground state |
| $1^+$, right | $b+3$ | $-1<b<0$ | Coleman–Thun diagram; internal zero reduces order two to one |
| $1^+$, left | $3-b$ | $0<b<1$ | Coleman–Thun diagram; same zero reduces order two to one |
| $1^+$, left | $3-b$ | $1<b<2$ | creation of the already identified $2^+$ |
| $2^+$, right | $b-1$ | $1<b<2$ | return process through ground state |
| $2^+$, left | $3-b$ | $1<b<2$ | return process involving $1^+$ |
| $2^+$, right | $b+1$, double | $1<b<2$ | two Coleman–Thun channels through the existing states |

The relevant internal zero for the order-reduced diagrams is $T_b(i\pi(1-b)/6)=0$. The angle (1-b) is positive only for $b<1$, explaining why the $1^+$ left pole changes interpretation at $b=1$. For the $2^+$ double pole, the internal ground-state transmissions at heights (3-b) and (b-1) are generically nonzero; the order is not reduced. At $b=1.6$, evaluating amplitudes at distances $\epsilon$ and $\epsilon/2$ from the poles gives orders (1.0000000,1.0000000,2.0000000), respectively.

These pole counts, the explicit internal zero, and the energy/transmission triangle reproduce Bajnok–Simon's closure on $0,1^+,2^+$ at the level claimed in the chapter. They do not establish uniqueness against arbitrary extra CDD factors, nor do they constitute an independent TCSA verification. Endpoint/coincident poles must be interpreted by limits.

### Independent finite-size derivation and the phase error

Define the crossed-channel insertion $D(\theta)=\bar T_b(i\pi/2-\theta)$. A one-particle mirror state has weight $D(\theta)e^{-mL\cosh\theta}$. Moving this multiplicative weight into a pseudoenergy gives

$$
\widetilde\epsilon(\theta)=mL\cosh\theta-\log D(\theta)
-\varphi*\log(1+e^{-\widetilde\epsilon}).
$$

Equivalently, with $\epsilon=\widetilde\epsilon+\log D$, the logarithm in the convolution is $\log(1+D e^{-\epsilon})$, as in BHW (95). The free energy is minus the log partition function, so its first term is

$$
\Delta E_0(L)=-m\int\frac{d\theta}{2\pi}\cosh\theta\,
D(\theta)e^{-mL\cosh\theta}.
$$

Both the negative logarithm in the driving term and the overall negative Lüscher sign in chapter 36 are correct. The expression is the finite-size contribution after the extensive and defect-energy terms have been separated; subtracting the defect-free finite-size contribution would instead replace $D$ by (D-1).

The next sentence in chapter 36, line 265, is false: BHW's real-rapidity phase condition has been changed there into a condition on $D$. At $b=-3$,

$$
D(\theta)=S(3i\pi/2-\theta)
=\frac{\cosh\theta-\sqrt3/2}{\cosh\theta+\sqrt3/2},
$$

which is positive and strictly less than one, **not a phase**. At $\theta=0.43$, $|D|=0.116259448$, whereas the real-rapidity $T_{-3},\bar T_{-3}$ are phases and obey Hermitian analyticity. The original BHW statement after (96) refers to $T_+$ on the real rapidity axis, not to its crossed-channel value.

The independently integrated Lüscher corrections, with $b=-3,m=1$, are

```text
L=6:  -4.80102270643e-05
L=10: -5.72060375327e-07
```

These provide an explicit sign check with a positive mirror weight.

### Ready-to-paste corrections

Chapter 36, line 265, replace the final sentence of the first bullet with:

```markdown
Bajnok, Holló and Watts regard the equations as reliable at least for $b=-3+i\alpha$, with real $\alpha$, where the transmission amplitudes at real rapidity are phases [@bhw2014, after eq. (96)]. This is not a phase condition on the crossed-channel insertion $\bar{\mathbb T}(i\pi/2-\theta)$; at $b=-3$ that insertion is positive and smaller than one for real $\theta$.
```

Chapter 36, line 295, replace its opening two sentences with:

```markdown
For generic real $b$ the factors are not phases. Hermitian analyticity, $\bar{\mathbb T}_b(\theta)=\mathbb T_b(-\theta^*)^*$, holds exactly when $\operatorname{Re}b\equiv3\pmod6$, or on the two representative lines $\operatorname{Re}b=\pm3$ modulo the period $12$ (@exr-ly-defect). Equivalently, $\theta_b$ is real modulo $i\pi$. The real endpoints $b=\pm3$ are therefore exceptions to the generic non-phase behaviour.
```

Chapter 36, line 293, qualify “the fundamental range” as “the fundamental ground-state range for real $b$, allowing reversal of the defect orientation.” Its excited-state and closure statements otherwise survive this check.

## Other observations — unverified

* Chapter 36, line 179, says that the doubled theory is not parity invariant even when the defect theory is. A parity operation combined with interchange of the two copies may require a more careful formulation; this was not part of the crossing test.
* The general claims equating pure transmission, conservation of asymptotic momentum, deformability of the defect, and an unperturbed categorical topological line deserve a separate review. The present calculations establish scattering identities, not all of those operator-level equivalences.
* The unrestricted extension in chapter 37 from integer topological-charge sectors to arbitrary $s$, including its physical Hilbert-space interpretation, was not audited here. The algebraic $q^{-s/2}\Phi$ normalization follows from the block calculation.
* The broader higher-spin anomaly claims in Penati–Polvara and the later-chapter claims about their classical counterparts were not independently recomputed. This report compares their stated momentum renormalization with Bajnok–Simon and does not settle the anomaly question.
