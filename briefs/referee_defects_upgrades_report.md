# Referee report: Part VIII proof upgrades

Date: 2026-10-13. Model: Codex (GPT-6).
Commit inspected: `4c78d75e56beabda1f0d2c38710e3fcb54f5419f`.
The six integration commits specified in the brief were also inspected with
`git show`; their relevant diffs are saved under
`referee_scratch/defects_upgrades/`.

This is a review of the eleven listed upgrades, not a second review of the
earlier chapters in full. No pre-existing repository file was edited, no
index-changing git command was run, and Quarto was not run. The computations
below import no book or previous-referee code. Finite-rank checks are identified
as such; they do not serve as proofs for arbitrary rank.

The principal mathematical upgrades survive this review. Corrections are needed
to the scalar normalization in chapter 42's limit statement, the unqualified
comparison with CZ2011's printed conventions, and the CDD discussion and its
attribution. There is also a sign error in a script docstring, without an error
in the executed spin-two calculation. A short all-rank proof of the second
module's conjugate-vector half is supplied in item 10.

Independent executable evidence:

| File, under `referee_scratch/defects_upgrades/` | Purpose | Recorded result |
|---|---|---|
| `algebra.py`, `algebra_output.txt` | Exact lattice propagation, Omega, displayed intertwiners, second module, cyclic inverse | PASS, $h=3,\ldots,7$; symbolic cyclic identities for $h=3,\ldots,6$ |
| `scalars.py`, `scalars_output.txt` | Symbolic linear elimination, both particle fusings, CDD examples, Gamma products | PASS; arbitrary symbolic Fourier/fusing phases |
| `spin_euler.py`, `spin_euler_output.txt` | Independent Laurent-polynomial Euler test of the $a_3$ flux | All 24 Euler derivatives vanish exactly; two controls fail |
| `doubling.py`, `doubling_output.txt` | Symbolic difference-operator module equations and arbitrary-operator obstructions | PASS; symbolic $q,q^n,q^s,y,y_d$ |

The source PDFs inspected visually include CZ2011 printed pp. 4 and 7 and
Kolb printed p. 15. The Gamma formulas were also checked against the local text
extractions. Other source checks and exact locations are specified below.

## 1. All-rank uniqueness

**Verdict: correct**, for the algebraic lattice module, nonzero spectral/module
parameters, $q-q^{-1}\ne0$, and the stated charge-conservation and
invertibility qualifications. Infinite dimension does not invalidate the
argument.

**Independent computation.** Put $n=N+s$, and use cyclic indices. For distinct
adjacent nodes $j,k=j+1$, expansion of the four coproduct products gives

\[
\begin{aligned}
\Omega^-_j={}&(f_jf_k-q^{-1}f_kf_j)\otimes k_j^{-1}k_k^{-1}
 +(q-q^{-1})f_k\otimes k_k^{-1}f_j,\\
\Omega^+_j={}&(e_je_k-qe_ke_j)\otimes1
 +(q^{-1}-q)e_k k_j\otimes e_j.
\end{aligned}
\]

The omitted terms are exactly zero on the second tensor factor, since

\[
f_jf_k=q^{-1}f_kf_j,\qquad e_je_k=qe_ke_j
\]

there. In the first formula the other mixed term cancels using
$k_j^{-1}f_k=q^{-1}f_k k_j^{-1}$. In the second it cancels using
$k_k e_j=q^{-1}e_j k_k$. Thus the signs and the ordering of the Cartan factor
in the displayed chapter formula are right.

On $V$, $f_jf_k=0$ and $f_kf_j=y^{-2}E_{j+2,j}$. The other term of
$\Omega^-_j$ is supported on $E_{j+2,j+1}$. Its image is therefore exactly
the $j+2$ line: the first coefficient is nonzero regardless of possible
cancellation in the second coefficient. On $\bar V$, the two supports are
$E_{j,j+2}$ and $E_{j+1,j+2}$, and the former is nonzero. The kernel is
exactly the coordinate hyperplane excluding $j+2$. For $\Omega^+_j$ these
two conclusions are exchanged: the kernel on $V$ excludes $j+2$, and the
image on $\bar V$ is the $j+2$ line. This works also for the cyclic edge
and for $h=3$; it is not an induction inferred from sample ranks.

Strictly, $\Omega_j$ maps one total-charge sector to another. This is no
problem: from
$E_{Q'}\Omega_{j,Q}=\Omega_{j,Q}E_Q$, the target image and source kernel
are invariant. Letting $Q$ run over the charge coset proves that every
$E_Q$ is diagonal. The nonzero diagonal entries of $1\otimes f_j$, or
$k_j\otimes e_j$, imply
$E_{Q\mp\alpha_j,ii}=E_{Q,ii}$. The single particle-changing entry then
equates the $j$ and $j+1$ diagonal values. The root shifts connect the
entire charge coset, so the scalar is global, not one scalar per sector.

Every *nonempty* sector has dimension $h$. For example its labels on $V$
are $Q-\bar e_i-s\in\Lambda_R$, where
$\bar e_i=e_i-h^{-1}\sum e_l$. If one is in the root lattice, all are,
since their differences are roots. The coset is
$s+\bar e_1+\Lambda_R$, or $s-\bar e_1+\Lambda_R$ for $\bar V$.
There is no restriction on the real offset beyond its zero sum.

Here is an independent determinant calculation. On a $V$ sector the
diagonal product is $q^{h-1}$, and the only other determinant term is the
cyclic permutation, giving

\[
\det M_V=q^{h-1}+(-1)^{h-1}X^h.
\]

For $\bar V$, define on the full lattice space

\[
A|i,n\rangle=q^{n_i}|i,n\rangle,\qquad
U|i,n\rangle=q^{n_{i+1}}|i+1,n-\alpha_i\rangle.
\]

Then $U^h=q^{-1}$, and direct multiplication of coefficients gives

\[
M_{\bar V}=\sum_{d=0}^{h-1}(XU)^d A,\qquad
(1-XU)M_{\bar V}=(1-X^h/q)A.
\]

On a sector, $\det A=q^{h-1}$ and
$\det(1-XU)=1-X^h/q$. Consequently
$\det M_{\bar V}=(q-X^h)^{h-1}$. The other determinants, by their two
cyclic determinant terms, are

\[
\begin{aligned}
\det M^{\rm rev}_{\bar V}&=q^{h-1}[1+(-1)^{h-1}qX^h],\\
\det M^+_{\bar V}&=q^{-(h-1)}[1+(-q)^{h-1}X_+^h],\\
\det M^+_V&=\frac{q^{-(h-1)}}{1-qX_+^h}.
\end{aligned}
\]

All are sector- and offset-independent. The last matrix must also be defined:
its inverse construction has a pole at $1-qX_+^h=0$. Thus “avoid determinant
zeros” is read together with the requirement that the displayed matrix exists.

The space is the algebraic direct sum of its finite charge sectors. A vector
has finite support, so the direct sum of the block inverses of $T^0$ is an
actual linear inverse. Multiplying the intertwining identity by it proves
that $(T^0)^{-1}T$ commutes with the Borel. No boundedness or Hilbert-space
inverse is required or claimed.

For real $\omega$, $q^m=1$ is equivalent to
$m(1-\omega)\in2\mathbb Z$; existence of such an integer $m>0$ is
equivalent to rational $\omega$. At $\omega=3/5$, $q=e^{2\pi i/5}$.
The shift by $5\alpha_1$ leaves all the generator coefficients unchanged:
its effect on each exponent is an integer multiple of 5. It commutes with
the module action and changes total charge. The counterexample is valid.
For non-root-of-unity $q$, differences of total weights in one coset are
integral roots, so the Cartan eigenvalues separate them, even at nonintegral
$s$.

`algebra.py` independently gives, for every $h=3,\ldots,7$:

```text
Omega image/kernel; four displayed intertwiners; second module on V and Vbar: exact PASS
```

No mathematical correction is required. For precision, chapter 41 line 72
could say “where the matrices are defined and invertible on their nonempty
total-charge sectors.”

## 2. Type-I+ vector intertwiner and reversed orientation

**Verdict: correct.** An all-rank verification of the reversed matrix can be
supplied in two lines of coefficient arithmetic; the finite-rank checks are
not essential. The chapter should include that verification to make the
proof self-contained.

**Independent computation.** Write $c=\mu/y_d$, and let $L^{\rm rev}$
denote the unflipped two-entry operator. For reversed orientation,
$f_j|n\rangle=cq^{-n_{j+1}}|n-\alpha_j\rangle$. Apply
$L^{\rm rev}\Delta(f_j)=\Delta^{\rm op}(f_j)L^{\rm rev}$ on
$\bar v_{j+1}\otimes|n\rangle$. The coefficient of
$\bar v_{j+1}\otimes|n-\alpha_j\rangle$ is

\[
X/y+cq=cq^{-1}.
\]

This is exactly $X=-(q-q^{-1})\mu y/y_d$. The other two possible outputs
on this input have, on both sides, coefficients
$y^{-1}q^{n_{j+1}}$ and $cXq^{1-n_{j+2}}$. On the other inputs the
equations are just the corresponding equalities of powers of $q$; an
input can be affected only if it is adjacent to the particle matrix entry.
The Cartan equations follow from the root shift. This verifies the formula
for arbitrary cyclic $h\ge3$, without a rank cutoff. The exact propagation
in `algebra.py` separately checks every node for $h=3,\ldots,7$.

The Chevalley identity is verified generator by generator; for example

\[
(\varpi\otimes\varpi)\Delta(e_j)
=f_j\otimes1+k_j^{-1}\otimes f_j
=\Delta^{\rm op}(f_j).
\]

Composition with $\varpi$ changes $V(y)$ to $\bar V(1/y)$.
Under $m=-n$, the positive-Borel action
$\mu y_dq^{n_{j+1}}|n+\alpha_j\rangle$ becomes
$\mu y_dq^{-m_{j+1}}|m-\alpha_j\rangle$: the reversed negative-Borel
module. Its spectral coefficient is
$-(q-q^{-1})\mu y_d/y=-X_+$. The reversed coproduct reverses the
intertwiner's direction, so one takes the inverse block, with $Q,s$ changed
to $-Q,-s$. Both reversals in the book are necessary and are correct.

Finally, set $b_n=q^{-|n|^2/2}$. Changing the defect basis to $b_n|n\rangle$
multiplies an entry with output $n-\alpha_i$ by

\[
\frac{b_n}{b_{n-\alpha_i}}=q^{-n_i+n_{i+1}+1}.
\]

Its old coefficient $Xq^{n_i-n_{i+1}}$ becomes $qX$. Thus keeping the
same $\mu$ gives $\hat x=qX$; preserving $\hat x=X$ requires
$\mu/q$. This is a basis computation, not a fitted numerical normalization.

Suggested addition at `part8-defects/41-soliton-transmission.qmd:80`, before
the Chevalley argument:

```markdown
The reversed matrix (-@eq-dt-reversed) intertwines for every $h\ge3$.
Writing $c=\mu/y_d$, its only nontrivial coefficient relation is
$X/y+cq=cq^{-1}$, obtained on $\bar v_{j+1}\otimes|N\rangle$;
the remaining coefficients are equal powers of $q$. This relation is
$X=-(q-q^{-1})\mu y/y_d$, as required.
```

## 3. Selection rule

**Verdict: correct for the relative selection rule.** The last sentence should
not be read as imposing the type-I− *absolute normalization* on the inverse
matrix of type I+.

**Independent computation.** At real $n,\theta,\theta_d$, $q$-powers have
modulus one. For $t=|X|>0$, the largest entry of $M_{\bar V}$ has modulus
$\max(1,t^{h-1})$. Hence the exact relative size of entry $j$ is

\[
\frac{t^j}{\max(1,t^{h-1})}
=\min(t^j,t^{j-h+1}).
\]

It is an equality, so the claimed upper bound holds. With
$t=e^{\omega(\theta-\theta_d)}$, every interior $j$ is exponentially
suppressed at a fixed nonzero $\theta-\theta_d$ as $\omega\to+\infty$.
At $\theta=\theta_d$ all the entries of this $M_{\bar V}$ have modulus
one. The two endpoints are the identity and $\pi(-e_i)=-e_{i-1}$.

For type I+, the inverse of the reversed cyclic matrix has a finite geometric
numerator, with entries proportional in modulus to
$1,t_+,\ldots,t_+^{h-1}$, $t_+=|X_+|$. Its common denominator is
$1-qX_+^h$. Thus the identical *relative* estimate applies with $V$ and
$\bar V$ exchanged and with the fast/slow sides reversed. At $t_+=1$,
the full-entry inverse matrix has equal moduli
$1/|1-qX_+^h|$, generally not one. This does not change the selection rule
or its scalar-independent suppression estimate.

Optional clarification at chapter 41 line 74, replacing its last sentence:

```markdown
For type I$_+$ the roles of $V$ and $\bar V$ are exchanged and the same
relative bound holds with $X_+$ in place of $X$. The inverse matrix for $V$
has a common scalar denominator, so its entries have equal, not necessarily
unit, modulus at $\theta=\theta_d$.
```

## 4. Classical linear transmission at every rank

**Verdict: correct.** The elimination in the book is genuinely symbolic in
the Fourier phase and rapidity variables. It is not a check at sampled
ranks or sampled couplings.

**Independent computation.** Let $c=e^{i\gamma}$, $p=e^\theta$,
$d=e^{\theta_d}$, and $\pi=c^2$ on one Fourier mode. Then

\[
D=\frac2{1-c^2},\quad D^T=\frac2{1-c^{-2}},\quad E=1-D,
\quad m_k=\frac{m(c-c^{-1})}{i}.
\]

The quadratic potential, differentiated before inserting waves, gives

\[
B_\phi=m[\sigma(\phi-c^{-2}\psi)+\sigma^{-1}(\phi-\psi)],\qquad
B_\psi=m[\sigma(\psi-c^2\phi)-\sigma^{-1}(\phi-\psi)].
\]

These follow also from $DD^T M^2/4=m^2$ on the mode, with the crossed
Hessian factors supplied by $D^2$ and $(D^T)^2$. Thus no rank-specific
root coordinates enter the reduction. Inserting incident/reflected waves
with $E_k=m_k(p+p^{-1})/2$, $p_k=m_k(p-p^{-1})/2$, solving the two
linear equations, and reversing $D,E$ while replacing $\sigma$ by
$\sigma^{-1}$ for type I+, gives

\[
R_-=R_+=0,\qquad
T_-={d/p+c\over d/p+c^{-1}},\qquad
T_+={p/d-c\over p/d-c^{-1}}.
\]

For incidence from the right the same calculation gives
$\bar R=0$ and $\bar T(\theta)=1/T(-\theta)$. In `scalars.py` the
two-by-two system is constructed afresh and inverted over the rational
function field. Its output is

```text
linear type=+1, left: R=0; transmission difference=0
linear type=+1, right: R=0; transmission difference=0
linear type=-1, left: R=0; transmission difference=0
linear type=-1, right: R=0; transmission difference=0
```

Here script `type=+1` denotes I−. The identities apply to every
$c=e^{i\pi k/h}$, $1\le k<h$, by specialization and meromorphic
continuation. The usual zero-momentum threshold is understood by that
continuation. The phases, poles, limits, and class shift follow algebraically
from these formulas, as stated. The book's elimination at lines 721–744 uses
the same general symbolic variables, correctly allowing $R=0$.

No correction to the upgrade is needed.

## 5. CZ2011 and CZ2009a scalar-product corrections

**Verdict: wrong as an unqualified assertion about the printed CZ2011
conventions; correct after an additional explicit spectral-variable
correction.** The Gamma endpoint, the factor-of-two error in the pole
positions, and the CZ2009a normalization at zero are correct.

**Independent computation of the products.** Set $a=\pi/(3\xi)$ and
$A_{k,l}=\tfrac12+(6k+l)a$. The equal-index CZ2011 product, normalized at
zero, is

\[
H(z)=\prod_{k\ge0}\frac{F_k(z)}{F_k(0)},\qquad
F_k(z)=\prod_{l=1,2}\Gamma(A_{k,l}+z)\Gamma(A_{k,l}-z)
\bigg/\prod_{l=4,5}\Gamma(A_{k,l}+z)\Gamma(A_{k,l}-z).
\]

Reading the printed formula from the PDF, its finite Gamma product obeys
exactly

\[
\frac{P_K(z)}{P_K(0)}
=\prod_{k=0}^{K-1}\frac{F_k(z)}{F_k(0)}
\prod_{l=1,2}\frac{\Gamma(A_{K,l}-z)}{\Gamma(A_{K,l})}.
\]

The endpoint behaves as $(6aK)^{-2z}$; the normalized equal-index product
has a summable $O(K^{-2})$ factor. Saying simply that the *unnormalized*
printed product drifts as $K^{-2z}$ suppresses an additional
$z$-independent normalization. The rapidity-dependent drift is as claimed.

For a hand check of the functional equations, write
$G(z)=\Gamma(\tfrac12+z)\Gamma(\tfrac12-z)=\pi/\cos\pi z$.
In $H(z-a)H(z+a)/H(z)$, all interior Gamma arguments cancel: the sequence
of multiplicities $+1,+1,0,-1,-1,0$, repeated in $l$, satisfies
$m_{l-1}+m_{l+1}-m_l=0$. Only the $l=0$ factor and a
constant endpoint remain. Normalizing at $z=0$ fixes the result:

\[
H(z-a)H(z+a)=B G(z)H(z),\qquad B=H(a)^2/\pi.
\]

Applying this at $z+a,z+2a$ and cancelling gives the second identity

\[
H(z)H(z+3a)=B^2G(z+a)G(z+2a).
\]

Consequently a *single* normalization

\[
\rho(z)=C e^{i\pi z}H(z),\qquad C=[2H(a)^2]^{-1},
\]

satisfies both (3.10) and (3.8), with

\[
\tau^2x=e^{-2\pi iz},\qquad q=i e^{i\pi a}.
\]

For example, $(1+e^{-2\pi iz})e^{i\pi z}=2\cos\pi z$ proves the
bootstrap, and the crossing polynomial becomes
$4e^{-i\pi(2z+3a)}\cos\pi(z+a)\cos\pi(z+2a)$.
Thus the two constants are compatible; checking only variation with $z$
would not have established that point.

`scalars.py` checks the *fully normalized* equations, not just their ratios
at different rapidities. It obtains $C=0.0241346292$ at $\xi=0.93$ and
$C=0.121702325$ at $\xi=1.73$; maximum crossing/bootstrap errors are
$1.6\times10^{-8}$ with double-precision log-Gamma summation. The larger
cutoff reaches the cancellation-error floor; the exact argument above is
the proof.

The outside Gamma poles are $z=\tfrac12+a$ and
$z=\tfrac12+2a$, hence

\[
\theta=\nu-i(\xi/2+\pi/3),\qquad
\theta=\nu-i(\xi/2+2\pi/3).
\]

At the corresponding *printed* positions the relevant argument is
$-1/2$, not a Gamma pole. With $\nu=\eta+i\pi/2$, the imaginary
rapidities are $\pi/6-\xi/2$ and $-\pi/6-\xi/2$. Requiring them
in $(-\pi/2,0)$ gives the doubled windows printed in the book.

**The missing qualification.** CZ2011 (3.2), printed p. 4, actually says
$x_i=q^{2\pi\theta_i/\xi}$, not $e^{2\pi\theta_i/\xi}$.
This is visible in the PDF, not an extraction error. With that printed
definition and $\tau=e^{-\pi\nu/\xi}$, the identity
$\tau^2x=e^{-2\pi iz}$ used by the chapter does not follow. Its logarithmic
derivative would be $(2\pi/\xi)\log q$, whereas the required one is
$2\pi/\xi$; no logarithm of a unit-modulus $q$ equals 1. The functional
equation verification therefore also repairs (3.2). The pole correction is
still justified in the repaired convention, but “so the slip is in the pole
positions” is incomplete unless this is disclosed.
Source: [CZ2011, equations (3.2), (3.8), (3.10), (3.12)](https://arxiv.org/abs/1012.4186).

**Paste-ready correction**, `appendices/D-corrections.qmd:153`: replace the
functional-equation sentence with

```markdown
For this comparison we read the spectral variable in (3.2) as
$x_i=e^{2\pi\theta_i/\xi}$; the PDF prints $q^{2\pi\theta_i/\xi}$.
Then $\tau^2x=e^{-2\pi iz}$. After normalization at $z=0$, removal of the
endpoint Gamma factors, and fixing one overall constant, the equal-$k$
product satisfies both (3.8) and (3.10). In this convention the displayed
Gamma arguments require the corrected pole positions above.
```

Make the same qualification at
`part8-defects/41-soliton-transmission.qmd:189`, replacing “With $z=\cdots$,
so that $\cdots$” by

```markdown
For this check we use $x=e^{2\pi\theta/\xi}$, correcting the printed
$q^{2\pi\theta/\xi}$ in (3.2), so that
$z=i(\theta-\nu)/\xi$ gives $\tau^2x=e^{-2\pi iz}$.
```

For CZ2009a, put $A_k=\tfrac12+(hk+(h-1)/2)\omega$. The printed
numerators are $\Gamma(A_k-z)\Gamma(A_{k-1}+z)$, the two denominators
have $A_{k-1}+\omega\pm z$, and the outside factor is
$\Gamma(A_0-z)$. Reindexing immediately gives

\[
\frac{P_K(z)}{P_K(0)}
=\prod_{k<K}\frac{F_k(z)}{F_k(0)}
\frac{\Gamma(A_K-z)}{\Gamma(A_K)}.
\]

The script compares the exact multisets of Gamma arguments for independent
choices of $h,K$, with zero discrepancy. The $z=0$ qualification now in
Appendix D is right. This part needs no correction.

## 6. Universal-R proof of the all-rank triangle relation

**Verdict: correct.** This proof does not need hypothesis (ii) of
`thm-t-triangle`. Its use of a bilateral Borel module is legitimate as a
formal principal-degree construction. Pinpoint references and a sentence
about the Cartan lift would improve it.

**Independent first-order computation.** Put $\kappa=q-q^{-1}$,
$n=N+s$. For the coproduct actually used in Part VI, take

\[
\mathcal R=q^t(1+\kappa\sum f_j\otimes e_j+\cdots),\qquad
\mathcal R_{21}^{-1}=(1-\kappa\sum e_j\otimes f_j+\cdots)q^{-t}.
\]

The following table is obtained by acting first with the rightmost factor,
including the change in the defect weight before the final Cartan factor.
The coefficient in the last column includes $X$ or $X_+$.

| Type and particle | Order zero | Order-one particle output and coefficient |
|---|---|---|
| I−, $V$ | $q^{-n_i}$ | $i\to i-1: X$ |
| I−, $\bar V$ | $q^{n_i}$ | $i\to i+1: Xq^{n_i+n_{i+1}}$ |
| I+, $\bar V$ | $q^{-n_i}$ | $i\to i-1: X_+q^{n_i-n_{i-1}-1}$ |
| I+, $V$ | $q^{n_i}$ | $i\to i+1: X_+q^{2n_{i+1}-1}$ |

Here $X=-\kappa\mu y/y_d$ and $X_+=\kappa\mu y_d/y$.
For example the first row is
$-\kappa(y/y_d)\mu q^{n_i}q^{-n_i}=X$. In the third row the raising
operator contributes $\mu(y_d/y)q^{n_i}$; the final Cartan factor is
$q^{-n_{i-1}-1}$. The fourth row also follows by expanding the inverse
from item 2, $(A+XB)^{-1}=A^{-1}-XA^{-1}BA^{-1}+O(X^2)$, with
$n\mapsto-n, X\mapsto-X_+$. It agrees, including the factor $q^{-1}$.
All four rows therefore reproduce the displayed matrices, not just their
gauge invariants. This calculation is valid for arbitrary $h\ge3$.

**Formal evaluation and the third leg.** In principal degree $m$, there
are finitely many nonnegative affine-root combinations of height $m$, and
the corresponding nilpotent homogeneous spaces are finite-dimensional
(they are quotients of the finite set of words of length $m$). Each term
acts by a finite number of lattice shifts. Local nilpotence on $\mathcal W$
is unnecessary. A series is interpreted coefficientwise on each finite
charge sector, with finite-support vectors extended to the Laurent-series
field. It is not asserted to converge analytically or in operator norm.

For I+, all third-leg factors of the evaluated universal Yang–Baxter
identity are in $U_q(\mathfrak b^+)$, including its Cartan algebra.
For I− they are in $U_q(\mathfrak b^-)$. The first two modules carry the
whole algebra, so both intermediate tensor products and both sides of the
identity are defined. The same observation applies to the intertwining
identity restricted to the appropriate Borel. No extension of
$\mathcal W$ to the other Borel is being assumed.

For I+, set $u=y_2/y_1$, $v=y_d/y_2$. A product of terms of degrees
$a,b,c\ge0$ from the three factors contributes $u^{a+b}v^{b+c}$.
At $u^mv^n$, $0\le b\le\min(m,n)$, so only finitely many terms
contribute. For I− the same count works with $u=y_1/y_2$, $v=y_2/y_d$.
This verifies the formal multiplication explicitly.

The proof in item 1 uses identities and divisions by nonzero elements only,
so it works unchanged over $\mathbb C((x))$. For the finite-dimensional
bulk pair, extending the rational function field to its Laurent-series
field preserves the rank of the intertwining equations: the same nonzero
minor stays nonzero, and no additional minor becomes nonzero. Consequently
the evaluated matrices are scalar series times the rational matrices in the
book. The four leading terms above give $\phi_a=1+O(x^2)$.

The product $\psi_{ab}(u)\phi_a(uv)\phi_b(v)$ is well-defined in
$\mathbb C((u))((v))$, has a nonzero lowest coefficient, and is invertible
there. It is the same on both sides. Cancelling it and using the injective
expansion of rational functions into that field proves the rational
triangle identity. Both choices of universal R give the same bulk
intertwiner up to a scalar by its generic uniqueness.

**Cartan and specialization.** The definition in chapter 7 uses the derived
algebra generated by $e_i,f_i,k_i^{\pm1}$, adjoining the derivation only
for the universal R. One may first affinize the modules, let the degree
operator act on the spectral variables, and then impose level zero. The
factor $q^{c\otimes d+d\otimes c}$ becomes 1. It does not leave an
unaccounted spectral shift. With a fixed logarithm of $q$, the explicit
weights $n=N+s$ give a lift of the torus action, and $q^t$ acts as
$q^{\lambda\cdot n}$. Integral $s$ is not required; all commutation
checks involve integral weight *differences*. This is the extra Cartan
information already fixed by the chapter's $q^{s_i}$ convention.

Finally regard $X,q,q^{s_i}$ as the algebraic parameters (with the
zero-sum constraint). Every matrix entry is rational in these variables.
An identity proved generically continues to roots of unity wherever its
displayed denominators are defined. This is a specialization of the
rational identity, not an evaluation of a singular universal-R product at
a root of unity.

**Sources.** Tolstoy–Khoroshkin 1992, §4, Theorem 2, printed Russian p. 88,
is the relevant affine factorization theorem; its §2 specifies the
coproduct and the added Cartan degree generator.
[Original paper](https://www.mathnet.ru/eng/faa771).
The 1991 paper concerns finite-dimensional contragredient algebras and is
not, by itself, the affine theorem.
[Author-uploaded primary paper](https://www.researchgate.net/publication/225799830_UniversalR-matrix_for_quantized_superalgebras).
Frenkel–Hernandez §5.1, especially (5.16) and the following L-operator
paragraph, and §7.1 give the completed Borel factorization and YBE, in the
opposite Borel ordering from the book's chosen $\mathcal R$.
Their transfer-matrix category hypotheses should not be applied to this
bilateral module: the finite principal-degree argument above is what
justifies the present evaluation.
[Frenkel–Hernandez](https://arxiv.org/html/1308.3444v4#S5.SS1).

Suggested source clarification at chapter 41 lines 88–90:

```markdown
Use the affine factorization of [@tolstoy_khoroshkin1992, Sec. 4, Thm. 2],
translated to our coproduct and Borel ordering. For completed Borel
L-operators and the universal Yang--Baxter identity see also
[@frenkel_hernandez2015, Secs. 5.1, 7.1]. Our bilateral module need not lie
in their category $\mathcal O$: finiteness in each principal degree gives
the coefficientwise evaluation used here. We fix a logarithm of $q$ and
the weight lift $N+s$ to define the Cartan factor.
```

## 7. Exact $a_3$ spin-two flux

**Verdict: correct.** The theorem upgrade for $a_3$, both types and both
chiralities, is justified. The executed current has the right relative
normalization. The `boundary_G` docstring has the wrong sign for the
Theta term and omits the chirality factor relating its result to the
derivative of the usual bulk charge.

**Independent derivation of the flux.** Write
$d_\epsilon=\partial_t+\epsilon\partial_x$,
$u=d_\epsilon\phi$, $w=d_\epsilon^2\phi$. For
$\phi_{tt}-\phi_{xx}=-\nabla V$, and the defect condition
$\phi_x=F_\phi$,

\[
w=2\phi_{tt}+2\epsilon D_tF_\phi+\nabla V(\phi).
\]

The same formula holds for $\psi$. This is precisely the expression
executed by `boundary_G`; its second derivative substitution is correct.

Use

\[
T=A(u,u,u)+w^TBu,\qquad
\Theta=\sum_\alpha(-B\alpha)\cdot u\,e^{\alpha\cdot\phi},
\qquad B=3(1-D).
\]

Differentiating on shell gives

\[
d_{-\epsilon}T-d_\epsilon\Theta
=-\sum_\alpha e^{\alpha\cdot\phi}
\{3A(u,u,\alpha)+2(\alpha^TBu)(\alpha\cdot u)\}.
\]

For $A=2\sum\mu_k^{\otimes3}$, this vanishes with $B=3(1-D)$.
The current is exactly three times the BCZ normalization
$T=\tfrac13A(u,u,u)+w^T(1-D)u$. A common overall factor of three is
harmless. Keeping the cubic coefficient used by the script but setting
$B=1-D$ would fail; the independent calculation explicitly tests that
control.

The conservation law implies

\[
\partial_t(T-\Theta)=\epsilon\partial_x(T+\Theta).
\]

Thus, for the sum of the left and right bulk integrals of $T-\Theta$,

\[
\dot Q_{\rm bulk}=\epsilon G,\qquad
G=(T+\Theta)_\phi-(T+\Theta)_\psi.
\]

This is the $G$ actually returned by the code. If the code constructs
$D_tF=G$, the conserved modified charge is $Q_{\rm bulk}-\epsilon F$.
The spin is $-2\epsilon$, consistently with the book's boost convention.
The existence of a local contribution is unaffected by this conventional
sign, but the docstring's $T-\Theta$ expression is not the flux.

**Independent exact computation.** `spin_euler.py` uses a rational
orthonormal Hadamard basis of the zero-sum subspace of $\mathbb R^4$.
It constructs $D=2(1-\pi)^{-1}$ directly, then uses independent formal
variables $z_i=e^{q_i/2}$ for the six defect field coordinates. All
exponentials become Laurent monomials. It tests

\[
\mathcal E_i(G)=\frac{\partial G}{\partial q_i}
-D_t\frac{\partial G}{\partial\dot q_i}
+D_t^2\frac{\partial G}{\partial\ddot q_i}=0
\]

without constructing a homotopy primitive. The output is

```text
bulk coefficient identities: 4 exact zeros; B=E control fails
type=+1, epsilon=-1, spin=+2: six Euler derivatives = 0; flux has 156 Laurent monomials
type=+1, epsilon=+1, spin=-2: six Euler derivatives = 0; flux has 156 Laurent monomials
type=-1, epsilon=-1, spin=+2: six Euler derivatives = 0; flux has 156 Laurent monomials
type=-1, epsilon=+1, spin=-2: six Euler derivatives = 0; flux has 156 Laurent monomials
nonintegrable D: Euler derivative 0 is nonzero; value at rational point = -1745615841528721/30240000000
PASS: independent exact a3 flux test, both types/chiralities, symbolic nonzero sigma.
```

The nonintegrable perturbation preserves $D+D^T=2$, so the failing
control is not merely a failure of antisymmetry. The exact Laurent
identities extend to nonzero complex $\sigma$, and by analytic
continuation to the imaginary-coupling convention with the appropriate
field rescaling. The local Euler criterion, on these field-coordinate
domains, proves local total-derivative exactness. This is independent
evidence that the homotopy construction has not manufactured a zero by
definition: an arbitrary flux fails the Euler test, as the control does.

**Paste-ready docstring correction** at
`code/defect_spin_two_a3.py:175`:

```python
    """G = (T + Theta)(phi) - (T + Theta)(psi), after eliminating x derivatives.
    For Q_bulk = integral(T - Theta), dQ_bulk/dt = eps * G.
    A primitive Dt(F) = G gives the defect contribution -eps * F.
    The charge spin is -eps * spin, with spin=2 or spin=1.
    """
```

No downgrade of `prp-spin-two` is warranted. The all-rank hierarchy remains
a separate sketch, as labelled.

## 8. Chapter 42: second module and diagonal obstructions

**Verdict: wrong as literally normalized in item (5); correct for the
claimed Borel symmetries and obstructions after explicitly removing the
scalar transmission factor.** The all-operator arguments do not depend on
the finite-window ansatz.

**Independent arbitrary-operator computation.** Let
$O=\left(\begin{smallmatrix}0&B\\C&0\end{smallmatrix}\right)$,
$D_0=\operatorname{diag}(A,D)$. Keep the torus scalar $t\ne0$ explicit.
The following uncancellable blocks follow by multiplying two-by-two
matrices over an arbitrary noncommutative operator algebra:

| Type and failed equation | Power | Residual block independent of unknowns |
|---|---|---|
| I−, untwisted $e_1$ | $y^2$ | lower-right $C$ |
| I−, swapped $f_1$ | $y^{-1}$ | upper-right $-t^{-1}D$ |
| I+, untwisted $f_1$ | $y^{-2}$ | lower-right $-t^{-1}B$ |
| I+, swapped $e_1$ | $y$ | upper-right $A$ |

For example the first two full residuals are

\[
\begin{pmatrix}-t\kappa C&0\\0&C\end{pmatrix},\qquad
\begin{pmatrix}0&-t^{-1}D\\D\kappa&0\end{pmatrix}.
\]

They cannot vanish since $A,B,C,D$ are invertible. Nothing restricts the
unknown $Y,\kappa$ to one lattice shift in this calculation.
`doubling.py` verifies all four formulas with noncommuting symbols.

The algebra defined in chapter 7, before adjoining the derivation to discuss
the universal R, is precisely the derived algebra with generators
$e_i,f_i,k_i^{\pm1}$. Level zero is a property of its modules, not a
replacement of this algebra by an unstated completion. Kolb's §3.1 defines
that same $U_q(\mathfrak g')$, and §3.2, Theorem 3.2 states the torus
semidirect diagram-automorphism classification, attributing it to
Twietmeyer Theorem 2.1. It is over $K(q)$. Accordingly the chapter's
generic-$q$ use is justified; it should not be promoted to a classification
of Hopf automorphisms of a specialized root-of-unity algebra or an enlarged
algebra containing the derivation.
[Kolb, §§3.1–3.2](https://arxiv.org/abs/1207.6036).
For $A_1^{(1)}$ the only diagram choices are identity and node exchange,
so these four failures cover every twist in the stated setting.

**The second module.** With $z=q^n$, $s_0=q^s$, the new module operators
are weighted shifts by $+1,-1,+2,-2$. Define the input-weighted shift by
$A_{a,k}|n\rangle=a(q^n)|n+k\rangle$. Their compositions are computed
exactly in `doubling.py` using

\[
(A_{a,k}A_{b,l})|n\rangle
=a(zq^l)b(z)|n+k+l\rangle,\qquad z=q^n.
\]

The output checks the full $e_1,e_0$ equations, both Cartan equations,
both inverse relations, and all four Cartan conjugation relations, with
symbolic $q,z,s_0,y,y_d$. In particular,

\[
\widetilde e_1\widetilde e_0
=\widetilde e_0\widetilde e_1
=-\frac{y_d^2q^2}{(q^2-1)^2}.
\]

The two cubic Serre expressions are therefore a common product multiplied
by $1-[3]_q+[3]_q-1=0$. This completes every generator and relation
needed for (2), not just the $e_1$ recurrence. A homogeneous constant in
that recurrence violates the own-node Cartan relation by a factor
$1-q^4$, and is excluded at generic $q$.

**Combining the halves.** For type I+, conjugation by $\Delta(k_1)$ gives
weights $q^{-2}$ and $q^2$ on the two terms of
$f_1\otimes k_0^{-1}+1\otimes f_0$. If $q^4\ne1$, they can be
separated by a linear combination. The algebra then contains $1\otimes f_0$.
Its intertwining equation has the $y^{-2}$ term
$-E_+O_+=\operatorname{diag}(-C,0)$, independent of its proposed module
action. It is nonzero. The claimed type-I+ proof is complete.

**Limits and their inverse-Cartan obstruction.** For the diagonal matrix
limit and torus twist, the coefficients of $y$ and $y^{-1}$ give

\[
\kappa=t^{-1}AD^{-1},\qquad
\kappa^{-1}=t^{-1}D^{-1}A.
\]

Both are $t^{-1}q^{-2n-s}$. Their product is
$t^{-2}q^{-4n-2s}$, and cannot be 1 for every $n$ at generic $q$.
The swapped twist already fails the particle-support test. For the
off-diagonal limit the swapped equations instead give

\[
\kappa=t^{-1}CB^{-1},\qquad
\kappa^{-1}=t^{-1}B^{-1}C.
\]

Both shift $n\to n+2$; their product shifts by 4 and cannot be the
identity. The unswapped twist fails the particle-support test. This proves
the statement for all Hopf twists classified above. Type I+ exchanges the
two limits and yields the same contradictions (with the corresponding
diagonal powers reversed). These are algebraic limits, sectorwise or on
finite-support vectors, not an operator-norm convergence claim for an
arbitrary nonunitary $q$.

**The normalization error.** Chapter 37 defines the *full* matrix as
$\mathbb T=\Phi M$. Therefore $P\mathbb T\ne D_0+yO$; that polynomial
is $PM$. For real parameters,

\[
|\Phi|=(1+Y^2)^{-1/2},\qquad Y=y/y_d.
\]

The converting entries of the full matrix divided by $Y$ have modulus
$(1+Y^2)^{-1/2}\to0$, while its diagonal entries divided by $Y$ also
tend to zero. Thus the limit claimed in (5) for $\mathbb L/Y$, with
$\mathbb L=P\mathbb T$ as defined, is the zero matrix, not the displayed
nonzero shifted conjugation. Removing the scalar fixes the statement and
does not affect any intertwining obstruction.

**Paste-ready corrections.** At
`part8-defects/42-defects-and-boundaries.qmd:91`, start item (5) with

```markdown
5. These are limits of the matrix part: write $\mathbb T=\Phi M$ and
$L_M=PM$. For type I$_-$, as $\theta_d\to+\infty$, $L_M$ tends to
the direct sum over $n$ of ${\rm diag}(q^{-n},q^{n+s})$.
```

Use $L_M$, rather than $\mathbb L=P\mathbb T$, in the subsequent
division by $y/y_d$. The rest of item (5) is unchanged. At line 95 replace
the opening with

```markdown
Remove the common scalar factor by writing $\mathbb T=\Phi M$, and set
$L_M=PM$. The scalar cancels from each intertwining equation. In this
proof write $L_M=\mathcal D+y\mathcal O$ for type I$_-$ and
$L_M=\mathcal D_++y^{-1}\mathcal O_+$ for type I$_+$.
```

Use this matrix-part notation consistently in the remainder of the proof.
The theorem labels for the algebraic claims need not be downgraded.

## 9. All-rank particle transmission identities

**Verdict: correct**, with the real-coupling interacting range $0<B<2$
for the assertion about actual bulk S-matrix poles, as its proof states.

**Independent computation.** Let
$F_-(x,u)=(x+u)/(x+u^{-1})$ and
$F_+(x,u)=(1-xu)/(1-x/u)$. The latter expresses the type-I+ factor in
the *same* $x=e^{-(\theta-\theta_d)}$ variable, eliminating a possible
confusion about the sign of the rapidity shifts. For either sign,

\[
\begin{aligned}
F_\pm(x,u)F_\pm(-x,-u^{-1})&=1,\\
F_\pm(x/v,u)F_\pm(xu,v)&=F_\pm(x,uv),\\
F_\pm(-xv,u)F_\pm(-x/u,v)&=F_\pm(x,-uv).
\end{aligned}
\]

`scalars.py` reduces all six differences to exact zero for symbolic
$x,u,v$. In particular the second fusing uses $-xv,-x/u$, not the
first fusing's shifts. An explicit type-I+ cancellation for that case is

\[
\frac{1+xuv}{1+xv/u}\,
\frac{1+xv/u}{1+x/(uv)}
=F_+(x,-uv).
\]

For the first fusing the cancelling factors are $1-xu/v$. This supplies
the details omitted from the chapter's outline of the positive type.

For real $x>0$, $|u|=1$, numerator and denominator are conjugates in
both formulas (using the real variable $1/x$ for the second). The crossed
factor is the conjugate of the opposite-rapidity factor, giving Hermitian
analyticity. Their poles are respectively
$-i(\pi-\gamma_a)$ and $-i\gamma_a$ modulo $2\pi i$, so neither
has a pole in $0<\operatorname{Im}(\theta-\theta_d)<\pi$.
The right-incidence statements follow from crossing. Item (5) follows from
the independent arbitrary-phase elimination in item 4 of this report.

The proof is therefore all-rank, not a test of $h\le8$. No correction to
the upgrade is needed. At the free endpoints the usual statement about
bulk *poles* must of course be replaced by its limiting interpretation;
the rational transmission identities themselves remain valid.

## 10. The explicit second module in chapter 43

**Verdict: correct.** The vector half has an all-rank proof. The
conjugate-vector half also admits the short algebraic proof below, so its
numerical-only qualification is conservative.

**Independent vector and module computation.** Put
$c=t/(q-q^{-1})$, $t=e^{\omega(\theta_d+i\pi/2)}$, and
$\delta_j=\alpha_j-\alpha_{j-1}$. Separating the vector equation by
powers of $y$ gives the following exhaustive conditions:

\[
\kappa_N=1,\quad
Y_{N-\alpha_j}-q^{-1}Y_N=tq^{-n_j},\quad
qY_N-Y_{N-\alpha_{j-1}}=tq^{-n_j},
\]

plus invariance under $-\alpha_{i-1}$ for $i\ne j,j+1$. Subtracting
the particular solution $cq^{-n_j}$, the two remaining recurrences and
those invariances say that $q^{-n_j}$ times the homogeneous solution is
constant on every root-lattice direction. Thus the complete solution in
this charge grade is

\[
Y_N=cq^{-n_j}+c_jq^{n_j}.
\]

The Cartan relation with $\widetilde k_j=S_{\delta_j}$ changes the two
terms by $q^2,q^{-2}$, respectively. It removes $c_j$ when $q^4\ne1$.
For different nodes,
$(\delta_i)_j=2\delta_{ij}-\delta_{i+1,j}-\delta_{i-1,j}=a_{ij}$,
where the two neighbours coincide for $h=2$. Hence

\[
\widetilde k_i\widetilde e_j\widetilde k_i^{-1}
=q^{a_{ij}}\widetilde e_j,\qquad
\prod_i\widetilde k_i=1.
\]

The adjacent $\widetilde e_j$'s $q^{-1}$-commute and the others
commute, yielding each quadratic Serre relation by the polynomial
$(c_0-q)(c_0-q^{-1})$. At $h=2$ the product of the two $e$'s is
central, and the cubic relation vanishes as in item 8. The Cartan-element
intertwining equation follows directly by tracking the two shifts. Thus
both the module relations and its compatibility are established for every
rank, not merely postulated from the recurrence.

**A short proof for $\bar V$.** The cyclic identity in item 1 gives an
especially simple sparse inverse:

\[
M_{\bar V}^{-1}=\frac{H}{1-X^h/q},\qquad
H|i,n\rangle=q^{-n_i}|i,n\rangle
-\frac Xq|i+1,n-\alpha_i\rangle.
\]

It suffices to verify

\[
H(Y+yE_{j,j-1}\kappa)
=(yE_{j+1,j}+\pi_{\bar V}(k_j)Y)H.
\]

This is the original equation multiplied by the inverse, with its sides
reversed; its scalar denominator cancels. Write $b=-1/(qt)$, so the
off-diagonal part of $H$ is $yb$ times the plain forward shift.
At order $y^0$, both coefficients are
$cq^{-n_j-n_i-(\alpha_j)_i}$. At order $y^2$, only the input $i=j-1$
contributes; both sides give coefficient $b$, output $j+1$, and label
$n-\alpha_{j-1}$.

At order $y$, the difference between the two terms involving $Y$ is

\[
bcq^{-n_j}{1-q^{(\alpha_i)_j-(\alpha_j)_{i+1}}}.
\]

It vanishes unless $i=j$ or $i=j-1$. In those cases the exponent is
$+2$ or $-2$, and

\[
bc(1-q^2)=1,\qquad bc(1-q^{-2})=-q^{-2}.
\]

These are precisely the two remaining diagonal/particle terms. This proves
the equation for every $h$. The Cartan equation follows by the same
diagonal-weight count, or directly from
$\delta_j=\alpha_j-\alpha_{j-1}$. Multiplying back proves the result
away from the inverse's pole; polynomial continuation gives the matrix
identity everywhere it is defined. No fusion embedding or numerical
bootstrap assumption is needed. Independently, `algebra.py` checks the
full equations on every node of $V,\bar V$ for $h=3,\ldots,7$, exactly.

**Rank-one reduction.** For
$d_n=\tau^{-n}q^{n(n-1)/2+ns/2}$, $\tau^2=q^{-1}$,

\[
\frac{d_n}{d_{n+1}}=\tau q^{-n-s/2},\qquad
\frac{d_n}{d_{n+2}}=q^{-2n-2-s}.
\]

The first converts $cq^{-n-s/2}S_{+1}$ to
$iy_dq^{-2n-s}S_{+1}/(q^2-1)$; the second converts the plain $S_{+2}$
to the stated $\widetilde k_1$. For the negative shifts the ratios give
$iy_dq^{2n+s}S_{-1}/(q^2-1)$ and
$q^{2n-2+s}S_{-2}$. Multiplying the vector matrix by $q^{s/2}$ gives
converting coefficients $-i(y/y_d)q^{-n}$ and
$-i(y/y_d)q^{n+s}$, exactly chapter 37. The basis change and offset are
right, including the factor of $q^{-1}$.

**Failed twists and combined algebra.** The vector matrix has an
off-diagonal part $X$ times $i\mapsto i-1$. Thus its $y^2$ composition
with $e_j$, on input $j+1$, exits at $j-1$. The corresponding right
side vanishes unless the twist sends $j\mapsto j-1$. Similarly, the
$y^{-1}$ term from $f_{\psi(j)}$ on the right acts on input $\psi(j)$,
where the left $f_j$ is zero unless $\psi(j)=j$. These support
obstructions allow arbitrary $Y,\kappa$ and establish all the listed
failures for every $r\ge2$.

For $r\ge2$, conjugating the two summands of
$e_j\otimes1+k_j\otimes e_{j-1}$ by $\Delta(k_j)$ gives weights
$q^2,q^{-1}$. They separate when $q^3\ne1$. The resulting
$e_j\otimes1$ has an uncancellable $y^2$ term, proving the coideal
obstruction. The $h=2$ version instead uses $q^4\ne1$, as chapter 42
correctly does.

There is no required correction. The sparse-inverse argument above is
ready to insert after `eq-dpt-second-module` at chapter 43 line 229 if an
all-rank theorem label is desired. That upgrade should explicitly retain
generic $q$, the fixed charge lattice, and nonzero finite parameters.

## 11. Chapter 43 claims inherited from the handoff

**Verdict: wrong in the blanket CDD comparison and insufficiently qualified
in the novelty claim; the residue, imaginary-parameter exceptions,
Lee–Yang direction, TBA algebra, and vacuum estimate are correct.** Each
component is checked separately below.

### Residue rank and scalar pole

On a $\bar V$ sector write $n^{(i)}=Q+\bar e_i$, and set
$z_i=q^{Q_i-1/h}$, so $\prod_i z_i=q^{-1}$. The row $k$, column $i$
coefficient is

\[
(M_{\bar V})_{ki}=qX^d\prod_{l=i}^{i+d}z_l,\qquad d=k-i\pmod h,quad0\le d<h.
\]

At $X^h=q$, define, using indices $0,\ldots,h-1$,

\[
u_k=X^k\prod_{l=0}^kz_l,\qquad
w_i=X^{-i}\left(\prod_{l=0}^{i-1}z_l\right)^{-1}.
\]

For a path that wraps, the extra factor is
$X^h\prod_lz_l=1$. Therefore the equality
$M_{\bar V}=q\,u w^T$ holds for *every* entry. Both vectors have nonzero
entries, so its rank is exactly 1 on every sector, for arbitrary offset.
This is stronger than inferring rank from its determinant.

At the species-$r$ pole, the constituent Gamma variable in the $m$-th
factor of the scalar bootstrap is $z_m=a_0-m\omega$,
$0\le m\le r-1$, where $a_0=\tfrac12+r\omega/2$. The two numerator
arguments in factor $k=0$ are

\[
m\omega,\qquad 1+(r-m)\omega,
\]

and the denominator arguments are
$(m+1)\omega, 1+(r+1-m)\omega$. For $\omega>0$, only the first,
with $m=0$, is zero; every other argument is positive. The remaining
convergent product is finite and nonzero. Thus $\bar g$ has one
uncancelled simple pole. The nonzero diagonal entries of the rank-one
matrix prevent its scalar residue from being cancelled as an operator.
The statement is correct *per charge sector*; on the whole infinite
direct sum the residue does not have finite rank.

### Soliton CDD condition and the even-$h$ exception

Let

\[
\bar f(D)=\prod_{m=0}^{h-2}
f\bigl(D+i\pi(2m-h+2)/h\bigr).
\]

The second crossing equation is exactly
$f(D-i\pi)\bar f(D)=1$, which is the displayed $h$-string condition.
The first crossing equation, $f(D)\bar f(D-i\pi)=1$, is the same
condition with its central argument shifted by
$-i\pi+2\pi i/h$. Hence the condition is necessary and sufficient for
preserving the two scalar crossing equations *given the bootstrap
definition of $\bar f$*. It is not by itself a minimality, analyticity,
or physical-spectrum condition. It also implies $2\pi i$-periodicity:
divide the condition at $D+2\pi i/h$ by that at $D$.

For $f=T^{\rm cl}_1$, let $u=e^{i\pi/h}$. With
$x=e^{-D}$ the string product is

\[
\prod_{k=0}^{h-1}\frac{u-xe^{-2\pi ik/h}}{u^{-1}-xe^{-2\pi ik/h}}
=\frac{u^h-x^h}{u^{-h}-x^h}=1.
\]

The example is right. Multiplying over the two breather strings specified
in the chapter, `scalars.py` independently finds

```text
omega=sqrt(2): -0.2199993795 - 0.9755000118 i
omega=sqrt(5):  0.2004425131 - 0.9797054654 i
```

These agree with the stated $a_3, B^{(2)}_1, D=0.41$ values and show
the claimed coupling dependence.

However, for *even* $h$, the simpler equation
$f(D)f(D+i\pi)=1$ is **sufficient**. Pair the string entry $k$ with
entry $k+h/2$; every pair is 1. It is not necessary. For odd $h>2$,
it is neither necessary nor sufficient. The sinh example gives string
values $0.5920002778-0.8059377588i$ at $h=3,D=0.41$, but exactly 1
at every even $h$; the script checks $h=4,6,8$ to roundoff. The current
sentence “For $h>2$ that simpler condition is neither necessary nor
sufficient” is false.

**Paste-ready replacement**, `part8-defects/43-particle-transmission.qmd:120`:

```markdown
For $h=2$ it is $f(D)f(D+i\pi)=1$. For $h>2$ that simpler condition is
not necessary. It is sufficient when $h$ is even, since the $h$-string
pairs arguments separated by $i\pi$, but is not sufficient in general
when $h$ is odd. The factor $\mathbb T^{\rm cl}_1$ satisfies
(-@eq-dpt-soliton-cdd) for every $h$, while generally failing the simpler
condition; $(\sinh D+0.3i)/(\sinh D-0.3i)$ satisfies the simpler condition
but fails the string condition at $h=3$.
```

### Imaginary defect rapidity and Lee–Yang poles

For $\theta_d=iv$, $x=e^{-\theta}e^{iv}$,

\[
|x+e^{i\gamma_a}|^2-|x+e^{-i\gamma_a}|^2
=4e^{-\theta}\sin v\sin\gamma_a.
\]

Since $0<\gamma_a<\pi$, the factor is a phase on the real axis exactly
when $v\in\pi\mathbb Z$. This is also exactly when
$\sigma=e^{-iv}$ is real and the stated $\Theta$ symmetry is
unbroken. At $v=\pi$, the particle pole has $u=\gamma_a$; it lies in
the positive-energy bound-state window $0<u<\pi/2$ exactly for
$a<h/2$. Its residue is

\[
\frac{e^{i\gamma_a}-e^{-i\gamma_a}}{e^{-i\gamma_a}}
=e^{2i\gamma_a}-1,
\]

as in the exercise. For sine-Gordon the same window becomes
$0<\pi-\pi/(2\lambda)<\pi/2$, equivalent to
$1/2<\lambda<1$. All the corrected exceptions are right. A pole at
$\sigma=-1$ does not cure the unbounded energy established below.

For Lee–Yang, $f_{2/3}(\theta)$ has poles at $i\pi/3,2i\pi/3$.
In $\bar{\mathbb T}_b=f_{2/3}(\theta+i\pi(3-b)/6)$, these become
$i\pi(b-1)/6,i\pi(b+1)/6$. In $\mathbb T_b$ they instead become
$i\pi(5-b)/6,i\pi(7-b)/6$. In the indicated fundamental ranges, the
two bound-state poles are therefore in the barred factor, as the corrected
caption says.

### TBA and Jiang

Let $x=e^{-(\theta-\theta_d)}>0$ and $u=e^{i\gamma_a}$. Direct
substitution into the proposed mirror fugacity gives

\[
F_a=\frac{-ix-u^{-1}}{-ix-u},\qquad
F_{\bar a}=\overline{F_a},\qquad
|F_a|^2=\frac{x^2+1-2x\sin\gamma_a}{x^2+1+2x\sin\gamma_a}<1.
\]

For the self-conjugate species this is
$F_{h/2}=(x-1)/(x+1)=-\tanh((\theta-\theta_d)/2)$. All three
identities are independently simplified in `scalars.py`. The strict
inequality holds at finite real rapidity; the limiting modulus at infinity
can equal 1. A logarithmic driving term needs a branch at the sign change,
whereas the fugacity form handles the zero directly. With conjugation-
symmetric kernels, paired pseudoenergies and compatible logarithm
branches give paired conjugate integrands in the energy. The chapter
correctly asserts reality only for a solution with that symmetry, not
existence or uniqueness of such a solution.

The change between the two forms is also algebraic: if the driven
pseudoenergy is $\epsilon^{\rm d}$, put
$\epsilon=\epsilon^{\rm d}+\log F$. Then
$e^{-\epsilon^{\rm d}}=F e^{-\epsilon}$, moving the defect from the
driving term into the logarithm. This agrees with Bajnok–Simon's §5.1;
their §3.1 (18) is the defect-energy formula, and §5.3 recovers it in the
ultraviolet. Their §7 discusses the extension to diagonal scattering.
[Bajnok–Simon](https://arxiv.org/abs/0712.4292).

Jiang's (5.5) and (5.6) give, respectively, the two-particle correction and
its defect contribution with integrands proportional to
$\widehat T(\theta)^2e^{-2mr\cosh\theta}$ and
$(\widehat T(\theta)^2-1)e^{-2mr\cosh\theta}$. Setting
$\widehat T=1$ removes the defect contribution exactly. His §§4–6 and
§7 support the revised bullet's UV/IR distinction and its stated examples.
It does not imply an all-$a_r$ entanglement calculation, and the chapter
does not claim one.
[Jiang](https://arxiv.org/abs/1703.03562).

### Unbounded energy and layer estimate

For $\phi=\psi=tv\ne0$, $g_j=1$ and
$f_j=e^{t\beta\alpha_j\cdot v}$. Since the affine roots span the
root space and sum to zero, a nonzero $v$ has
$A=\max_j\beta\alpha_j\cdot v>0$. Thus for negative $\sigma$
the defect potential is bounded above by $-c e^{At}+O(1)$, $c>0$,
for both types (with $\sigma$ exchanged with $1/\sigma$).

To check that bulk energy does not invalidate this conclusion, use on
each side a linear ramp from 0 to $tv$, of width $\ell=e^{-At/2}$.
The total gradient energy is $t^2|v|^2/\ell$. Every bulk exponential
along the ramp is at most $e^{At}$, so the total potential energy is
at most $C\ell e^{At}+C'\ell$. Therefore

\[
E_{\rm total}\le-c e^{At}+C''(t^2+1)e^{At/2}+O(1)\longrightarrow-\infty.
\]

These ramps are admissible finite-energy $H^1$ configurations. They need
not solve the static Euler–Lagrange equations to test boundedness of the
energy functional. Smoothing their outer corners preserves the estimates.
The explicit direction and the stated big-O estimate are correct.

### Citations and novelty

The attribution to CZ2007 (7.11) and its comparison with (7.5) is right.
The sentence after (7.11), printed p. 23 in the local copy, explicitly
uses that agreement to justify the choice of $\Delta$ in (6.12).
The new wording accurately says that, rather than attributing a broader
claim to them.

There is nevertheless a direct antecedent for the CDD condition: immediately
after CZ2007 (6.9), the paper writes the three-factor condition

\[
h(\theta)h(\theta+2\pi i/3)h(\theta+4\pi i/3)=1.
\]

Changing the argument by $-i\pi$ makes it exactly the chapter's condition
for $h=3$. Thus “New here: ... the soliton CDD condition” needs the
qualification *all-rank form/derivation*, and a citation to this antecedent.
This is a specific source comparison, not a claim to have searched all
prior literature.
[CZ2007, after (6.9) and after (7.11)](https://arxiv.org/abs/0705.1066).

**Paste-ready correction**, `part8-defects/43-particle-transmission.qmd:288`:
replace the CDD item in the novelty list with

```markdown
the all-rank form and derivation of the soliton CDD condition
(-@eq-dpt-soliton-cdd), whose $h=3$ case is already stated by
Corrigan and Zambon [@cz2007, after eq. (6.9)];
```

Also add this antecedent after `eq-dpt-soliton-cdd` in the main discussion.
For the residue sentence in the spectrum table, adding “on each fixed-charge
sector” would prevent the rank-one statement from being misread as a
finite rank for the whole infinite-dimensional operator.

## Other observations — unverified

- The new all-rank second-module proof does not address the type-I+ doubled
  module, intermediate fused multiplets, or a classification of larger
  symmetry coideals. Those are separate questions.
- No attempt was made to determine the physical Tzitzéica scalar factor or
  remove its CDD ambiguity. The repaired Gamma functional equations do not
  settle that physical identification.
- No exhaustive priority search was made for the all-rank particle ansatz,
  the second module, or the spin-two preservation result. The specific
  CZ2007 CDD antecedent above is verified; broader novelty remains unverified.
- No existence, uniqueness, or ultraviolet solution of the proposed
  multi-species defect TBA is established here.

## Integration (2026-10-13)

The referee's scripts in `briefs/referee_defects_upgrades_scripts/` (`algebra.py`, `scalars.py`, `spin_euler.py`, `doubling.py`) were not used as verification. The integrating session (Claude Opus 5.5) re-derived every item by hand and with its own scripts, written from the chapters' definitions and importing at most the book helper `defect_common.sg_transmission_scalar`: `integration_item1_3_sectors.py` (numpy: sector determinants and moduli of the type-I$_+$ $V$ matrix), `integration_item2_reversed.py` (sympy, exact, symbolic labels, $h=3,4,5$: `eq-dt-reversed`), `integration_item5_cz2011.py` (mpmath: CZ2011 (3.12) against (3.8), (3.10) in both readings of (3.2); the finite identity; pole arguments), `integration_item8_scalar_limits.py` (the scalar factor in the limits of `prp-diagonal-borel` (5)), `integration_item9_typeplus.py` (sympy: the type-I$_+$ particle identities in $x$), `integration_item10_second_module.py` (sympy, exact, symbolic labels, $h=3,4,5$: the sparse inverse, the $\bar V$ half, the Cartan equations, the module relations, with controls), `integration_item11_cdd_rank.py` (the simpler CDD condition for even and odd $h$; rank one of the $\bar V$ residue), `integration_refcheck.py`; all with `*_output.txt`, and the final book-script outputs as `integration_book_*_output.txt`. Sources read: CZ2011 (1012.4186) PDF p. 4–5, (3.2), (3.5)–(3.13); CZ2007 (0705.1066) (6.5)–(6.9) and the sentence after (6.9). Tolstoy–Khoroshkin 1992 and Frenkel–Hernandez 2015 are not available locally, so the suggested pinpoint citations were not added.

- **Item 1. Verified; optional precision added.** The determinant list was re-derived (the type-I$_+$ $V$ matrix, an inverse, has $\det=q^{-r}/(1-qX_+^h)$: poles, not zeros; checked numerically for $h=3,4,5$). Change: ch. 41, `prp-an-transmission` (3) now assumes $\theta$ such that the matrices are "defined and invertible on their nonempty sectors of fixed total charge", and says that the type-I$_+$ $V$ matrix has poles rather than zeros.
- **Item 2. Verified; coefficient argument added.** By hand and exactly ($h=3,4,5$): every residual of the $f_j$-equation for `eq-dt-reversed` is a multiple of $X/y+cq-c/q$, $c=\mu/y_d$, on the inputs $\bar v_{j+1}\otimes|N\rangle$ and $\bar v_j\otimes|N\rangle$; with $X=-(q-q^{-1})\mu y/y_d$ none remain. Change: the proof of item 2 now begins with this all-rank check (the Chevalley step is unchanged).
- **Item 3. Verified; clarification added.** At $\theta=\theta_d$ the type-I$_+$ $V$ matrix has entries of equal modulus $1/|1-qX_+^h|$ (0.637 for $h=3$, $\omega=1.37$), not 1; the relative sizes are $\min(t^d,t^{d-h+1})$, $t=|X_+|$. Change: the last sentence of item 5 (same relative bound with $X_+$; equal, not necessarily unit, modulus).
- **Item 4. Verified (as at the ch. 39 integration); no change.**
- **Item 5. Verified.** The PDF prints $x_i=q^{2\pi\theta_i/\xi}$ in (3.2). With $x=e^{2\pi\theta/\xi}$ the equal-$k$ product with one constant, $C=1/(2H(a)^2)$, solves (3.10) and (3.8) to $10^{-13}$ at $\xi=0.93$, $1.73$ ($C=0.0241346292$, $0.1217023253$); with the printed $x$ both fail (0.08–2.8). The finite identity for the printed grouping holds to $10^{-28}$; the unnormalized printed products diverge also at $z=0$ (like $K^{3a}$ numerically), so "drift like $K^{-2z}$" was qualified by "after division by their value at $z=0$". Changes: ch. 41 `sec-a22-transmission` (we read $x=e^{2\pi\theta/\xi}$; the printed variable fails; one constant normalization; their $\tau$), ch. 41 Notes (new (3.2) bullet; "in this convention the slip is in the pole positions"; normalization), Appendix D (new item "(3.2), the spectral variable", *Misprint*; the (3.12) item qualified). Script check (j) gains the control with the printed (3.2) (fails by 11). The CZ2009a entry is unchanged.
- **Item 6. Verified.** The order-0 and order-1 terms of $\mathcal R_{21}^{-1}$ on $V$, $\bar V$ and of $\mathcal R$ on $\bar V$ were recomputed by hand and reproduce the matrix parts exactly. Changes: ch. 41 proof of item 4: the Cartan lift (a fixed logarithm of $q$ defines $q^{(N+s)_i}$, $q^t$ acts through $N+s$, only integral weight differences enter) and the remark that $\mathcal W$ need not lie in Frenkel and Hernandez's category $\mathcal O$, the finiteness of each coefficient being all that the evaluation needs. *Not done:* the pinpoint citations (sources not available locally).
- **Item 7. Verified.** By hand: the code returns $(T+\Theta)(\phi)-(T+\Theta)(\psi)$; $d_{-\epsilon}T=d_\epsilon\Theta$ is $\partial_t(T-\Theta)=\epsilon\partial_x(T+\Theta)$, so $\dot Q_{\rm bulk}=\epsilon G$, the defect contribution is $-\epsilon F$ for $D_tF=G$, and the spin is $-\epsilon\cdot$spin. Change: the `boundary_G` docstring of `code/defect_spin_two_a3.py` (`defect_classical_an.py` has no such docstring). Rerun: 9 PASS, 151 s. No chapter change.
- **Item 8. Verified.** $|\Phi|^2(1+Y^2)=1$, $\Phi\to1$ as $\theta_d\to+\infty$, and with $\Phi$ the entries of $\mathbb L/(y/y_d)$ tend to zero as $\theta_d\to-\infty$; the scalar cancels from every intertwining equation. Since $M(\theta)$ already denotes the doubled multiplet in ch. 42, the matrix part is written $\mathbb T=\Phi\hat{\mathbb T}$, $\hat{\mathbb L}=P\hat{\mathbb T}$ rather than $L_M=PM$. Changes: `prp-diagonal-borel` (5) (limits of $\hat{\mathbb L}$; the full $\mathbb L$ has the same $\theta_d\to+\infty$ limit; $\mathbb L/(y/y_d)\to0$ and $\mathbb L/(\Phi Y)$ has the conjugation limit with $|\Phi Y|\to1$; type-I$_+$ twist matrices ${\rm diag}(q^n,q^{-n-s})$, beyond the report), the proof (opening sentence; $\hat{\mathbb L}$ throughout), the caption of `fig-dbd-symmetries` ("becomes proportional to" in the $-\infty$ limit), the docstrings of `defect_doubling.py` and of the figure script, which repeated the slip. No label change. Rerun: 20 PASS; figure script PASS, regenerated SVG restored (ids only).
- **Item 9. Verified; optional rewrite done.** The type-I$_+$ outline of `prp-an-particle-transmission` is replaced by $\mathbb T^+_a=(1-xu)/(1-x/u)$ in the same variable, with both fusings displayed (cancelling factors $1-xu/v$ and $1+xv/u$); exact sympy check.
- **Item 10. Verified; label upgraded.** The sparse inverse was re-derived by hand ($M_{\bar V}=\sum_{d<h}(XU)^dA$, $U^h=q^{-1}$, $H=A^{-1}(1-XU)$) and exactly for $h=3,4,5$, together with the $\bar V$ equation (directly and in the inverted form), its controls, the Cartan equations for $V$ and $\bar V$, the general solution for $V$, and the module relations; the failed-twist supports and the coideal obstruction ($q^2$, $q^{-1}$; $q^3\ne1$) were re-derived by hand. **Label change:** `prp-dpt-doubled` [Numerical] ($a_2$, $a_3$) → **[Theorem]** for type I$_-$, every $r\ge2$, generic $q$ ($q^3,q^4\ne1$ suffice), finite $\theta_d$, any offset; the statement now includes the failure of the other twists for arbitrary $\theta$-independent $Y$, $\kappa$ and the coideal obstruction. The proof block contains every step: the $f$-type half, the $V$ half, the module relations, the $\bar V$ half, the $h=2$ reduction, the failed twists and the coideal obstruction. Other changes: problem 12 (points to the theorem), Notes novelty ("the proof, for every $r$"), ch. 42 "Beyond sine-Gordon" ("proves this for type I$_-$ of every $a_r^{(1)}$"); script check (h) gains $HM_{\bar V}=1-X^h/q$ on exact labels. The proof comes from this report, so the upgrade has had a referee pass; the integration text itself has not.
- **Item 11. Verified, with the two corrections and the qualifier.** The sinh-type factor passes the string condition for $h=2,4,6,8$ and fails it for $h=3,5,7$; $\mathbb T^{\rm cl}_1$ passes it for every $h$ and fails the simpler condition for $h>2$; CZ2007 states the $h=3$ condition right after (6.9). Rank one of the $\bar V$ residue: all $2\times2$ minors vanish exactly ($h=3,4,5$). Changes: the sentence after `eq-dpt-soliton-cdd` (not necessary; sufficient for even $h$; not sufficient in general for odd $h$; CZ2007 cited), the Notes novelty item ("the all-rank form and derivation … whose $h=3$ case is stated by Corrigan and Zambon"), "on each fixed-charge sector" in `tbl-dpt-spectrum` and problem 10, check (i) of `defect_particle_transmission.py` (docstring, report text, the sinh-type factor passes at $h=4,6$ and fails at $h=3,5$), the status file and the conventions brief. The residue, imaginary-parameter, Lee–Yang, TBA and vacuum statements were verified at the ch. 43 integration and are unchanged.

**Rejected:** nothing. **Modified:** item 8 (notation $\hat{\mathbb L}$ instead of $L_M$, because of the clash with $M(\theta)$; type-I$_+$ twist matrices added); item 6 (no pinpoints).

**Other observations of the report.** Type-I$_+$ doubled module, fused multiplets and larger coideals: not addressed (problem 12 lists the type-I$_+$ case). The physical Tzitzéica scalar factor: unchanged, still open. Priority searches: no new novelty claims beyond the qualified CDD item. Defect TBA existence: unchanged.

**Book scripts after integration:** `defect_spin_two_a3.py` 9 PASS (151 s); `defect_an_transmission.py` 28 PASS (192 s); `defect_doubling.py` 20 PASS (23 s); `defect_particle_transmission.py` 11 PASS (about 1 min). `defect_common.py` unchanged. Reference check over the 65 book `.qmd` files: 1164 ids, no new ids, no duplicates, every cross-reference resolves, every citation key is in `references.bib`.

**Labels changed.** `prp-dpt-doubled`: [Numerical] → [Theorem] (type I$_-$, every $r\ge2$). No other label changed, none downgraded. The upgrades refereed by this report (`prp-an-transmission` (2)–(5), `prp-an-linear-transmission`, `prp-spin-two` (3) for $a_3$, `prp-diagonal-borel` (2)–(5), `prp-an-particle-transmission`) now have had an independent referee pass.

**Pending (outside the files this integration could edit).** Ch. 40 `tbl-boundary-defect-algebra`, doubled-boundary status cell: "$a_2$, $a_3$ [Numerical]" → "$a_r$, type I$_-$ [Theorem]". Optionally `index.qmd` l. 51, "or for $a_2^{(1)}$ and $a_3^{(1)}$ by the rotation …" → "or for $a_r^{(1)}$ …" (not wrong as it stands).
