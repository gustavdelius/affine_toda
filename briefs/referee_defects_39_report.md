# Deep referee report on chapter 39

Date: 2026-10-12. Model: GPT-6 (Codex). Starting commit: `23dff2d0f4b25bcb8ffbd6859f1f69a26c9feb08`. Commit at the manuscript/reference check: `418d1924251b67e0c0f4bfa38002fea38c555e59`. Other agents advanced HEAD during this review.

This report follows `briefs/referee_defects_39.md`. I edited no pre-existing file, ran no Git command that changes the index or working tree, and did not run Quarto. The correction target is `part8-defects/39-defect-solitons.qmd`, unless another file is specified. Its audited SHA-256 is `b378c80b1a3fc9c8c507e2180aa58ef84e8a56ebce3523273673f4e2d1b6db49`; line numbers below refer to that version.

The evidence is independent: none of the new scripts imports or runs the book's verification scripts. The numerical computations construct the roots, weights, matrices, fields, potential gradients and scattering systems anew. In particular, they evaluate the exponentials of the fields, not just identities between the claimed delay formulas. The symbolic computation verifies rational identities for arbitrary soliton amplitude and rapidity, covering all species and classes of both requested algebras and, in fact, general rank.

Reproduce with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_39/exact_checks.py
PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_39/independent_checks.py
PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_39/supplemental_checks.py
PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_39/metadata_checks.py
```

The main formulas survive this review. The substantive corrections concern the scope of the absorption theorem, the distinction between real energy and a real energy–momentum pair, and the direction assumed in the selection rule. The frozen-soliton factorization, its obstruction for two defects of one type, and the proposed correction of CZ2007 (7.11) pass independent checks.

## 1. Type-I+ delays and jump classes

**Verdict: correct as meromorphic solution identities.** The signs of both class shifts and the time reversal are correct. “Any complex rapidity and position” needs the usual exclusion of zeros of the tau-functions and a limiting interpretation at poles of the delay; it does not establish globally regular solutions.

### Exact independent computation

Use the book's roots, but construct their vectors and fundamental weights independently. If both vacua are shifted by $2\pi\lambda/\beta$, then $f_j$ acquires $e^{2\pi i\alpha_j\cdot\lambda}=1$ and $g_j$ is unchanged. Thus only the difference of the vacua matters. For a right-hand shift by $2\pi\lambda_p/\beta$,

\[
\alpha_j\cdot D\lambda_p
=2(\lambda_j-\lambda_{j-1})\cdot\lambda_p
\equiv-2p/h\pmod2.
\]

Consequently $f_j\mapsto\omega^{-p}f_j$, $g_j\mapsto\omega^pg_j$. The coefficients of these two families are reciprocal, so the effective coefficient of $f_j$ is

\[
S=\sigma\omega^{-p}\quad(\mathrm I_-),\qquad
S=\sigma^{-1}\omega^{-p}\quad(\mathrm I_+).
\]

This proves the two effective rapidities $\theta_d+2\pi ip/h$ and $\theta_d-2\pi ip/h$, including their different signs.

Here is an exact sewing calculation independent of first-order matching. Put $v=e^{i\pi a/h}$, $q=v^2$, $K=e^\theta$, and let $s=1$ for I− and $s=-1$ for I+. Use a formal variable $w=\omega^{ja}E$, and set

\[
L=\frac w{1+w},\quad L_+=\frac{qw}{1+qw},\quad
R=\frac{zw}{1+zw},\quad R_+=\frac{zqw}{1+zqw},
\]
\[
F=\frac{(1+qw)(1+zw/q)}{(1+w)(1+zw)},\quad
F_+=\frac{(1+q^2w)(1+zw)}{(1+qw)(1+zqw)},\quad
G_+=\frac{(1+w)(1+zqw)}{(1+qw)(1+zw)}.
\]

Set $m=\beta=1$, $k_x=m_a\cosh\theta$, $k_t=-m_a\sinh\theta$, $m_a=(v-v^{-1})/i$. Multiply each vector sewing residual by $1-\pi$, which is invertible on the root space, and use $D^\mathsf T=-D\pi$. Its coefficients in the affine roots, divided by $i$, are

\[
C_1=k_x(L-L_+)+s k_t(L+L_+-2R)-SF+S^{-1}G_+,
\]
\[
C_2=k_x(R-R_+)+s k_t(2L_+-R-R_+)-SF_++S^{-1}G_+.
\]

For

\[
z=\frac{K^{-s}/S-iv}{K^{-s}/S-i/v},
\]

**both rational functions are identically $C_1=C_2=-S+S^{-1}$**. Their contribution to the vector residual vanishes because $\sum_j\alpha_j=0$. This verifies all orders in $E$, with arbitrary $K,S,v$; it is stronger than checking a few times or the linearized equations. Substitution of the two $S$'s gives exactly lines 34–36. Restoring $m,\beta$ multiplies both sides of the sewing equations by the same dimensional factors.

`exact_checks.py` performs this cancellation symbolically and separately verifies the weight congruence with exact rational matrices for $a_2,a_3$:

```text
type I-: exact sewing root-coefficient residuals = [0, 0]
  a_2: all 2 species x 3 classes, 12 identities are specializations of zero
  a_3: all 3 species x 4 classes, 24 identities are specializations of zero
type I+: exact sewing root-coefficient residuals = [0, 0]
  a_2: all 2 species x 3 classes, 12 identities are specializations of zero
  a_3: all 3 species x 4 classes, 24 identities are specializations of zero
a_2: exact alpha_j D lambda_p + 2p/h is even for every j,p
a_3: exact alpha_j D lambda_p + 2p/h is even for every j,p
```

The separate direct-field computation in `independent_checks.py` uses $\phi=i\sum\alpha_j\log\tau_j$, differentiates the tau-functions analytically, and evaluates the potential gradients from the original exponentials:

```text
1 direct exponentials, exact one-soliton profiles: 216 samples, residual 1.315e-13
1 negative control correct charged delay: sewing residual=4.14813554e-15
1 negative control wrong class-zero delay: sewing residual=6.67418564
```

It covers both types, all classes and species of $a_2,a_3$, real and complex rapidities, and points with $x\ne0$. The nonzero negative control excludes a test that merely agrees with its own delay formula.

Time reversal changes $k_t$ and both kinetic sewing matrices in sign. The I− potential at $\sigma$ is the I+ potential at $1/\sigma$. The constant jump does not change. Therefore $\theta\mapsto-\theta$, $\theta_d\mapsto-\theta_d$, $p\mapsto p$, which reproduces the displayed I+ formula without changing the sign of $p$ by hand.

**Paste-ready precision edit, line 42:** replace “with any complex $\theta$ and $\xi$” by:

```markdown
with complex $\theta$ and $\xi$, interpreted on domains where the displayed $\tau$-functions are nonzero and by a limiting configuration at a zero or pole of $z$
```

The exact calculation above also supplies the all-orders justification missing from “as they must for a Bäcklund pair” at line 51. The Bäcklund interpretation alone would not prove that this particular finite-amplitude ansatz solves its equations.

## 2. Absorption and emission: `prp-ds-absorption`

**Verdict: wrong as an unrestricted [Theorem]; correct after excluding the zero-envelope case and stating the meromorphic/localized-solution scope.** The pole locations, real-rapidity congruence, charged $a_5$ example, and energy–momentum difference are correct. The proof does not justify an absorption interpretation for every class at every pole.

### Poles, limits and energy independently

Let $\epsilon=1$ for I− and $\epsilon=-1$ for I+. A denominator zero gives

\[
\theta_* =\theta_d+i\epsilon\left(\frac{2\pi p}h+\gamma_a-\frac\pi2\right)
\pmod{2\pi i},
\]

which is line 95. There are two well-defined rational limits: keep $E$ fixed and send $z\to\infty$, giving a soliton on the left and its right vacuum on the right; or keep $zE$ fixed and send $E\to0$, giving a vacuum on the left and a soliton on the right. Direct substitution of these two configurations, rather than large finite delays, gives:

```text
2 pole limits: 768 configurations, residual 5.380e-13
```

These tests cover $h=3,4,6$, every species and class, both defect types and $\theta_d=\pm0.6$. They verify the sewing identity even in degenerate cases; that does not make those cases localized solitons.

Writing $Q=2\pi p/h$, the class energy difference is

\[
\Delta B
=\hat M[\cosh(\theta_d+i\epsilon Q)
-\cosh(\theta_d+i\epsilon(Q+2\gamma_a))]
=-2i\epsilon\hat M\sin\gamma_a\,
\sinh(\theta_d+i\epsilon(Q+\gamma_a))
=M_a\cosh\theta_*.
\]

The analogous difference-of-sinh identity gives $\Delta P=M_a\sinh\theta_*$. This proves the energy claim exactly, for both types and all $p$, as an algebraic identity. The numerical sweep through $h=2,\ldots,9$ has maximum residual $4.632\times10^{-14}$.

At $h=6,a=1,p=1$, $\theta_*=\theta_d$ and $M_1=12$ when $m=\beta=1$. At $\theta_d=0.6$:

```text
type=-1: stored EP=[14.22558262+0j, 7.63984299-3.55e-15j], error=3.662e-15
type=+1: stored EP=[14.22558262+0j, 7.63984299+3.55e-15j], error=3.662e-15
```

Thus the example is correct. In I− the two class energies are

\[
B_1=\hat M\left(\tfrac12\cosh\theta_d
-i\tfrac{\sqrt3}2\sinh\theta_d\right),\qquad
B_2=\hat M\left(\tfrac32\cosh\theta_d
-i\tfrac{\sqrt3}2\sinh\theta_d\right).
\]

Their difference is the real $M_1\cosh\theta_d$, although each energy is complex when $\theta_d\ne0$. I+ gives their conjugates.

The real-pole condition is $a+2p\equiv h/2\pmod{2h}$, with $p$ taken modulo $h$. It has no solution for odd $h$. Independent enumeration gives:

```text
h=4: [(2, 0)]
h=6: [(1, 1), (3, 0), (5, 5)]
h=8: [(2, 1), (4, 0), (6, 7)]
```

Here each pair is (species, initial absorption class). The parity condition in the chapter is necessary and sufficient for some class to exist, not a list of all classes for a fixed species.

### What the theorem misses

For $\theta_*=\theta_d+i\chi$,

\[
|E|=\exp\{m_a\cos\chi\cosh\theta_d
(x-t\tanh\theta_d)+\operatorname{Re}\xi\}.
\]

When $\cos\chi\ne0$, the envelope velocity is indeed $\tanh\theta_d$. Its sign establishes the stated left/right orientation. When $\cos\chi<0$, factor 

\[
1+\omega^{ja}E=\omega^{ja}E
\left(1+\omega^{-ja}E^{-1}\right).
\]

This gives the reversed species $\bar a$, rapidity $\theta_*+i\pi$ modulo $2\pi i$, and reverses the class transition. This part of item 2 is correct.

But **$\cos\chi=0$ is omitted**. It occurs precisely when $a+2p\equiv0\pmod h$. For $h=4,a=2,p=1$, I− has $\theta_*=0.6+i\pi/2$. Then $\operatorname{Re}\cosh\theta_*=\operatorname{Re}\sinh\theta_*=0$, so $\lvert E\rvert=e^{\operatorname{Re}\xi}$ everywhere. There is no localized envelope, no vacuum limit along real $x$, and no passage that can be called absorption or emission. The script outputs 

```text
x=0,t=-1: |E|=1.0000000000000000
x=0,t=2:  |E|=0.9999999999999999
x=.7,t=.3: |E|=1.0000000000000000
```

Also $\theta_d=0$ makes any nondegenerate envelope stationary. It is a stationary limiting configuration, not an incoming or outgoing event.

Finally, a genuinely complex rapidity does not generically give a globally regular real-spacetime soliton. On its envelope centre, for $\xi=0$,

\[
x=t\tanh\theta_d,\qquad
\operatorname{Im}\log E=-m_a\sin\chi\,t/\cosh\theta_d.
\]

If $\sin\chi\ne0$, this crosses the phases for which some $\tau_j=0$. At $h=3,a=1,p=0,\theta_d=0.6$, the absorption pole has $\chi=-\pi/6$. `supplemental_checks.py` finds

```text
x=-2.309523725, t=-4.300392118, |tau_0|=5.666e-16
```

The singularity is inside the left half-line. The book may consistently discuss these as meromorphic configurations, as chapter 11 does; it must not silently promote them to globally smooth absorbed particles.

### Paste-ready corrections

**Replace line 101 and the opening of item 1 at line 103 with:**

```markdown
**[Theorem for the meromorphic one-soliton configurations]** Let $\theta_d\in\mathbb R$, and consider a type-I$_\mp$ defect of class $p$ and species $a$. Set $\theta_* = \theta_{\rm abs}$ and $\chi=\operatorname{Im}(\theta_*-\theta_d)$. The statements about localized absorption and emission below assume $\cos\chi\ne0$ and $\theta_d\ne0$; the field identities hold wherever the limiting $\tau$-functions are nonzero.

1. At $\theta=\theta_*$ the two pole limits give a soliton on $x<0$ with $\vec\psi$ constant, and $\vec\phi$ constant with a soliton on $x>0$. They satisfy the defect conditions as meromorphic identities. Their localized envelopes have velocity $\tanh\theta_d$. For $\theta_d>0$ they describe absorption from the left and emission to the right; for $\theta_d<0$, emission to the left and absorption from the right. At $\theta_d=0$ the envelope is stationary. If $a+2p\equiv0$ mod $h$, then $\cos\chi=0$: the limiting fields still obey the sewing identities, but have no localized envelope and are not absorption or emission configurations.
```

Keep items 2–3 with the new assumptions. Add at the end of the proof, line 113:

```markdown
At genuinely complex rapidity these identities do not assert global regularity on real spacetime: the phase on the envelope centre varies with time and can produce zeros of the $\tau$-functions. The absorption interpretation is an analytic continuation of the localized pole configuration.
```

**Replace the first two sentences of line 116 with:**

```markdown
For type I$_\mp$, a jump from class 0 to class $a$ stores the formal energy and momentum $M_a(\cosh\theta_*,\sinh\theta_*)$, with $\theta_* = \theta_d\mp i(\pi/2-\gamma_a)$. At real $\theta_d\ne0$ its energy is real only for $a=h/2$. At $\theta_d=0$ the energy is real for every $a$, while the momentum is real only for $a=h/2$.
```

The existing proof establishes the corrected local/meromorphic theorem; its current unrestricted physical wording is not justified.

## 3. Stored charge, mass shell and selection rule

**Verdict: correct core formulas; wrong energy-only reality claim, and incomplete incidence and self-conjugate qualifications.** The sharp step bound and the assertion about missing charges hold within the one-soliton scattering family under discussion.

### Class energies and reality

At a constant jump, $f_j=\omega^{-p}$, $g_j=\omega^p$. Substitution into the original potential and momentum gives line 61–62 directly. With $Q=2\pi p/h$, 

\[
B_p=\hat M[\cosh d-\cosh(d+i\epsilon Q)],\qquad
P_p=\hat M[\sinh d-\sinh(d+i\epsilon Q)].
\]

After subtracting $\hat M(\cosh d,\sinh d)$, the energy–momentum pair is 

\[
(-\hat M\cosh(d+i\epsilon Q),-\hat M\sinh(d+i\epsilon Q)),
\]

whose Minkowski square is exactly $\hat M^2$. It is an algebraic mass shell; complex classes are not real timelike four-vectors. Independent numerical residual: $1.449\times10^{-13}$ through $h=9$.

For real $d$,

\[
\operatorname{Im}B_p=-\epsilon\hat M\sinh d\sin Q,\qquad
\operatorname{Im}P_p=-\epsilon\hat M\cosh d\sin Q.
\]

Thus **both** are real iff $2p\equiv0\pmod h$; the energy alone is real also at $d=0$, for every class. A counterexample to line 65 is

```text
h=3, p=1, theta_d=0, m=beta=1:
EP=[9+0j, 0-5.19615242j]
```

Conjugation $p\leftrightarrow h-p$ is correct, including this special case.

**Replace the last two sentences of line 65 with:**

```markdown
The energy and momentum are both real exactly when $2p\equiv0$ mod $h$. The energy alone is real when either this condition holds or $\theta_d=0$. Classes $p$ and $h-p$ have conjugate energy–momentum pairs.
```

### Charge bookkeeping and the actual cycle

For positive rapidity, as $t$ goes from $-\infty$ to $+\infty$, the defect sees the right vacuum and then the left vacuum of each profile. Hence

\[
\Delta(\psi-\phi)=\frac{2\pi}{\beta}(\lambda_L-\lambda_R).
\]

Both charges belong to the same representation class, so this difference lies in the root lattice; class and static energy are unchanged. For negative rapidity the temporal order is reversed. With fixed left/right labels the general formula is

\[
\Delta(\psi-\phi)=\frac{2\pi}{\beta}\operatorname{sign}(\theta)
(\lambda_L-\lambda_R).
\]

If the labels genuinely mean incoming and outgoing, charge conservation still says $\lambda_{\rm in}-\lambda_{\rm out}$, but the outgoing phase shift for right incidence is $-\arg z$, since $\phi$ is then the outgoing field. Lines 76 and 85 allow all real rapidities while using the positive-rapidity proof and phase convention.

Independent continuation at $h=5,a=2$, with equal relative rapidities but opposite velocities:

```text
theta=+0.8: jump=[0,-1,1,-1,1], lambda_left-lambda_right=[0,-1,1,-1,1]
theta=-0.8: jump=[0,1,-1,1,-1], lambda_left-lambda_right=[0,-1,1,-1,1]
residual with sign(theta)=1.665e-16 in both cases
```

For the cycle, derive the charge directly from continuation along real $x$:

\[
\lambda_a(\eta)=-\frac1{2\pi}\sum_j\alpha_j
\left[\eta+\frac{2\pi ja}h\right]_{(-\pi,\pi)}.
\]

The singular boundaries are $\eta=\pi-2\pi ja/h$ modulo $2\pi$. There are $n=h/g$ distinct boundaries, spaced by $\delta=2\pi g/h$, where $g=\gcd(a,h)$. A phase shift $2\pi a/h=(a/g)\delta$ gives $\lambda_a(\eta+2\pi a/h)=\pi\lambda_a(\eta)$. Since $a/g$ is coprime to $h/g$, this phase cycle is exactly one $\pi$-orbit. Adjacent phase sectors need not be adjacent applications of $\pi$; the chapter correctly records this distinction.

Multiplying the delay numerator by the conjugate denominator gives a positive real denominator and a numerator $x^2+e^{2i\gamma_a}$. Therefore its continuous phase lies between $0$ and $2\gamma_a$ for $a<h/2$, and between $0$ and $-2\gamma_{\bar a}$ for $a>h/2$. An interval of length less than $K\delta$ crosses at most $K$ equally spaced boundaries, where

\[
K=\frac{\min(a,h-a)}{g}.
\]

The upper bound is attained at a **finite** relative rapidity: choose the phase shift between $(K-1)\delta$ and $K\delta$, then choose the incoming phase so that it crosses $K$ boundaries. The lower bound is attained by making the shift small and staying inside one sector. I+ exchanges the large- and small-rapidity limits. When $a=h/2$, the phase is $0$ or $\pi$, giving $k=0$ or $1$ away from the absorption pole; the printed item 3 omits this case.

There is a further incidence qualification: this sharpness argument lets the relative rapidity range over the whole real line. It does not establish that both endpoints are accessible for a fixed defect rapidity with incidence restricted to the left. For example, take $h=5,a=2$, type I+, $\theta_d=0$, and $\theta>0$. Its delay has $x=e^\theta>1$, so
\[
0<\arg(x^2+e^{4\pi i/5})<2\pi/5=\delta.
\]
Only zero or one sector boundary can be crossed; the formal bound $a/g=2$ is never attained in that left-incidence experiment. Conversely, for I− with $\theta_d=-3$ and $\theta>0$, the phase exceeds $\delta$, so the zero-step endpoint is inaccessible. Both bounds are sharp over the full parameter family, or for the algebraic left/right profile relation with unrestricted rapidities. The physical incidence statement must say which meaning is intended.

Independent enumeration:

```text
(h,a)=(3,1): sectors=3, spacing/pi=.666667, attainable |steps|=[0,1]
(3,2): sectors=3, spacing/pi=.666667, attainable |steps|=[0,1]
(4,2): sectors=2, spacing/pi=1, attainable |steps|=[0,1]
(5,2): sectors=5, spacing/pi=.4, attainable |steps|=[0,1,2]
(6,2): sectors=3, spacing/pi=.666667, attainable |steps|=[0,1]
(8,3): sectors=8, spacing/pi=.25, attainable |steps|=[0,1,2,3]
```

Multiplying $E$ by a scalar cannot leave this orbit. Thus no missing weight is generated by this one-soliton transmission ansatz. This is not a no-go theorem for arbitrary multi-soliton or singular defect processes.

**Paste-ready edits:**

At line 76, replace “and real rapidity $\theta$” by:

```markdown
and positive real rapidity $\theta$, incident from the left, with incoming and outgoing phases off the singular sector boundaries
```

At line 80 replace its last sentence, beginning “Both bounds are attained”, with:

```markdown
The bounds are sharp as the relative rapidity and incoming phase vary over their full ranges. The upper endpoint corresponds to the limit $\vec\lambda_{\rm out}=\pi\vec\lambda_{\rm in}$, approached at large $\theta-\theta_d$ for type I$_-$ and large $\theta_d-\theta$ for type I$_+$. For a fixed $\theta_d$ and incidence restricted to the left, $\theta>0$ can restrict the accessible steps. For $a=h/2$, the outgoing charge either stays fixed or moves one step, away from the absorption pole. For incidence from the right, exchange the incoming and outgoing fields and use the phase shift $-\arg z$.
```

At line 81 replace the last sentence by:

```markdown
One-soliton transmission of this form does not create the missing charges of @sec-missing-charges.
```

## 4. Frozen-soliton comparison: `prp-defect-frozen-soliton`

**Verdict: correct, with real $\theta_c$, finite nonzero defect parameters, and the explicitly stated CZ energy normalization.** The theorem is a factorization and energy identity, not a theorem that the pair is dynamically identical to a soliton in every respect.

### Independent factorization against chapters 10 and 16

Let $u=\theta-\theta_c$, $y=e^{-u}$, $C=\gamma_c$, $B=\gamma_b$, and $\vartheta=\theta_c+i(\pi/2+C)$. Direct simplification gives

\[
z_b(\theta-\vartheta)
=\frac{y-e^{i(B-C)}}{y-e^{-i(B+C)}},
\]
\[
z_b(\bar\vartheta-\theta)
=e^{2iB}\frac{y-e^{i(C-B)}}{y-e^{i(C+B)}}.
\]

Their product is

\[
e^{2iB}\frac{y^2-2y\cos(C-B)+1}{y^2-2y\cos(C+B)+1}
=\omega^b\frac{\cosh u-\cos(C-B)}{\cosh u-\cos(C+B)}
=\omega^b A_{cb}(u),
\]

using exactly `eq-interaction-coefficient` in chapter 10. Continuing the incident soliton rapidity to $\theta+i\pi/2$ in **both** defect factors gives $\cosh(u+i\pi/2)=i\sinh u$. Chapter 16's independent expression

\[
X_{cb}(z)=\frac{z-\cos(C-B)}{z-\cos(C+B)}
\]

then proves the particle statement, with precisely the same $\omega^b$. The I+ factor must have argument $\bar\vartheta-\theta-i\pi/2$; the calculation uses that sign.

The CZ energy–momentum sum is

\[
-\hat M\bigl(\cosh\vartheta+\cosh\bar\vartheta,
\sinh\vartheta+\sinh\bar\vartheta\bigr)
=2\hat M\sin C(\cosh\theta_c,\sinh\theta_c)
=M_c(\cosh\theta_c,\sinh\theta_c).
\]

This is exact. In the book's normalization the vacuum energy of a class-zero defect is zero even at complex parameter. One must actually use the CZ constants to make the stated comparison; the chapter does say so. For the real-parameter charged realization the book's pair has the extra $2\hat M(\cosh\theta_c,\sinh\theta_c)$, as line 191 correctly explains.

### Same-type obstruction without a zero–pole search

For any finite nonzero parameter or class shift, an I− factor tends to $1$ at $u\to-\infty$ and $\omega^b$ at $u\to+\infty$. A product of two I− factors consequently has limiting ratio $\omega^{2b}$. I+ gives the reciprocal ratio. But $\kappa A_{cb}^{\pm1}$ tends to the **same** constant $\kappa$ at both ends. Equality requires

\[
\omega^{2b}=1\quad\Longleftrightarrow\quad b=h/2,
\qquad 1\le b<h.
\]

This also covers repeated zeros and poles without a genericity assumption about their locations. For $b=h/2$, $z_b(-u)=-z_b(u)$; replacing the I+ factor in the displayed mixed pair by I− changes the sign, which cancels $\omega^b=-1$. Thus two I− factors do reproduce $A_{cb}$ in that channel. Shifting the parameters by $i\pi$ inverts them and gives $A_{cb}^{-1}$. There is no same-type pair that works for every species when $h\ge3$.

### Real-parameter realization

The two effective rapidities of real positive $\sigma=e^{-\theta_c}$, both of class $p$, are $\theta_c\pm2\pi ip/h$. They equal the required pair modulo $2\pi i$ iff

\[
4p\equiv h+2c\pmod{4h}.
\]

An integral class exists exactly when $h+2c\equiv0\pmod4$. Allowing negative real $\sigma$ does not enlarge the possible $(h,c)$: for even $h$ its extra $i\pi$ is already the class shift $h/2$, and for odd $h$ the required congruence still has incompatible parity.

Independent output:

```text
h=2..10: soliton factor residual 8.465e-16
particle factor residual 1.643e-14; EP residual 1.599e-14
real-sigma realizations h<=8:
(2,1,1), (4,2,2), (6,1,2), (6,3,3), (6,5,4),
(8,2,3), (8,4,4), (8,6,5)
```

The triples are $(h,c,p)$.

**Paste-ready hypothesis clarification, line 176:**

```markdown
**[Theorem]** Let $1\le b,c<h$, $\theta_c\in\mathbb R$, and $\omega^b=e^{2i\gamma_b}$. Defect parameters are finite and nonzero; in item 2 the displayed rapidities are the effective parameters after any class shifts.
```

With these hypotheses the theorem label is justified by the algebra, including its same-type obstruction. It does not require the reported grid search as evidence.

## 5. Linear transmission, crossing and the CZ2007 misprint

**Verdict: correct.** This includes both types, $R=0$, the analytic continuation of $z$, physical unimodularity at real uncharged parameters, algebraic unitarity, crossing, the class shifts, and the corrected source denominator.

### Independent scattering system

Differentiate the original exponential potential twice at the vacuum. On the Fourier mode $k$, put $d=1+i\cot\gamma_k$, $e=1-d$, and $w=e^{2i\gamma_k}$. For I− the gradients are

\[
B_\phi=m[(\sigma+\sigma^{-1})\phi
-(\sigma w^{-1}+\sigma^{-1})\psi],
\]
\[
B_\psi=m[(\sigma+\sigma^{-1})\psi
-(\sigma w+\sigma^{-1})\phi].
\]

With incident amplitude one, the two sewing equations are

\[
ip(1-R)=-i\nu[e(1+R)+dT]
-m[(\sigma+\sigma^{-1})(1+R)-(\sigma w^{-1}+\sigma^{-1})T],
\]
\[
ipT=-i\nu[d^*(1+R)-eT]
+m[(\sigma+\sigma^{-1})T-(\sigma w+\sigma^{-1})(1+R)].
\]

Here $d^*=1-i\cot\gamma_k$, $\nu=m_k\cosh\theta$, $p=m_k\sinh\theta$. Elimination gives

\[
R=0,\qquad
T_k=\frac{e^{-(\theta-\theta_d)}+e^{i\gamma_k}}
{e^{-(\theta-\theta_d)}+e^{-i\gamma_k}}
=z_k(\theta-\theta_d+i\pi/2).
\]

Changing the kinetic signs and swapping $\sigma$ and $\sigma^{-1}$ gives the printed I+ factor. For incidence from the right, use $\psi=e^{-ipx}+\bar R e^{ipx}$, $\phi=\bar T e^{-ipx}$; solving the same equations yields

\[
\bar R=0,\qquad \bar T_k(\theta)=1/T_k(-\theta)
=T_{h-k}(i\pi-\theta).
\]

Consequently $T_k(\theta)T_{h-k}(\theta+i\pi)=1$. For real $\theta,\theta_d$, the numerator and denominator of each vacuum factor are conjugate, so $|T|=1$. About a charged vacuum the phases in the Hessian implement the effective parameter shift. Algebraic crossing/unitarity survives, but $|T|=1$ requires $\sin(2\pi p/h)=0$.

`independent_checks.py` forms the full ambient-space Hessian from root outer products, solves the overdetermined vector equations for both incidence directions, and checks the result against the formulas. It does not assume the exercise's scalar gradients when forming this system:

```text
h=2..7, both types, all modes and classes:
R / sewing residual 3.238e-13
T formula 1.045e-14
crossing 3.878e-15
algebraic unitarity 1.724e-15
```

### The source's own conventions

I read equations (7.5)–(7.11) in the supplied text of [CZ2007](https://arxiv.org/abs/0705.1066), including the original linearized field equations (7.7). Set $D_{\rm CZ}=D^\mathsf T$, $E_{\rm CZ}=1-D_{\rm CZ}$, and $\rho_k=-\sum_j\omega^{jk}\alpha_j$, and differentiate their potential. Solving gives

\[
T^{\rm CZ}_k=
\frac{ie^{-\theta}m_k-m\sigma(1-\omega^{-k})}
{ie^{-\theta}m_k+m\sigma(1-\omega^k)}.
\]

This derivation confirms both the denominator correction and the missing explicit $m$. Using

\[
1-e^{-2i\gamma}=2ie^{-i\gamma}\sin\gamma,\qquad
1-e^{2i\gamma}=-2ie^{i\gamma}\sin\gamma,
\]

and $m_k=2m\sin\gamma$, it becomes

\[
\frac{x-e^{-i\gamma_k}}{x-e^{i\gamma_k}}
=T^{\rm book}_{h-k},\qquad x=e^{-(\theta-\eta)}.
\]

For $h=3$, converting these two rational functions to sinh ratios gives the two displayed factors (7.5), with their respective phases. This is also checked directly against the source's sinh expressions:

```text
source linearized system vs corrected formula: 2.483e-15
theta=.83, eta=.27, a2 k=1: corrected-(7.5)=1.665e-16; printed |T|=.327943516
theta=.83, eta=.27, a2 k=2: corrected-(7.5)=5.933e-16; printed |T|=.381619961
```

The printed denominator fails substantially. The chapter's correction is justified, not merely inferred from a desired modulus. No substantive replacement is required in this section. Its [Numerical] status is conservative: the scalar elimination above proves the result for every $a_r$.

## 6. Other proposed corrections to the sources

**Verdict: correct mathematical extensions, with a further energy-versus-momentum qualification.** The CZ2009a statement after (4.8) is a sufficient condition, not an asserted “if and only if”; the chapter should call its addition a clarification rather than a false theorem in that source.

### Sector spacing

The explicit boundary set $\{\pi-2\pi ja/h\}$ derived in item 3 has spacing $2\pi g/h$. Equality with $2\gamma_a=2\pi a/h$ holds iff $g=a$, equivalently $a\mid h$. This exactly identifies the scope of the older spacing statement. For $h=5,a=2$, the boundaries are

\[
\{\pi/5,3\pi/5,\pi,7\pi/5,9\pi/5\},
\]

not separated by $4\pi/5$. The independent charge enumeration exhibits the two-step transitions that the older one-step argument misses. CZ2009a's $2\pi/\widetilde n_a$, with $\widetilde n_a=h/g$, agrees. I checked the passages following (2.22) in the 2007 source and (4.10) in the 2009 source directly.

### Real absorption

For class zero, the denominator-zero calculation of item 2 reduces to $a=h/2$ at real rapidity. For charged classes it gives the displayed congruence instead. The explicitly constructed $h=6,a=1,p=1$ pair satisfies the full nonlinear conditions and has the correct real energy increment. Thus the proposed extension beyond root-lattice jumps is justified. It extends the backgrounds treated by the source formula; it does not disprove that formula in its stated class-zero setting.

### Static energies

For the source's (4.8), restore its opposite momentum convention, but take the imaginary parts directly. Simultaneous reality requires

\[
\sin\bigl(2\pi(a-b)/h\bigr)=0.
\]

Thus equal classes and, when $h$ is even, classes differing by $h/2$, both give real energy **and** momentum. For energy alone, $\eta=0$ is an additional possibility for every class difference. The source text says equal-representation weights give real energy and momentum; it does not deny the second possibility. See [CZ2009a](https://arxiv.org/abs/0902.1307), (4.8) and its following paragraph.

**Replace the last bullet, line 236, with:**

```markdown
- [@cz2009a, Sec. 4], after (4.8): in addition to equal classes, classes differing by $h/2$ when $h$ is even also give real energy and momentum. The source states the equal-class sufficient condition without claiming necessity. Energy alone is real for every class when $\theta_d=0$.
```

At line 234 the quoted bound $a/g$ should be stated with its species range, to avoid losing the sharper conjugate-species result already present in the proposition. Replace its last sentence by:

```markdown
A soliton can therefore move up to $\min(a,h-a)/\gcd(a,h)$ phase-sector steps, in the direction specified in @prp-stored-charge; species 2 of $a_4$ can move two steps.
```

## 7. The limit $\theta_d\to-\infty$

**Verdict: correct for the jump-free branch as a classical sewing limit.** In a charged class, the local gluing also contains the vacuum translation and its class-dependent constant energy generally diverges. A statement about convergence of the full action or quantum line operator would require more than the soliton delay limit.

### Independent field-level derivation

For I−, $\sigma\to\infty$. Bounded derivatives in the sewing conditions force the leading gradient of the $f$-part to vanish. Since $D$ is invertible on the root space,

\[
\sum_j\alpha_j f_j=0.
\]

For affine $a_r$, the only linear relation among these $h$ roots has all coefficients equal; therefore all $f_j$ are equal. Also $\prod_jf_j=1$. On the branch connected to the uncharged vacuum, $f_j\to1$, and locally on that branch

\[
D^\mathsf T\phi+D\psi\to0.
\]

The independently verified matrix identity $D^\mathsf T=-D\pi$ then gives $\psi=\pi\phi$. The derivative gluing follows too: writing $A$ for the finite limiting $f$-gradient in the first equation, its contribution to the second is $-\pi A$. Eliminating $A$, using $\psi_t=\pi\phi_t$, and commuting $D,E,\pi$, yields $\psi_x=\pi\phi_x$. Hence the limiting local conditions are those of the invertible diagram-rotation line.

On tau-functions, $z\to\omega^a$ sends $\tau_j$ to $\tau_{j+1}$, and

\[
\frac i\beta\sum_j\alpha_j\log\tau_{j+1}
=\frac i\beta\sum_j\alpha_{j-1}\log\tau_j
=\pi\phi
\]

up to the vacuum shift associated with continued logarithms. On particles the same calculation gives the eigenvalue $\omega^k$. This verifies the orientation of the rotation, not just the existence of some cyclic symmetry.

```text
D*pi=-D^T residual, h=3,4,6: 2.859e-15
z_limit=omega^a algebraically
```

For a fixed class $p$, the branch has $f_j\to\omega^{-p}$, so instead $\psi=\pi\phi+2\pi\lambda/\beta$ modulo root-lattice translations. Its vacuum energy has leading term

\[
B_p\sim\frac{\hat M}{2}e^{-\theta_d}
(1-e^{-2\pi ip/h}),
\]

which diverges unless $p=0$. This does not affect the local sewing equations, but prevents identifying a charged limit with the normalized finite-energy rotation line without a subtraction.

**Replace the sentence at line 54 beginning “Limit 3 says more” with:**

```markdown
On the branch connected to the jump-free vacuum, the limit $\theta_d\to-\infty$ of the type-I$_-$ sewing conditions is the invertible $\mathbb Z_h$ symmetry line, $\vec\psi=\pi\vec\phi$ with the corresponding derivative gluing. In a class-$p$ branch there is also a vacuum translation of class $p$; comparing its energy with that of the rotation line requires subtraction of the class-dependent constant vacuum energy, which generally diverges in this limit.
```

This qualification preserves the intended classical claim without inferring a general operator limit from one family of solutions.

## 8. Reality and advances

**Verdict: correct for the conjugation of delays and the delay/advance criterion; wrong blanket statements about complex energy and imaginary defect rapidity.**

### Independent computation

Conjugation of the fundamental delay gives

\[
\overline{z_a(u)}=z_{h-a}(\bar u),
\]

because $e^{i\gamma_{h-a}}=-e^{-i\gamma_a}$. Applying this to the two class-shifted arguments gives exactly line 137, including $p\mapsto-p$ and $\theta_d\mapsto\bar\theta_d$. A real constant weight jump changes sign under $\Theta$, establishing the same class action directly on fields.

For real incident and defect rapidities put $x=re^{iQ}$, $r>0$, $Q=2\pi p/h$. The numerator and denominator of the soliton delay obey

\[
|x-ie^{i\gamma_a}|^2-|x-ie^{-i\gamma_a}|^2
=4r\sin\gamma_a\cos Q.
\]

The exact script obtains this polynomial identity, not a fit to sampled moduli. Since $\sin\gamma_a>0$, the chapter's criterion follows: delay for $\cos Q>0$, advance for $\cos Q<0$. When $\cos Q=0$, $|z|=1$, so there is a phase/charge effect but no envelope displacement. The real part of the I+ argument changes $r$, not this sign criterion.

```text
Theta delay residual through h=9: 3.202e-15
advanced samples with cos(2*pi*p/h)<0: 216, all |z|<1
```

The odd-sector identity follows by substituting $e^{-(u+i\pi)}=-e^{-u}$: $z_a(u+i\pi)=1/z_{h-a}(u)$. The class $h/2$, when it exists, is fixed by $p\mapsto-p$ and gives real energy and momentum; for odd $h$ the phase $-1$ is not an allowed class multiplier.

The first sentence at line 142 repeats the energy-only error of item 3. Non-fixed classes have a non-real energy–momentum pair, but at $\theta_d=0$ their energy is real. The last sentence at line 144 is also too strong: $\theta_d=i\pi$ is imaginary but $\sigma=-1$ is real, and the defect equations are $\Theta$-symmetric. More generally, directly conjugating the original potential shows that the criterion for the full theory is

\[
\sigma=\bar\sigma
\quad\Longleftrightarrow\quad
\theta_d\in\mathbb R+i\pi\mathbb Z,
\]

not “$\theta_d$ has no imaginary part.” This agrees with the precise proposition in chapter 38. A class is additionally fixed modulo the invisible root-lattice translations iff $2p\equiv0\pmod h$.

**Replace the opening of line 142 with:**

```markdown
Classes $2p\not\equiv0$ mod $h$ are not fixed by $\Theta$. For real $\theta_d$ their energy–momentum pairs are non-real and occur in conjugate classes $p,h-p$; at $\theta_d=0$ their energies alone are real. A soliton of real rapidity is delayed if $\cos(2\pi p/h)>0$, advanced if it is negative, and has no envelope displacement if it vanishes.
```

**Replace the last sentence of line 144 with:**

```markdown
The defect theory is $\Theta$-symmetric precisely when $\sigma$ is real, equivalently $\theta_d\in\mathbb R+i\pi\mathbb Z$. For other complex parameters $\Theta$ maps it to the conjugate-parameter theory. At real $\sigma$, only classes with $2p\equiv0$ mod $h$ are fixed as classes by $\Theta$.
```

For a conjugate breather pair in class zero, the two delays are conjugate, so their product in the centre-of-mass exponential is $|z|^2$. With $\theta_{1,2}=\theta_0\pm iu$, its envelope displacement is

\[
\Delta x=-\frac{\log|z|}{m_a\cos u\cosh\theta_0}.
\]

This follows by applying the multiplier to $E_1E_2$, whose spatial exponent is $2m_a\cos u\cosh\theta_0$. It is real. A change in $\arg z$ also shifts the internal phase; a real centre displacement does not mean the complete breather is changed only by a spatial translation.

## 9. Status labels, exercises and bibliography

**Verdict: mostly correct; the absorption status needs the restrictions in item 2, the selection proposition needs its omitted cases, and the energy exercise must distinguish energy from energy–momentum.** No unresolved chapter references or bibliography keys were found.

### Independent exercise spot checks

1. **`exr-ds-class`, line 196.** The exact weight congruence and the direct potential substitution in items 1 and 3 solve the exercise. Its answer to “For which classes is the energy real?” must include **all classes at $\theta_d=0$**. To prevent the surrounding text from suggesting a false answer, replace that last question with:

   ```markdown
   For which classes are both energy and momentum real? When is the energy alone real, including $\theta_d=0$?
   ```

2. **`exr-ds-linear`, line 200.** The independently differentiated Hessian gives exactly its two gradients. The two scalar sewing equations and their solution are written in item 5. Continuing the incident rapidity to $\theta-i\pi/2$ returns the soliton delay. This exercise is correct.

3. **`exr-ds-time-reversal`, line 204.** Negating the time derivatives, the kinetic matrices and $\theta$, while replacing the I+ parameter by $1/\sigma$, proves the assertion directly. For $h=2$, $z_1(-u)=+\coth(u/2)$. Type I+ reaches the rotation limit when $\theta-\theta_d\to-\infty$. Correct.

4. **`exr-ds-selection`, line 208.** For $h=5,a=2$, the boundary set is the five odd multiples of $\pi/5$. In the ambient basis $\alpha_j=e_j-e_{j+1}$, the successive charges, starting just above $\pi/5$, are

   \[
   \frac15(-2,3,-2,-2,3),\quad
   \frac15(-2,3,-2,3,-2),\quad
   \frac15(3,-2,-2,3,-2),\quad
   \frac15(3,-2,3,-2,-2),\quad
   \frac15(-2,-2,3,-2,3).
   \]

   At I− relative rapidity $u=3$, an incoming phase $0.55$ moves from the last of these sectors to the second, crossing the boundaries $\pi/5$ and $3\pi/5$. The script independently calculates the two charge vectors. Spacing $2\gamma_a$ is equivalent to $a\mid h$; the bound is at most one iff $\min(a,h-a)=g$, equivalently $\min(a,h-a)\mid h$. Correct for the stated one-soliton, left-incidence phase convention.

5. **`exr-ds-absorption`, line 212.** The energy identity and complete real-pole lists for $h=4,6,8$ are given in item 2. Its final question has an important answer: $a+2p\equiv0\pmod h$ gives no localized envelope. This is already a good exercise, but its answer contradicts the unqualified item 1 of the current theorem. The theorem must incorporate it.

6. **`exr-ds-bootstrap`, line 216.** Write $A=\pi a/h$, $B=\pi b/h$, $x=e^{-u}$. Multiplying the two rational delays at $u+iB$ and $u-iA$ cancels an inner factor and leaves

   \[
   \frac{x-ie^{i(A+B)}}{x-ie^{-i(A+B)}}=z_{a+b}(u).
   \]

   The same identity holds with the two shifts reversed; this is useful for I+. Analytic continuation by $i\pi/2$ proves the linear-factor version. Independent maximum residual through $h=8$: $2.809\times10^{-15}$ for either choice of the two opposite shifts. Correct.

7. **`exr-ds-breather`, line 220.** The odd-sector identity and the conjugate-product displacement are derived in item 8. The class-zero hypothesis of the breather assertion should be explicit: real $\sigma$ alone does not fix a general charged class. Replace “a defect with real $\sigma$ preserves” by “a class-zero defect with real $\sigma$ preserves”. More generally, a class fixed by $\Theta$, including $h/2$ for even $h$, preserves the corresponding conjugation relation.

8. **`exr-ds-frozen`, line 224.** At $h=2,c=1$, the mixed pair gives $-A_{11}$. Since $z^{\rm I_+}_{\rm odd}=-z^{\rm I_-}_{\rm odd}$, this is equivalent to $(z^{\rm I_-}_{\rm odd})^2=A_{11}$, exactly as the exercise requests after cancelling the common sign. At $h=4,c=2$, $p=2$, $\hat M=8m/\beta^2$, $M_2=16m/\beta^2$. The pair has CZ energy $M_2\cosh\theta_c$, and book energy $(2\hat M+M_2)\cosh\theta_c$. The $a_2$ obstruction follows from the unequal limiting constants of two I− factors. Correct.

The class-zero breather qualification also belongs at line 140. Replace “passes a defect with real $\sigma$” by “passes a class-zero defect with real $\sigma$”.

### Status labels and several solitons

The delay proposition's “Established; numerically checked” status is supported by the exact rational identity in item 1; numerical agreement alone was not needed. The frozen-soliton theorem has a complete algebraic proof once its real-$\theta_c$ hypothesis is explicit. The absorption theorem requires the corrections above. The linear-transmission [Numerical] label is conservative, since the two-equation derivation works in every rank.

I independently built the two-soliton tau-functions using chapter 10's $A_{ab}$, evaluated the original exponential gradients, and multiplied each incident exponential by its respective delay. For $h=3,4$, both types and classes 0 and 1, the maximum sewing residual was $1.013\times10^{-13}$. This supports line 122. **The general $N$-soliton claim is unclear as a fully supplied proof in this chapter:** the text acknowledges that its $a_r$ Bianchi argument is not written out. My two-soliton computation does not promote it to an independently proved all-$N$ theorem.

### Bibliography and references

I checked the relevant formulas and passages in the supplied local source texts: CZ2007 (2.20)–(2.23), (7.5)–(7.11); CZ2009a (4.5)–(4.10); BCZ2004b (6.2)–(6.3); and BCZ2003 (3.1), (3.5)–(3.6). In the last source the soliton exponent has $+t\sinh\theta$, so the required rapidity reversal really is necessary. Direct substitution into BCZ2004b (6.2) after $\theta\mapsto-\theta$ gives the chapter's delay. CZ2007 (2.20) gives the same formula after $a\mapsto h-a$.

The arXiv records agree with the principal bibliography entries: [CZ2007](https://arxiv.org/abs/0705.1066), [CZ2009a](https://arxiv.org/abs/0902.1307), [BCZ2004b](https://arxiv.org/abs/hep-th/0401020), and [BCZ2004 / 2003 preprint](https://arxiv.org/abs/hep-th/0305022). The chapter's six citation keys exist in `references.bib`. This is a spot check of bibliographic identity and the cited classical passages, not a rederivation of all quantum historical claims.

```text
distinct references/citations=62; unresolved=[]
citation keys=bcz2004b, bcz2005, cz2004, cz2007, cz2009a, cz2010
chapter sha256=b378c80b1a3fc9c8c507e2180aa58ef84e8a56ebce3523273673f4e2d1b6db49
```

Use `metadata_checks.py` for that output. The deliberately simple reference scan at the end of `supplemental_checks.py` retains sentence-final periods and reports nine apparent misses; the punctuation-aware check resolves all of them. They are not broken references in the manuscript.

## Other observations — unverified

- **Unverified:** the identification of every charged complex class with a quantum sector excluded by $\Theta$. The classical conjugation properties checked above do not themselves construct the quantum state space or its physical inner product.
- **Unverified:** global regularity, admissible integration contours, and conserved integrated energies through the real-spacetime singularities of complex-rapidity pole configurations. This report verifies the meromorphic identities and identifies explicit singularities; it does not supply a prescription across them.
- **Unverified:** whether the mixed pair has the auxiliary degrees of freedom and complete dynamics of a particular fused type-II defect. Agreement of delay, linear scattering and a chosen constant energy normalization is insufficient to establish that equivalence.
- **Unverified:** the historical priority claim that all the class-shifted or time-reversed observations are new. The cited local papers were checked where relevant, but this was not an exhaustive search of the later literature.

## Integration (2026-10-13)

The referee's four check scripts were copied to `briefs/referee_defects_39_scripts/` with their file paths made relative to the repository root (`ROOT=Path(__file__).resolve().parents[2]`), so they run from any directory; `format_report.py`, a one-off formatter that rewrites this report, was not copied. All four were rerun from another directory and reproduce every number quoted above (`exact_checks.py` 44 s, `independent_checks.py` 1.7 s, `supplemental_checks.py` 0.2 s, `metadata_checks.py` 0.05 s); their outputs at the start of the integration are saved as `*_output.txt`. The book cites only `textbook_code/`. Independent checks were made in a scratch script of the integrating session that uses only `defect_common.py` and its own formulas (no referee imports): the sewing residual at $2h+4$ arbitrary complex values of $E$ for $a_2$, $a_3$, all species, classes and types at complex $\theta$ ($6\times10^{-14}$; wrong-class control $\ge0.05$); $\cos\chi=0\iff a+2p\equiv0$ mod $h$ for $h\le10$; $|E|\equiv1$ at $h=4,a=2,p=1$; the zero of $\tau_0$ at $(-2.3095,-4.3004)$ for $h=3,a=1$ ($|1+E|=10^{-15}$); energy at $\theta_d=0$; step sets at fixed $\theta_d$ ($h=5,a=2$: $\{0,1\}$ for I$_+$ at $\theta_d=0$, $\{1,2\}$ for I$_-$ at $\theta_d=-3$); $a=h/2$ steps; the limits $1$, $\omega^b$ of each defect factor; the real-$\sigma$ realizations; an exact sympy elimination of the two scalar linear equations in the variables $c=e^{i\gamma_k}$, $e^\theta$, $e^{\theta_d}$ for both types and both directions ($R=\bar R=0$, $T$ and $\bar T(\theta)T(-\theta)=1$ identically); the $\theta_d\to-\infty$ class energy; $|z|=1$ at $\cos\frac{2\pi p}h=0$; the class-1 breather control. The source passages were read: CZ2007 after (2.22) and (2.23), CZ2009a around (4.8) and after (6.15).

- **Item 1. Accepted.** The exact identity (`exact_checks.py`) and the direct-exponential test with its negative control were rerun. Independently: the conditions see $(x,t)$ only through $E$, and the residual times $\prod_j\tau_j\tau'_j$ is a polynomial of degree at most $2h$ in $E$, so its vanishing at $2h+4$ arbitrary complex $E$ (scratch and now the book script) proves the all-orders statement and item 4 without appeal to the Bäcklund property. Changes: ch. 39 `prp-an-defect-delay` statement ("complex $\theta$ and $\xi$"; meromorphic identity, nonzero $\tau$-functions, limiting configurations at a zero or pole of $z$); proof (the degree argument replaces "as they must for a Bäcklund pair"). Script check 1: new degree-argument test.
- **Item 2. Accepted** (in the book's label vocabulary: **[Theorem]** kept, scope written into the statement). All counterexamples reproduced independently. In addition, the book script's check 3 included $(a_2,a=1,p=1)$, which is the zero-envelope case $a+2p\equiv0$ mod $h$: its stored-energy test passed only because $\cos\chi=6\times10^{-17}$ in floating point gave $T=5\times10^{17}$. Changes: ch. 39 `prp-ds-absorption` (real $\theta_d$; meromorphic solutions; velocity only for $a+2p\not\equiv0$; stationary at $\theta_d=0$; the case $a+2p\equiv0$ has no envelope and is neither absorption nor emission; exactly one class per species in item 3), its proof ($\cos\chi=\sin\frac{\pi(a+2p)}h$; the $\tau$-zeros at complex $\theta_{\rm abs}$ with the $a_2$ example, $m=1$; the uniqueness of $p$), and the paragraph after it (energy and momentum of the stored soliton; reality for $\theta_d\ne0$ only for $a=h/2$, energy real for every $a$ at $\theta_d=0$). Script check 3: the zero-envelope case excluded from the energy test and checked separately ($|E|$ constant for $a_2$, $a=1$, $p=1$ and $a_3$, $a=2$, $p=1$; the configurations still solve the conditions); $\cos\chi=0\iff a+2p\equiv0$ mod $h$; the $\tau_0$ zero; one class per species.
- **Item 3. Accepted.** The $\theta_d=0$ counterexample, the sign$(\theta)$ bookkeeping, the incidence restriction and the $a=h/2$ case were reproduced independently. Changes: ch. 39 stored-energy paragraph (algebraic mass shell; imaginary parts; energy and momentum both real iff $2p\equiv0$, energy alone also at $\theta_d=0$; conjugate pairs); "Charge bookkeeping" ($\theta>0$, incidence from the left); `prp-stored-charge` (hypotheses: $\theta>0$ from the left, phases off the boundaries; item 3 with $a=h/2$, attainability over the full family and the fixed-$\theta_d$ restriction; item 4 "one-soliton transmission"); its proof (the bound $K=\min(a,h-a)/g$, attainment at finite rapidity, $a=h/2$, the $a_4$ example); a new paragraph for $\theta<0$ (phase shift $-\arg z$, sign$(\theta)(\vec\lambda_L-\vec\lambda_R)$). Script check 2: energy reality at $\theta_d=0$, the $\theta<0$ jump, step sets at fixed $\theta_d$, $a=h/2$.
- **Item 4. Accepted.** The asymptotic argument was checked by hand and numerically (limits $1$ and $\omega^b$ for finite complex parameters; two type-I$_-$ factors give $A_{c,h/2}$). Changes: ch. 39 `prp-defect-frozen-soliton` hypotheses ($1\le b,c<h$, real $\theta_c$, finite nonzero effective parameters), item 1 (two type-I$_-$ defects do work for $b=h/2$), proof of item 1 replaced by the limit argument (covers coincident zeros and poles), the zero–pole matching kept for item 2, and the grid search described as corroboration. Script check 7: the limits and the $b=h/2$ pair.
- **Item 5. Correct; label strengthened.** The referee calls the [Numerical] label conservative. The exact sympy elimination (scratch, and now in the book script) proves $R=0$, the closed forms and the right-incidence factor for every $\gamma_k$, so for every $r$; the reduction to $\vec\rho_k$ is exact because all matrices commute with $\pi$, and item 4 follows from `eq-ds-class-shift`. Changes: ch. 39 `prp-an-linear-transmission` **[Numerical]**, $r\le4$ → **[Theorem]** for every $r$; proof rewritten (elimination; the vector solver as a check of the reduction); the redundant "by hand" sentence removed. Script check 6: sympy elimination.
- **Item 6. Accepted.** CZ2009a after (4.8) reads "when the two weights … belong to the same representation the energy and momentum will be real", a sufficient condition; by the rules of Appendix D (misprint, error, incomplete list, disputed, clarification) the addition of the $h/2$ classes is not a correction. CZ2007 after (2.23) concerns the delay (2.20), so the root-lattice restriction is a hypothesis that makes the statement right: a *Clarification*. Changes: ch. 39 Notes (spacing bullet: $\min(a,h-a)/\gcd(a,h)$ steps in the direction of `prp-stored-charge`; absorption bullet: "the setting of their delay"; the (4.8) bullet removed from the corrections and replaced by a sentence in the Notes paragraph, including $\theta_d=0$); Appendix D (CZ2007 (2.23) entry relabelled *Clarification*; CZ2009a (4.8) entry deleted; the known-correction entry on CZ2007 (2.22) now says $\min(a,h-a)/\gcd(a,h)$).
- **Item 7. Accepted.** The class-$p$ vacuum energy was checked independently ($\frac{\hat M}2e^{-\theta_d}(1-e^{-2\pi ip/h})$; the book script now checks it to $4\times10^{-11}$ relative at $\theta_d=-12$). Changes: ch. 39 paragraph after `prp-an-defect-delay` (limit on the class-0 branch, the $f_j$ argument with $\mathsf D^{\mathsf T}=-\mathsf D\pi$, the class-$p$ divergence and subtraction, classical statement only). Script check 3: the divergence.
- **Item 8. Accepted.** Changes: ch. 39 `sec-defect-reality` (breather: class 0 or $h/2$, no other class; "a defect with real $\sigma$ in a class fixed by $\Theta$"; charged classes: non-real energy–momentum pairs, energy alone real at $\theta_d=0$, $|z|=1$ when $\cos\frac{2\pi p}h=0$ with the modulus identity; $\Theta$-symmetry iff $\theta_d\in\mathbb R+i\pi\mathbb Z$, classes fixed iff $2p\equiv0$). Minimal consistency edits: ch. 38 (paragraph after `prp-defect-theta`: imaginary $\theta_d\notin i\pi\mathbb Z$ breaks $\Theta$); ch. 43 `sec-defect-spectrum` (two sentences) (soliton-pole energy complex for real $\theta_d\ne0$ except at real absorptions; imaginary $\theta_d\notin i\pi\mathbb Z$ makes $\sigma$ a non-real phase). Script check 5: $|z|=1$ at $\cos=0$ (previously skipped), the class-1 breather control, $\Theta$-symmetry at $\theta_d=i\pi$.
- **Item 9. Accepted.** Changes: `exr-ds-class` (both real; energy alone, including $\theta_d=0$), `exr-ds-breather` (class-zero defect; class 1), `sec-defect-multi` (beyond two solitons the statement is now called unproved). `exr-ds-absorption` needed no change; its final question is now answered by item 1 of the theorem.

**Other observations — unverified.** None was resolved here; all four are recorded as open items in `briefs/part8_defects_status.md` (ch. 39). The book now calls the complex-rapidity configurations meromorphic and states no prescription across their singularities; the quantum reading of the charged classes stays a sentence of analogy; the type-II question was already open in ch. 43. While checking item 8 the integrator noticed that for $\theta_d\in\mathbb R+i\frac\pi h\mathbb Z$ a coweight translation composed with $\Theta$ maps the $\sigma$-theory to itself; this is recorded as an open question, not stated in the book.

**Book script after integration:** `defect_solitons_an.py` 9 PASS (48 s on a shared machine), with the additions listed above. Final check over all 64 book `.qmd` files: 1149 ids, no duplicates, every `@sec/eq/prp/thm/cnj/tbl/fig/exr` reference resolves, every citation key is in `references.bib`.
