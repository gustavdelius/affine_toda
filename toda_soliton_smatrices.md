# Exact S-matrices for affine Toda solitons: $c_n^{(1)}$, $a_{2n-1}^{(2)}$ and $e_6^{(2)}$

*7 October 2026*

## Abstract

We construct exact soliton S-matrices for the imaginary-coupling affine Toda theories based on $c_n^{(1)}$, $a_{2n-1}^{(2)}$ and $e_6^{(2)}$. Among the non-simply-laced theories, these and $f_4^{(1)}$ were the cases still missing from the literature.

Each amplitude has the form $S=c\,f\,\check R$. Here $\check R$ is the R-matrix of the dual quantum affine algebra, $U_q(d_{n+1}^{(2)})$, $U_q(b_n^{(1)})$ or $U_q(f_4^{(1)})$, in its smallest representation, and $c\,f$ is a Gamma-function scalar factor whose zeros sit at the R-matrix poles. Crossing and unitarity fix those zeros only modulo 1. We fix the remaining integer lifts by requiring the soliton masses to be exactly proportional to the particle masses.

The resulting amplitudes are unitary and crossing symmetric, and every bound-state pole carries the predicted multiplet. In each case the lowest breathers scatter exactly as the conjectured real-coupling particle S-matrices of Delius–Grisaru–Zanon and Corrigan–Dorey–Sasaki, continued to imaginary coupling, with no free parameter. The soliton masses follow in closed form: $M_a=2M_n\sin(a\pi/H)$, with the same floating Coxeter number $H$ as the particles. The construction also reproduces the conjectured coupling dependence $B(\beta)$ of the particle S-matrices exactly, which those authors checked only to one loop.

For $c_3^{(1)}$ we give all six soliton–soliton amplitudes and the lowest-breather amplitudes in closed form. For $f_4^{(1)}$ we construct the $U_q(e_6^{(2)})$ R-matrix and a soliton amplitude on the $\mathbf{27}$ that satisfies every axiom, including its self-fusion bootstrap. A one-loop computation from the Lagrangian rules this amplitude out: perturbation theory confirms the Corrigan–Dorey–Sasaki relation $H=12+3B$, whereas the $\mathbf{27}$ amplitude implies $H=12+\tfrac92 B$. The same computation confirms our $e_6^{(2)}$ result. The one-loop data also fix the integer lifts of the $\mathbf{27}$ amplitude and the first-order form of the correction its scalar factor still needs. As a check, the same method reproduces Takács's $d_4^{(3)}$ and $g_2^{(1)}$ results.

## 1 Introduction

Affine Toda field theory with imaginary coupling has topologically charged solitons with real masses [1]. For a Kac–Moody algebra $\hat{\mathfrak g}$, the theory carries non-local charges that generate the quantum affine algebra $U_q(\hat{\mathfrak g}^\vee)$ of the dual algebra. Solitons therefore come in finite-dimensional $U_q(\hat{\mathfrak g}^\vee)$ multiplets, and their S-matrix is an R-matrix of $U_q(\hat{\mathfrak g}^\vee)$ multiplied by a scalar factor [2, 3].

This programme has been carried out for $a_n^{(1)}$ [1], $d_{n+1}^{(2)}$ [4], $b_n^{(1)}$ [5], and for the pair $d_4^{(3)}$, $g_2^{(1)}$ [6, 7]. It had not been done for $c_n^{(1)}$, $a_{2n-1}^{(2)}$, $f_4^{(1)}$ and $e_6^{(2)}$. Two obstacles stood in the way.

- **Multiplicities.** The tensor-product-graph method [8] needs multiplicity-free tensor products, and the relevant products fail this for $c_n^{(1)}$ with $n\ge3$.
- **Unusual representations.** The twisted and exceptional algebras need representations that had to be built from scratch.

We construct representations directly from the defining relations. We obtain R-matrices from the intertwining equations, reducing to highest-weight vectors when the space is large, and use exact block arithmetic for everything downstream.

The main results are:

1. **$c_n^{(1)}$.** The spinor-soliton S-matrix for all $n$, built on a closed-form $U_q(d_{n+1}^{(2)})$ spinor R-matrix and verified for $n=2,3,4$. For $c_3^{(1)}$, all six soliton amplitudes are obtained by fusion (Sections 3–4).
2. **$a_{2n-1}^{(2)}$.** The spinor-soliton S-matrix from a closed-form $U_q(b_n^{(1)})$ spinor R-matrix, verified for $n=3,4$ (Section 5).
3. **$e_6^{(2)}$.** The S-matrix of the 26-dimensional soliton, from $U_q(f_4^{(1)})$ (Section 6).
4. **Breather–particle identification.** In every case above, the lowest-breather amplitudes coincide exactly with the conjectured particle S-matrices [9, 10] continued to imaginary coupling. The soliton masses are $M_a=2M_n\sin(a\pi/H)$, with the particles' floating Coxeter number $H$.
5. **$f_4^{(1)}$.** The $U_q(e_6^{(2)})$ R-matrix on $\mathbf{27}=\mathbf{26}\oplus\mathbf 1$ and a soliton amplitude that satisfies every axiom, but which one-loop perturbation theory rules out (Section 7).
6. **Check.** The same method independently reproduces Takács's $d_4^{(3)}$ and $g_2^{(1)}$ results (Section 8).

Two companion notes bear on these results. The first resolves the one-loop disagreement for $c_n^{(1)}$ soliton masses between Delius–Grisaru [11] and MacKay–Watts [12]. It traces the disagreement to a Jordan block in the fluctuation operator, and finds the Delius–Grisaru value correct. Our exact $c_3^{(1)}$ masses reproduce that value. The second note establishes the existence and Krein-space structure of the lattice theory, the setting in which a non-Hermitian S-matrix of this kind is meaningful.

## 2 The construction

Every amplitude in this paper is built in the same five steps. The case-specific inputs are the representation, the R-matrix and the integer lifts of the scalar factor's zeros.

### 2.1 Amplitude, gradation and rapidity variable

A soliton multiplet is a finite-dimensional representation $V(x)$ of $U_q(\hat{\mathfrak g}^\vee)$, with $x$ the spectral parameter. The two-soliton amplitude is

$$
S_{ab}(\theta)=\Phi_{ab}(\theta)\,\check R_{ab}\big(x(\theta)\big),\qquad x(\theta)=e^{2T\theta}.
$$

Here $\check R_{ab}$ maps $V_a(x)\otimes V_b(1)$ to $V_b(1)\otimes V_a(x)$ and is normalised to 1 on the highest-weight channel. Gandenberger, MacKay and Watts [5] fix the gradation from the Lorentz spins of the non-local charges:

$$
2T=\frac{4\pi h}{\beta^2}-h^\vee .
$$

In this formula $h$ and $h^\vee$ are the Coxeter and dual Coxeter numbers of $\hat{\mathfrak g}$. We write $q=e^{-i\pi\omega}$ and use the rapidity variable $t=T\theta/(i\pi)$, so that $x=e^{2\pi i t}$ and the physical strip is $0<t<T$. A pole of $\check R$ at $x=\varepsilon q^{-l}$ then corresponds to positions $t\equiv l\omega/2+\arg(\varepsilon)/2\pi \pmod 1$.

### 2.2 Unitarity and crossing

With the normalisation above, $\check R(x)\check R(1/x)=1$, so unitarity requires $\Phi(\theta)\Phi(-\theta)=1$. The dual representation satisfies $V^*(x)\cong V(x\,x_c)$ for a unique crossing point $x_c$, which we compute numerically in each case. The universal R-matrix then gives

$$
\big(R(x)^{-1}\big)^{t_1}=c_u(x)\,(C^{-1}\otimes 1)\,R(x\,x_c)\,(C\otimes 1),
$$

with $c_u$ rational in $x$. Crossing requires $x(i\pi)=x_c$, which links $T$ to $\omega$, and $\Phi(i\pi-\theta)=\Phi(\theta)\,c_u\big(1/x(\theta)\big)$.

### 2.3 The scalar factor

Following [5], we take $\Phi=F=c\,f$, where $c$ vanishes at every pole of $\check R$:

$$
c(t)=\prod_j\sin\pi(t-a_j),\qquad f^{(1)}(t)=\pi^{-N}\prod_j\Gamma(t-a_j)\,\Gamma(1+t+a_j),
$$

$$
f(t)=\prod_{k\ge1}\frac{f^{(1)}\big(t+2T(k-1)\big)\,f^{(1)}\big(-t+2T(k-\tfrac12)\big)}{f^{(1)}\big(t+2T(k-\tfrac12)\big)\,f^{(1)}\big(-t+2Tk\big)} .
$$

The product makes $f(T-t)=f(t)$ and $f(t)f(-t)=1/\big(c(t)c(-t)\big)$, so $F$ is unitary and carries the crossing factor. Then $c\,\check R$ has no poles at all. The only physical-strip poles of $S$ come from the Gamma functions, at $t=a_j-k$ with $k=0,1,2,\dots$, and their crossed images. At $t=a_j$ the residue is $\check R$'s residue at the corresponding R-matrix pole, so it projects onto the bound-state multiplet.

### 2.4 The integer lifts

Crossing fixes each $a_j$ only modulo 1. Shifting $a_j$ to $a_j+1$ multiplies $F$ by the crossing-symmetric unitary CDD factor $(a_j+1)(T-a_j-1)$, in the block notation of Section 2.5. This moves the top of the $j$-th pole tower by one step, which changes the corresponding bound-state mass at order $\beta^2$.

We fix the lifts by one principle: the fusion angles must make the soliton mass ratios equal, exactly, to the particle mass ratios of the same theory. For $c_n^{(1)}$ and $a_{2n-1}^{(2)}$ this is equivalent to $T-a_j$ being an integer multiple of $\omega+\tfrac12$.

### 2.5 The breather test

The singlet channel of the soliton $n$ with itself has a tower of breather poles, the lowest at $t=T-1$. Its residue is rank one. Because a breather is a singlet, its scattering with a soliton is a scalar times the identity. Breather–breather amplitudes then follow by multiplying two shifted breather–soliton scalars, with no further numerics. All of these are finite products of the blocks

$$
(y)=\frac{\sinh\!\left(\frac{\theta}{2}+\frac{i\pi y}{2T}\right)}{\sinh\!\left(\frac{\theta}{2}-\frac{i\pi y}{2T}\right)},\qquad (y)(-y)=1,\quad (y+2T)=(y),\quad (T)=-1 .
$$

We compare these products, in exact arithmetic in $\omega$, with the Delius–Grisaru–Zanon (DGZ) [9] and Corrigan–Dorey–Sasaki (CDS) [10] particle S-matrices. A DGZ/CDS block $(x)_H$ corresponds to our block $(y)$ with $y/T=x/H$.

### 2.6 The crossing sign

Our numerical crossing test compares $S(i\pi-\theta)$ with $(C\otimes 1)\,S_{21}(\theta)^{t_1}\,(C^{-1}\otimes 1)$. For the $c_n^{(1)}$ spinor it returns $(-1)^n$. For $n=1$ (sine-Gordon, whose S-matrix is known exactly) and for $n=2$ (where the S-matrix of [5] is known), the physical relation holds with the same $c\,f$. So the test form differs from physical crossing by the convention factor $(-1)^n$, and no extra CDD factor is needed. The exact breather identities below confirm this for $n=3$ and $n=4$. A derivation of the convention factor from the charge-conjugation structure is still owed.

## 3 $c_n^{(1)}$ solitons

The spinor soliton of $c_n^{(1)}$ has an S-matrix with a closed-form R-matrix and exact masses for all $n$. Its lowest breather reproduces the DGZ particle amplitude exactly for $n=2,3,4$.

### 3.1 The spinor representation and its R-matrix

The symmetry algebra is the twisted $U_q(d_{n+1}^{(2)})$, with horizontal subalgebra $U_q(b_n)$. On the $2^n$-dimensional spinor we use spin operators on $n$ sites, with long simple roots ($q_i=q^2$) for $i<n$ and short roots ($q_i=q$) for $i=n$ and $i=0$:

$$
e_i=\sigma_i^+\sigma_{i+1}^-\ (i<n),\qquad e_n=\sigma_n^+,\qquad e_0=x\,\sigma_1^-,\qquad k_i=q^{2(\alpha_i,\,\mathrm{wt})}.
$$

Under $U_q(b_n)$, spinor $\otimes$ spinor decomposes into the exterior powers $\Lambda^0,\dots,\Lambda^n$ of the vector representation, without multiplicities. Write $\langle l\rangle_s=(x-s\,q^l)/(1-s\,x\,q^l)$, which satisfies $\langle l\rangle_s(x)\,\langle l\rangle_s(1/x)=1$. Then

$$
\check R_{nn}(x)=\sum_{k=0}^{n}\Big[\prod_{j=1}^{k}\langle 2j\rangle_{(-1)^{j+1}}\Big]P_{\Lambda^{n-k}} .
$$

We verified this against direct solutions of the intertwining equations for $n=2,3,4$, to $10^{-13}$. The crossing point is $x_c=(-1)^{n+1}q^{-2n}$.

### 3.2 Parameters and the scalar factor

For $c_n^{(1)}$, $h=2n$ and $h^\vee=n+1$. Crossing then fixes

$$
q=-e^{-4\pi^2 i/\beta^2},\qquad \omega=\frac{4\pi}{\beta^2}-1,\qquad T=n\omega+\frac{n-1}{2} .
$$

The zeros of $c$ sit at $a_j=j\omega+(j-1)/2$ for $j=1,\dots,n$, so $a_n=T$. These are the natural lifts of the R-matrix poles at $x=(-1)^{j-1}q^{-2j}$. Section 3.4 shows that the lifts matter.

### 3.3 Bound states and exact masses

The pole of $S$ at $t=a_j$ is the bound state $n+n\to n-j$. Its residue projects onto the Kirillov–Reshetikhin module formed by the channels that contain $\langle 2j\rangle$. For $n=4$ the residue ranks are $130=84+36+9+1$, $46=36+9+1$ and $10=9+1$. The breathers sit at $t=T-k$, with rank-one residues.

Since $T-a_{n-a}=a(\omega+\tfrac12)$, the fusion angles give the soliton masses in closed form:

$$
\frac{M_a}{M_n}=2\sin\frac{a\pi}{H},\qquad H=\frac{2T}{\omega+\frac12}=2n+B,\qquad B=-\frac{1}{\omega+\frac12}=-\frac{2\beta^2}{8\pi-\beta^2},\qquad a=1,\dots,n-1 .
$$

This is the classical $c_n^{(1)}$ soliton spectrum with $h$ replaced by $H$. For $n=3$ it reproduces the Delius–Grisaru one-loop corrections [11], since $H=6+\beta^2/2\pi+O(\beta^4)$. This confirms their value in the dispute with [12].

### 3.4 Breather–particle identification

The lowest breathers have mass ratios $m_a/m_n=\sin(a\pi/H)/\sin(n\pi/H)$, exactly the DGZ $c_n^{(1)}$ particle masses [9]. The lowest-breather self-amplitude $S[B^n,B^n]$ is identical, block by block and including the sign, to the DGZ amplitude

$$
S_{nn}=\prod_{p=1}^{n}\{2p-1\}_H\,\{H-2p+1\}_H,\qquad \{x\}_H=\frac{(x-1)_H\,(x+1)_H}{(x-1+B)_H\,(x+1-B)_H},
$$

continued to imaginary coupling. We checked this for $n=2,3,4$ at two couplings, $\omega=2.37$ and $1.61$. The value of $B$ is DGZ's $B(\beta)=\dfrac{\beta^2/2\pi}{1+\beta^2/4\pi}$ evaluated at $\beta_{\rm DGZ}^2=-\beta^2/2$: the sign is the continuation, and the factor $\tfrac12$ is the relative root normalisation.

The lifts are physical. At $n=4$, an alternating choice with $a_3=3\omega$ also gives a unitary, crossing-symmetric amplitude. Its breather amplitude, however, differs from DGZ's by the CDD-type factor

$$
\frac{(2-2B)_H\,(2+2B)_H\,(H-2-2B)_H\,(H-2+2B)_H}{(2)_H^2\,(H-2)_H^2}=1+O(B^2),
$$

so the two choices disagree first at one loop, in the breather masses as well as the amplitude. The rule $a_j=j\omega+(j-1)/2$ removes the discrepancy; for $n\le3$ the two choices coincide.

## 4 The complete $c_3^{(1)}$ S-matrix

For $c_3^{(1)}$ we obtain all six soliton–soliton amplitudes by fusion from $S_{33}$. The lowest-breather amplitudes are identical to all six DGZ particle amplitudes. Here $\omega=4\pi/\beta^2-1$, $T=3\omega+1$, and we write $u_1=\pi(2\omega+\tfrac12)/T$ and $u_2=\pi\omega/T$.

### 4.1 Soliton–soliton amplitudes

Solitons 1 and 2 are the bound states $3+3\to1$ and $3+3\to2$. Every amplitude takes the form

$$
S_{ab}(\theta)=k_{ab}(x)\prod_{\alpha}F_{33}(\theta+i\alpha)\;\check R_{ab}(x),
$$

where $\check R_{ab}$ is the unique intertwiner normalised on the top channel, and $k_{ab}$ is a rational fusion factor, verified exactly with constant $\pm1$.

| $ab$ | shifts $\alpha$ | $k_{ab}(x)$ |
| --- | --- | --- |
| 33 | $0$ | $1$ |
| 31 | $\pm u_1/2$ | $\langle0\rangle_i\,\langle2\rangle_{-i}$ |
| 32 | $\pm u_2/2$ | $-\langle1\rangle_+$ |
| 11 | $\pm u_1,\ 0,\ 0$ | $-\langle2\rangle_+\,\langle4\rangle_-/\langle2\rangle_-$ |
| 22 | $\pm u_2,\ 0,\ 0$ | $-\langle2\rangle_+$ |
| 12 | $\pm(u_1+u_2)/2,\ \pm(u_1-u_2)/2$ | $\langle1\rangle_i\,\langle3\rangle_{-i}$ |

The reversed amplitudes follow from $k_{ba}(x)=1/k_{ab}(1/x)$. The fused amplitudes close on their multiplets to $10^{-6}$ and are unitary to the truncation accuracy of the infinite product.

### 4.2 Spectrum

We tested every candidate pole, taken from the shifted Gamma-function poles of $F_{33}$. Each fundamental fusion is a simple pole with the correct multiplet, and occurs in exactly the channels the fusion rules allow. All were confirmed at $\omega=2.37$ and $\omega=1.61$.

| Fusion | $t$ | Residue |
| --- | --- | --- |
| $3+3\to2$ | $\omega$ | 29-plet |
| $3+3\to1$ | $2\omega+\tfrac12$ | 8-plet |
| $1+1\to2$ | $\omega+\tfrac12$ | 29-plet |
| $1+2\to1$ | $\tfrac52\omega+\tfrac34$ | 8-plet |
| $1+3\to3$ | $2\omega+\tfrac34$ | 8-plet |
| $2+3\to3$ | $\tfrac52\omega+1$ | 8-plet |

There is no pole for $2+2\to2$. Two excited solitons appear, each classically degenerate with its ground state.

- **$2^*$.** It appears in $1+1$ at $t=\omega$ as a simple pole and in $2+2$ as a triple pole, with $M_{2^*}=2M_1\cos(\pi\omega/2T)$. Its excitation energy is close to $m/2$, the one-loop bound-mode frequency.
- **$3^*$.** It is a spinor 8-plet in both $3+1$ and $3+2$, of mass $1.1687\,M_3$ at $\omega=2.37$.

The breathers form complete towers: integer $p$ for soliton 3, and half-integer $p$ for solitons 1 and 2. Even-order poles follow the usual affine Toda pattern of anomalous thresholds. Some single-channel tower poles, and three light 29-plet poles in $S_{22}$, still need a Coleman–Thun analysis.

### 4.3 Breather amplitudes

Breather–soliton amplitudes are finite block products, with blocks in the notation of Section 2.5. The three amplitudes $S[B^a,3]$ were identified numerically. All the others follow from them by exact bootstrap arithmetic, and five of the derived forms were confirmed against independent fused numerics. Every residual scales exactly as $1/J$ under truncation of the infinite product, extrapolating to about $10^{-5}$.

| Amplitude | Blocks |
| --- | --- |
| $S[B^3,3]$ | $(\tfrac{\omega}{2}+\tfrac12)(\tfrac{3\omega}{2})(\tfrac{3\omega}{2}+1)(\tfrac{5\omega}{2}+\tfrac12)(-\tfrac{\omega}{2})(-\tfrac{5\omega}{2}-1)$ |
| $S[B^1,3]$ | $-(\tfrac{\omega}{2}+\tfrac12)(\tfrac{5\omega}{2}+\tfrac12)$ |
| $S[B^2,3]$ | $(\tfrac14)(\omega+\tfrac34)(2\omega+\tfrac14)(3\omega+\tfrac34)$ |
| $S[B^1,2]=S[B^2,1]$ | $(\tfrac12)(\omega+\tfrac12)(2\omega+\tfrac12)(3\omega+\tfrac12)$ |
| $S[B^1,B^1]$ | $-(\tfrac12)(\omega+\tfrac12)(2\omega+\tfrac12)(3\omega+\tfrac12)(-\omega-1)(-2\omega)$ |
| $S[B^1,B^3]$ | $(\omega+\tfrac12)(2\omega+\tfrac12)(-\omega+\tfrac12)(-2\omega-\tfrac32)$ |
| $S[B^3,B^3]$ | $-(1)(\omega)(\omega+\tfrac12)(2\omega+\tfrac12)(2\omega+1)(3\omega)(-\omega+\tfrac12)(-\omega-1)(-2\omega)(-2\omega-\tfrac32)$ |

All six breather–breather amplitudes, including $S[B^1,B^2]$, $S[B^2,B^2]$ and $S[B^2,B^3]$ not listed here, are identical to the DGZ $c_3^{(1)}$ particle amplitudes. The identification is $H=2T/(\omega+\tfrac12)=6+B$ with $B=-1/(\omega+\tfrac12)$, and the comparison is done in exact arithmetic in $\omega$. Among the lowest breathers, which are the particles, we find the fusions $1+1\to2$, $1+2\to1$, $1+2\to3$, $1+3\to2$ and $2+3\to1$. The remaining poles of the breather amplitudes are DGZ's own displaced and split double poles. These are inherited from the particle S-matrix, not artefacts of the soliton construction.

## 5 $a_{2n-1}^{(2)}$ solitons

The spinor soliton of $a_{2n-1}^{(2)}$ is governed by $U_q(b_n^{(1)})$. Its S-matrix has the same structure as for $c_n^{(1)}$, and its lowest breather reproduces the DGZ $a_{2n-1}^{(2)}$ particle amplitude exactly for $n=3,4$.

### 5.1 Representation and R-matrix

The finite part is the $U_q(b_n)$ spinor of Section 3.1. Only node 0 changes: it is now the long root $-(\varepsilon_1+\varepsilon_2)$.

$$
e_0=x\,\sigma_1^-\sigma_2^-,\qquad f_0=x^{-1}\sigma_1^+\sigma_2^+,\qquad k_0=q^{-2(w_1+w_2)},\qquad q_0=q^2 .
$$

All relations, including the $q$-Serre relations, hold to $10^{-16}$. The R-matrix is again multiplicity-free, with $\check R=\sum_k\rho_k P_{\Lambda^{n-k}}$ and

$$
\rho_k(x)=\prod_{i=1}^{\lceil k/2\rceil}\langle 4k-8i+6\rangle_+ ,
$$

so that for $n=3$ the four eigenvalues are $1$, $\langle2\rangle$, $\langle6\rangle$ and $\langle2\rangle\langle10\rangle$. We found this formula at $n=3$ and confirmed it at $n=4$ against an independent solution, to $10^{-12}$. The crossing point is $x_c=q^{-2(2n-1)}$, where $2n-1$ is the dual Coxeter number of $b_n$.

### 5.2 Parameters, bound states and masses

For $a_{2n-1}^{(2)}$, $h=2n-1$ and $h^\vee=2n$. The scalar factor has zeros at the lifts

$$
a_j=(2j-1)\,\omega+(j-1),\qquad T=a_n=(2n-1)\,\omega+n-1,\qquad \omega=\frac{2\pi}{\beta^2}-1,\qquad q=-e^{-2\pi^2 i/\beta^2} .
$$

This $T$ coincides exactly with the gradation of Section 2.1. Again $T-a_j=2(n-j)(\omega+\tfrac12)$, and the other lift choices we tested fail the breather test.

The bound state $n+n\to n-j$ sits at $t=a_j$. At $n=4$ the residue ranks are $93=84+9$ (the Kirillov–Reshetikhin module $W^{(3)}$), $37=36+1$ and $9$ (the vector). The masses are

$$
\frac{M_a}{M_n}=2\sin\frac{a\pi}{H},\qquad H=\frac{T}{\omega+\frac12}=2n-1+\frac{B}{2},\qquad B=-\frac{2\beta^2}{4\pi-\beta^2} .
$$

For this family the soliton masses are exactly proportional to the DGZ particle masses, $m_a/m_n=2\sin(a\pi/H)$.

### 5.3 Breather–particle identification

The lowest-breather amplitude $S[B^n,B^n]$ is identical to the DGZ $a_{2n-1}^{(2)}$ particle amplitude

$$
S_{nn}=\prod_{a=0}^{n-1}\frac{(2a)_H\,(H-2a)_H}{(2a+B)_H\,(H-2a-B)_H},
$$

for $n=3$ ($a_5^{(2)}$) and $n=4$ ($a_7^{(2)}$). Here $B$ is DGZ's $B(\beta)$ at $\beta_{\rm DGZ}^2=-\beta^2$, with no normalisation factor for this family. The crossing test gives $-1$ at $n=3$, as for $c_3^{(1)}$; the exact identity confirms that no extra CDD factor is needed.

## 6 $e_6^{(2)}$ solitons

The 26-dimensional $e_6^{(2)}$ soliton transforms under $U_q(f_4^{(1)})$. Its scalar factor is fixed in advance by the CDS mass pattern, and its lowest breather then reproduces the CDS particle amplitude $S_{11}$ exactly.

### 6.1 Representation and R-matrix

The $\mathbf{26}$ of $f_4$ is quasi-minuscule: its weights are the 24 short roots, together with a 2-dimensional zero-weight space. We solved for the generators $e_i$, $f_i$ from all the defining relations at unimodular $q$, including the $q$-Serre relations and the affine node $\alpha_0=-\theta$ attached to node 1. The residual is $2\times10^{-16}$. No extra singlet is needed, as for the $b_n^{(1)}$ vector.

The full intertwiner system on 676 dimensions exceeds memory. We therefore work with highest-weight vectors: $\check R$ acts on each component of

$$
\mathbf{26}\otimes\mathbf{26}=\mathbf{324}\oplus\mathbf{273}\oplus\mathbf{52}\oplus\mathbf{26}\oplus\mathbf 1
$$

as a scalar, and the $e_0$ intertwining relation becomes a linear system for the five eigenvalues, with residual $10^{-15}$. The eigenvalues are exact:

| Channel | $\mathbf{324}$ | $\mathbf{273}$ | $\mathbf{52}$ | $\mathbf{26}$ | $\mathbf 1$ |
| --- | --- | --- | --- | --- | --- |
| Eigenvalue | $1$ | $\langle2\rangle$ | $\langle8\rangle$ | $\langle2\rangle\langle12\rangle$ | $\langle8\rangle\langle18\rangle$ |

The crossing point is $x_c=q^{-18}$, where 18 is twice the dual Coxeter number of $f_4$. The bound states $\mathbf{26}+\mathbf{26}\to\mathbf{299}\ (=\mathbf{273}+\mathbf{26})$, $\mathbf{53}\ (=\mathbf{52}+\mathbf1)$ and $\mathbf{26}$ sit at $x=q^{-2}$, $q^{-8}$ and $q^{-12}$, and the breathers sit at $q^{-18}$.

### 6.2 Fixing the zeros by the mass pattern

The CDS masses [10] are $m_1=\sin(\pi/H)\sin(2\pi/H')$, $m_2=\sin(3\pi/H)\sin(\pi/H')$ and $m_3=\sin(2\pi/H)\sin(2\pi/H')$, with $1/H+1/H'=1/6$. Classically the soliton mass ratios $2\cos(\pi/18)$ and $2\cos(4\pi/18)$ equal $m_3/m_1$ and $m_2/m_1$ at $H=18$. The constraint on $H'$ gives the identity

$$
\frac{m_2}{m_1}=2\cos\Big(\frac{\pi}{H}+\frac{\pi}{6}\Big),\qquad \frac{m_3}{m_1}=2\cos\frac{\pi}{H} .
$$

Requiring the solitons $\mathbf{26}$, $\mathbf{53}$ and $\mathbf{299}$ to carry the masses of particles 1, 2 and 3 exactly places the fusions at $t=2T/H+T/3$ for $\mathbf{53}$, $t=2T/H$ for $\mathbf{299}$, and $t=2T/3$ for the self-fusion. The last is a coupling-independent angle $2\pi/3$, matching CDS's self-coupling 111. Combined with the gradation ($h=9$, $h^\vee=12$, $\omega=2\pi/\beta^2-1$), this gives

$$
a\in\{\omega,\ 4\omega+1,\ 6\omega+2\},\qquad T=9\omega+3,\qquad H=\frac{2T}{\omega}=18+\frac{6}{\omega} .
$$

This $H$ lies beyond the $e_6^{(2)}$ end of the CDS floating range $12\le H\le18$, as an imaginary-coupling continuation should.

### 6.3 Result

At $\omega=2.37$ the amplitude is unitary, and its residue ranks are 299, 53, 26 and 1 at the predicted positions. The lowest-breather amplitude is identical to the CDS $e_6^{(2)}$ amplitude

$$
S_{11}=[1]_0\,[H/3+1]_0,\qquad [x]_0=\{x\}_0\{H-x\}_0,\qquad B=\frac{H-12}{3},
$$

at the predicted $H=18+6/\omega$, with no adjustment. The lowest breather of the $\mathbf{26}$ is therefore CDS particle 1. A one-loop computation from the $e_6^{(2)}$ Lagrangian gives $\delta H/\tilde B_{\rm tree}=-3.000000$, exactly CDS's $H=18-3\tilde B$ at this end of their range (Appendix A), so this S-matrix is also consistent with perturbation theory. The same holds when the breather amplitude is compared directly with the Lagrangian, without reference to CDS. Its first-order term has exactly the tree-level $\theta$-dependence at eight rapidities, and that term fixes the coupling map $1/\omega=\beta^2/2\pi$, which is the gradation assumed. With that map the amplitude predicts the one-loop shifts of all three mass ratios exactly.

## 7 The $f_4^{(1)}$ solitons: an open case

For $f_4^{(1)}$ we construct the $U_q(e_6^{(2)})$ R-matrix on $\mathbf{27}=\mathbf{26}\oplus\mathbf1$ and a soliton amplitude that passes every structural check. Its lowest breather, however, does not reproduce any CDS $f_4^{(1)}$ particle amplitude.

### 7.1 Representation and R-matrix

The twisted algebra $U_q(e_6^{(2)})$ has horizontal subalgebra $U_q(f_4)$. Its affine node is $\alpha_0=-\theta_s=-\varepsilon_1$, the highest short root, attached to the short end of the diagram. We extend the $f_4$ generators on the $\mathbf{26}$ by a singlet and solve for $e_0$ and $f_0$ from all relations. The residual is $4\times10^{-16}$, and $e_0$ genuinely couples the singlet, so the $\mathbf{27}$ is indecomposable. Gandenberger, MacKay and Watts constructed this R-matrix earlier [5]; our construction is independent.

The product $\mathbf{27}\otimes\mathbf{27}$ contains the $\mathbf{26}$ three times and the singlet twice. $\check R$ therefore acts on those isotypic blocks as a $3\times3$ and a $2\times2$ matrix, sixteen unknown functions in all. We align the copies of each component by building them from the same words in the lowering operators $\Delta(f_i)$, and impose the $e_0$ and $f_0$ intertwining relations on the highest-weight vectors. The resulting linear system closes to $10^{-15}$, and the full R-matrix commutes with $U_q(f_4)$ and intertwines $e_0$ to $10^{-13}$.

| $x$ | $q^{-2}$ | $q^{-8}$ | $-q^{-6}$ | $-q^{-12}$ |
| --- | --- | --- | --- | --- |
| Residue rank | 378 | 27 | 79 | 1 |
| Meaning | $\mathbf{27}+\mathbf{27}\to\mathbf{378}$ | self-fusion | $\mathbf{27}+\mathbf{27}\to\mathbf{79}$ | singlet (crossing point) |

The classical mass ratios $2\cos(\pi/12)$ and $2\cos(\pi/4)$ are again $m_3/m_1$ and $m_2/m_1$ for the CDS particles at the $f_4^{(1)}$ end, $H=12$.

### 7.2 The amplitude

The mass pattern and the self-fusion at exactly $2\pi/3$ give $T=6\omega+\tfrac32$ and zeros at $\omega+1$, $3\omega+\tfrac32$, $4\omega+1$ and $T$, with $\omega=4\pi/\beta^2-1$. The amplitude is unitary, crossing holds with factor exactly $+1$, and the residue ranks 378, 79, 27 and 1 appear at the predicted positions.

### 7.3 The obstruction

The lowest-breather amplitude matches no CDS amplitude $S_{aa}$. We searched:

- all four particles;
- every $H$ obtained by pairing one of our poles with a CDS block position, each linear in $H$;
- integer shifts of the zeros in $\{-2,\dots,2\}^3$, applied exactly through the CDD shift factors of Section 2.4;
- with and without a sign-fixing CDD factor at $\theta=i\pi/2$;
- four values of $T$: $6\omega+\tfrac12$, $6\omega+\tfrac32$, $6\omega+\tfrac52$ and $6\omega-\tfrac32$.

The best near-miss shares 6 of 14 blocks. The obstruction is specific. Our amplitude has the pole of $B_1+B_1\to B_2$ at $y=1$, fixed by the breather tower. CDS's $S_{11}$ has its coupling-dependent counterpart at $y=-BT/H=\tfrac32$. Moving our pole there requires a value of $T$ incompatible with the self-fusion $\mathbf{27}+\mathbf{27}\to\mathbf{27}$ at exactly $2\pi/3$, which equal masses enforce.

### 7.4 Diagnosis: ruled out at one loop

Five further tests show that the $\mathbf{27}$ amplitude is a consistent solution of the axioms, but not the $f_4^{(1)}$ soliton S-matrix.

- **The method works on the exact analogue.** The $c_3^{(1)}$ soliton 1 is the $\mathbf8=\mathbf7\oplus\mathbf1$ of $U_q(d_4^{(2)})$, structurally identical to the $\mathbf{27}$. Both are twisted Kirillov–Reshetikhin modules whose square contains the quasi-minuscule representation three times and the singlet twice. Its R-matrix poles come in sign pairs, $\pm q^{-2}$ and $\pm q^{-6}$, and the minimal construction with zeros at $\omega$, $\omega+\tfrac12$, $3\omega+\tfrac12$ and $T$ reproduces the known, fused $S_{11}$ exactly, up to overall sign.
- **No pole is missing.** A scan over all 12th roots of unity times $q^{-l}$, for $l\le30$, finds no poles of the $\mathbf{27}$ R-matrix beyond the four listed in Section 7.1.
- **The amplitude is self-consistent.** The bound state $\mathbf{27}+\mathbf{27}\to\mathbf{27}$ at $2\pi/3$ scatters exactly like the $\mathbf{27}$, so the self-fusion bootstrap holds. The breather sector is bootstrap-consistent too: the 111 self-coupling holds, and $S_{12}$ built through $1+1\to2$ equals $S_{12}$ built through $1+3\to2$.
- **It is a different particle S-matrix.** The breather amplitude has CDS's fusion poles for particle 1 (111, 112, 113) at $H=2T/(\omega+1)$, but different coupling-dependent blocks. At weak coupling $\log S_{11}$ equals $\tfrac23$ of CDS's, with exactly the same $\theta$-dependence. So the two imply different relations between $H$ and the tree-level strength: $H-12=\tfrac92 B_{\rm tree}$ for ours, $3B_{\rm tree}$ for CDS.
- **Perturbation theory decides.** From the $f_4^{(1)}$ Lagrangian we computed the one-loop mass shifts and the tree-level $1+1\to1+1$ amplitude. Their ratio $\delta H/B_{\rm tree}$ is independent of normalisations, and the same method gives exactly 1 for $c_3^{(1)}$, DGZ's $H=h+B$. For $f_4^{(1)}$ it gives $3.000000$, confirming CDS and ruling out $\tfrac92$. All three mass ratios give the same $\delta H$, which also confirms CDS's floating-mass formula at one loop.

**The verdict does not depend on CDS.** Compared directly with the $f_4^{(1)}$ Lagrangian, the breather amplitude's first-order term has exactly the tree-level $\theta$-dependence at eight rapidities. That term fixes the coupling map, and with it the amplitude predicts one-loop shifts of $(m_2/m_1)^2$ and $(m_3/m_1)^2$ that are exactly $1.5$ times the Lagrangian's. The identical comparison reproduces the $c_3^{(1)}$ and $e_6^{(2)}$ amplitudes exactly. The tree term also gives $1/\omega=\beta^2/8\pi$, half the value the gradation implies; for $c_3^{(1)}$ and $e_6^{(2)}$ the two agree.

No CDD factor built from up to two crossing-symmetric block pairs converts our breather sector into CDS's.

### 7.5 What this implies

The obstruction is structural. Four facts fix the $\mathbf{27}$ amplitude:

- its R-matrix fixes the fusion positions linearly in $\omega$;
- its single singlet pole gives a breather tower of unit spacing;
- its crossing point $-q^{-12}$ fixes $T$ modulo 1;
- the $\mathbf{27}+\mathbf{27}\to\mathbf{27}$ fusion must sit at exactly $2\pi/3$.

No $T$ compatible with all four yields $H-12=3B_{\rm tree}$.

**Topology fixes the multiplet.** The vacua sit at $\vec\phi=\tfrac{2\pi}{\beta}\vec\lambda$, with $\vec\lambda$ in the co-weight lattice of $f_4$. The charges of the species attached to the affine node therefore form the Weyl orbit of the fundamental co-weight $\omega_1^\vee$, 24 vectors, which are exactly the nonzero weights of the $\mathbf{27}$. Folding $f_4^{(1)}$ from $e_6^{(1)}$, under which the $\mathbf{27}$ of $e_6$ branches as $\mathbf{26}\oplus\mathbf1$, gives the three zero weights as well, and closure under $U_q(e_6^{(2)})$ forces the singlet. The $\mathbf{27}$ is therefore the multiplet of the lightest $f_4^{(1)}$ soliton, and the obstruction must lie in the scalar factor. Topology alone fixes only weights: it cannot count zero-charge states or separate modules with the same weights, which is why closure and folding are needed.

**No finite CDD factor removes the obstruction.** A CDD factor $Z$ multiplies the $\mathbf{27}$ amplitude by a product of factors $\sinh(\theta/2+i\pi c/2T)$ with integer exponents $n_c$, at positions $c=a\omega+b$. Converting the breather sector into CDS's, and preserving the self-fusion bootstrap, require

$$
Z(\theta+i\sigma)\,Z(\theta)^2\,Z(\theta-i\sigma)=\frac{S_{11}^{\rm CDS}}{S_{11}},\qquad \sigma=T-1,\qquad Z(\theta)=Z\big(\theta+\tfrac{i\pi}{3}\big)\,Z\big(\theta-\tfrac{i\pi}{3}\big).
$$

Both conditions are linear in the exponents. We solved them exactly over the position lattice, with $a$ taken modulo 12 and $b$ in quarter steps up to $|b|\le24$, so that arbitrarily many block pairs inside this window are allowed. There is no solution in any window, with or without the bootstrap condition. The least-squares residual falls only as the inverse of the window size, and the best-fit exponents are non-integer and grow without bound. Any solution would therefore be an infinite product, not a CDD factor.

**Conclusion.** The representation, its R-matrix and the gradation are all fixed. If the $\mathbf{27}$ is the right multiplet, the physical scalar factor must differ from $c\,f$ in its Gamma-function structure itself, not merely by a finite CDD factor. No other case in this paper needed this, including the structurally identical $c_3^{(1)}$ vector. Finding that scalar factor is the open problem. Section 7.6 derives what the one-loop data fix and tests what the missing factor could be.

### 7.6 What the one-loop data fix

The one-loop data, meaning the tree amplitude and the bubble-diagram mass shifts computed from the $f_4^{(1)}$ Lagrangian, fix the lifts and the first-order form of the missing correction. They do not fix the scalar factor to all orders.

1. **The tree strength does not depend on the lifts.** In every case tested (the choice of Section 7.2 and candidates (a) and (b) below), the breather amplitude's first-order term has exactly the Lagrangian's shape and implies $1/\omega=\beta^2/8\pi$.
2. **The masses select the lifts.** In these amplitudes the particle mass ratios equal soliton fusion ratios, for example $m_3/m_1=2\cos(\pi a_{378}/2T)$, so they depend only on $T$ and the lifts. Write $T=6\omega+d$, $a_{378}=\omega+c$ and $a_{79}=3\omega+c'$, and use the gradation's coupling map $1/\omega=\beta^2/4\pi$, which is correct for $c_3^{(1)}$ and $e_6^{(2)}$. The one-loop shifts then require
$$
c-\frac d6=\frac14,\qquad \frac{c'}{3}-\frac d6=\frac1{12},
$$
and the self-fusion at $2\pi/3$ requires $d=\tfrac32+3k$. With integer lifts there are exactly two solutions, listed in the table below. Both reproduce the one-loop shifts of $(m_2/m_1)^2$ and $(m_3/m_1)^2$ exactly, whereas the choice $T=6\omega+\tfrac32$ of Section 7.2 gives three times the Lagrangian value. Under the amplitude's own tree map no integer lifts work: that would need quarter-integer offsets.
3. **What remains is a factor of 2 in $B$.** For (a) and (b), every particle fusion pole coincides exactly with CDS's $S_{11}$, at $H=2T/(\omega+1)$ and $H=2T/\omega$ respectively. The only difference is the coupling-dependent blocks: ours are displaced from the particle poles by 1, CDS's by $\tfrac12$. In DGZ/CDS terms our breather amplitude is CDS's $S_{11}$ with $B$ doubled at the same $H$, which is exactly the factor 2 in the tree strength.
4. **The first-order correction is fixed.** A correction $Z$ to the soliton scalar factor must satisfy
$$
Z(\theta+i\pi)\,Z(\theta)^2\,Z(\theta-i\pi)=S_{BB}^{-1/2}\qquad\text{up to }O(\beta^4).
$$
In words, the coupling-dependent shifts must be halved.

| Candidate | $T$ | Zeros of $c$ |
| --- | --- | --- |
| (a) | $6\omega+\tfrac92$ | $\omega+1,\ 3\omega+\tfrac52,\ 4\omega+3$ |
| (b) | $6\omega-\tfrac32$ | $\omega,\ 3\omega-\tfrac12,\ 4\omega-1$ |

Two simple explanations fail. A missing half-integer breather tower is impossible, because the $\mathbf{27}$'s R-matrix has no singlet point where it could sit. Nor is our breather the bound state of two CDS particle-1's at their coupling-dependent pole: the bootstrap amplitude of that bound state differs from ours.

**The half-integer conjecture is false.** A natural guess was that halving the coupling-dependent shifts needs Gamma-function towers of half-integer spacing. By the duplication formula $\Gamma(2z)\propto 2^{2z}\,\Gamma(z)\,\Gamma(z+\tfrac12)$, that is the construction of Section 2.3 with extra zeros at half-shifted positions. Crossing requires them in pairs $\{e,T-e\}$, and each pair multiplies the scalar factor by an infinite-product factor that adds zeros but no poles. Its effect on the breather amplitude has a closed form. Crossing and unitarity turn the breather shifts $\pm(T-1)$ into steps of $\mp1$, and each zero $a\in\{e,T-e\}$ contributes

$$
W_a(t)=\tan\frac{\pi(t-1-a)}{2T}\,\tan\frac{\pi(t+a)}{2T}\,\cot\frac{\pi(t-a)}{2T}\,\cot\frac{\pi(t+1+a)}{2T},
$$

which we verified against the infinite product. We then solved exactly for any real combination of such pairs, on a quarter-step grid of positions (173 pairs for (a), 111 for (b)), against the CDS target. There is no solution: the overdetermined systems have full rank and residuals $0.135$ and $0.178$.

**What a correction would have to be.** The same step-1 structure holds for any crossing-symmetric, unitary correction $Y$: it changes the breather amplitude by

$$
\frac{Y(t)^2}{Y(t+1)\,Y(t-1)} .
$$

In the exponents of the factors $\sinh(\theta/2+i\pi c/2T)$ this is a second difference along lines of unit step. On each line, a finite CDD factor requires the zeroth and first moments of the target exponents to vanish, and a Gamma-type tower requires only the zeroth. For the CDS target the zeroth moment is nonzero on 16 of the 24 lines, for both candidates. Any correction reproducing CDS would therefore need exponents growing linearly along those towers, a Barnes double-Gamma structure in the step direction. No known soliton scalar factor has that structure.

Two possibilities remain. Either the $f_4^{(1)}$ soliton amplitude really carries such a factor, or CDS's $S_{11}$ differs from the true breather sector beyond one loop, which our one-loop data cannot test. A second-order check from the Lagrangian, for instance the one-loop residue of $S_{11}$ at a bound-state pole, would decide between them.

## 8 Check on known cases: $d_4^{(3)}$ and $g_2^{(1)}$

Takács constructed the $d_4^{(3)}$ soliton S-matrix from $U_q(g_2^{(1)})$ [6] and the $g_2^{(1)}$ soliton S-matrix from the 8-dimensional $U_q(d_4^{(3)})$ R-matrix [7]. He found consistency with the real-coupling particle S-matrices through breather–particle correspondence. Applied blind, our method reproduces this, which tests the zero-fixing principle on a pair where the answer is known.

### 8.1 $d_4^{(3)}$ solitons

We solve $U_q(g_2^{(1)})$ on the $\mathbf7$ (the six short roots of $g_2$ plus a zero weight) from the relations, with residual $10^{-16}$. The R-matrix lives in a 4-dimensional space, since $\mathbf7\otimes\mathbf7=\mathbf1\oplus\mathbf7\oplus\mathbf{14}\oplus\mathbf{27}$. Its poles are at $x=q^{-2}$ (rank 15), $x=q^{-8}$ (rank 7, the self-fusion) and $x=q^{-12}$ (rank 1, the singlet and crossing point).

The self-fusion at exactly $2\pi/3$ forces $T=6\omega+3k$. The mass pattern puts the 15-plet zero at $2T/H$, and the breather pole at $y=1$ gives $H=12T/(T-3)$. The unique solution is

$$
a\in\{\omega,\ 4\omega+2\},\qquad T=6\omega+3,\qquad H=12+\frac{6}{\omega},
$$

fixed before any comparison. The lowest-breather amplitude is then identical to the CDS $(g_2^{(1)},d_4^{(3)})$ amplitude $S_{11}=\{1\}_0\,\{H/2\}_{1/2}\,\{H-1\}_0$ at this $H$.

### 8.2 $g_2^{(1)}$ solitons

On the $\mathbf8=\mathbf7\oplus\mathbf1$ of $U_q(d_4^{(3)})$, the special points carry cube-root-of-unity phases: $x=q^{-2}$ (rank 29), $x=e^{\pm2\pi i/3}q^{-4}$ (rank 8 each), and $x=q^{-6}$ (singlet and crossing point). With $T=3\omega+1$ and zeros at $\omega$, $2\omega+\tfrac13$, $2\omega+\tfrac23$ and $T$, the lowest-breather amplitude is identical to the CDS amplitude

$$
S_{22}=\{H/3-1\}_1\,\{H/3+1\}_1\,\{2H/3-1\}_1\,\{2H/3+1\}_1,\qquad H=\frac{6T}{T+1} .
$$

Here the lowest breather is particle 2, not particle 1. For a twisted multiplet the lowest breather can therefore be a heavier particle, which matters for the $f_4^{(1)}$ search.

## 9 Discussion and open problems

A single construction covers every case treated here. The amplitude is $S=c\,f\,\check R$, where $c$ vanishes at the R-matrix poles, and the integer lifts make the soliton–particle mass relations exact. The breathers then reproduce the conjectured particle S-matrices exactly, with $H$ continued past the end of the floating range belonging to the theory.

| Theory | Symmetry and multiplet | $T$ | Floating Coxeter number $H$ | Lowest breather |
| --- | --- | --- | --- | --- |
| $c_n^{(1)}$ | $U_q(d_{n+1}^{(2)})$, spinor | $n\omega+\frac{n-1}{2}$ | $2n-\frac{1}{\omega+1/2}$ | particle $n$ |
| $a_{2n-1}^{(2)}$ | $U_q(b_n^{(1)})$, spinor | $(2n-1)\omega+n-1$ | $2n-1-\frac{1/2}{\omega+1/2}$ | particle $n$ |
| $e_6^{(2)}$ | $U_q(f_4^{(1)})$, $\mathbf{26}$ | $9\omega+3$ | $18+\frac{6}{\omega}$ | particle 1 |
| $d_4^{(3)}$ [6] | $U_q(g_2^{(1)})$, $\mathbf7$ | $6\omega+3$ | $12+\frac{6}{\omega}$ | particle 1 |
| $g_2^{(1)}$ [7] | $U_q(d_4^{(3)})$, $\mathbf8$ | $3\omega+1$ | $\frac{6T}{T+1}$ | particle 2 |
| $f_4^{(1)}$ | $U_q(e_6^{(2)})$, $\mathbf{27}$ | $6\omega+\frac92$ or $6\omega-\frac32$ (fixed by one-loop masses) | minimal $\mathbf{27}$ amplitude fails at one loop; no finite CDD fix | open |

The identifications fix the particle S-matrices' coupling dependence completely. For $c_n^{(1)}$ and $a_{2n-1}^{(2)}$ they give exactly DGZ's $B(\beta)=\dfrac{\beta^2/2\pi}{1+\beta^2/4\pi}$, continued to imaginary coupling. DGZ had checked this form only to $O(\beta^4)$. For $e_6^{(2)}$ they give a definite $H(\beta)$, which CDS left unspecified.

The soliton masses are given exactly by the classical formula with $h$ replaced by $H$. For $c_n^{(1)}$, their one-loop expansion agrees with Delius–Grisaru [11]. The companion analysis traces the disagreement with MacKay–Watts [12] to a Jordan block in the fluctuation operator, which must be counted with its algebraic multiplicity.

Open problems:

1. **$f_4^{(1)}$.** Find the scalar factor of the $\mathbf{27}$ amplitude. Topology fixes the multiplet, the minimal amplitude fails at one loop, and no finite CDD factor repairs it, so the Gamma-function structure itself must change (Section 7.5). The one-loop data fix its lifts and its first-order form (Section 7.6).
2. **Other species.** Each theory has several fundamental soliton species, obtainable by fusion; we carried this out only for $c_3^{(1)}$. The S-matrices of the vector and Kirillov–Reshetikhin solitons of $a_{2n-1}^{(2)}$ and $e_6^{(2)}$ follow mechanically but are not written down here.
3. **Coleman–Thun analysis.** Explain the soliton poles not accounted for here: single-channel tower poles, and three light 29-plet poles in the $c_3^{(1)}$ amplitude $S_{22}$.
4. **Crossing convention.** Derive the factor $(-1)^n$ relating our numerical crossing test to physical crossing for the $c_n^{(1)}$ spinor.
5. **Larger $n$.** Checks for $n\ge5$ need spinor R-matrix projectors of dimension $2^{2n}$, at increasing cost. The closed forms of Sections 3 and 5 are conjectured for all $n$.
6. **Unitarity.** These S-matrices are not unitary in the Hilbert-space sense. The companion foundations note shows they are unitary with respect to an indefinite (Krein) inner product. A unitary theory requires a restriction of RSOS type, which remains to be carried out for these algebras.

## Appendix A: Numerical methods and verification

Every closed form in this paper was identified numerically and then confirmed exactly, by exact rational arithmetic in $\omega$ or by matching at a second coupling.

**Representations.** The generators $e_i$, $f_i$ are unknown complex matrices supported on the allowed weight pairs. We solve the full set of defining relations by Levenberg–Marquardt least squares at unimodular $q=e^{-i\pi\omega}$, with $\omega=2.37$ or $1.61$. These are $[e_i,f_j]=\delta_{ij}[k_i]_{q_i}$ and the $q$-Serre relations. A residual at machine precision ($10^{-16}$) proves that the representation exists.

**R-matrices.** For dimension up to 256 we solve the linear intertwining equations directly. For $\mathbf{26}\otimes\mathbf{26}$ we use highest-weight vectors of $U_q(f_4)$ and the $e_0$ relation, five unknowns per spectral parameter. For $\mathbf{27}\otimes\mathbf{27}$ we use aligned multiplicity copies and the $e_0$ and $f_0$ relations, sixteen unknowns. For the small $(g_2,d_4^{(3)})$ cases we fix a basis of the span of $\check R(x)$ first. Eigenvalues are fitted as rational functions of $x$, and their zeros and poles are identified as $\pm q^l$, $\pm i q^l$, or cube-root phases times $q^l$.

**Scalar factor.** The infinite product $f$ is evaluated as a sum of differences of log-Gamma functions at nearby arguments, which keeps it stable for large truncation $J$. It is normalised by unitarity at a regular point. The truncation error scales exactly as $1/J$: quadrupling $J$ from 4000 to 16000 reduces it by a factor of 3.99–4.00.

**Poles.** Every pole of a fused amplitude is a shifted Gamma-function pole of one factor, so the candidate list is complete. At each candidate we read the pole order from the norm ratio at two distances, and the residue rank from the singular values of the leading coefficient.

**One-loop test.** From the Toda Lagrangian with normal-ordered exponentials we compute the masses and the cubic and quartic couplings. The one-loop mass shifts come from the two-dimensional bubble diagram; tadpoles vanish under normal ordering. The tree-level amplitude $1+1\to1+1$ comes from $s$-, $t$- and $u$-channel exchange plus the contact term. The ratio $\delta H/B_{\rm tree}$ does not depend on how roots or coupling are normalised. It equals $1.000000$ for $c_3^{(1)}$, $3.000000$ for $f_4^{(1)}$, and $-3.000000$ for $e_6^{(2)}$ measured from its own end of the CDS range, in exact agreement with DGZ and CDS. In each case the tree amplitude has exactly the $\theta$-dependence of the exact S-matrix's first-order term.

**Direct comparison.** For the breather amplitudes we also bypass the particle S-matrices entirely. The weak-coupling term $f(\theta)=\lim_{\omega\to\infty}\omega\,\log S/i$ of the amplitude is compared with the Lagrangian tree amplitude at eight rapidities; a constant ratio means identical shape, and that constant fixes the coupling map between $1/\omega$ and $\beta^2$. The amplitude's bound-state poles give the mass ratios as exact functions of $\omega$, so with that map it predicts the one-loop mass shifts, which we compare with the bubble-diagram values. The results are exact agreement for $c_3^{(1)}$ and $e_6^{(2)}$, and a factor $1.5$ for $f_4^{(1)}$.

**Exact checks.**

| Check | Accuracy |
| --- | --- |
| Representation relations (all cases) | $10^{-16}$ |
| Closed-form R-matrices vs direct solution | $10^{-12}$ to $10^{-14}$ |
| Crossing of $S_{33}$ ($c_3^{(1)}$) | $10^{-9}$ |
| Unitarity of fused amplitudes | $10^{-4}$ to $10^{-2}$ (truncation) |
| Breather block forms at a second coupling | $2\times10^{-3}$, extrapolating to $10^{-5}$ |
| Breather–particle identities | exact (rational arithmetic in $\omega$) |

The Python code for every step accompanies these notes.

## References

1. T.J. Hollowood, Solitons in affine Toda field theory, *Nucl. Phys.* B384 (1992) 523.
2. D. Bernard, A. LeClair, Quantum group symmetries and non-local currents in 2D QFT, *Commun. Math. Phys.* 142 (1991) 99.
3. G.W. Delius, Exact S-matrices with affine quantum group symmetry, [hep-th/9503079](https://arxiv.org/abs/hep-th/9503079).
4. G.M. Gandenberger, N.J. MacKay, Exact S-matrices for $d_{n+1}^{(2)}$ affine Toda solitons and their bound states, [hep-th/9506169](https://arxiv.org/abs/hep-th/9506169).
5. G.M. Gandenberger, N.J. MacKay, G.M.T. Watts, Twisted algebra R-matrices and S-matrices for $b_n^{(1)}$ affine Toda solitons and their bound states, [hep-th/9509007](https://arxiv.org/abs/hep-th/9509007).
6. G. Takács, Quantum affine symmetry and scattering amplitudes of the imaginary coupled $d_4^{(3)}$ affine Toda field theory, *Nucl. Phys.* B502 (1997) 629, [hep-th/9701118](https://arxiv.org/abs/hep-th/9701118).
7. G. Takács, The R-matrix of the $U_q(d_4^{(3)})$ algebra and $g_2^{(1)}$ affine Toda field theory, [hep-th/9702196](https://arxiv.org/abs/hep-th/9702196).
8. G.W. Delius, M.D. Gould, Y.-Z. Zhang, Solving the quantum Yang–Baxter equation with the tensor product graph, *Nucl. Phys.* B432 (1994) 377.
9. G.W. Delius, M.T. Grisaru, D. Zanon, Exact S-matrices for nonsimply-laced affine Toda theories, *Nucl. Phys.* B382 (1992) 365, [hep-th/9201067](https://arxiv.org/abs/hep-th/9201067).
10. E. Corrigan, P.E. Dorey, R. Sasaki, On a generalised bootstrap principle, [hep-th/9304065](https://arxiv.org/abs/hep-th/9304065).
11. G.W. Delius, M.T. Grisaru, Toda soliton mass corrections and the particle–soliton duality conjecture, [hep-th/9411176](https://arxiv.org/abs/hep-th/9411176).
12. N.J. MacKay, G.M.T. Watts, Quantum mass corrections for affine Toda solitons, [hep-th/9411169](https://arxiv.org/abs/hep-th/9411169).
13. A.B. Zamolodchikov, Al.B. Zamolodchikov, Factorized S-matrices in two dimensions as the exact solutions of certain relativistic quantum field theory models, *Ann. Phys.* 120 (1979) 253.
14. Companion note: Exceptional points and the one-loop masses of $c_n^{(1)}$ Toda solitons (this project, 2026).
15. Companion note: Foundations of imaginary-coupling affine Toda field theory (this project, 2026).
