# Appendix — The escape criterion (★): statement and self-contained proof

*Draft for a paper appendix, prepared 2026-10-12. It has not yet been reviewed externally or verified formally (e.g. in a proof assistant).*

## A.0 Scope and attribution

The criterion proved here is the case $d = m_1m_2/K$ of **Theorem 1 of Yoshida (1972)** [Y72, p. 400]. It is also the case $d = m_1m_2/K$ of **criterion A′ of Yoshida (1974)** [Y74, p. 369]; see Remark A.4.2.

The corollary to Theorem 1 in [Y72, p. 402] does not assume $r(0)\le d$. **Yoshida himself corrected this, first in an erratum [Y73] and then in [Y74, p. 369, Note].** The erratum [Y73] replaces condition (26) of the corollary by $d\ge r(0)$; it also corrects inequality (40) on p. 405 of [Y72]. The later Note [Y74, p. 369] calls conditions (25)–(26) of [Y72] "incorrect" and replaces them by $\rho(0)\ge ad\ge a\,r(0)$. The statement below assumes $r(t_0)<d_1$ explicitly, so that issue does not affect it.

Lemmas A and B below are the case $n=1$ of Lemma 1 of [Y72, p. 393] and Lemma 2 of [Y72, p. 396]. Yoshida proves them by induction on $n$ and by analysing the $\theta$-derivative, respectively. We give different, self-contained proofs based on the Legendre (multipole) expansion. This is a standard technique, and we do not claim originality for it.

Everything below is proved in full except the classical facts about collisions and their regularization. Section A.2 lists those facts and states exactly how we use them.

---

## A.1 Setting and notation

Three point masses $m_0,m_1,m_2>0$ move in $\mathbb R^2$ or $\mathbb R^3$ under Newtonian gravity with $G=1$. We work in the barycentric frame, where the total momentum is zero. The body $m_0$ is the candidate escaper. The other two are labelled so that
$$m_1\le m_2 .$$
Define
$$
M=m_0+m_1+m_2,\qquad \lambda=\frac{m_2}{m_1+m_2},\qquad \nu=\frac{m_1}{m_1+m_2},\qquad
m=\frac{m_1m_2}{m_1+m_2},\qquad \mu=\frac{m_0(m_1+m_2)}{M}.
$$
Hence $\lambda\ge\nu>0$, $\lambda+\nu=1$, $\mu M\lambda=m_0m_2$ and $\mu M\nu=m_0m_1$.

**Jacobi vectors.**
- $\mathbf r=\mathbf x_2-\mathbf x_1$.
- $\boldsymbol\rho=\mathbf x_0-\mathbf g$, where $\mathbf g=(m_1\mathbf x_1+m_2\mathbf x_2)/(m_1+m_2)$.
- Then $\mathbf x_2-\mathbf g=\nu\mathbf r$ and $\mathbf x_1-\mathbf g=-\lambda\mathbf r$. Thus the lighter body $m_1$ is the one farther from $\mathbf g$.
- Write
$$
r=|\mathbf r|,\quad \rho=|\boldsymbol\rho|,\quad
r_1=|\mathbf x_0-\mathbf x_2|=|\boldsymbol\rho-\nu\mathbf r|,\quad
r_2=|\mathbf x_0-\mathbf x_1|=|\boldsymbol\rho+\lambda\mathbf r| .
$$
- For $r>0$ let $c=\cos\theta=\boldsymbol\rho\cdot\mathbf r/(\rho r)$. Then
$$
r_1^2=\rho^2+\nu^2r^2-2\nu\rho r c,\qquad r_2^2=\rho^2+\lambda^2r^2+2\lambda\rho r c .
$$

**Radial velocity.** $\dot\rho:=\dfrac{d}{dt}|\boldsymbol\rho|=\dfrac{\boldsymbol\rho\cdot\dot{\boldsymbol\rho}}{\rho}$. This is the radial component, *not* the speed $|\dot{\boldsymbol\rho}|\ge|\dot\rho|$.

**Energy.** In the barycentric frame, and with **vector** velocities,
$$
H=\tfrac12 m|\dot{\mathbf r}|^2+\tfrac12\mu|\dot{\boldsymbol\rho}|^2-\frac{m_1m_2}{r}-\mu M\Big(\frac{\lambda}{r_1}+\frac{\nu}{r_2}\Big). \tag{E}
$$
We assume $H<0$ and put $K=-H>0$ and $d_1=m_1m_2/K$.

**Equation of motion of $\boldsymbol\rho$.**
$$
\ddot{\boldsymbol\rho}=-M\Big(\frac{\lambda}{r_1^3}+\frac{\nu}{r_2^3}\Big)\boldsymbol\rho+\lambda\nu M\Big(\frac1{r_1^3}-\frac1{r_2^3}\Big)\mathbf r=:\Phi(\boldsymbol\rho,\mathbf r). \tag{EM}
$$

**Radial acceleration.** Where $\rho>0$, Lagrange's identity gives
$$
\ddot\rho=\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}+\frac{|\boldsymbol\rho|^2|\dot{\boldsymbol\rho}|^2-(\boldsymbol\rho\cdot\dot{\boldsymbol\rho})^2}{\rho^3}\;\ge\;\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}. \tag{RA}
$$
The dropped term equals $|\boldsymbol\rho\times\dot{\boldsymbol\rho}|^2/\rho^3\ge0$. Nothing is assumed about the angular momentum; in particular $L=0$ is allowed.

**Auxiliary functions.** For $\rho>0$, $r\ge0$ and $c\in[-1,1]$ with $r_1,r_2>0$, set
$$
P(\rho,r,c)=\frac{\lambda}{r_1}+\frac{\nu}{r_2},\qquad
F(\rho,r,c)=\frac{\lambda(\rho-\nu rc)}{r_1^3}+\frac{\nu(\rho+\lambda rc)}{r_2^3}=-\partial_\rho P .
$$
- The potential term in (E) is $\mu M P$.
- By (EM), $\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}/\rho=-M\,F(\rho,r,c)$.
- For $r=0$ these reduce to $P=1/\rho$ and $F=1/\rho^2$ whatever the value of $c$.

For $d\ge 0$ and $\rho>\lambda d$, set
$$
W(\rho,d)=\frac{\lambda}{\rho+\nu d}+\frac{\nu}{\rho-\lambda d},\qquad
\hat F(\rho,d)=\frac{\lambda}{(\rho+\nu d)^2}+\frac{\nu}{(\rho-\lambda d)^2}=-\partial_\rho W(\rho,d).
$$
Then $W(\rho,d)=P(\rho,d,-1)$ and $\hat F(\rho,d)=F(\rho,d,-1)$. Geometrically, this is the collinear configuration in which the lighter component $m_1$, at distance $\lambda d$ from $\mathbf g$, lies between $\mathbf g$ and $m_0$.

---

## A.2 Classical facts on collisions (cited, not proved here)

**(S1) Painlevé** [P1897]. In the three-body problem every singularity is a collision. Precisely: if a solution on $[t_0,T)$ with $T<\infty$ cannot be continued beyond $T$, then the positions converge as $t\to T$, and at least two of them coincide in the limit.

**(S2) Regularization of binary collisions** [Su1912; LC1920 (planar); KS1965 (spatial)]. Suppose a solution tends to a collision of $m_1$ and $m_2$ at time $t_c$, while $r_1$ and $r_2$ stay bounded away from $0$ near $t_c$. Then:
- it has a unique continuation beyond $t_c$ as a solution of the regularized equations (physically, an elastic bounce);
- along this continuation, $\boldsymbol\rho$ and $\dot{\boldsymbol\rho}$ are continuous in $t$;
- (EM) holds for $t\neq t_c$, and $H$ is conserved.

**(S3) Non-accumulation** [Su1912; SM1971, Ch. 1]. Binary collision times of a regularized solution cannot accumulate at a finite time, as long as $m_0$ stays at a distance bounded below by some $\delta>0$ from both $m_1$ and $m_2$.

**Convention.** "Solution" means the solution continued through binary collisions of $m_1,m_2$ by (S2). Let $[t_0,T^*)$ be its maximal interval of existence.

**What we assume.** Facts (S1)–(S3) enter the proof only through the following statement:

> **(R)** If $\inf_{[t_0,T^*)}\min(r_1,r_2)>0$, then $T^*=\infty$.

*(R) follows from (S1)–(S3).* Suppose $T^*<\infty$ and $\min(r_1,r_2)\ge\delta>0$ on $[t_0,T^*)$.
- By (S3) there are only finitely many binary collisions in $[t_0,T^*)$. So the solution is classical on $[t_\ell,T^*)$ for some $t_\ell$.
- By (S1) the solution ends in a collision at $T^*$. This collision must be between $m_1$ and $m_2$, because $r_1,r_2\ge\delta$.
- By (S2) the solution continues beyond $T^*$. This contradicts maximality.

**Regularity used below.** Let $I\subset[t_0,T^*)$ be an interval on which $\min(r_1,r_2)$ is bounded below by a positive constant. Then $t\mapsto\boldsymbol\rho(t)$ is $C^2$ on $I$, with $\ddot{\boldsymbol\rho}=\Phi(\boldsymbol\rho,\mathbf r)$. The reasons:
- $\Phi$ is continuous on $\{r_1>0,\,r_2>0\}$, including the points with $\mathbf r=0$.
- By (S2), $t\mapsto(\boldsymbol\rho(t),\mathbf r(t))$ is continuous. Hence $t\mapsto\Phi(\boldsymbol\rho(t),\mathbf r(t))$ is continuous.
- $\dot{\boldsymbol\rho}$ is continuous, and on each collision-free interval it is a primitive of this function. Collision times are isolated by (S3), so it is a primitive everywhere.

Consequently, wherever $\rho>0$, $\rho$ is $C^2$ and (RA) holds.

**Caveat on the citations.** We checked the bibliographic data of [P1897, Su1912, LC1920, KS1965, SM1971]. We did **not** check the exact wording of the cited statements in these texts.
- Sundman's global results assume a nonzero total angular momentum. As we read it, that assumption serves to exclude triple collision and to keep the third body away from a colliding pair.
- In our application the total angular momentum is $L=0$. Separation is guaranteed instead by the a priori bound (SEP) proved in Step 2 of A.4.
- We have not found a verbatim statement of (S3) for $L=0$ under such an a priori bound. Remark A.4.5 sketches why it holds; that sketch is not part of the proof.

---

## A.3 Two inequalities

> **Lemma A.** Let $\lambda\ge\nu>0$, $\lambda+\nu=1$, $d\ge0$ and $\rho>\lambda d$. Then for all $0\le r\le d$ and all $c\in[-1,1]$,
> $$P(\rho,r,c)\le W(\rho,d).$$
>
> **Lemma B.** Under the same hypotheses,
> $$F(\rho,r,c)\le \hat F(\rho,d).$$

Both inequalities become equalities at $(r,c)=(d,-1)$, so they are sharp.

Combining Lemma B with (EM) gives the following. Along a solution, if $\rho>\lambda d$ and $r\le d$, then
$$
\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}=-M F\ \ge\ -M\hat F(\rho,d)=M\,\partial_\rho W(\rho,d). \tag{B′}
$$

**Facts on Legendre polynomials $P_n$.**
- (L1) Generating function [DLMF, Eq. 18.12.11]: $(1-2xz+z^2)^{-1/2}=\sum_{n\ge0}P_n(x)z^n$ for $x\in[-1,1]$ and $|z|<1$.
- (L2) Laplace's integral [DLMF, Eq. 18.10.5]: $P_n(\cos\vartheta)=\frac1\pi\int_0^\pi(\cos\vartheta+i\sin\vartheta\cos\phi)^n\,d\phi$. Since $|\cos\vartheta+i\sin\vartheta\cos\phi|^2=\cos^2\vartheta+\sin^2\vartheta\cos^2\phi\le1$, we get $|P_n(x)|\le1$ on $[-1,1]$.
- (L3) In (L1), replacing $(x,z)$ by $(-x,-z)$ gives $P_n(-x)=(-1)^nP_n(x)$. Setting $x=1$ gives $\sum P_n(1)z^n=(1-z)^{-1}$, so $P_n(1)=1$. Hence $P_n(-1)=(-1)^n$.

**Proof of Lemmas A and B.**

*Step 1 (normalization).*
- Put $s=r/\rho$, $a=\nu s$ and $b=\lambda s$. Since $\lambda r\le\lambda d<\rho$, we have $0\le a\le b<1$.
- Now $r_1=\rho\sqrt{1-2ac+a^2}$ and $\rho-\nu rc=\rho(1-ac)$.
- Similarly $r_2=\rho\sqrt{1-2b(-c)+b^2}$ and $\rho+\lambda rc=\rho(1-b(-c))$.
- Therefore
$$
\rho P=\lambda G(a,c)+\nu G(b,-c),\qquad \rho^2F=\lambda g(a,c)+\nu g(b,-c),
$$
where
$$
G(a,x)=(1-2ax+a^2)^{-1/2},\qquad g(a,x)=(1-ax)(1-2ax+a^2)^{-3/2}.
$$

*Step 2 (expansions).*
- By (L1), $G(a,x)=\sum_{n}a^nP_n(x)$.
- One checks directly that $g=G+a\,\partial_aG$. Hence $g(a,x)=\sum_n(n+1)a^nP_n(x)$.
- For $0\le a<1$ and $|x|\le1$ both series converge absolutely, since $\sum(n+1)a^n|P_n(x)|\le(1-a)^{-2}$ by (L2). So they may be added and regrouped term by term.

*Step 3 (combination).* By Steps 1 and 2 together with (L3),
$$
\rho^2F=\sum_{n\ge0}(n+1)s^n\beta_nP_n(c),\qquad \rho P=\sum_{n\ge0}s^n\beta_nP_n(c),\qquad
\beta_n:=\lambda\nu^n+\nu(-\lambda)^n .
$$
- $\beta_0=1$, and $\beta_1=\lambda\nu-\nu\lambda=0$ (no dipole term about $\mathbf g$).
- **For every $n$, $(-1)^n\beta_n\ge0$.**
  - For even $n$, $(-1)^n\beta_n=\lambda\nu^n+\nu\lambda^n>0$.
  - For odd $n$, $(-1)^n\beta_n=\nu\lambda^n-\lambda\nu^n=\lambda\nu(\lambda^{n-1}-\nu^{n-1})\ge0$, because $\lambda\ge\nu$.
  - This is the only place where $m_1\le m_2$ is used.
- Consequently $|\beta_n|=(-1)^n\beta_n$.

*Step 4 (maximum over the angle).* By (L2) and Step 3,
$$
\rho^2F\le\sum_n(n+1)s^n|\beta_n|=\lambda\sum_n(n+1)(-\nu s)^n+\nu\sum_n(n+1)(\lambda s)^n
=\frac{\lambda}{(1+\nu s)^2}+\frac{\nu}{(1-\lambda s)^2}.
$$
In the same way,
$$
\rho P\le\sum_n s^n|\beta_n|=\frac{\lambda}{1+\nu s}+\frac{\nu}{1-\lambda s}.
$$
The right-hand sides equal $\rho^2F(\rho,r,-1)$ and $\rho P(\rho,r,-1)$, respectively.

*Step 5 (maximum over the separation).*
- Each term $(n+1)s^n|\beta_n|$ and $s^n|\beta_n|$ is nondecreasing in $s\ge0$.
- Since $r\le d$, we have $s\le s_d:=d/\rho<1/\lambda$, so all series still converge at $s_d$.
- Therefore
$$
\rho^2F\le\frac{\lambda}{(1+\nu s_d)^2}+\frac{\nu}{(1-\lambda s_d)^2}=\rho^2\hat F(\rho,d),\qquad
\rho P\le\frac{\lambda}{1+\nu s_d}+\frac{\nu}{1-\lambda s_d}=\rho\,W(\rho,d).
$$
Dividing by $\rho^2$ and by $\rho$ gives the lemmas. $\square$

**Remarks on the range.**
- (i) No upper bound on $\rho$ is needed.
- (ii) The bounds blow up as $\rho\to\lambda d^+$, but they hold at every point with $\rho>\lambda d$. No uniformity is used.
- (iii) $r=0$ is included; this covers the instant of a binary collision.
- (iv) Only $c=\cos\theta$ enters, so the lemmas hold in the plane and in space.
- (v) If $\lambda=\nu$, all odd $\beta_n$ vanish and equality also holds at $c=+1$.

**Remark (second proof of Lemma A).** Fix $r$. Then:
- $\partial_cP=\lambda\nu\rho r\,(r_1^{-3}-r_2^{-3})$, and $r_2^2-r_1^2=r\big((\lambda-\nu)r+2\rho c\big)$.
- So $P$ decreases in $c$ up to $c^*=-(\lambda-\nu)r/(2\rho)$ and increases after it. Its maximum over $[-1,1]$ is therefore attained at $c=\pm1$.
- $P(\rho,r,-1)-P(\rho,r,1)=2\lambda\nu r\big[(\rho^2-\lambda^2r^2)^{-1}-(\rho^2-\nu^2r^2)^{-1}\big]\ge0$.
- Finally, $P(\rho,r,-1)=W(\rho,r)$ and $\partial_dW=\lambda\nu\big[(\rho-\lambda d)^{-2}-(\rho+\nu d)^{-2}\big]\ge0$.

---

## A.4 The criterion

> **Theorem (★).** Let $m_1\le m_2$ and $H<0$, in the barycentric frame. Suppose that at some time $t_0$
>
> - **(C1)** $r(t_0)<d_1=m_1m_2/K$,
> - **(C2)** $\rho(t_0)>\lambda d_1$,
> - **(C3)** $\dot\rho(t_0)>\sqrt{2M\,W(\rho(t_0),d_1)}$.
>
> Let $C:=\sqrt{\dot\rho(t_0)^2-2M\,W(\rho(t_0),d_1)}>0$. Then the solution, continued through binary collisions as in A.2, exists for all $t\ge t_0$. Moreover, for all $t\ge t_0$:
> 1. $r(t)\le \dfrac{m_1m_2}{K+\frac12\mu C^2}<d_1$;
> 2. $\tfrac12\dot\rho(t)^2-M\,W(\rho(t),d_1)\ge\tfrac12C^2$ and $\dot\rho(t)>C$;
> 3. $\rho(t)\ge\rho(t_0)+C\,(t-t_0)$.
>
> In particular $\rho\to\infty$ while the binary separation stays bounded; that is, the motion is hyperbolic–elliptic.

**Proof.**

*Step 1 (bootstrap set).* Let
$$
\mathcal G=\big\{t\in[t_0,T^*):\ r(s)\le d_1\ \text{and}\ \dot\rho(s)\ge0\ \text{for all } s\in[t_0,t]\big\},\qquad \tau=\sup\mathcal G .
$$
$\mathcal G$ is an interval containing $t_0$. By (C1) and (C3), $r(t_0)<d_1$ and $\dot\rho(t_0)>0$, and both $r$ and $\dot\rho$ are continuous. Hence $\tau>t_0$.

*Step 2 (estimates on $\mathcal G$).* Fix $t\in\mathcal G$.

(a) Separation.
- $\rho$ is nondecreasing on $[t_0,t]$. So for $s\in[t_0,t]$ we have $\rho(s)\ge\rho(t_0)>\lambda d_1\ge\lambda r(s)$.
- Put $\delta:=\rho(t_0)-\lambda d_1>0$. Then on $[t_0,t]$
$$
r_2\ge\rho-\lambda r\ge\delta,\qquad r_1\ge\rho-\nu r\ge\rho(t_0)-\nu d_1\ge\delta. \tag{SEP}
$$

(b) Radial velocity.
- For $s\in[t_0,t]$, (RA) and (B′) with $d=d_1$ give $\ddot\rho(s)\ge M\,\partial_\rho W(\rho(s),d_1)$.
- Multiplying by $\dot\rho(s)\ge0$ gives $\dot\rho\ddot\rho\ge\frac{d}{ds}\big[M\,W(\rho(s),d_1)\big]$.
- Integrating (this uses $\rho\in C^2$, see A.2):
$$
\tfrac12\dot\rho(t)^2-M\,W(\rho(t),d_1)\ \ge\ \tfrac12\dot\rho(t_0)^2-M\,W(\rho(t_0),d_1)=\tfrac12 C^2 .
$$
- Since $W>0$ and $\dot\rho(t)\ge0$, it follows that $\dot\rho(t)>C$.

(c) Binary separation.
- If $r(t)=0$ there is nothing to prove.
- Otherwise use (E), $\tfrac12 m|\dot{\mathbf r}|^2\ge0$, $|\dot{\boldsymbol\rho}|\ge|\dot\rho|$, Lemma A (valid since $r(t)\le d_1$ and $\rho(t)>\lambda d_1$), and (b):
$$
\frac{m_1m_2}{r(t)}=K+\tfrac12 m|\dot{\mathbf r}|^2+\tfrac12\mu|\dot{\boldsymbol\rho}|^2-\mu M P
\ \ge\ K+\mu\Big(\tfrac12\dot\rho^2-M\,W(\rho,d_1)\Big)\ \ge\ K+\tfrac12\mu C^2 .
$$
- Hence $r(t)\le m_1m_2/(K+\tfrac12\mu C^2)<d_1$.

*Step 3 ($\tau=T^*$).* Suppose $\tau<T^*$.
- By continuity, $r(\tau)\le d_1$ and $\dot\rho(\tau)\ge0$. Hence $\tau\in\mathcal G$.
- By Step 2, $r(\tau)<d_1$ and $\dot\rho(\tau)>C>0$.
- By continuity these strict inequalities persist on $[\tau,\tau+\varepsilon)$ for some $\varepsilon>0$. So $\tau+\varepsilon/2\in\mathcal G$, which contradicts $\tau=\sup\mathcal G$.
- Therefore $\mathcal G=[t_0,T^*)$.

*Step 4 (global existence).* By (SEP), $\min(r_1,r_2)\ge\delta>0$ on $[t_0,T^*)$. So $T^*=\infty$ by (R).

*Step 5 (conclusion).* Items 1 and 2 hold for every $t\ge t_0$ by Steps 2–4. Item 3 follows by integrating $\dot\rho>C$. $\square$

### Remarks

**A.4.1 (angular momentum).** The only place where information is discarded is (RA), and the discarded term is nonnegative. The theorem therefore applies to $L=0$, which is the case of the Pythagorean problem.

**A.4.2 (relation to Yoshida's statements).** These relations are based on the statements of [Y72], [Y73] and [Y74] as transcribed from the original scans.

*Theorem 1 of [Y72].* Its hypotheses are:
- (13) $\rho(0)\ge ad>a\,r(0)$ for some $a>\lambda$;
- (14) $\dot\rho(0)>\sqrt{2MW(\rho(0),d)}$;
- (15) $\dot\rho(0)^2-2MW(\rho(0),d)>2H/\mu+2m_1m_2/(\mu d)$.

For $d=d_1$:
- The right side of (15) is $0$, so (15) follows from (14).
- (13) is equivalent to (C1) together with (C2): take $a\in(\lambda,\rho(t_0)/d_1]$.

So (★) is Theorem 1 of [Y72] with $d=m_1m_2/K$.

*Criterion A′ of [Y74].* With $d=d_1$ it requires $\rho(0)\ge ad_1\ge a\,r(0)$, which allows the boundary case $r(0)=d_1$. That case cannot occur together with (C2) and (C3): Step 2(c) at $t=t_0$ would force $r(t_0)<d_1$.

*Corollary of [Y72].* It has the gap described in A.0, which was corrected in [Y73] and again in [Y74, p. 369]. We do not use it.

**A.4.3 (ordering of the masses).** If $m_1>m_2$, relabel the binary components. With $\lambda$ and $\nu$ interchanged in $W$, the function becomes smaller, and Lemmas A–B fail, for example at $(r,c)=(d,-1)$. The ordering is essential.

**A.4.4 (checkability).** (C1)–(C3) are strict inequalities in $(r,\rho,\dot\rho)$ with exact constants. It suffices to check them for every state in a rigorous enclosure of the true state at a single time $t_0$. The value of $K$ must be the exact energy of the true solution, not the energy recomputed from an enclosure.

**A.4.5 (why (S3) holds here — a sketch, not part of the proof).** Assume (SEP) holds on a finite interval.

The perturbed binary.
- The relative motion satisfies $\ddot{\mathbf r}=-k\,\mathbf r/r^3+\mathbf f$, with $k=m_1+m_2$ and $\mathbf f=m_0\big[(\boldsymbol\rho-\nu\mathbf r)/r_1^3-(\boldsymbol\rho+\lambda\mathbf r)/r_2^3\big]$.
- By (SEP), $|\mathbf f|\le 2m_0/\delta^2$.
- Also $|\ddot{\boldsymbol\rho}|\le M/\delta^2$. So $|\dot{\boldsymbol\rho}|$ is bounded on finite intervals.

Bounds on the Kepler energy.
- Let $h=\tfrac12|\dot{\mathbf r}|^2-k/r$ be the Kepler energy per unit reduced mass, so $mh=H-\tfrac12\mu|\dot{\boldsymbol\rho}|^2+\mu MP$.
- Upper bound: $mh\le-K-\tfrac12\mu C^2$, by Step 2(b).
- Lower bound: $mh\ge-K-\tfrac12\mu|\dot{\boldsymbol\rho}|^2$, which is bounded below on finite intervals.
- So $h\in[-A,-a]$ with $0<a\le A<\infty$.

Levi-Civita variables (planar case).
- Write $\mathbf r\simeq z=u^2\in\mathbb C$ and $dt=|u|^2ds$.
- A direct computation (checked by computer algebra) gives
$$
u''=\tfrac h2u+\tfrac12|u|^2\bar u f,\qquad |u'|^2=\tfrac12\big(k+h|u|^2\big).
$$
- So $u$ is a perturbed harmonic oscillator. Its frequency $\sqrt{-h/2}$ is bounded above and below, and $u$ passes through $0$ with $|u'|=\sqrt{k/2}$.
- Each revolution therefore takes a fictitious time bounded above and contributes a physical time $\int|u|^2ds$ bounded below. Collisions cannot accumulate at finite $t$.

Making this rigorous requires the standard perturbation estimates, which we do not reproduce.

---

## A.5 Application to the Pythagorean problem

**Initial data** [Bu1913; SP1967]. Masses $3,4,5$ start at rest at $(1,3)$, $(-2,-1)$ and $(1,-1)$, respectively.
- The barycentre is at the origin, and the total momentum and angular momentum vanish.
- The mutual distances are $|3\text{–}4|=5$, $|3\text{–}5|=4$ and $|4\text{–}5|=3$.
- Hence, exactly,
$$
H=-\Big(\frac{3\cdot4}{5}+\frac{3\cdot5}{4}+\frac{4\cdot5}{3}\Big)=-\frac{769}{60},\qquad K=\frac{769}{60}.
$$
- The energy is conserved along the (regularized) solution, so this $K$ is the value to use at any later time.

**Parameters.** The escaper is $m_0=3$, and the binary is $m_1=4\le m_2=5$:
$$
M=12,\quad \lambda=\tfrac59,\quad \nu=\tfrac49,\quad \mu=\tfrac94,\quad m=\tfrac{20}9,\quad
d_1=\frac{m_1m_2}{K}=\frac{1200}{769}\approx1.560468140442,
$$
$$
\lambda d_1=\frac{2000}{2307}\approx0.866926744690,\qquad \nu d_1=\frac{1600}{2307}\approx0.693541395752 .
$$

**Conditions.** Let $\mathbf x_3,\mathbf x_4,\mathbf x_5$ be the positions of the masses $3,4,5$. Then
$$
r=|\mathbf x_5-\mathbf x_4|,\qquad \boldsymbol\rho=\mathbf x_3-\frac{4\mathbf x_4+5\mathbf x_5}{9},\qquad \dot\rho=\frac{\boldsymbol\rho\cdot\dot{\boldsymbol\rho}}{\rho},
$$
and the conditions of Theorem (★) read
$$
\text{(C1)}\ r<\frac{1200}{769},\qquad \text{(C2)}\ \rho>\frac{2000}{2307},\qquad
\text{(C3)}\ \dot\rho>\sqrt{24\,W(\rho)},\quad W(\rho)=\frac{5/9}{\rho+1600/2307}+\frac{4/9}{\rho-2000/2307}.
$$
In $W$, the coefficient $5/9$ goes with the **larger** distance $\rho+\nu d_1$, and $4/9$ with the **smaller** distance $\rho-\lambda d_1$.

**Conclusion.** If (C1)–(C3) hold for the true state at some $t_0$, then for all $t\ge t_0$:
- the $4$–$5$ separation stays below $m_1m_2/(K+\tfrac12\mu C^2)<1200/769$;
- the body of mass $3$ recedes at least linearly: $\rho(t)\ge\rho(t_0)+C(t-t_0)$.

Here $\mu=9/4$ and $C$ is defined from the state at $t_0$.

Verifying (C1)–(C3) for the actual trajectory is a separate, computer-assisted step described in the main text. It is not part of this appendix.

**Historical note.** Yoshida [Y72, p. 406] already applied the corollary to this problem at $t=59.72$, using an earlier non-rigorous numerical solution. At that point $r\approx0.42<d_1$, so the corrected condition was also satisfied.

---

## A.6 Independent cross-checks (not part of the proof)

These checks look for errors in the proof above; they do not replace it.
- **Computer algebra (SymPy).** We checked every algebraic identity used in A.1 and A.3: (E) and (EM) derived from the Cartesian Newton equations, (RA) in 2D and 3D, the expansions of Steps 1–3 up to finite order, the closed forms of Step 4, the second proof of Lemma A, and the Levi-Civita relations of A.4.5. Deliberately wrong variants were rejected.
- **Adversarial numerical search** (50-digit arithmetic, 40 000 points concentrated near the boundaries $\rho\to\lambda d^+$, $r\in\{0,d\}$, $c=\pm1$, for several $\lambda\in[\tfrac12,1)$). We found no violation of Lemmas A–B. The maximum ratio was $1+O(10^{-50})$, attained at $(r,c)=(d,-1)$.
- **Dynamical falsification test.** We took 800 states with masses $3,4,5$ and various energies satisfying (C1)–(C3), and integrated them to $t\le40$. No conclusion of Theorem (★) was violated. Variants with the factor $2$ in (C3) omitted, with (C1) dropped, or with $\lambda,\nu$ interchanged in $W$ did produce violations, including the escaper returning.

---

## References

- **[Bu1913]** C. Burrau, *Numerische Berechnung eines Spezialfalles des Dreikörperproblems*, Astron. Nachr. **195** (1913) 113–118. doi:10.1002/asna.19131950602.
- **[DLMF]** NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov/, Eqs. 18.10.5 and 18.12.11 (accessed 2026-10-12).
- **[KS1965]** P. Kustaanheimo, E. Stiefel, *Perturbation theory of Kepler motion based on spinor regularization*, J. Reine Angew. Math. **218** (1965) 204–219.
- **[LC1920]** T. Levi-Civita, *Sur la régularisation du problème des trois corps*, Acta Math. **42** (1920) 99–144. doi:10.1007/BF02404404.
- **[P1897]** P. Painlevé, *Leçons sur la théorie analytique des équations différentielles, professées à Stockholm (1895)*, Paris, 1897.
- **[SM1971]** C. L. Siegel, J. K. Moser, *Lectures on Celestial Mechanics*, Grundlehren der mathematischen Wissenschaften **187**, Springer, Berlin, 1971 (reprinted in Classics in Mathematics, 1995).
- **[Su1912]** K. F. Sundman, *Mémoire sur le problème des trois corps*, Acta Math. **36** (1912) 105–179. doi:10.1007/BF02422379.
- **[SP1967]** V. Szebehely, C. F. Peters, *Complete solution of a general problem of three bodies*, Astron. J. **72** (1967) 876. doi:10.1086/110355.
- **[Y72]** J. Yoshida, *Improved criteria for hyperbolic-elliptic motion in the general three-body problem*, Publ. Astron. Soc. Japan **24** (1972) 391–408. doi:10.1093/pasj/24.3.391.
- **[Y73]** J. Yoshida, *Erratum: Improved criteria for hyperbolic-elliptic motion in the general three-body problem*, Publ. Astron. Soc. Japan **25** (1973) 285. doi:10.1093/pasj/25.2.285.
- **[Y74]** J. Yoshida, *Improved criteria for hyperbolic-elliptic motion in the general three-body problem. II*, Publ. Astron. Soc. Japan **26** (1974) 367–377.
