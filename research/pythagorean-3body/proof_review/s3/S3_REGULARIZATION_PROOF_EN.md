# Global existence through binary collisions under a separation bound (planar case)

*Draft, 2026-10-13. Intended to replace the cited facts (S1)–(S3) of Appendix A.2 for the planar case.
Not externally reviewed. Algebraic identities checked by computer algebra (`s3_algebra_check.py`).*

## 0. What this replaces

Appendix A.2 of the paper cites three classical facts. Statement (R) is the only place where the proof of (★)
uses them:

> **(R)** If $\inf_{[t_0,T^*)}\min(r_1,r_2)>0$, then $T^*=\infty$.

The paper also uses a regularity property: $\boldsymbol\rho\in C^2$ in $t$ with $\ddot{\boldsymbol\rho}=\Phi$.
The weak point was (S3), the non-accumulation of collisions. No verbatim statement of it was found for $L=0$
under a separation bound.

This note proves (R) and the regularity property directly for the **planar** problem. That case covers the
Pythagorean application. The proof uses:

- the Levi-Civita change of variables, used as a *definition* of the continued solution;
- the standard existence, uniqueness and extension theorem for ODEs with locally Lipschitz right-hand side
  [Hartman, *Ordinary Differential Equations*, 2nd ed., SIAM 2002, Ch. II; exact theorem numbers not checked].

It cites no fact about collisions. In particular it does not use Painlevé's theorem (S1), Sundman's
regularization (S2) or the non-accumulation fact (S3). The spatial case would need the
Kustaanheimo–Stiefel analogue. It is not treated here, and the application does not need it.

Notation is that of A.1: $k=m_1+m_2$, $\mathbf r=\mathbf x_2-\mathbf x_1$, and
$r_1=|\boldsymbol\rho-\nu\mathbf r|$, $r_2=|\boldsymbol\rho+\lambda\mathbf r|$.

Plane vectors are identified with complex numbers. For complex $a,b$ we write
$\langle a,b\rangle=\operatorname{Re}(\bar a b)$, which is the dot product.

## 1. The regularized system

The relative motion satisfies
$$\ddot{\mathbf r}=-k\,\mathbf r/r^3+\mathbf f,\qquad
\mathbf f=m_0\Big[\frac{\boldsymbol\rho-\nu\mathbf r}{r_1^3}-\frac{\boldsymbol\rho+\lambda\mathbf r}{r_2^3}\Big],$$
and $\ddot{\boldsymbol\rho}=\Phi(\boldsymbol\rho,\mathbf r)$ by (EM).

The state is $Y=(t,u,w,h,\boldsymbol\rho,\boldsymbol\sigma)\in\mathbb R\times\mathbb C\times\mathbb C\times\mathbb R\times\mathbb C\times\mathbb C$. Write $z=u^2$; it plays the role of $\mathbf r$. On
$$D=\{Y:\ |\boldsymbol\rho-\nu z|>0,\ |\boldsymbol\rho+\lambda z|>0\}$$
consider the system, with $'=d/ds$:
$$
\begin{aligned}
t'&=|u|^2, & u'&=w, & w'&=\tfrac h2\,u+\tfrac12|u|^2\,\bar u\,\mathbf f(\boldsymbol\rho,z),\\
h'&=\langle 2uw,\mathbf f(\boldsymbol\rho,z)\rangle, & \boldsymbol\rho'&=|u|^2\boldsymbol\sigma, & \boldsymbol\sigma'&=|u|^2\,\Phi(\boldsymbol\rho,z).
\end{aligned}\tag{LC}
$$

$D$ is open. On $D$ the right-hand side is real-analytic, because its only denominators are $r_1^3$ and
$r_2^3$. The point $u=0$, a binary collision, is an ordinary point of (LC).

**Lemma 1 (identities).** Along every solution of (LC) in $D$:
1. $I:=2|w|^2-h|u|^2$ is constant.
2. $\mathcal H:=m\,h+\tfrac12\mu|\boldsymbol\sigma|^2-\mu M\,P(\boldsymbol\rho,z)$ is constant. Here $P=\lambda/r_1+\nu/r_2$ and $m=m_1m_2/k$.

Suppose in addition that $I=k$ and $u(s)\neq0$ on an $s$-interval $J$. Put $\dot{}=d/dt=|u|^{-2}\,d/ds$. Then on $J$:

3. $\dot z=2w/\bar u$ and $\ddot z=-k\,z/|z|^3+\mathbf f$;
4. $h=\tfrac12|\dot z|^2-k/|z|$, the Kepler energy;
5. $\dot{\boldsymbol\rho}=\boldsymbol\sigma$ and $\ddot{\boldsymbol\rho}=\Phi$;
6. $\mathcal H$ equals the Newtonian energy (E) of the corresponding three-body state in the barycentric frame.

The relations $\ddot{\mathbf r}$ and $\ddot{\boldsymbol\rho}$ above agree with Newton's equations under the Jacobi
change of variables.

*Proof.* All items are polynomial or rational identities. `s3_algebra_check.py` checks them in sympy. It
treats $r_1,r_2$ as symbols, applies the chain rule, and reduces modulo
$r_1^2=|\boldsymbol\rho-\nu z|^2$ and $r_2^2=|\boldsymbol\rho+\lambda z|^2$.

By hand, item 3 is short:
- $\dot z=2uw/|u|^2=2w/\bar u$.
- $\frac{d}{ds}\dot z=2w'/\bar u-2|w|^2/\bar u^2$.
- Insert $w'$ and $2|w|^2=k+h|u|^2$:
$$\ddot z=|u|^{-2}\big(hu/\bar u+|u|^2\mathbf f-(k+h|u|^2)/\bar u^2\big)=\mathbf f-\frac{k}{|u|^2\bar u^2}=\mathbf f-\frac{k\,z}{|z|^3}.$$

Item 1 follows from $\frac d{ds}2|w|^2=2h\langle w,u\rangle+2|u|^2\langle w,\bar u\mathbf f\rangle$ and $\frac d{ds}(h|u|^2)=|u|^2\langle 2uw,\mathbf f\rangle+2h\langle u,w\rangle$, together with $\langle w,\bar u\mathbf f\rangle=\operatorname{Re}(\bar w\bar u\mathbf f)=\langle uw,\mathbf f\rangle$. $\square$

**Definition (continued solution).** Let a classical solution be given near $t_0$ with $r(t_0)>0$. Choose $u_0$ with $u_0^2=\mathbf r(t_0)$, and put $w_0=\bar u_0\dot{\mathbf r}(t_0)/2$ and $h_0=\tfrac12|\dot{\mathbf r}(t_0)|^2-k/r(t_0)$. Then $I=k$.

The *continued solution* is the maximal solution $Y(s)$, $s\in[0,S^*)$, of (LC) in $D$ with $Y(0)=(t_0,u_0,w_0,h_0,\boldsymbol\rho(t_0),\dot{\boldsymbol\rho}(t_0))$. We read it in physical time through $t(s)$, and put $T^*:=\lim_{s\to S^*}t(s)\in(t_0,\infty]$.

The other choice $-u_0$ gives the solution $(t,-u,-w,h,\boldsymbol\rho,\boldsymbol\sigma)$. That is the same physical motion.

This is Levi-Civita's regularization. At a collision $u=0$ the solution passes smoothly, and physically this is the elastic bounce.

**Lemma 2 (properties of the continued solution).**
- (i) $t(s)$ is $C^1$ and strictly increasing. The zeros of $u$, the binary collisions, are isolated in $s$. So on $[0,S^*)$ they form a discrete set, and in $t$ they have no accumulation point in $[t_0,T^*)$.
- (ii) Let $s(t)$ be the inverse of $t(s)$. Then $\boldsymbol\rho(t):=\boldsymbol\rho(s(t))$ is $C^2$ on $[t_0,T^*)$, with $\dot{\boldsymbol\rho}=\boldsymbol\sigma(s(t))$ and $\ddot{\boldsymbol\rho}=\Phi(\boldsymbol\rho,\mathbf r)$, including at collision instants. Also $\mathbf r(t)=u(s(t))^2$ is continuous.
- (iii) On every interval free of collisions, the functions $\mathbf r,\boldsymbol\rho$ solve Newton's equations, in Jacobi form. The energy (E) holds there with the constant value $\mathcal H=H$.

*Proof.*
(i) $t'=|u|^2\ge0$. By Lemma 1(1), at a zero of $u$ we have $|u'|^2=|w|^2=k/2>0$. So each zero of $u$ is a simple zero of the $C^1$ map $u:\mathbb R\to\mathbb C$, and is isolated. Hence $t'>0$ except on a discrete set, and $t$ is strictly increasing.

(ii) $t(s)$ is $C^1$ and increasing, and $\boldsymbol\sigma$ is continuous. The substitution rule gives
$$\boldsymbol\rho(s(t))-\boldsymbol\rho(t_0)=\int_0^{s(t)}\boldsymbol\sigma\,t'\,ds=\int_{t_0}^{t}\boldsymbol\sigma(s(\tau))\,d\tau .$$
The integrand is continuous in $\tau$, so $\boldsymbol\rho\in C^1$ with $\dot{\boldsymbol\rho}=\boldsymbol\sigma\circ s$. The same argument with $\boldsymbol\sigma'=|u|^2\Phi$ gives $\dot{\boldsymbol\sigma}=\Phi\circ s$, which is continuous.

(iii) This is Lemma 1(3)–(6). $\mathcal H$ is constant by Lemma 1(2). It equals $H$ at $s=0$. $\square$

## 2. Proposition (R)

> **Proposition.** Let $Y$ be the continued solution. Suppose $\min(r_1,r_2)\ge\delta>0$ along it, that is, on $[t_0,T^*)$. Then $T^*=\infty$.

*Proof.* Suppose $T^*<\infty$. We first derive bounds that hold for all $s\in[0,S^*)$. Write $\Delta=T^*-t_0$.

**(B1) Forces.** For any $(\boldsymbol\rho,z)$ with $r_1,r_2\ge\delta$:
- $|\mathbf f|\le m_0(r_1^{-2}+r_2^{-2})\le 2m_0/\delta^2=:F$;
- $P\le(\lambda+\nu)/\delta=1/\delta$;
- $|\Phi|\le M/\delta^2$.

The last bound holds because $\Phi=\tfrac{M}{k}\ddot{\mathbf x}_0$ and $|\ddot{\mathbf x}_0|\le m_1/r_2^2+m_2/r_1^2\le k/\delta^2$. This is an algebraic identity in $(\boldsymbol\rho,z)$ (Lemma 1, item e2 of the script), so it holds also at $z=0$.

**(B2) $\boldsymbol\sigma$, $\boldsymbol\rho$.**
- $|\boldsymbol\sigma(s)-\boldsymbol\sigma(0)|\le\int_0^s|u|^2|\Phi|\,ds\le (M/\delta^2)\,(t(s)-t_0)\le M\Delta/\delta^2$. So $|\boldsymbol\sigma|\le\Sigma:=|\dot{\boldsymbol\rho}(t_0)|+M\Delta/\delta^2$.
- Likewise $|\boldsymbol\rho(s)-\boldsymbol\rho(0)|\le\Sigma\Delta$.

**(B3) $h$.** By Lemma 1(2), $m h=H-\tfrac12\mu|\boldsymbol\sigma|^2+\mu MP$. Hence $|h|\le A:=\big(|H|+\tfrac12\mu\Sigma^2+\mu M/\delta\big)/m$.

**(B4) $u$, $w$.**
- Let $R:=|u|^2$, which is $r$. At points with $R\ge1$ there is no collision. There, by Lemma 1(3)–(4),
$$\Big|\frac{dR}{dt}\Big|\le|\dot z|=\sqrt{2h+2k/R}\le\sqrt{2A+2k}=:V .$$
- So $R(s)\le R_0+V\Delta$ for all $s$, where $R_0=\max(1,r(t_0))$. To see this, consider the last time before $s$ at which $R\le R_0$; if there is none, start from $s=0$. On the interval after it, $R>R_0\ge1$, so the speed bound applies.
- Put $U:=\sqrt{R_0+V\Delta}$. Then $|u|\le U$.
- By Lemma 1(1), $|w|^2=\tfrac12(k+h|u|^2)\le\tfrac12(k+AU^2)=:W_b^2$.

**Case A: $S^*<\infty$.**
- By (B2)–(B4) and the hypothesis, $Y(s)$ stays for all $s\in[0,S^*)$ in
$$\mathcal K=\{t\in[t_0,T^*],\ |u|\le U,\ |w|\le W_b,\ |h|\le A,\ |\boldsymbol\rho|\le|\boldsymbol\rho(t_0)|+\Sigma\Delta,\ |\boldsymbol\sigma|\le\Sigma,\ r_1\ge\delta,\ r_2\ge\delta\}.$$
- $\mathcal K$ is closed and bounded, so it is compact, and $\mathcal K\subset D$.
- This contradicts the extension theorem: a maximal solution with $S^*<\infty$ leaves every compact subset of $D$.

**Case B: $S^*=\infty$.** Then $\int_0^\infty|u|^2ds=\Delta<\infty$. We show this integral diverges.
- Choose $\eta\in(0,1]$ with:
  - $A\eta^2\le k/2$;
  - $2\eta+\varepsilon<\sqrt k/2$, where $\varepsilon:=\tfrac A2\eta+\tfrac F2\eta^3$.

  Such $\eta$ exists, since the left-hand sides tend to $0$ as $\eta\to0$.

- *Claim: on every interval $[a,a+1]$ there is a point with $|u|>\eta$.* Suppose instead that $|u|\le\eta$ on all of $[a,a+1]$.
  - There, $|w'|\le\tfrac A2|u|+\tfrac12|u|^3F\le\varepsilon$, so $|w(s)-w(a)|\le\varepsilon$.
  - By Lemma 1(1), $|w(a)|^2=\tfrac12(k+h|u|^2)\ge\tfrac12(k-A\eta^2)\ge k/4$.
  - Hence $|u(a+1)-u(a)|=\big|w(a)+\int_a^{a+1}(w-w(a))\big|\ge\sqrt k/2-\varepsilon>2\eta$.
  - This contradicts $|u(a)|,|u(a+1)|\le\eta$.
- Now let $n\ge1$, and let $s^*\in[n,n+1]$ have $|u(s^*)|>\eta$.
  - Since $|u'|=|w|\le W_b$, we have $|u|\ge\eta/2$ on $[s^*-\ell,s^*+\ell]$, where $\ell=\min(\eta/(2W_b),1)$.
  - This interval lies in $[n-1,n+2]$. So $\int_{n-1}^{n+2}|u|^2ds\ge\eta^2\ell/2=:c>0$.
- The intervals $[n-1,n+2]$ for $n=1,4,7,\dots$ do not overlap. So $\int_0^\infty|u|^2ds=\infty$, a contradiction.

Both cases are impossible. Hence $T^*=\infty$. $\square$

## 3. Effect on the proof of (★)

In Appendix A, replace A.2 for the planar case by the following:
- the Definition above (the continued solution);
- Lemma 2 (regularity, isolated collisions, energy);
- the Proposition, which is (R).

Then:
- Steps 1–3 of A.4 use only the continuity of $r$ and $\dot\rho$, the $C^2$ regularity of $\boldsymbol\rho$, and the energy identity (E) at times with $r>0$. These are Lemma 2(ii)–(iii).
- Step 4 uses (R), the Proposition, with $\delta=\rho(t_0)-\lambda d_1$ from (SEP).
- The Pythagorean problem starts at rest in the plane, so it stays planar.
- The R4 enclosure on $[0,61]$ is a collision-free classical solution, so it is the continued solution there: there $u\neq0$, and the map to (LC) is one to one up to the sign of $u$.

The remaining citations are the ODE theorems and, for attribution only, Levi-Civita (1920).

## 4. Adversarial self-check (where could this break?)

| Attack | Answer |
|---|---|
| Is the "continued solution" the right object? | It is the standard Levi-Civita continuation, adopted as a definition. The theorem's statement already says "continued through binary collisions". A different convention would be a different theorem. |
| Triple collision? | Excluded: $r_1,r_2\ge\delta$ by (SEP). |
| $h\to0$ or $h>0$ (binary becomes unbound)? | No sign of $h$ is used, only $\lvert h\rvert\le A$ on finite time. |
| $\lvert\dot{\boldsymbol\rho}\rvert$ unbounded? | Bounded on finite time by $\lvert\Phi\rvert\le M/\delta^2$. The bound is integrated in $s$ with weight $t'$, so it holds through collisions. |
| Is $\Phi$ (LC) the same function as Newton's $\ddot{\boldsymbol\rho}$ at $z=0$? | $\Phi$ is an explicit function of $(\boldsymbol\rho,z)$, continuous on $D$. The bound in (B1) is algebraic and does not need $z\ne0$. |
| Hidden use of $L\ne0$? | None. $L$ never appears. |
| Circularity of bounds (B3) and (B4)? | (B2) uses only (B1) and $t(s)\le T^*$. (B3) uses (B2). (B4) uses (B3). There is no loop. |
| Is $\mathcal K\subset D$ compact? | It is closed. All coordinates are bounded, including $\boldsymbol\rho$ by (B2). $r_1,r_2\ge\delta>0$. |
| Spatial problem? | Not covered. KS regularization would be needed. Not needed for the application. |
| Algebra errors? | Identities are checked by sympy with exact symbolic reduction. Item 3 is also derived by hand above. |

**What remains cited:** standard ODE theory (existence, uniqueness, extension), and nothing else.
**Status:** written and self-checked; computer-algebra check of the identities (see the script output);
not externally reviewed.
