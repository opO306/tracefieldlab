# 부록 — 탈출 조건 (★): 정리와 자기완결 증명 (한국어 대조본)

*논문 부록 초안. 2026-10-12 작성. 영어본 `PROOF_STAR_FINAL_EN.md`와 절 번호·식·결론이 같다. 외부 검토와 형식 검증(증명 보조기)은 아직 받지 않았다.*

---

## 처음 읽는 분을 위한 쉬운 설명

**비유 (실제 의미와 구분해서 읽을 것).**
- 별 3이 쌍성(별 4·5)이 만든 언덕을 굴러 올라가고 있다. 쌍성이 돌면서 언덕 모양은 계속 바뀐다.
- 그래서 증명은 **"쌍성이 아무리 바뀌어도 생길 수 있는 가장 높은 언덕"**을 미리 계산해 둔다.
- 별 3의 속도가 그 최악의 언덕도 넘을 만큼 크면, 실제 언덕은 반드시 넘는다.

**실제 의미.**
- **보조정리 A·B:** 최악의 언덕은 일렬 배치에서 생긴다. 가벼운 별 4가 별 3 쪽을 향하고 쌍성 간격이 한계 d에 꽉 찬 배치다. 이 배치에서 퍼텐셜(높이)도, 끌어당기는 힘(기울기)도 가장 크다.
  - 증명 방법은 끌어당김을 "기본 모양 블록"(르장드르 다항식)들의 합으로 쪼개는 것이다.
  - 각 블록은 각도에 따라 −1배에서 1배 사이로만 변한다.
  - 최악 배치에서는 모든 블록이 같은 방향으로 더해진다.
- **정리 (★):** 어느 순간 세 조건 (C1) 쌍성 간격 < d₁, (C2) 거리 ρ > λd₁, (C3) 바깥 속도가 최악 언덕을 넘을 만큼 큼이 모두 성립한다고 하자.
  - 그러면 이 조건들이 서로를 지켜 주며 영원히 유지된다.
  - 속도가 유지되면 에너지 보존 때문에 쌍성이 묶여 있고, 쌍성이 묶여 있으면 최악 언덕 계산이 유효해서 속도가 유지된다.
  - 결국 별 3은 적어도 일정한 속도 C로 영원히 멀어진다.
- **우리가 증명하지 않고 빌려 쓰는 것:** "쌍성 두 별이 충돌해도 계산을 이어 갈 수 있고, 그런 충돌이 유한 시간 안에 무한히 쌓이지 않는다"는 고전 결과다(정칙화). §A.2에 정확히 적었다.

**출처 표기에 관한 주의.**
- 이 정리는 **Yoshida(1972)의 Theorem 1**에서 d=m₁m₂/K로 둔 경우다.
- 1972년 논문의 Corollary가 r(0)≤d 조건을 빠뜨린 문제는 **Yoshida 본인이 먼저 1973년 Erratum(PASJ 25, 285)에서, 이어서 1974년 논문(PASJ 26, p.369, Note)에서 정정했다.** 우리 발견이 아니다.
  - 근거: `LITERATURE_RECHECK_KO.md`. 원문·문헌 재확인 세션이 1972·1974 원문 스캔을 직접 판독했다.
  - 1973 Erratum은 조율 세션이 Crossref로 확인했고, 논문 세션이 ADS 스캔으로 판독했다. 이 세션은 Crossref 서지 정보(PASJ 25(2), 285, 1973)만 직접 확인했다.
  - 이 문서의 Yoshida 원문 쪽수·식 번호 인용은 모두 그 판독에 기댄다. 이 세션은 원문 스캔을 직접 보지 않았다.
- 보조정리의 르장드르 전개 증명은 표준 기법이다. 새로움을 주장하지 않는다.

---

## A.0 범위와 출처

여기서 증명하는 조건은 **Yoshida(1972) Theorem 1** [Y72, p.400]에서 $d=m_1m_2/K$로 둔 경우다. 또한 **Yoshida(1974) criterion A′** [Y74, p.369]에서 $d=m_1m_2/K$로 둔 경우이기도 하다(비고 A.4.2).

[Y72, p.402]의 Theorem 1 Corollary는 $r(0)\le d$를 가정하지 않는다. **이 점은 Yoshida 본인이 먼저 Erratum [Y73]에서, 이어서 [Y74, p.369, Note]에서 정정했다.** Erratum [Y73]은 Corollary의 조건 (26)을 $d\ge r(0)$로 바꾸고, [Y72] p.405의 부등식 (40)도 고친다. 뒤의 Note [Y74, p.369]는 [Y72]의 조건 (25)–(26)을 "incorrect"라고 부르고 $\rho(0)\ge ad\ge a\,r(0)$로 바꾼다. 아래 진술은 $r(t_0)<d_1$을 명시적으로 가정하므로 이 문제의 영향을 받지 않는다.

아래 보조정리 A·B는 [Y72] Lemma 1(p.393)의 $n=1$인 경우와 Lemma 2(p.396)에 해당한다. Yoshida는 이를 각각 $n$에 대한 귀납법과 $\theta$-미분 분석으로 증명한다. 우리는 르장드르(다중극) 전개에 기반한 다른 자기완결 증명을 쓴다. 이것은 표준 기법이며, 새로움을 주장하지 않는다.

충돌과 정칙화에 관한 고전 사실(§A.2)만 빼고는 모두 여기서 증명한다. §A.2는 그 사실들을 적고, 증명에서 정확히 어떻게 쓰는지 밝힌다.

---

## A.1 설정과 기호

질점 $m_0,m_1,m_2>0$가 $\mathbb R^2$ 또는 $\mathbb R^3$에서 뉴턴 중력($G=1$)으로 움직인다. 무게중심 좌표계(전체 운동량 0)에서 다룬다. $m_0$는 탈출 후보이고, 나머지 둘은 다음이 되도록 이름을 붙인다.
$$m_1\le m_2 .$$
다음과 같이 정의한다.
$$
M=m_0+m_1+m_2,\qquad \lambda=\frac{m_2}{m_1+m_2},\qquad \nu=\frac{m_1}{m_1+m_2},\qquad
m=\frac{m_1m_2}{m_1+m_2},\qquad \mu=\frac{m_0(m_1+m_2)}{M}.
$$
따라서 $\lambda\ge\nu>0$, $\lambda+\nu=1$, $\mu M\lambda=m_0m_2$, $\mu M\nu=m_0m_1$이다.

**야코비 벡터.**
- $\mathbf r=\mathbf x_2-\mathbf x_1$.
- $\boldsymbol\rho=\mathbf x_0-\mathbf g$이고, $\mathbf g=(m_1\mathbf x_1+m_2\mathbf x_2)/(m_1+m_2)$이다.
- 그러면 $\mathbf x_2-\mathbf g=\nu\mathbf r$, $\mathbf x_1-\mathbf g=-\lambda\mathbf r$이다. 즉 **가벼운 $m_1$이 $\mathbf g$에서 더 멀다.**
- 다음과 같이 쓴다.
$$
r=|\mathbf r|,\quad \rho=|\boldsymbol\rho|,\quad
r_1=|\mathbf x_0-\mathbf x_2|=|\boldsymbol\rho-\nu\mathbf r|,\quad
r_2=|\mathbf x_0-\mathbf x_1|=|\boldsymbol\rho+\lambda\mathbf r| .
$$
- $r>0$이면 $c=\cos\theta=\boldsymbol\rho\cdot\mathbf r/(\rho r)$로 두고, 다음이 성립한다.
$$
r_1^2=\rho^2+\nu^2r^2-2\nu\rho r c,\qquad r_2^2=\rho^2+\lambda^2r^2+2\lambda\rho r c .
$$

**반지름 방향 속도.** $\dot\rho:=\dfrac{d}{dt}|\boldsymbol\rho|=\dfrac{\boldsymbol\rho\cdot\dot{\boldsymbol\rho}}{\rho}$. 반지름 방향 성분이며 속력 $|\dot{\boldsymbol\rho}|\ge|\dot\rho|$가 **아니다**.

**에너지.** 무게중심 좌표계에서, **벡터 속도**로 쓰면 다음과 같다.
$$
H=\tfrac12 m|\dot{\mathbf r}|^2+\tfrac12\mu|\dot{\boldsymbol\rho}|^2-\frac{m_1m_2}{r}-\mu M\Big(\frac{\lambda}{r_1}+\frac{\nu}{r_2}\Big). \tag{E}
$$
$H<0$을 가정하고, $K=-H>0$, $d_1=m_1m_2/K$로 둔다.

**$\boldsymbol\rho$의 운동방정식.**
$$
\ddot{\boldsymbol\rho}=-M\Big(\frac{\lambda}{r_1^3}+\frac{\nu}{r_2^3}\Big)\boldsymbol\rho+\lambda\nu M\Big(\frac1{r_1^3}-\frac1{r_2^3}\Big)\mathbf r=:\Phi(\boldsymbol\rho,\mathbf r). \tag{EM}
$$

**반지름 방향 가속도.** $\rho>0$인 곳에서 라그랑주 항등식에 의해 다음이 성립한다.
$$
\ddot\rho=\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}+\frac{|\boldsymbol\rho|^2|\dot{\boldsymbol\rho}|^2-(\boldsymbol\rho\cdot\dot{\boldsymbol\rho})^2}{\rho^3}\;\ge\;\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}. \tag{RA}
$$
버린 항은 $|\boldsymbol\rho\times\dot{\boldsymbol\rho}|^2/\rho^3\ge0$이다. 각운동량에 대해서는 아무것도 가정하지 않으므로 $L=0$이어도 된다.

**보조 함수.** $\rho>0$, $r\ge0$, $c\in[-1,1]$이고 $r_1,r_2>0$일 때 다음과 같이 둔다.
$$
P(\rho,r,c)=\frac{\lambda}{r_1}+\frac{\nu}{r_2},\qquad
F(\rho,r,c)=\frac{\lambda(\rho-\nu rc)}{r_1^3}+\frac{\nu(\rho+\lambda rc)}{r_2^3}=-\partial_\rho P .
$$
- (E)의 퍼텐셜 항은 $\mu M P$이다.
- (EM)에서 $\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}/\rho=-M\,F(\rho,r,c)$이다.
- $r=0$이면 $c$와 상관없이 $P=1/\rho$, $F=1/\rho^2$이다.

$d\ge0$, $\rho>\lambda d$에 대해 다음과 같이 둔다.
$$
W(\rho,d)=\frac{\lambda}{\rho+\nu d}+\frac{\nu}{\rho-\lambda d},\qquad
\hat F(\rho,d)=\frac{\lambda}{(\rho+\nu d)^2}+\frac{\nu}{(\rho-\lambda d)^2}=-\partial_\rho W(\rho,d).
$$
$W(\rho,d)=P(\rho,d,-1)$, $\hat F(\rho,d)=F(\rho,d,-1)$이다. 기하학적으로는 가벼운 성분 $m_1$이 $\mathbf g$에서 거리 $\lambda d$에 있으면서 $\mathbf g$와 $m_0$ 사이에 놓인 일렬 배치다.

---

## A.2 충돌에 관한 고전 사실 (인용, 여기서 증명하지 않음)

**(S1) Painlevé** [P1897]. 3체 문제의 모든 특이점은 충돌이다. 정확히 말하면 다음과 같다. $[t_0,T)$, $T<\infty$에서 정의된 해를 $T$ 너머로 연장할 수 없다면, $t\to T$일 때 위치들이 수렴하고, 극한에서 적어도 두 위치가 일치한다.

**(S2) 쌍성 충돌의 정칙화** [Su1912; LC1920 (평면); KS1965 (공간)]. 해가 시각 $t_c$에 $m_1$–$m_2$ 충돌로 다가가고, $t_c$ 근처에서 $r_1,r_2$가 0에서 떨어져 있다고 하자. 그러면 다음이 성립한다.
- 해는 정칙화 방정식의 해로서 $t_c$ 너머로 유일하게 연장된다(물리적으로는 탄성 반사).
- 연장된 해에서 $\boldsymbol\rho$, $\dot{\boldsymbol\rho}$는 $t$에 대해 연속이다.
- $t\neq t_c$에서 (EM)이 성립하고, $H$는 보존된다.

**(S3) 충돌이 쌓이지 않음** [Su1912; SM1971, 1장]. $m_0$가 $m_1$, $m_2$ 모두로부터 어떤 $\delta>0$ 이상 떨어져 있는 동안, 정칙화된 해의 쌍성 충돌 시각은 유한 시각에 쌓일 수 없다.

**약속.** "해"는 (S2)로 쌍성 충돌을 넘어 연장한 해를 뜻한다. $[t_0,T^*)$를 그 최대 존재 구간으로 둔다.

**우리가 가정하는 것.** (S1)–(S3)은 증명에서 다음 형태로만 쓰인다.

> **(R)** $\inf_{[t_0,T^*)}\min(r_1,r_2)>0$이면 $T^*=\infty$이다.

*(R)은 (S1)–(S3)에서 나온다.* $T^*<\infty$이고 $[t_0,T^*)$에서 $\min(r_1,r_2)\ge\delta>0$이라고 하자.
- (S3)에 의해 $[t_0,T^*)$의 쌍성 충돌은 유한 개다. 따라서 어떤 $t_\ell$ 이후 $[t_\ell,T^*)$에서 해는 고전적이다.
- (S1)에 의해 $T^*$에서 충돌로 끝난다. $r_1,r_2\ge\delta$이므로 그 충돌은 $m_1$–$m_2$ 충돌일 수밖에 없다.
- (S2)에 의해 해는 $T^*$ 너머로 연장된다. 이는 최대성에 모순이다.

**아래에서 쓰는 정칙성.** $[t_0,T^*)$ 안의 구간 $I$에서 $\min(r_1,r_2)$가 양의 상수로 아래로 유계라고 하자. 그러면 $I$에서 $t\mapsto\boldsymbol\rho(t)$는 $C^2$이고 $\ddot{\boldsymbol\rho}=\Phi(\boldsymbol\rho,\mathbf r)$이다. 이유는 다음과 같다.
- $\Phi$는 $\{r_1>0,\,r_2>0\}$에서 연속이다. $\mathbf r=0$인 점에서도 그렇다.
- (S2)에 의해 $t\mapsto(\boldsymbol\rho(t),\mathbf r(t))$가 연속이다. 따라서 $t\mapsto\Phi(\boldsymbol\rho(t),\mathbf r(t))$도 연속이다.
- $\dot{\boldsymbol\rho}$는 연속이고, 충돌이 없는 각 구간에서 이 함수의 원시함수다. 충돌 시각은 (S3)에 의해 고립되어 있으므로 전체에서 원시함수다.

따라서 $\rho>0$인 곳에서 $\rho$는 $C^2$이고 (RA)가 성립한다.

**인용에 관한 주의.** [P1897, Su1912, LC1920, KS1965, SM1971]의 서지 정보는 확인했다. 그러나 인용한 진술이 원문에 어떤 문구로 있는지는 **확인하지 않았다**.
- Sundman의 전역 결과는 전체 각운동량이 0이 아니라고 가정한다. 우리가 읽기로는, 그 가정은 3중 충돌을 배제하고 충돌하는 쌍에서 셋째 천체를 떼어 놓는 역할을 한다.
- 우리 적용 대상은 $L=0$이다. 대신 A.4 2단계에서 증명하는 사전 하한 (SEP)이 그 분리를 보장한다.
- $L=0$에서 이런 사전 하한이 있을 때의 (S3)을 그대로 적은 진술은 찾지 못했다. 비고 A.4.5에 그 이유를 스케치했지만, 증명의 일부는 아니다.

---

## A.3 두 부등식

> **보조정리 A.** $\lambda\ge\nu>0$, $\lambda+\nu=1$, $d\ge0$, $\rho>\lambda d$라고 하자. 그러면 모든 $0\le r\le d$, $c\in[-1,1]$에 대해
> $$P(\rho,r,c)\le W(\rho,d).$$
>
> **보조정리 B.** 같은 가정에서
> $$F(\rho,r,c)\le \hat F(\rho,d).$$

두 부등식 모두 $(r,c)=(d,-1)$에서 등호이므로 더 좋은 상수는 없다.

보조정리 B와 (EM)을 합치면 다음을 얻는다. 해를 따라 $\rho>\lambda d$이고 $r\le d$이면
$$
\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}=-M F\ \ge\ -M\hat F(\rho,d)=M\,\partial_\rho W(\rho,d). \tag{B′}
$$

**르장드르 다항식 $P_n$에 관한 사실.**
- (L1) 생성함수 [DLMF, 식 18.12.11]: $x\in[-1,1]$, $|z|<1$에서 $(1-2xz+z^2)^{-1/2}=\sum_{n\ge0}P_n(x)z^n$.
- (L2) 라플라스 적분 [DLMF, 식 18.10.5]: $P_n(\cos\vartheta)=\frac1\pi\int_0^\pi(\cos\vartheta+i\sin\vartheta\cos\phi)^n\,d\phi$. $|\cos\vartheta+i\sin\vartheta\cos\phi|^2=\cos^2\vartheta+\sin^2\vartheta\cos^2\phi\le1$이므로 $[-1,1]$에서 $|P_n(x)|\le1$이다.
- (L3) (L1)에서 $(x,z)$를 $(-x,-z)$로 바꾸면 $P_n(-x)=(-1)^nP_n(x)$이다. $x=1$을 넣으면 $\sum P_n(1)z^n=(1-z)^{-1}$이므로 $P_n(1)=1$이다. 따라서 $P_n(-1)=(-1)^n$이다.

**보조정리 A·B의 증명.**

*1단계 (정규화).*
- $s=r/\rho$, $a=\nu s$, $b=\lambda s$로 둔다. $\lambda r\le\lambda d<\rho$이므로 $0\le a\le b<1$이다.
- $r_1=\rho\sqrt{1-2ac+a^2}$, $\rho-\nu rc=\rho(1-ac)$이다.
- $r_2=\rho\sqrt{1-2b(-c)+b^2}$, $\rho+\lambda rc=\rho(1-b(-c))$이다.
- 따라서
$$
\rho P=\lambda G(a,c)+\nu G(b,-c),\qquad \rho^2F=\lambda g(a,c)+\nu g(b,-c),
$$
여기서
$$
G(a,x)=(1-2ax+a^2)^{-1/2},\qquad g(a,x)=(1-ax)(1-2ax+a^2)^{-3/2}.
$$

*2단계 (전개).*
- (L1)에 의해 $G(a,x)=\sum_n a^nP_n(x)$이다.
- 직접 계산하면 $g=G+a\,\partial_aG$이다. 따라서 $g(a,x)=\sum_n(n+1)a^nP_n(x)$이다.
- $0\le a<1$, $|x|\le1$이면 (L2)에 의해 $\sum(n+1)a^n|P_n(x)|\le(1-a)^{-2}$이다. 그래서 두 급수는 절대수렴하고, 항별로 더하고 묶을 수 있다.

*3단계 (합치기).* 1·2단계와 (L3)에서
$$
\rho^2F=\sum_{n\ge0}(n+1)s^n\beta_nP_n(c),\qquad \rho P=\sum_{n\ge0}s^n\beta_nP_n(c),\qquad
\beta_n:=\lambda\nu^n+\nu(-\lambda)^n .
$$
- $\beta_0=1$이고, $\beta_1=\lambda\nu-\nu\lambda=0$이다($\mathbf g$ 기준이라 쌍극자 항이 없다).
- **모든 $n$에서 $(-1)^n\beta_n\ge0$이다.**
  - $n$이 짝수면 $(-1)^n\beta_n=\lambda\nu^n+\nu\lambda^n>0$이다.
  - $n$이 홀수면 $(-1)^n\beta_n=\nu\lambda^n-\lambda\nu^n=\lambda\nu(\lambda^{n-1}-\nu^{n-1})\ge0$이다. $\lambda\ge\nu$이기 때문이다.
  - $m_1\le m_2$를 쓰는 곳은 **여기뿐**이다.
- 따라서 $|\beta_n|=(-1)^n\beta_n$이다.

*4단계 (각도에 대한 최댓값).* (L2)와 3단계에서
$$
\rho^2F\le\sum_n(n+1)s^n|\beta_n|=\lambda\sum_n(n+1)(-\nu s)^n+\nu\sum_n(n+1)(\lambda s)^n
=\frac{\lambda}{(1+\nu s)^2}+\frac{\nu}{(1-\lambda s)^2}.
$$
같은 방식으로
$$
\rho P\le\sum_n s^n|\beta_n|=\frac{\lambda}{1+\nu s}+\frac{\nu}{1-\lambda s}.
$$
우변은 각각 $\rho^2F(\rho,r,-1)$, $\rho P(\rho,r,-1)$와 같다.

*5단계 (간격에 대한 최댓값).*
- 각 항 $(n+1)s^n|\beta_n|$, $s^n|\beta_n|$은 $s\ge0$에서 줄어들지 않는다.
- $r\le d$이므로 $s\le s_d:=d/\rho<1/\lambda$이다. 그래서 $s_d$에서도 급수가 수렴한다.
- 따라서
$$
\rho^2F\le\frac{\lambda}{(1+\nu s_d)^2}+\frac{\nu}{(1-\lambda s_d)^2}=\rho^2\hat F(\rho,d),\qquad
\rho P\le\frac{\lambda}{1+\nu s_d}+\frac{\nu}{1-\lambda s_d}=\rho\,W(\rho,d).
$$
각각 $\rho^2$, $\rho$로 나누면 보조정리가 나온다. $\square$

**범위에 관한 비고.**
- (i) $\rho$의 상한은 필요 없다.
- (ii) $\rho\to\lambda d^+$이면 상계가 발산하지만, $\rho>\lambda d$인 모든 점에서 부등식은 성립한다. 균등성은 쓰지 않는다.
- (iii) $r=0$(쌍성 충돌 순간)도 포함한다.
- (iv) $c=\cos\theta$만 들어가므로 평면과 공간 모두에서 성립한다.
- (v) $\lambda=\nu$이면 홀수 $\beta_n$이 모두 0이어서 $c=+1$에서도 등호가 성립한다.

**비고 (보조정리 A의 두 번째 증명).** $r$을 고정하자.
- $\partial_cP=\lambda\nu\rho r\,(r_1^{-3}-r_2^{-3})$이고, $r_2^2-r_1^2=r\big((\lambda-\nu)r+2\rho c\big)$이다.
- 그래서 $P$는 $c^*=-(\lambda-\nu)r/(2\rho)$까지 감소하고 그 뒤로 증가한다. $[-1,1]$에서의 최댓값은 $c=\pm1$에서 나온다.
- $P(\rho,r,-1)-P(\rho,r,1)=2\lambda\nu r\big[(\rho^2-\lambda^2r^2)^{-1}-(\rho^2-\nu^2r^2)^{-1}\big]\ge0$이다.
- 끝으로 $P(\rho,r,-1)=W(\rho,r)$이고 $\partial_dW=\lambda\nu\big[(\rho-\lambda d)^{-2}-(\rho+\nu d)^{-2}\big]\ge0$이다.

---

## A.4 탈출 조건

> **정리 (★).** $m_1\le m_2$, $H<0$이고 무게중심 좌표계라고 하자. 어느 시각 $t_0$에서 다음이 성립한다고 하자.
>
> - **(C1)** $r(t_0)<d_1=m_1m_2/K$
> - **(C2)** $\rho(t_0)>\lambda d_1$
> - **(C3)** $\dot\rho(t_0)>\sqrt{2M\,W(\rho(t_0),d_1)}$
>
> $C:=\sqrt{\dot\rho(t_0)^2-2M\,W(\rho(t_0),d_1)}>0$로 두자. 그러면 해(A.2처럼 쌍성 충돌을 넘어 연장한 해)는 모든 $t\ge t_0$에서 존재한다. 또한 모든 $t\ge t_0$에서 다음이 성립한다.
> 1. $r(t)\le \dfrac{m_1m_2}{K+\frac12\mu C^2}<d_1$
> 2. $\tfrac12\dot\rho(t)^2-M\,W(\rho(t),d_1)\ge\tfrac12C^2$이고 $\dot\rho(t)>C$
> 3. $\rho(t)\ge\rho(t_0)+C\,(t-t_0)$
>
> 특히 $\rho\to\infty$이고 쌍성 간격은 유계다. 즉 쌍곡–타원 운동이다.

**증명.**

*1단계 (부트스트랩 집합).* 다음과 같이 둔다.
$$
\mathcal G=\big\{t\in[t_0,T^*):\ \text{모든 } s\in[t_0,t]\text{에서 } r(s)\le d_1,\ \dot\rho(s)\ge0\big\},\qquad \tau=\sup\mathcal G .
$$
$\mathcal G$는 $t_0$를 포함하는 구간이다. (C1)·(C3)에서 $r(t_0)<d_1$, $\dot\rho(t_0)>0$이고 $r$, $\dot\rho$는 연속이다. 그래서 $\tau>t_0$이다.

*2단계 ($\mathcal G$ 위의 추정).* $t\in\mathcal G$를 하나 고르자.

(a) 분리.
- $[t_0,t]$에서 $\rho$는 줄어들지 않는다. 그래서 $s\in[t_0,t]$에 대해 $\rho(s)\ge\rho(t_0)>\lambda d_1\ge\lambda r(s)$이다.
- $\delta:=\rho(t_0)-\lambda d_1>0$로 두면, $[t_0,t]$에서
$$
r_2\ge\rho-\lambda r\ge\delta,\qquad r_1\ge\rho-\nu r\ge\rho(t_0)-\nu d_1\ge\delta. \tag{SEP}
$$

(b) 반지름 방향 속도.
- $s\in[t_0,t]$에서 (RA)와 (B′)($d=d_1$)에 의해 $\ddot\rho(s)\ge M\,\partial_\rho W(\rho(s),d_1)$이다.
- 양변에 $\dot\rho(s)\ge0$을 곱하면 $\dot\rho\ddot\rho\ge\frac{d}{ds}\big[M\,W(\rho(s),d_1)\big]$이다.
- 적분하면($\rho\in C^2$ 사용, A.2 참조)
$$
\tfrac12\dot\rho(t)^2-M\,W(\rho(t),d_1)\ \ge\ \tfrac12\dot\rho(t_0)^2-M\,W(\rho(t_0),d_1)=\tfrac12 C^2 .
$$
- $W>0$이고 $\dot\rho(t)\ge0$이므로 $\dot\rho(t)>C$이다.

(c) 쌍성 간격.
- $r(t)=0$이면 할 일이 없다.
- 아니면 다음을 쓴다: (E), $\tfrac12 m|\dot{\mathbf r}|^2\ge0$, $|\dot{\boldsymbol\rho}|\ge|\dot\rho|$, 보조정리 A($r(t)\le d_1$, $\rho(t)>\lambda d_1$이므로 쓸 수 있다), (b).
$$
\frac{m_1m_2}{r(t)}=K+\tfrac12 m|\dot{\mathbf r}|^2+\tfrac12\mu|\dot{\boldsymbol\rho}|^2-\mu M P
\ \ge\ K+\mu\Big(\tfrac12\dot\rho^2-M\,W(\rho,d_1)\Big)\ \ge\ K+\tfrac12\mu C^2 .
$$
- 따라서 $r(t)\le m_1m_2/(K+\tfrac12\mu C^2)<d_1$이다.

*3단계 ($\tau=T^*$).* $\tau<T^*$라고 가정하자.
- 연속성에 의해 $r(\tau)\le d_1$, $\dot\rho(\tau)\ge0$이다. 따라서 $\tau\in\mathcal G$이다.
- 2단계에 의해 $r(\tau)<d_1$, $\dot\rho(\tau)>C>0$이다.
- 연속성에 의해 어떤 $\varepsilon>0$에 대해 $[\tau,\tau+\varepsilon)$에서도 이 엄격한 부등식이 유지된다. 그러면 $\tau+\varepsilon/2\in\mathcal G$가 되어 $\tau=\sup\mathcal G$에 모순이다.
- 따라서 $\mathcal G=[t_0,T^*)$이다.

*4단계 (전역 존재).* (SEP)에 의해 $[t_0,T^*)$에서 $\min(r_1,r_2)\ge\delta>0$이다. (R)에 의해 $T^*=\infty$이다.

*5단계 (결론).* 2–4단계에 의해 결론 1, 2가 모든 $t\ge t_0$에서 성립한다. 결론 3은 $\dot\rho>C$를 적분하면 나온다. $\square$

### 비고

**A.4.1 (각운동량).** 정보를 버리는 곳은 (RA) 한 곳뿐이고, 버린 항은 0 이상이다. 그래서 $L=0$(피타고라스 문제)에도 그대로 적용된다.

**A.4.2 (Yoshida 진술과의 관계).** 아래 관계는 [Y72], [Y73], [Y74]의 진술을 원문 스캔에서 옮긴 것(`LITERATURE_RECHECK_KO.md` §2)에 기댄다.

*[Y72] Theorem 1.* 가정은 다음과 같다.
- (13) 어떤 $a>\lambda$에 대해 $\rho(0)\ge ad>a\,r(0)$
- (14) $\dot\rho(0)>\sqrt{2MW(\rho(0),d)}$
- (15) $\dot\rho(0)^2-2MW(\rho(0),d)>2H/\mu+2m_1m_2/(\mu d)$

$d=d_1$로 두면 다음이 성립한다.
- (15)의 우변이 0이 되어(`final_doc_checks.py` F3), (15)는 (14)에서 나온다.
- (13)은 (C1)과 (C2)를 합친 것과 동치다. $a\in(\lambda,\rho(t_0)/d_1]$을 고르면 된다.

따라서 (★)는 [Y72] Theorem 1에서 $d=m_1m_2/K$로 둔 것이다.

*[Y74] criterion A′.* $d=d_1$로 두면 $\rho(0)\ge ad_1\ge a\,r(0)$를 요구하므로 경계 경우 $r(0)=d_1$을 허용한다. 그러나 이 경우는 (C2)·(C3)과 함께 일어날 수 없다. $t=t_0$에서 2단계 (c)를 적용하면 $r(t_0)<d_1$이 되기 때문이다.

*[Y72] Corollary.* A.0에서 말한 빈틈이 있고, [Y73]에서, 이어서 [Y74, p.369]에서 정정됐다. 우리는 이것을 쓰지 않는다.

**A.4.3 (질량 순서).** $m_1>m_2$이면 쌍성의 이름을 바꿔 붙인다. $W$에서 $\lambda$와 $\nu$를 바꾸면 함수가 작아지고, 보조정리 A·B가 깨진다(예: $(r,c)=(d,-1)$). 순서는 본질적이다.

**A.4.4 (검사 가능성).** (C1)–(C3)은 $(r,\rho,\dot\rho)$에 대한 엄격한 부등식이고 상수는 정확하다. 한 시각 $t_0$에서 참 상태를 감싸는 엄밀한 상자의 모든 상태가 이 조건을 만족하면 충분하다. $K$는 참 해의 정확한 에너지여야 한다. 상자에서 다시 계산한 에너지를 쓰면 안 된다.

**A.4.5 ((S3)이 여기서 성립하는 이유 — 스케치, 증명의 일부 아님).** 유한 구간에서 (SEP)가 성립한다고 하자.

섭동받는 쌍성.
- 상대 운동은 $\ddot{\mathbf r}=-k\,\mathbf r/r^3+\mathbf f$이다. 여기서 $k=m_1+m_2$, $\mathbf f=m_0\big[(\boldsymbol\rho-\nu\mathbf r)/r_1^3-(\boldsymbol\rho+\lambda\mathbf r)/r_2^3\big]$이다.
- (SEP)에 의해 $|\mathbf f|\le 2m_0/\delta^2$이다.
- 또 $|\ddot{\boldsymbol\rho}|\le M/\delta^2$이므로, 유한 구간에서 $|\dot{\boldsymbol\rho}|$는 유계다.

케플러 에너지의 범위.
- $h=\tfrac12|\dot{\mathbf r}|^2-k/r$(단위 환산질량당 케플러 에너지)로 두면 $mh=H-\tfrac12\mu|\dot{\boldsymbol\rho}|^2+\mu MP$이다.
- 위쪽: 2단계 (b)에 의해 $mh\le-K-\tfrac12\mu C^2$이다.
- 아래쪽: $mh\ge-K-\tfrac12\mu|\dot{\boldsymbol\rho}|^2$이고, 유한 구간에서 아래로 유계다.
- 따라서 $0<a\le A<\infty$인 상수로 $h\in[-A,-a]$이다.

레비치비타 변수(평면).
- $\mathbf r\simeq z=u^2\in\mathbb C$, $dt=|u|^2ds$로 둔다.
- 직접 계산하면(`final_doc_checks.py` F1–F2로 확인) 다음이 나온다.
$$
u''=\tfrac h2u+\tfrac12|u|^2\bar u f,\qquad |u'|^2=\tfrac12\big(k+h|u|^2\big).
$$
- 그래서 $u$는 섭동받는 조화진동자다. 진동수 $\sqrt{-h/2}$는 위아래로 유계이고, $u$는 $|u'|=\sqrt{k/2}$로 0을 지난다.
- 따라서 한 바퀴의 가상 시간은 위로 유계이고, 한 바퀴가 기여하는 물리 시간 $\int|u|^2ds$는 아래로 유계다. 그래서 유한 시각에 충돌이 쌓일 수 없다.

엄밀하게 하려면 표준 섭동 추정이 필요한데, 여기서는 옮기지 않았다.

---

## A.5 피타고라스 문제에 적용

**초기 조건** [Bu1913; SP1967]. 질량 3, 4, 5가 각각 $(1,3)$, $(-2,-1)$, $(1,-1)$에서 정지 상태로 출발한다.
- 무게중심은 원점이고, 전체 운동량과 각운동량은 0이다.
- 서로의 거리는 $|3\text{–}4|=5$, $|3\text{–}5|=4$, $|4\text{–}5|=3$이다.
- 따라서 정확히
$$
H=-\Big(\frac{3\cdot4}{5}+\frac{3\cdot5}{4}+\frac{4\cdot5}{3}\Big)=-\frac{769}{60},\qquad K=\frac{769}{60}.
$$
- 에너지는 (정칙화된) 해를 따라 보존되므로, 이후 어느 시각에서든 이 $K$를 쓴다.

**매개변수.** 탈출체는 $m_0=3$, 쌍성은 $m_1=4\le m_2=5$이다.
$$
M=12,\quad \lambda=\tfrac59,\quad \nu=\tfrac49,\quad \mu=\tfrac94,\quad m=\tfrac{20}9,\quad
d_1=\frac{m_1m_2}{K}=\frac{1200}{769}\approx1.560468140442,
$$
$$
\lambda d_1=\frac{2000}{2307}\approx0.866926744690,\qquad \nu d_1=\frac{1600}{2307}\approx0.693541395752 .
$$
(정확값과 소수는 `final_doc_checks.py` F4로 확인했다.)

**조건.** 질량 3, 4, 5의 위치를 $\mathbf x_3,\mathbf x_4,\mathbf x_5$로 쓰면
$$
r=|\mathbf x_5-\mathbf x_4|,\qquad \boldsymbol\rho=\mathbf x_3-\frac{4\mathbf x_4+5\mathbf x_5}{9},\qquad \dot\rho=\frac{\boldsymbol\rho\cdot\dot{\boldsymbol\rho}}{\rho}
$$
이고, 정리 (★)의 조건은 다음과 같다.
$$
\text{(C1)}\ r<\frac{1200}{769},\qquad \text{(C2)}\ \rho>\frac{2000}{2307},\qquad
\text{(C3)}\ \dot\rho>\sqrt{24\,W(\rho)},\quad W(\rho)=\frac{5/9}{\rho+1600/2307}+\frac{4/9}{\rho-2000/2307}.
$$
$W$에서 $5/9$는 **더 먼** 거리 $\rho+\nu d_1$에, $4/9$는 **더 가까운** 거리 $\rho-\lambda d_1$에 붙는다.

**결론.** 어느 $t_0$에서 참 상태가 (C1)–(C3)을 만족하면, 모든 $t\ge t_0$에서 다음이 성립한다.
- 별 4–5의 간격은 $m_1m_2/(K+\tfrac12\mu C^2)<1200/769$ 아래에 머문다.
- 질량 3인 별은 적어도 선형으로 멀어진다: $\rho(t)\ge\rho(t_0)+C(t-t_0)$.

여기서 $\mu=9/4$이고, $C$는 $t_0$의 상태로 정한다.

실제 궤적에서 (C1)–(C3)을 확인하는 일은 본문에서 다루는 별도의 컴퓨터 보조 단계다. 이 부록에 포함되지 않는다.

**역사 비고.** Yoshida [Y72, p.406]는 이미 이전의 비엄밀 수치해를 써서 $t=59.72$에서 Corollary를 이 문제에 적용했다. 그 지점에서 $r\approx0.42<d_1$이었으므로 정정된 조건도 만족했다(`LITERATURE_RECHECK_KO.md` §2).

**프로젝트 메모 (논문 본문에는 넣지 않음).**
- 판정기 Y1이 위 정의를 따르는지는 이 세션이 직접 보지 않았다. 확인한 쪽은 조율 세션이다.
  - 조율 세션의 전달(2026-10-12)에 따르면 Y1 보고서의 식은 $\mathbf g=(m_1q_1+m_2q_2)/(m_1+m_2)$, $W=\lambda/(\rho+\nu d)+\nu/(\rho-\lambda d)$, $\dot\rho=\boldsymbol\rho\cdot\dot{\boldsymbol\rho}/\rho$이다.
  - 재유도 값과 $3\times10^{-12}$ 이내로 일치했다고 한다.
- 이 문서는 `PROOF_REVIEW_KO.md` §7에서 지적한 ①②③을 반영한 정리본이다.
  - ① 에너지식을 벡터형 (E)로 썼다.
  - ② 표준 결과 (S1)–(S3)과 가정 (R)을 명시했다.
  - ③ 부트스트랩 집합 $\mathcal G$와 $\tau=\sup\mathcal G$를 정의했다.

---

## A.6 독립 교차 점검 (증명의 일부 아님)

아래 점검은 위 증명의 오류를 찾기 위한 것이며, 증명을 대신하지 않는다.
- **기호 계산(SymPy).** A.1과 A.3에서 쓴 대수 항등식을 모두 확인했다. 직교좌표 뉴턴 방정식에서 (E)·(EM) 유도, 2D·3D에서 (RA), 1–3단계 전개(유한 차수까지), 4단계 닫힌꼴, A의 두 번째 증명, A.4.5의 레비치비타 관계가 포함된다. 일부러 틀린 식은 모두 거부됐다.
  - 파일: `symbolic_checks.py`(33개, 음성 대조 4개), `final_doc_checks.py`(7개, 음성 대조 1개).
- **적대적 수치 탐색.** 50자리로 40,000점을 검사했다. 경계 근처($\rho\to\lambda d^+$, $r\in\{0,d\}$, $c=\pm1$)에 점을 모았고, 여러 $\lambda\in[\tfrac12,1)$를 썼다.
  - 보조정리 A·B 위반은 없었다. 최대 비율은 $1+O(10^{-50})$이었고, $(r,c)=(d,-1)$에서 나왔다.
  - 파일: `lemma_adversarial.py`.
- **동역학 반증 시험.** 질량 3, 4, 5이고 에너지가 여러 값인 상태 중 (C1)–(C3)을 만족하는 800개를 $t\le40$까지 적분했다.
  - 정리 (★)의 결론이 깨진 경우는 없었다.
  - (C3)에서 2배를 빼거나, (C1)을 빼거나, $W$에서 $\lambda,\nu$를 바꾼 변형에서는 결론이 깨졌다. 탈출체가 되돌아오는 경우도 있었다.
  - 파일: `star_dynamic_falsify.py`, `dynamic_summary.txt`. 자세한 내용은 `PROOF_REVIEW_KO.md` §6.

---

## 참고문헌

- **[Bu1913]** C. Burrau, *Numerische Berechnung eines Spezialfalles des Dreikörperproblems*, Astron. Nachr. **195** (1913) 113–118. doi:10.1002/asna.19131950602.
- **[DLMF]** NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov/, 식 18.10.5, 18.12.11 (2026-10-12 접속).
- **[KS1965]** P. Kustaanheimo, E. Stiefel, *Perturbation theory of Kepler motion based on spinor regularization*, J. Reine Angew. Math. **218** (1965) 204–219.
- **[LC1920]** T. Levi-Civita, *Sur la régularisation du problème des trois corps*, Acta Math. **42** (1920) 99–144. doi:10.1007/BF02404404.
- **[P1897]** P. Painlevé, *Leçons sur la théorie analytique des équations différentielles, professées à Stockholm (1895)*, Paris, 1897.
- **[SM1971]** C. L. Siegel, J. K. Moser, *Lectures on Celestial Mechanics*, Grundlehren der mathematischen Wissenschaften **187**, Springer, Berlin, 1971 (Classics in Mathematics 재간행, 1995).
- **[Su1912]** K. F. Sundman, *Mémoire sur le problème des trois corps*, Acta Math. **36** (1912) 105–179. doi:10.1007/BF02422379.
- **[SP1967]** V. Szebehely, C. F. Peters, *Complete solution of a general problem of three bodies*, Astron. J. **72** (1967) 876. doi:10.1086/110355.
- **[Y72]** J. Yoshida, *Improved criteria for hyperbolic-elliptic motion in the general three-body problem*, Publ. Astron. Soc. Japan **24** (1972) 391–408. doi:10.1093/pasj/24.3.391.
- **[Y73]** J. Yoshida, *Erratum: Improved criteria for hyperbolic-elliptic motion in the general three-body problem*, Publ. Astron. Soc. Japan **25** (1973) 285. doi:10.1093/pasj/25.2.285.
- **[Y74]** J. Yoshida, *Improved criteria for hyperbolic-elliptic motion in the general three-body problem. II*, Publ. Astron. Soc. Japan **26** (1974) 367–377.

**참고문헌 확인 수준.**
- 이 세션은 웹 검색으로 서지 정보를 확인했다: 권, 쪽, DOI.
- 원문 본문은 읽지 않았다.
- [LC1920]의 DOI는 처음에 1906년 다른 논문(Acta Math. 30, 305–327)의 것으로 잘못 적었다. 조율 세션의 지적을 받아, 이 세션이 Crossref로 10.1007/BF02404404(Acta Math. 42, 99–144, 1920)를 확인하고 고쳤다.
- 예외는 [Y72], [Y73], [Y74]의 쪽·식 번호다. 이것은 원문·문헌 재확인 세션이 원문 스캔을 직접 판독한 결과(`LITERATURE_RECHECK_KO.md`)에서 옮겼다.
- 불확실한 서지 항목:
  - [KS1965] 끝쪽 219는 출처 하나에서만 확인했다.
  - [Bu1913] 시작쪽은 출처에 따라 113 또는 114로 다르다. ADS 서지코드가 113이라 113으로 적었다.
  - [P1897]의 출판사는 확인하지 못해서 적지 않았다.
