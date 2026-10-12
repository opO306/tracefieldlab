# 판정 표 초안 (구현자 보고서·검사기 열람 전에 작성, 2026-10-12 KST 15:4x 경)

> 이 파일은 독립성 기록용이다. `R3C_RESULT_KO.md`, `check_r3c.py`, `R3C_CHECK_REPORT_KO.md`, 공유 노트의
> r3c-impl / r3c-check 인계 노트를 **읽기 전에** 코드(lohner_step.py, dyad.py, taylor_step.py, r3c_run.py,
> r4_run.py)와 명세(GATE_R3C §3, R0, R1)만 보고 쓴 판정이다. 최종 보고서의 표는 이 초안에서 시험 결과만 덧붙인다.

| # | 단계 | 코드 위치 | 근거(요지) | 판정 |
|---|---|---|---|---|
| a1 | 정확한 역행렬 | lohner_step.py:513-539 | Fraction 가우스–조르당. 피벗은 절댓값 최대(엄격히 >0)만 선택 → 0 피벗 없음. 열의 남은 성분이 모두 0이면 None(=특이) — 이 경우 실제로 특이. 반환은 [A|I]→[I|A⁻¹] 의 오른쪽 절반. 모든 연산 정확 | 엄밀 |
| a2 | 역행렬을 구간으로 감싸기 | lohner_step.py:686, dyad.py:193-216, 243-257 | DI(q): lo = p비트 내림, hi = p비트 올림(_q_floor_ceil: |x| ∈ [q2^k,(q+1)2^k), r≠0 이면 올림쪽 q+1) → q ∈ [lo,hi] | 엄밀 |
| a3 | 역행렬 대상 = 저장되는 틀 | lohner_step.py:679-686, 695 | B_new 하나의 객체를 역행렬 계산과 state_new.B 양쪽에 씀. prec≥53: Fraction(float) 정확, prec<53: 내린 값으로 B_new를 먼저 확정한 뒤 역행렬 | 엄밀 |
| b | 이름표 갱신 R' | lohner_step.py:688-693 | 실수 항등식 r' = B'⁻¹J*B r + B'⁻¹(φ(c)−c'). 각 실수 인자가 대응 구간 안(B'⁻¹∈[B'⁻¹], J*∈[J], B 정확, r∈R, φ(c)−c'∈TZ−c') → 구간 행렬곱·합의 포함 단조성. 괄호 순서는 실수에서 결합법칙이 성립하므로 무관 | 엄밀 |
| c1 | [Dx_[m]] (X̂ 위) | lohner_step.py:152-306, 580-583 | _GD 전진 자동미분: 합·차·곱·몫·역수·제곱·제곱근 규칙 모두 정확한 미분 공식 + 값 구간이 상자 전체를 감쌈. 상수의 기울기 0(None)은 정확. 재귀는 taylor_step 원본 그대로 호출 | 엄밀 |
| c2 | [J] 조립 | lohner_step.py:653-662 | Σ_{m≤p}[Dx_[m]](X̂)h^m + h^{p+1}([Dx_[p+1]](P)·W). (3a) V의 τ^m 계수 = Dx_[m](y) (혼합편미분 교환), (3b) 성분별 라그랑주 나머지, 나머지 시각 ξ에서 V(s+τ)=Dφ_τ(x(s))V(s) 로 (1/(p+1)!)V^{(p+1)}(ξ)=Dx_[p+1](x(ξ))V(ξ) | 엄밀 |
| c3 | 변분 사전 상자 W | lohner_step.py:442-459, 636-641 | 선형 방정식 V'=A(τ)V, A(τ)=Df(x(τ;y))∈[Df](P)(=D x_[1] on P). 피카르 연산자가 닫힌 집합 C([0,h],Wc)를 자기 자신으로(평균은 볼록 구간 안) + 축약 → 진짜 V∈Wc, 다듬기로 V∈Z | 엄밀 |
| c4 | 평균값 정리의 범위 | lohner_step.py:573-574 | X̂ = hull(외접 상자, c): 볼록, c와 X_k 포함 → 선분 c+s(x−c) ⊂ X̂ 에서 [J] 유효. R이 0을 안 포함해도(c가 집합 밖) 성립 | 엄밀 |
| d1 | 경로 상자 P (Picard) | taylor_step.py:267-323, lohner_step.py:629 | X0+[0,h]F(Bc) ⊆ Bc (끝점 정확 비교, 경계 일치 허용). F(Bc) 계산 성공 = Bc에 충돌 없음 → f가 Bc 근방에서 C¹ → 가중 sup 노름 축약, 닫힌 집합에서 바나흐 고정점. 반환은 다듬은 Z | 엄밀 |
| d2 | h_used 일관성 | lohner_step.py:623-646 | P, W, Hk, info, t 모두 마지막 성공 반복의 h_used. W 실패 → h_used/2 로 P부터 다시 | 엄밀 |
| e | 중심 반올림 | lohner_step.py:665, 692 | c' = 정확한 중점의 p비트 내림. E = TZ − DI(c') (c'는 p비트라 DI 정확)가 φ(c)−c' 를 포함 — 반올림 오차가 E에 그대로 들어감. c'∈TZ 는 엄밀성에 필요 없음(조임만) | 엄밀 |
| f | enclose_at (걸음 안 시각 상자, R4 AUX) | lohner_step.py:706-717 | 0<τ≤h_used 검사. T_τ(c)+τ^{p+1}[x_[p+1]](P) + ([J_τ]·B)·R, [J_τ]에 같은 W·P (구간 [0,h_used] ⊇ [0,τ]에서 유효). 클로저가 잡는 값(prep, gP, DX, RemJ, BD, R, h_used)은 모두 그 걸음 것 | 엄밀 |
| g | 피벗 선택 무관성 | lohner_step.py:466-510, 667-686 | (6) 논증은 B'이 가역이기만 하면 성립. Q는 피벗과 무관하게 Householder 곱(부동소수 근사 직교). 가역성은 정확 역행렬 계산이 판정(특이하면 오류) | 엄밀(주장 맞음) |
| h1 | 정밀도 하나 규칙 | lohner_step.py:318-363, 569-571, 614-616, 696 | 입력 c,B,R 끝점 가수 ≤ prec 확인(아니면 ValueError), prep 정밀도 불일치 ValueError, 출력도 같은 검사. 내부 상수·DI는 호출 시 정밀도 | 엄밀(규칙 지켜짐) |
| h2 | 혼합 정밀도 입력 | 위와 같음 | 높은 정밀도 DI, 비이진 분수, float → 거부. 엄밀성 자체는 정밀도와 무관(모든 연산이 현재 정밀도로 바깥 반올림) | 엄밀 |
| i | 테일러 나머지 Z | lohner_step.py:636, 650 | gP = P 위 p+1차까지. Z_i = h^{p+1}[x_[p+1],i](P), 성분별 ξ_i, x(ξ_i;c)∈P (c∈X̂) | 엄밀 |
| j | 초기화 | lohner_step.py:370-383 | c = 내림, R = DI(x0−c) 바깥 반올림, B=I → x0 ∈ c+R | 엄밀 |
| k | 외접 상자 | lohner_step.py:386-403 | DI(c)+Σ r_j·DI(B_ij), 0인 B 성분만 생략(정확) | 엄밀 |
| l | dyad 산술 | dyad.py 전체 | _dn(>> = floor), _up(ceil), _qdn/_qup(정수 몫 floor / 나머지 있으면 +1), 곱(네 끝점 최소·최대), 나눗셈(부호별 끝점 선택), 짝·홀 거듭제곱, sqrt(isqrt 내림/+1 올림) 방향 모두 맞음 | 엄밀 |
| m | 계수 재귀 | taylor_step.py:122-238 | s^{-3/2} 공식 손유도 일치(α=−3/2: (α+1)i−k = −(i+2k)/2), 코시 곱·대칭 제곱, 중력 부호 | 엄밀 |
| n | r3c_run 격자·AUX, r4_run 저장 | r3c_run.py:347-365, r4_run.py:49-67 | AUX τ = t_k − t_prev ∈ (0,h_used), X1은 새 상태 외접 상자, P는 걸음 경로 상자 — 셋 다 그 시각 참 상태를 감쌈. r4_run 의 k 표기는 id(행)으로, 행 객체가 실행 내내 살아 있어 유일 | 엄밀 |
| o | 견고성: QR에 NaN/inf | lohner_step.py:671-682 | A가 유한이어도 |A|≳1e154이면 제곱합 넘침 → Q에 NaN → Fraction(nan) ValueError(LohnerStepError 아님). 잘못된 상자를 내지 않고 멈춤(r3c_run은 EXCEPTION 처리). R4 범위(폭 ≤1)에서 도달 불가로 추정 | 의심(견고성, 엄밀성 무관) |
| p | 주석 불일치 | lohner_step.py:411 | `_imat_mul` 주석 "0인 점 원소는 건너뛴다" — 코드는 건너뛰지 않음. 값은 같으므로 무해 | 의심(문서) |
| q | 인용 미대조 | taylor_step.py:49-50, 65-66 | Nedialkov–Jackson–Corliss 1999, Jorba–Zou 2005 "기억에 의한 인용". 논증 자체가 자기완결이라 엄밀성 근거로 쓰이지 않음 | 의심(문서) |
