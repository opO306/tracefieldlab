# 피타고라스 삼체 문제 — 기존 연구 조사 v0 (2026-10-11)

## 문제
질량 3, 4, 5인 세 별이 변 길이 3:4:5 직각삼각형의 꼭짓점에 **정지 상태**로 놓인다. 각 별은 자기 질량과 같은 숫자의 변 맞은편 꼭짓점에 있다. 이후 뉴턴 중력만으로 어떻게 움직이는가?
- Burrau(1913)가 제기했다.
- 근접 조우가 반복되고 혼돈이 심해 손 계산으로는 끝까지 풀 수 없었다.

## 조사 방법과 한계
- 학술 검색
  - **scite**: `"Pythagorean" AND "three-body problem"` 검색에 88건이 나왔고, 상위 15건의 제목·초록·인용 분류(지지/반박/언급)를 읽었다. 본문 발췌는 읽지 않았다.
  - **반박 인용**: 0건.
  - **Consensus**: 이달 무료 검색 한도(30회)를 다 써서 **사용하지 못함**.
- 웹 검색 4회: 1967년 해법, Brutus의 신뢰성 비판, Liao의 CNS, 엄밀한(컴퓨터 보조) 증명 여부.
- **원 논문 본문은 아직 읽지 않았다.** 아래 내용은 초록·인용·2차 자료 수준이다.

## 알려진 것

| 연도 | 연구 | 방법 | 주장 | 근거 수준 |
|---|---|---|---|---|
| 1967 | Szebehely·Peters, *Astronomical Journal* 72, 876 | 근접 조우 정칙화(regularization) + 수치 적분 | 결국 두 별이 영구 쌍성이 되고 하나는 무한히 멀어진다. 탈출 시각은 자료마다 t≈60–67(정확한 값 미확인) | 수치 계산(배정밀도 시대) |
| 1969/72 | (PASJ) 쌍곡-타원 운동 판정 기준 개선 | 해석적 판정 기준 | 어느 시점 이후 "영구 탈출"을 보장하는 기준이 이 문제에 효과적이라고 함 | 해석적 기준. 세부 미확인 |
| 1994 | "Global chaoticity in the Pythagorean three-body problem" (*Cel. Mech. Dyn. Astron.*) | 초기조건 주변 탐색 | 이 문제 주변이 전역적으로 혼돈 | 수치 계산 |
| 2014–15 | Boekholt·Portegies Zwart, Brutus 코드 ("On the reliability of N-body simulations") | Bulirsch–Stoer + **임의 정밀도** 산술. 정밀도를 올려도 답이 같으면 "수렴해"로 봄 | 피타고라스 문제의 수렴해를 제시 | 수치 수렴(증명 아님) |
| 2018–20 | 같은 그룹: 시간 역행성 검증, "Planck 길이" 논문 | Brutus | 일부 삼체계는 플랑크 길이보다 작은 정밀도가 필요해 근본적으로 예측 불가 | 수치 계산 |
| 2021 | 상대론적 피타고라스 삼체 문제 (*Phys. Rev. D*) | Brutus + 후뉴턴 보정 | 상대론 효과가 결과에 주는 영향 | 수치 계산 |
| 2013– | Liao의 CNS (고차 테일러 + 수천 자리 정밀도) | 정밀도를 올린 계산과 비교해 수렴 확인 | 주로 주기 궤도. 피타고라스 문제에 쓴 기록은 이번에 확인 못 함 | 대부분 Liao 그룹 자체 연구. 로렌츠 계는 다른 그룹이 재현 |

## 반론·한계 (찾은 것)
- Brutus의 "수렴해"는 **정한 자릿수와 정한 시간 범위 안에서만** 참이다. 수렴은 증명이 아니라 "정밀도를 높여도 안 바뀐다"는 관찰이다.
- 그림자(shadowing) 연구에 따르면 **탈출·포획 영역에서는 그림자 보정이 실패하기 쉽다**(Sitnikov 문제 연구). 피타고라스 문제의 끝(탈출)이 바로 그런 영역이다.
- Liao CNS의 삼체 결과는 외부 독립 검증이 적다.

## 빈자리 후보 (확인 안 됨)
- 구간 산술(CAPD, Lohner 방법) 기반 **엄밀한 컴퓨터 보조 증명**은 제한 삼체 문제와 주기 궤도에는 여러 편 있다.
- 그러나 **피타고라스 문제의 최종 결과(쌍성 + 탈출)를 엄밀하게 증명한 논문은 이번 검색에서 찾지 못했다.**
  - 검색 범위가 좁아 "없다"고 단정할 수 없다(UNKNOWN).
  - 다음 확인: CAPD 그룹(Zgliczyński, Kapela 등)의 논문 목록, "Burrau problem validated numerics".
- 이 빈자리가 실제라면, 의미 있는 목표는 다음과 같다.
  1. 출발부터 시각 T까지 **오차 상자로 궤적을 감싼다**(엄밀한 상한).
  2. 시각 T에서 해석적 "영구 탈출" 판정 기준을 만족함을 보인다.
  3. → 최종 결과의 **컴퓨터 보조 증명**.

## 이 Lab 계산기와의 관계 (추론)
- Lab의 V2 오차 상한 방식은 Lohner 계열이다. **CAPD도 같은 일을 할 수 있다** → 새 방법이라고 주장할 수 없다. 비교 대상이다.
- 어려운 점
  - 근접 조우에서 힘이 거의 무한대가 된다. 그래서 정칙화(좌표 변환)가 필요할 가능성이 크다.
  - 기존 Gate 기록: V1 방식은 약 1,500–1,800 step에서 인증이 멈췄다. 2C에서는 혼돈 궤도 13,906 step까지 갔다. 다만 이 수치는 Hénon–Heiles 결과다.
  - 피타고라스 문제에서 얼마나 갈지는 **UNKNOWN**이다.

## v1 추가 (같은 날): 원 논문 본문 확인

**읽은 것**
- Szebehely·Peters 1967 전체 8쪽. ADS 스캔본의 OCR 글자를 읽었다. 그림과 표의 숫자 일부는 OCR이 깨졌다.
- Boekholt·Portegies Zwart(arXiv 1411.6671)의 1–5쪽(방법, 피타고라스 문제 절).

### 고정한 사실 (1967 본문)

| 항목 | 값 | 출처 |
|---|---|---|
| 질량 | m1=3, m2=4, m3=5 | 1967 p.876 |
| 위치 | 각 별은 자기 질량과 같은 길이의 변 맞은편 꼭짓점. 초속도 0, 질량중심 원점, G=1 | 1967 p.876 |
| 좌표 | (1,3), (−2,−1), (1,−1) | **기억에서 가져온 표준 좌표.** 논문 그림 1은 OCR로 읽히지 않음. 아래 정확 계산으로 논문 서술과 일치 확인 |
| 총에너지 | −769/60 = −12.8166… | 1967 p.876. 위 좌표로 분수 정확 계산해 **일치 확인** |
| 가장 가까운 접근 | t=15.8299236에서 m2–m3 거리 약 4×10⁻⁴, 상대속도 약 171 | 1967 p.878–879 |
| 결말 | **m2·m3가 쌍성**(주기 약 0.9), **m1이 탈출**. t≈59.4–59.7에 m1이 쌍성 사이를 관통. Alexeev 탈출 조건은 t=69에 만족. t=102까지 적분 | 1967 p.879 |
| 방법 | Levi-Civita 정칙화(두 별 거리 < 2일 때 켬) + 5차 Runge–Kutta(Zonneveld) 가변 보폭 | 1967 p.881–882 |
| 정확도 점검 | t=32에서 거꾸로 적분하면 초기 위치가 **소수 10째 자리**까지 복원, t=62에서는 **소수 3째 자리**까지만. 에너지 오차 < 10⁻¹⁰ | 1967 p.883 |

### Brutus 논문 (2014 arXiv)
- t=100까지 적분했다.
- **소수 3자리 수렴**에 Bulirsch–Stoer 허용오차 10⁻¹⁴, 가수 88비트가 필요했다.
- 배정밀도 Hermite 적분기는 마지막 근접 조우(t≈60)까지는 맞다. 그 **이후에는 수렴해와 어긋난다**(탈출 방향이 조금 달라짐). 보폭을 더 줄이면 반올림 오차 때문에 오히려 나빠진다.
- 이 부분에서는 어느 별이 쌍성이 되는지를 따로 확인하지 못했다(1967과 같은지 미확인).

### 엄밀한 증명 선행 연구 (추가 검색 2회, 총 3회)
- 구간 산술·CAPD 기반 증명은 다음에 있다.
  - 주기 궤도(Kapela·Zgliczyński 2003 초-8자, Kapela·Simó 2007 8자 궤도)
  - 제한 삼체 문제의 충돌 궤도(arXiv 2205.03922)
- **피타고라스(Burrau) 문제의 결말을 엄밀히 증명한 논문은 3번 찾아도 나오지 않았다.**
- 여전히 "없다"고 단정하지 않는다(색인 누락 가능).

### 이 결과가 우리 계획에 주는 의미 (추론)
- 1967의 역적분 점검이 t=62에서 소수 3째 자리로 떨어지고, Brutus 결과에서 배정밀도가 t≈60 이후 어긋난다. 그래서 **최종 상태의 정량값(탈출 방향·속도)은 배정밀도로 믿을 수 없다**는 데 두 자료가 일치한다.
- 반면 **질적 결말(쌍성 + 탈출)**은 1967과 Brutus가 같은 방향이다. 다만 Brutus에서 쌍성을 이루는 별이 1967과 같은지는 미확인이다.
- 엄밀 증명의 난관
  - 최소 거리 4×10⁻⁴, 상대속도 171인 근접 조우를 지나야 한다. → 정칙화 없이는 오차 상자가 터질 가능성이 크다.
  - 1967도 정칙화를 썼다.

## 다음 단계 제안
1. 원 논문(1967)과 Brutus 논문 본문 읽기. 초기조건·단위·탈출 시각을 정확히 고정한다.
2. 엄밀한 증명 선행 연구를 더 찾는다(위 "다음 확인").
3. 작은 타당성 시험(Gate 계획 봉인 후)
   - 정칙화 없이 Lab 방식으로 **어디까지 인증되는지** 잰다(첫 근접 조우에서 멈출 것으로 예상).
   - 멈추면 정칙화가 필요하다는 근거가 된다.

## 출처
- Szebehely & Peters 1967 (ADS): https://articles.adsabs.harvard.edu/pdf/1967AJ.....72..876S
- Physics Today 1968, "Yale Astronomers Solve Old Three-Body Problem": https://doi.org/10.1063/1.3034832
- Global chaoticity in the Pythagorean three-body problem (1994): https://doi.org/10.1007/bf00692114
- Improved Criteria for Hyperbolic-Elliptic Motion (PASJ): https://doi.org/10.1093/pasj/24.3.391
- On the Reliability of N-body Simulations: https://arxiv.org/pdf/1411.6671
- On the minimal accuracy required…: https://ar5iv.arxiv.org/html/1402.6713
- Numerical verification of microscopic time reversibility: https://arxiv.org/pdf/1802.00970
- Gargantuan chaotic gravitational three-body systems…: https://arxiv.org/pdf/2002.04029
- Relativistic Pythagorean three-body problem (PRD 2021): https://doi.org/10.1103/physrevd.104.083020
- The three-body problem (Rep. Prog. Phys. 2014 review): https://doi.org/10.1088/0034-4885/77/6/065901
- Liao, Chaotic motion of three-body problem: https://arxiv.org/pdf/1304.2089
- Li & Liao, Clean Numerical Simulation: https://arxiv.org/pdf/1609.09344
- Shadowing unstable orbits of the Sitnikov problem: https://arxiv.org/pdf/1005.1225
- CAP drift orbits, restricted 3BP (CAPD): https://arxiv.org/pdf/2105.07204
- CAP collision orbits, restricted 3BP: https://arxiv.org/pdf/2205.03922
- Zgliczyński, Computer assisted proofs in dynamics (lecture): https://ww2.ii.uj.edu.pl/~zgliczyn/cap07/lect1.pdf
