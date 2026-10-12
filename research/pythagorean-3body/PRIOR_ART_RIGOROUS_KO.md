# 피타고라스 삼체 문제 — 엄밀 증명 선행 연구 + 영구 탈출 판정 기준 (2026-10-11)

작성 범위: 읽기 전용 조사. 코드·계산기·다른 문서는 건드리지 않았다. 내려받은 PDF는 세션 스크래치 폴더에만 두었고 저장소에는 복사하지 않았다.

---

## 0. 쉬운 설명

**질문 1 — "이미 누가 컴퓨터로 엄밀하게 증명했나?"**
- 이번 조사에서는 **찾지 못했다.** 다만 "없다"가 증명된 것은 아니다(검색 범위가 한정됨. 6절 참조).
- 피타고라스 문제의 결말(쌍성 + 탈출)은 지금까지 모두 "컴퓨터로 계산해 보니 그렇더라"(수치 계산)였다. 정밀도를 높여도 안 바뀐다는 확인(Brutus)은 있어도, 오차 상자로 증명한 글은 못 찾았다.
- 부품은 이미 있다. 같은 도구(구간 산술 + 정칙화)로 **제한 삼체 문제**의 궤도를 증명한 논문들이 있다. 그래서 "새로운 것은 방법이 아니라 결과"라는 로드맵 판단이 맞다. 결과가 새로운지는 계속 "찾지 못함" 상태다.
- 1974년 논문 하나가 우리 계획과 **방법이 닮았다.** 수치 궤적 + 해석적 판정 기준으로 탈출을 보였다. 저자가 직접 "수치 오차의 엄밀한 평가는 없다"고 적었다. 이것은 우리가 메우려는 바로 그 빈틈이다.

**질문 2 — "'영구 탈출'을 수학적으로 어떻게 판정하나?"**
- 비유: 쌍성(4·5)을 하나의 무거운 점으로 보고, 질량 3이 그 점에서 멀어지는 거리 ρ를 본다. 충분히 멀고, 충분히 빠르게 멀어지고 있고(ρ̇가 문턱값보다 큼), 쌍성 안쪽 거리 r이 작으면, 이후에는 **중력이 되돌릴 수 없다**는 것을 부등식으로 증명한 정리들이 있다. 구체적으로는 "어느 순간 하나의 부등식 묶음이 참이면 그 뒤로 영원히 참이다"라는 꼴이다.
- 실제 기술적 의미: 이 정리들은 **한 시각의 상태(위치·속도 12개 숫자)** 만 보고 결론을 낸다. 그래서 오차 상자로 그 시각 T까지 가둔 뒤, 상자 안의 모든 점이 부등식을 만족함을 확인하면 끝이다.
- 가장 쓸 만한 기준은 Yoshida(1972)의 정리다. 원문에서 그가 이 문제에 직접 적용한 예에서는 **t=59.72** 에서 만족했다. 1967년 논문이 쓴 Alexeev 조건(t=69)보다 약 9 시간 단위(무차원) 일찍 만족한다. 단, **여유가 아주 작다**(0.37%). 4절 "우리 계산에서 확인해야 할 양" 참조.

---

## 1. 도구와 검색 기록 (전부)

### 1-1. 학술 검색 (scite `search_literature`) — 8개 질의
각 질의에서 **상위 15~25건의 제목·연도·저자·인용 분류만** 읽었다. 본문은 읽지 않은 것이 대부분이다. 총 건수는 scite가 알려 준 값이다.

| # | 질의(요지) | 총 건수 | 결과 |
|---|---|---|---|
| 1 | "Burrau" AND three-body AND (rigorous / computer-assisted / interval arithmetic / validated) | 26 | 관련 엄밀 증명 없음. 제한 삼체 문제 학위논문, Szebehely 1972 등 |
| 2 | "three-body problem" AND computer-assisted proof AND (escape / capture / hyperbolic) | 85 | 전부 제한 삼체 문제·Arnold 확산·이종궤도(Wilczak–Zgliczyński, Capiński 등). 피타고라스 아님 |
| 3 | "hyperbolic-elliptic" AND "three-body" AND (criterion / criteria) | 19 | Yoshida 1972·1974, Georgakarakos 개관 등 |
| 4 | Taylor model/interval/validated numerics AND N-body AND (rigorous/verified) AND (escape/close encounter/regularization) | 5530 | 잡음이 많아 의미 있는 항목 없음 |
| 5 | (interval/Taylor model/Lohner/validated/rigorous numerics) AND (three-body/N-body) AND (collision/close encounter/regularization/triple collision) | 131 | 제한 삼체 문제의 충돌·근충돌 궤도(Capiński–Kepley), 안무(Kapela–Zgliczyński) 등 |
| 6 | "Pythagorean" AND three-body AND (Szebehely OR Burrau) | 32 | 수치 연구뿐(Brutus, 상대론, 기호적분기, 시간 재규격화 등) |
| 7 | computer-assisted AND (free fall / equal mass / general three-body) AND (proof / rigorous) | 203 | Moeckel 2004(Saari 추측), 중심배열 등. 피타고라스 아님 |
| 8 | Marchal / escape / sufficient conditions / Yoshida / "test of escape" | 매우 큼 | Marchal 1974, Marchal–Yoshida 1984, Standish 1971 등 서지 확인 |

- 추가로 DOI로 서지·초록을 확인: 22편(Yoshida 1974, Szebehely 1972, Marchal 1974 등).
- **인용 추적**(`citation_graph`, 방향=인용한 쪽, 시드 4개: Szebehely–Peters 1967(AJ), Szebehely 1967(PNAS), Yoshida 1972, Aarseth 등 1994 "Global chaoticity"): 190개 인용 간선, 181편. **제목만 훑었다.** 제목에 "rigorous/interval/verified/validated"가 들어간 항목은 **0건**이었다. 의심스러운 제목(Taylor 계열, 정칙화 계열)의 초록은 별도로 확인했다.
- 철회·정정 공지: 조회된 서지 기록에 철회 공지는 나타나지 않았다.

### 1-2. 웹 검색 — 14회
1. Pythagorean + rigorous + computer-assisted + interval arithmetic
2. Burrau + validated numerics + Taylor model + free-fall escape
3. (확장 검색) Pythagorean/Burrau + rigorous/interval/validated/verified + arXiv
4. Pythagorean 한계·비판·수치 아티팩트
5. Yoshida 판정 기준·Alexeev·Standish·Marchal 진술
6. Zgliczyński·Kapela·Wilczak CAPD + 삼체 + 정칙화 + Pythagorean/Burrau/free fall
7. Berz·Makino Taylor model + 천체역학 + 근접 조우 + Levi-Civita
8. (확장 검색) 특정 초기조건 탈출의 computer-assisted proof, 쌍곡-타원 최종 운동
9. Pythagorean "never been proven / no rigorous proof / mathematically rigorous"
10. 구간 산술 + Kepler + Levi-Civita/KS 정칙화 + verified ODE solver
11. Nedialkov VNODE-LP + N-body/three-body
12. Pythagorean 2024–2026 + rigorous error bounds + Taylor + verified
13. Capiński–Gidea "Arnold diffusion in the full three-body problem"
14. Marchal 1990 + test of escape + Yoshida 기준

### 1-3. 원문을 직접 읽은 것
- **Yoshida 1972** (PASJ 24, 391–408) 전체 18쪽: ADS의 스캔 PDF를 Windows 내장 PDF 렌더러로 이미지로 바꿔 쪽마다 눈으로 읽었다. (텍스트 층이 없는 스캔이라 OCR 없이 직접 읽음. 식을 옮기다 틀렸을 가능성이 있어, 필요한 식은 4절의 확인 계획과 6절 미확인 항목에 반영.)
- **Yoshida 1974** (PASJ 26, 367–377) 1쪽(초록·서론).
- 후속 문헌에서 인용문만 확인: Georgakarakos 2008 (arXiv 1408.5431), Li–Fu–Sun 2009(A&A, 웹 페이지 요약), Antoñana 등 2020 (arXiv 2001.01221), Janssens (arXiv 1103.3662), Boekholt–Portegies Zwart 등.

### 1-4. 사용하지 못한 도구
- **Consensus**: 이달 한도 소진(지시사항).
- zbMATH, MathSciNet, Google Scholar, ADS 전문 검색은 **사용하지 않았다**(접근 도구 없음 또는 미시도).
- 유료 원문(Alekseev 1961, Standish 1971, Marchal 1974, Marchal–Yoshida 1984, Marchal 1990 단행본)은 **읽지 못했다.**

---

## 2. Q1 결론

> **엄밀한(구간 산술·검증 수치·컴퓨터 보조) 피타고라스 문제 결말 증명은 이번 조사에서 찾지 못했다. (NOT FOUND — "없음"의 증명은 아님.) 부분 선행은 있다: ① 같은 도구로 제한 삼체 문제를 증명한 논문들(방법의 선행), ② 수치 궤적 + 해석 판정 기준으로 탈출을 보인 비엄밀 연구들(구조의 선행).**

### 2-1. "피타고라스 문제"를 직접 다룬 것 — 전부 수치 계산
| 연구 | 한 일 | 엄밀성 |
|---|---|---|
| Szebehely & Peters 1967 (AJ 72, 876; DOI 10.1086/110355), Szebehely 1967 (PNAS 58, 60) | Levi-Civita 정칙화 + 5차 RK로 t=102까지. 결말: 질량 4·5 쌍성 + 3 탈출. Alexeev 조건을 t=69에 확인 | 수치. 저자는 "임의의 시간에 유효한 완전 해"라고 주장하나 수치 계산에 근거 |
| Yoshida 1972 (PASJ 24, 391) | 해석 판정 기준을 개선하고 이 문제의 수치 궤적(Yuasa 1971 개인 통신 자료)에 적용 | 기준은 해석적. 입력 궤적은 비엄밀 |
| Szebehely 1972 (PNAS 69, 1077) | 삼체 계 다수를 수치 실험 | 수치 |
| Aarseth·Anosova·Orlov 1994 (CMDA 58, 1; Szebehely를 저자로 적은 자료도 있음) | 초기조건 주변이 전역적 혼돈 | 수치 |
| Boekholt & Portegies Zwart 2015 (Brutus, 임의 정밀도 + Bulirsch–Stoer) | 정밀도를 올려 수렴한 해를 얻음. 배정밀도 코드는 t≈60 이후 어긋남 | 수렴 관찰(증명 아님). 저자도 "참인 해는 지정한 자릿수까지"라고 한정 |
| Rose & Dullin 2013 (CMDA 117, 169) | 정칙화·대칭 축소한 기호적분기. 피타고라스 궤도를 시험 | 수치 |
| Antoñana·Chartier·Makazaga·Murua 2020 (SIAM, arXiv 2001.01221) | N체 전역 시간 재규격화: 해가 가상 시간에서 복소수 띠로 정칙 확장됨을 **증명**(일반 N체). 피타고라스 t∈[0,63] 수치 시험 | 이론 부분은 증명이나 피타고라스 결과는 수치. (초록·일부 발췌만 읽음) |
| Hirayama 2019 (Springer, 4배 정밀도 고차 Taylor) | 수치 | 수치 (유료, 초록 수준) |
| 상대론적 피타고라스(2021, PRD / arXiv 2109.07013), Chitan 등 2021, Parischewsky 등 2023, Portegies Zwart 2026(arXiv 2608.17032, 부동소수점 비교) | 수치 연구 | 수치 |
| Moeckel–Montgomery–Venturelli 2012, Montgomery 2022–23 | 영속도 출발("brake") 삼체 궤도의 위상수학·변분 연구 | 피타고라스 결말 증명 아님 |

### 2-2. 같은 도구로 증명된 근처 결과 (방법의 선행)
- 주기 궤도: Kapela & Zgliczyński 2003 (안무 궤도, CAPD), Kapela & Simó 등.
- 제한 삼체 문제: 이종궤도(Wilczak–Zgliczyński 2003), 충돌·근충돌 궤도(Capiński–Kepley–Mireles James, arXiv 2205.03922), 표류 궤도·Arnold 확산·Lyapunov 궤도(Capiński, Gidea 등), 무한대로 가는 진동 운동(Capiński–Guàrdia–Martín–Seara–Zgliczyński 2022).
- **Levi-Civita 정칙화를 구간 산술 검증에 쓴 사례**: Capiński–Pasiut 2026 (arXiv 2608.05400, 지구–달 제한 삼체 문제). *웹 검색 요약으로만 확인했고 본문은 읽지 않았다.*
- 완전한(일반) 삼체 문제의 컴퓨터 보조 증명: Moeckel 2004(Saari 추측, 대수·기하 계산), Capiński–Gidea 2025(arXiv 2504.09273, 비교적 작은 질량 한 개를 가진 평면 삼체 문제의 Arnold 확산).
- 검증 적분기 일반: Berz–Makino(Taylor model)의 소행성 1997 XF11·Apophis 지구 근접 통과 검증 적분(ESA 연구 보고, 도메인 분할 사용), Nedialkov VNODE-LP(Lohner 방법), CAPD.
- **위 어느 것도 "일반 삼체 문제에서 특정 초기조건의 쌍성 형성 + 탈출을 컴퓨터 보조로 증명"한 것은 아니다**(제목·초록·요약 수준에서 확인). 이 일반 질문(Q1 후반)에 대한 답: **그런 사례를 찾지 못했다.**

### 2-3. 구조가 닮은 비엄밀 선행 (조심할 부분)
- **Yoshida 1974**(PASJ 26, 367) 초록·서론(쪽 367): 수치 계산한 궤적(Becker 1920 예)에 해석 판정 기준을 적용해 "t→+∞에서 쌍곡-타원"이라고 결론. 원문은 "수치 오차의 엄밀한 평가를 무시하고"(In disregard of the rigorous estimations of numerical errors)라고 적고 서론에서 "수치 오차 평가가 완전하지 않다"고 인정한다. → **"수치 궤적 + 해석 기준"의 조합은 이미 1972–74년에 있었다. 새로 더할 것은 궤적을 오차 상자로 엄밀히 가두는 부분**이다.
- Yoshida 1972도 마찬가지로 Yuasa(1971) 개인 통신 수치 자료를 그대로 입력으로 썼다(쪽 406).

---

## 3. Q2 — 영구 탈출 판정 기준 (원문 기준)

### 3-0. 용어 (Yoshida 1972, 쪽 392)
- m₀ = 탈출하는 질량, (m₁, m₂) = 쌍성. 피타고라스 문제에서는 **m₀=3, m₁=4, m₂=5** (쪽 406). G=1.
- ρ = 쌍성 질량중심에서 m₀까지의 거리, ρ̇ = 그 시간 미분(m₀의 쌍성 질량중심에 대한 시선 방향 속도). r = 쌍성 내부 거리 |m₁m₂|.
- M = m₀+m₁+m₂ (=12), μ = m₀(m₁+m₂)/M, Σ = m₀m₁ + m₁m₂ + m₂m₀ (=47), λ = m₂/(m₁+m₂) (=5/9), ν = m₁/(m₁+m₂) (=4/9).
- H = 총에너지, K = −H > 0. 피타고라스: H = −769/60 = −12.81667 (쪽 406 수치와 일치).
- **"쌍곡-타원 운동"** = ρ→∞ (ρ̇의 하한 극한이 양수)이고 r이 유계. 한 시각 t₀의 조건에서 이후 모든 시각을 결론으로 얻는다(쪽 391, 400–402).

### 3-1. 주 기준: Yoshida 1972, 정리 1의 따름정리 (쪽 400–402, 식 (13)–(15), (24)–(27))
**가정**: K = −H > 0, m₁ ≤ m₂(아니면 역할을 바꾸어 아래 변형 사용). 한 시각(t₀)에 양수 d와 a(>λ)가 존재하여

- (24) d ≥ m₁m₂/K
- (25) ρ ≥ a·d
- (26) ρ/r ≥ a
- (27) ρ̇ > √( 2M · { λ/(ρ + νd) + ν/(ρ − λd) } )

**결론** (쪽 402): 운동은 t→+∞에서 쌍곡-타원 유형이며, 모든 t ≥ t₀에서 r < m₁m₂/K.
- m₁ ≥ m₂인 경우: a>ν로 바꾸고 (27)의 분모를 (ρ−νd), (ρ+λd)로 바꾼다.
- 쪽 403 Note 1: 이 정리에서 r은 **가장 작은 상호 거리일 필요가 없다**(쌍성 내부 거리면 됨).
- 쪽 403 Note 2: 이 따름정리는 Merman 1954(M3), Tevzadze 1962(T), Standish 1971(S)의 기준을 **포함**한다. 즉 더 약한 조건으로 같은 결론을 낸다.
- **내가 유도한 간단화(추론, 원문 명시 아님)**: (27)의 문턱값은 a에 의존하지 않고 d에 대해 증가하므로 최적은 d = m₁m₂/K. 이때 a가 존재할 조건은 **ρ > λd 이고 ρ > λr** 둘뿐이다. 이렇게 되면 우리가 확인할 부등식은 (i) ρ > λ·d₁, (ii) ρ > λ·r, (iii) ρ̇ > 문턱값, 세 개다. (d₁ = m₁m₂/K = 1200/769.)

**증명의 구조** (쪽 401–402 직접 읽음): 연속성으로 구간 [0,t₁]에서 ρ/r ≥ a, r ≤ d, ρ̇ ≥ 0을 확보 → Lemma 2(ρ̈의 하한)와 에너지 적분으로 ρ̇ > 0, r < d를 유지 → 최대 시각 τ의 모순 논법으로 모든 t에 확장. 쌍성 사이의 이중 충돌은 해석적 연장이 가능하다고 보고 허용한다("ρ/r>max(λ,ν)인 동안은 m₁–m₂ 충돌만 가능", 쪽 401). 삼중 충돌은 가정으로 배제된다.

**엄밀성**: 해석적 정리와 증명이 있다(수치 가정 없음). 다만 (a) 내가 증명을 한 줄씩 검산하지는 않았다. (b) 기계 검증(형식 증명)이나 독립 재검증 문헌은 찾지 못했다. → **우리 쪽에서 별도 세션이 증명을 독립 재도출하는 것을 권한다.**

### 3-2. 피타고라스 문제에 적용된 예 (쪽 406–407)
**① 따름정리, t=59.72** (Yuasa 1971 수치 자료, 소수 5자리):

| 양 | 값 |
|---|---|
| H=−K | −12.81667 |
| d₁ = m₁m₂/K | 1.56047 |
| ρ | 1.69338 |
| r | 0.42196 |
| ρ̇ | **4.31611** |
| (27)의 문턱값 √(2M{…}) | **4.30030** |
| a | 1.0 (> λ=0.55556) |

→ 조건 만족. 결론 r < 1.56047. (수치 해는 t≥59.4에서 r<1.1.)

**내가 한 재계산(독립이라고 부를 수 없음)**: 원문이 인쇄한 12개 위치·속도 숫자(질량 3,4,5 순서)로 ρ, ρ̇, r, 문턱값을 다시 계산하니 ρ=1.69337, ρ̇=4.31611, r=0.42195, 문턱값=4.30031로 **원문과 일치**했다. 같은 숫자에서 질량중심·총운동량≈0, H=−12.8172(정확값 −12.81667과 5×10⁻⁴ 차이, 5자리 반올림 범위). 이것은 "원문의 산술과 행 순서를 확인"한 것이지 **궤적 자체를 확인한 것은 아니다.**

> **여유가 작다.** ρ̇ − 문턱값 = 0.0158 (ρ̇의 약 0.37%). 이 시각의 오차 상자는 ρ̇를 최소 약 0.3% 이내로 가두어야 한다. 더 늦은 시각의 여유는 **우리가 직접 재야** 알 수 있다(UNKNOWN).

**② 정리 2, t=59.99**: a=1.0, d=1.65, ρ=2.73984, ρ̇=3.65062, 문턱 3.11280, 좌변 8.14230 > 우변 8.04156이라고 원문에 적혀 있다.
- **내 재계산은 이 두 값을 재현하지 못했다**(좌변 8.2157, 우변 7.2146; 부등식 방향은 같음). 원인(내 식 해석 오류, 원문 오기, 다른 반올림) **미해결.** 정리 2는 쓸 예정이 없으므로 따름정리와 별개로 둔다.

### 3-3. 오래된 기준 (Yoshida가 쪽 398–400에 재진술; 원문 PDF는 읽지 못함)
아래 식은 **Yoshida의 재진술 이미지를 옮겨 적은 것**이라 전사 오류 가능성이 있다. 쓰기 전에 원문(Standish 1971 등)과 대조해야 한다.

| 기준 | 가정 | 조건 (t₀ 시각) | 이 문제에서 만족한 시각 |
|---|---|---|---|
| Khil'mi 1951 (정리 K) | H>0 | — | 해당 없음(H<0) |
| Merman 1952·1953 (M1, M2) | 에너지 부호 무관 | 복잡(상수 A′, B′, …) | — |
| **Alekseev 1961 (정리 A, 쪽 399)** | K>0, 쌍성 에너지 음. a>max(λ,ν) | ρ̇>0, 4MA < ρ·\|ρ̇\|·ρ̇, 쌍성의 (m₁+m₂)/r − ṙ²/2 > [ (3/(2√2))·m₀(m₁+m₂)/(ρ²\|ρ̇\|)·B₀ + (a(m₁+m₂)/ρ)^{3/2} ]^{2/3}. 결론: ρ(t)의 하한, **r < ρ(0)/a** | Szebehely–Peters가 **t=69, ρ=26.2**에서 확인(쪽 406–407). Yoshida는 a=6으로 t=61.31, ρ=6.89에서 확인하여 r<1.15 얻음 |
| Merman 1954 (M3, 쪽 400) | K>0 | ρ>2Σ/K, ρ̇>√(8M/ρ) | **t=64.5, ρ=15.4** (Σ/K=3.667이므로 2Σ/K=7.33) |
| Tevzadze 1962 (T, 쪽 400) | K>0 | ρ>Σ/K, ρ̇ ≥ √(2M{λ/(ρ−νΣ/K)+ν/(ρ−λΣ/K)}) | **t=60.45, ρ=4.30** |
| **Standish 1971 (S, 쪽 400)** | K>0 | ρ>Σ/K, ρ̇ ≥ √(2M{1/ρ + λν(Σ/K)²/(ρ²(ρ−Σ/K))}) | **t=60.52, ρ=4.52** |

- 시간순으로 가장 빠른 것은 Yoshida 따름정리(59.72), 그 다음 T(60.45), S(60.52), Alekseev a=6(61.31), M3(64.5), Alekseev 원 조건(69).
- Alekseev 1961의 원문은 읽지 못했고(러시아어, Vestn. Mosk. Univ. Ser. Fiz. Astron. 1, 67), 조건은 Yoshida의 재진술(Lemma 3, 쪽 397의 a₁, a₂, B₀ 정의 포함)에 의존한다.
- Standish·Tevzadze·Merman의 정리는 **에너지 K>0만 필요하고 각운동량은 안 쓴다**(쪽 391). 각운동량이 0인 이 문제에 맞는다.

### 3-4. 쓸 수 없거나 다른 성격인 것
- **Yoshida 1974** (PASJ 26, 367; 초록·서론만 읽음): 각운동량이 **0이 아닐 것**을 요구해 완화된 조건을 준다. 피타고라스 문제는 정지 출발이라 총 각운동량이 0(원문이 아닌 내 추론)이므로 **적용 불가** 가능성이 높다. 본문 미확인.
- **Marchal 1974** (CMDA 9, 381), **Marchal–Yoshida–Sun 1984** (Cel. Mech. 33, 193; 34, 65; Acta Astronaut. 11, 415), Marchal 1990(단행본): "아주 작은 상호 거리에서도 유효한 탈출 시험". 서지·소개만 확인(Georgakarakos 2008, Li–Fu–Sun 2009). 원문 미확인. Li–Fu–Sun(A&A 504, 277)의 설명(웹 요약, 일부 식 잘림)에 따르면 Marchal의 시험은 1차원 투영 λ(거리 하한)에 대한 부등식 하나(dλ/dt ≥ √(2GMĀ/λ)이면 영원히 성립)이고 일반 N체에 적용되며 **"탈출(분리)"만 결론**한다. 쌍성의 유계성(r<d) 결론은 따름정리가 더 직접적이다(추론).
- **Sosnitskii 2012** (arXiv 1210.0321, "Lagrange stability … final evolutions"): 쌍성 Hill 안정 쌍 + 거리 운동을 이용한 결과. 쪽 단위로 읽지 않았으며 이 문제에 적용 가능한지 **UNKNOWN**.
- **Georgakarakos 2008** (CMDA; arXiv 1408.5431)의 개관(쪽 6–8): Standish 1971, Yoshida 1972, Griffith–North 1973, Marchal 1974의 기준을 "쌍곡-타원 운동의 **충분조건**"으로 묘사하고 Marchal 등(1984)의 기준을 "더 강한 기준"이라 소개. 2차 자료.

---

## 4. 우리 계산에서 확인해야 할 양 (Yoshida 따름정리 기준)

각 시각 T에서 **상태 상자(위치 6 + 속도 6, 질량 3·4·5)** 에 대해 구간 산술(바깥 반올림)로 계산한다.

1. **역할 식별**: m₀=3, 쌍성=(4,5). 어떤 별이 탈출하는지는 상자 자체가 보여야 하며(예: ρ가 이 라벨에서 크고 r이 작음), 라벨 선택은 수동이 아니라 검사로 증명한다.
2. 쌍성 질량중심 G = (4·x₁ + 5·x₂)/9, 속도 같은 식. ρ = |x₀ − G|, ρ̇ = (x₀−G)·(v₀−v_G)/ρ, r = |x₂−x₁|.
3. **정확 상수**: K = 769/60, d₁ = m₁m₂/K = 1200/769, λ=5/9, ν=4/9, M=12 (모두 유리수. 진 해의 에너지 보존을 이용한 정확값이므로 상자 폭과 무관).
4. 부등식 3개를 **상자 전체에서** 증명:
   - ρ > λ·d₁ 및 ρ > λ·r
   - ρ̇ > √(2M{ λ/(ρ+νd₁) + ν/(ρ−λd₁) })  (문턱값은 ρ에 대해 감소 → 최악 값은 ρ 하한에서 계산. 의존성 때문에 단순 구간 평가는 보수적임)
   - 즉 **하한(ρ̇) > 상한(문턱값)** 비교만 하면 증명.
5. 시각 T 선택은 보수적으로 여러 후보(59.72, 60.5, 61.3, 64.5, 69)에서 각 기준의 여유를 재서 정한다. 이 여유는 **지금 모르는 값(UNKNOWN)**.
6. 이중 확인 권장: Yoshida 따름정리와 별개의 저자 기준(Standish S, Tevzadze T, Alekseev A)을 **다른 시각**에 각각 증명해 서로를 교차 점검한다. 한쪽 정리의 증명 오류나 내 전사 오류에 대한 보험이다(추론).
7. 상태 상자를 진 에너지 곡면에 교집합하면 상자를 줄일 수 있는지는 별도 검토(추론, 미시험).

---

## 5. 반론·한계

- 위 기준은 **충분조건**이다. 만족하지 않아도 탈출이 아닌 것은 아니다.
- 원문 수치(t=59.72)는 Yuasa의 **개인 통신 자료**(정밀도 미상)이므로 인증 자료가 아니다. 이 값이 우리 상자 안에 들어가는지는 우리 계산이 판정한다.
- 따름정리의 **여유 0.37%**는 얇다. 1967년 논문의 역적분 점검이 t=62에서 소수 3째 자리까지만 복원했고, Brutus는 배정밀도가 t≈60 이후 어긋난다고 보고(이미 SURVEY에 기록). 마지막 근접 통과(t≈59.4–59.7) 직후라서 상자 폭에 가장 불리한 시각 근처다. 이 판정 시각이 우리 R2–R5 단계에서 실제로 도달 가능한지는 UNKNOWN.
- "hyperbolic-elliptic"은 ρ→∞이며 r 유계라는 **질적** 결론이다. 쌍성이 점근적으로 케플러 타원이 되는 정량 결론(탈출 속도·방향의 수치)은 이 기준이 주지 않는다. 우리 목표("쌍성 + 영구 탈출")에는 충분하나 "탈출 방향" 같은 정량은 아니다.
- 피타고라스 문제에 대한 수치 결말의 신뢰는 여러 독립 코드(1967, Brutus, 후속)가 질적으로 같다는 데 있다. 그러나 1967 논문의 "임의의 시간에 유효"라는 표현은 **수치에 근거한 주장**이다. 이 주장에 대한 직접적인 반박 논문은 찾지 못했다.
- 일반 삼체 문제에서 쌍곡-타원 운동은 영측도 부근의 예외적 집합이라는 취지의 서술을 웹 검색 요약에서 보았으나 원문 확인은 못 했다. 그래서 "특정 초기조건을 직접 증명해야 한다"는 우리 접근과는 일치하나 이 서술 자체는 **UNKNOWN**.
- 새로움 주장의 위험: 이 방법 부품(구간 산술, Taylor/Lohner, Levi-Civita)은 이미 CAPD·VNODE-LP·Berz–Makino 계열에 있다. 새로움은 **결과에만** 걸어야 한다. 그리고 그 결과가 정말 새로운지는 아래 6절 때문에 아직 확정 못 한다.

---

## 6. 미확인 (중요도 순)

1. **[가장 중요] "엄밀 증명이 없다"는 것은 확인된 사실이 아니다.** 사용하지 않은 곳: zbMATH, MathSciNet, Google Scholar, ADS 전문 검색, 박사 논문 저장소, 러시아·일본어 문헌, CAPD/Taylor model 그룹의 미출판 자료·발표 자료. 결과가 새롭다고 주장하기 전에 최소한 이 곳들과 CAPD 그룹에 직접 문의가 필요하다.
2. 원문을 읽지 못한 기준: Alekseev 1961, Standish 1971, Marchal 1974, Marchal–Yoshida 1984, Marchal 1990. 3-3절 식은 Yoshida의 재진술을 눈으로 옮긴 것이라 오류가 있을 수 있다.
3. Yoshida 따름정리 증명의 정확성(내가 줄 단위로 검산하지 않았다).
4. 3-2절 ②(정리 2, t=59.99)의 인쇄 수치와 내 재계산이 맞지 않는 원인.
5. 우리 궤적에서 각 기준의 **실제 여유와 가능한 최소 T**(계산 전에는 알 수 없음).
6. Capiński–Pasiut 2026, Antoñana 등 2020, Hirayama 2019, Sosnitskii 2012, Janssens의 본문 내용(초록·요약 수준만 확인).
7. Yoshida 1974 본문(각운동량 조건의 정확한 형태. 이 문제가 정말 적용 불가인지).
8. 이 보고서의 학술 검색은 상위 15~25건만 읽은 것으로, 목록 밖에 관련 논문이 있을 가능성이 있다.

---

## 7. 출처

원문을 읽은 것(전문 또는 일부):
- Yoshida, J. (1972) Improved Criteria for Hyperbolic-Elliptic Motion in the General Three-Body Problem. PASJ 24, 391–408. https://doi.org/10.1093/pasj/24.3.391 (ADS: 1972PASJ...24..391Y)
- Yoshida, J. (1974) … II. PASJ 26, 367–377. https://doi.org/10.1093/pasj/26.3.367 (1쪽만)
- Szebehely & Peters (1967) AJ 72, 876 (이전 세션에서 읽음. DOI 10.1086/110355)

초록·요약 수준으로 확인:
- Szebehely (1967) PNAS 58, 60, https://doi.org/10.1073/pnas.58.1.60 ; (1972) PNAS 69, 1077, https://doi.org/10.1073/pnas.69.5.1077
- Georgakarakos, Stability criteria for hierarchical triple systems, https://arxiv.org/abs/1408.5431
- Li, Fu & Sun (2009) A&A 504, 277, https://doi.org/10.1051/0004-6361/200809392
- Marchal (1974) CMDA 9, 381, https://doi.org/10.1007/bf01228579 ; Marchal–Yoshida 1984, https://doi.org/10.1007/bf01230502 ; Standish (1971), https://doi.org/10.1007/bf01230319
- Boekholt & Portegies Zwart (2015), https://doi.org/10.1186/s40668-014-0005-3
- Rose & Dullin (2013), https://doi.org/10.1007/s10569-013-9503-x
- Antoñana, Chartier, Makazaga, Murua (2020), https://arxiv.org/abs/2001.01221
- Janssens, On Free Fall in the Three Body Problem, https://arxiv.org/abs/1103.3662
- Sosnitskii (2012), https://arxiv.org/abs/1210.0321
- Moeckel (2004) Trans. AMS 357, 3105, https://doi.org/10.1090/s0002-9947-04-03527-5
- Kapela & Zgliczyński (2003), https://doi.org/10.1088/0951-7715/16/6/302
- Capiński, Kepley, Mireles James, https://arxiv.org/abs/2205.03922 ; Capiński & Pasiut, https://arxiv.org/abs/2608.05400 ; Capiński & Gidea, https://arxiv.org/abs/2504.09273
- Portegies Zwart (2026), https://arxiv.org/abs/2608.17032
- ESA/Berz–Makino Apophis 2029 Taylor model 보고서 (웹 요약만): https://www.esa.int/gsp/ACT/doc/ARI/ARI%20Study%20Report/ACT-RPT-MAD-08-4303-Encounter2029_TaylorModels-Milano.pdf
- Nedialkov, VNODE-LP: https://www.cas.mcmaster.ca/~nedialk/vnodelp
