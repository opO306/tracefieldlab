# Gate FC1 독립 검사기 보고 (`check_fc1.py`) — 2026-10-12

> **상태: 검사기와 자기시험까지 완료. 실제 `fastcore/`(C++ 구현)에는 실행하지 않았다(NOT_TESTED).**
> 이 문서의 수치는 **내가 만든 참조 대상 3개와 가짜 대상 44개(+ 동치 변형 3개)**에 대한 것이지, 실제 구현의 성적이 아니다.
> 최종 실행(검사기 × 실제 구현)과 판정은 조율 세션이 한다.

## 0. 한눈에 보기

| 항목 | 결과 | 근거 |
|---|---|---|
| 봉인 해시 | 계획서·`dyad.py`·`taylor_step.py`·`lohner_step.py` 모두 `GATE_FC1_SEAL.json`과 **일치**(노트북·서버, 모든 실행의 시작·끝에서 재확인) | §2, 파일 해시는 §9 |
| 실제 `fastcore/` 구현 | **NOT_TESTED** — 열지도 실행하지도 않았다. 이 보고서는 검사기 자체의 보고다 | §1 |
| 검사 내용 | `plan` 프로파일 **23,452 사례** = 계획 FC-1/2 격자 14,000(280조합 × 50상자, DI·GD 같은 상자) + 경계 입력 8,192 + 정밀도 초과 입력 448 + 계획 밖 812. 입력 범주 28종, 정밀도 18종(8~512), 차수 K 12종(0~41), 체 수 1~4, DI·GD 두 경로 | §3 |
| 추가 점검 | 입력 불변 / 호출 뒤 정밀도 불변 / `system.dim` / 시스템을 **다른 정밀도에서 만들어** 호출(결정 A) / 시스템 재사용 + 재호출 / 작업자 충돌 복구 / 오라클 해시 재확인 / API 경고(잘못된 질량의 예외 종류, 튜플·제너레이터 질량) | §2 |
| 단위 시험 | 비교기·생성기·불러오기 **43개 통과** (노트북·서버) | §5 |
| 오탐(참조 대상을 틀렸다고 함) | **0건** — 참조 3개(원본 감싸기, 새 객체로 재구성, 가짜 공장 무패치 사본) × 1,516 사례. 또 `plan` 프로파일(23,452 사례) 전체를 `ref_rebuilt`에 45분 동안 돌려 **9,808 사례(42%, 사례 순서가 섞여 있어 무작위 표본) 불일치 0, 충돌 0** | §5.1 |
| 가짜 적발 | **44개 중 44개** — 자기시험 프로파일(1,516 사례)에서 43개, 나머지 1개(`exp_wrap32`)는 `plan` 프로파일에서 사례 #252 로 적발. 동치 변형 3개(비트 동일)는 불일치 0 = 정상 | §5.2~5.4 |
| 드문 효과 가짜 | 606 사례 중 1~2%만 걸리는 가짜 3개(`pairs_reversed`, `consts_cached`, `sqrt_exact_ignored`) 포함 — 검사기가 작은 표본으로는 놓칠 수 있는 유형을 `plan` 규모에서 잡는다는 근거 | §5.5 |
| 속도 | 자기시험(50개 대상 동시) 약 8분 / 민감도 약 2분(서버 2코어). `plan` 프로파일 1개 대상: 오라클만 약 47분(추정) + 입력 생성·비교, 합계 약 1.1시간(추정) | §7 |
| 실행 명령 | `python3 -I research/pythagorean-3body/fastcore_check/check_fc1.py --target research/pythagorean-3body/fastcore/fast_system.py --profile plan --workers 2 --outdir <새 폴더>` | §8 |
| 종료 코드 | 0 모두 PASS / 1 불일치·충돌 / 2 미실행 / 3 검사기 오류·봉인 해시 불일치 | §4 |


## 1. 쉬운 설명

- **비유(실제와 구분):** 같은 주문서를 받은 두 공장이 만든 제품을 **한 알갱이까지** 맞대어 보는 검수원이다. 원래 공장(원본 파이썬)이 정답지이고, 새 공장(C++)은 같은 순서·같은 반올림으로 만들어야 한다. 마지막 자릿수 하나만 달라도 불합격이다. 검수원은 평범한 주문서만 내지 않는다. 일부러 까다로운 주문서(숫자 길이가 딱 경계인 것, 0에 걸친 것, 두 천체가 부딪히는 것, 지수가 엄청 큰 것, 피타고라스 3-4-5의 딱 떨어지는 거리 등)를 많이 낸다.
- **실제:** 입력은 구간 상자 목록 `X`(각 끝점이 `(가수, 지수)` 정수 쌍)이다. `taylor_step.make_nbody(masses).coeffs(X, K)`(오라클)와 `make_nbody_fast(masses).coeffs(X, K)`(대상)의 결과를 **모든 끝점의 정수 쌍, 리스트 모양, 객체 종류, `None` 위치, 예외 종류**까지 비교한다.
- **검수원이 믿을 만한지도 확인했다.** 일부러 틀리게 만든 가짜 구현 44개를 검수원이 전부 잡는지, 올바른 참조 구현은 틀렸다고 하지 않는지 시험했다(§5). 가짜는 원본 소스의 **한 군데만** 바꾼 사본이라, 적발은 그 한 군데 때문이라고 말할 수 있다.
- **독립성:** 나는 `fastcore/` 안의 C++ 소스·빌드 파일·구현자 시험·보고서를 **열지 않았다.** 폴더 안의 파일 목록도 보지 않았다(작업 시작 때 `fastcore/` 폴더가 아직 없음을 확인한 것이 전부). 사용한 것은 계획서(해시 확인), 원본 `dyad.py`·`taylor_step.py`·`lohner_step.py`(오라클), 이전 검사 보고서(`R1_CHECK_REPORT_KO.md`)·`check_r1.py`의 형식, 공유 메모 폴더의 이전 검사 세션 인계 메모 1건(`r3c-check_handover.md`)뿐이다.

## 2. 무엇을 검사하나 — 비트 동일의 정의

| 비교 대상 | 같아야 하는 것 | 다르면 부르는 이름 |
|---|---|---|
| `DI` 끝점 `_l`, `_h` | `tuple(int, int)` 이고 `(가수, 지수)`가 **정확히 같음**(표준형 포함: 같은 값이라도 `(6,-1)` 과 `(3,0)` 은 다름) | `endpoint_value` / `endpoint_repr`(값은 같고 표현만 다름) / `endpoint_form`(튜플·`int` 형이 아님, `bool`·`float`·`mpz` 포함) |
| 결과 객체의 클래스 | 정확히 `dyad.DI`, `lohner_step._GD`(하위 클래스·다른 복사본 불가) | `type` |
| `_GD` | 값 `v`, 기울기 목록 `g` 전체 — `None`(정확한 0)의 **자리**까지 | `none_position` |
| 리스트 모양 | 성분 수·차수 수·기울기 수, 바깥·안쪽이 `list` | `length`, `form`(튜플 등 모양만 다름) |
| 예외 | 정확히 같은 클래스와 같은 `args` (예: `taylor_step._Singular`) | `exception_type`, `exception_args`, `result_vs_exception` |
| 입력 `X` | 호출 뒤에도 그대로(리스트·`_GD.g` 포함) | `input_mutated` |
| 정밀도 | 호출 뒤 `dyad.get_prec()` 불변 | `prec_side_effect` |
| `system.dim` | `4N` | `api_dim` |
| 시스템 생성 | 올바른 질량이면 오류 없이 생성 | `make_error` |

- **호출 시점 정밀도(R0 결정 A):** 시스템을 **다른 정밀도에서 만든 뒤** 사례의 정밀도로 바꿔 부른다(시스템 생성 정밀도는 사례마다 무작위). 또 40%의 사례에서는 시스템 객체를 **여러 사례에 걸쳐 재사용**하고, 재사용 사례는 같은 system에 **새 입력으로 한 번 더** 불러 다시 비교한다(내부 캐시·정적 버퍼 오염 점검).
- **오라클이 못 푸는 입력은 비교하지 않는다:** 오라클 자신이 `MemoryError` 등 자원 한계로 죽는 사례는 `INVALID`로 따로 세고, 1%를 넘으면 판정을 `UNRESOLVED`로 낮춘다.
- **충돌(세그폴트 등):** 작업자를 별도 프로세스로 돌리고 결과를 줄마다 바로 기록한다. 작업자가 죽으면 죽은 (사례, 대상)을 `crash`로 기록하고(불일치로 센다) 나머지를 이어서 돌린다.
- **봉인 해시:** 시작 때와 끝에 `GATE_FC1_SEAL.json`의 해시와 오라클 파일(`dyad.py`, `taylor_step.py`, `lohner_step.py`, 계획서)의 해시를 비교한다. 다르면 종료 코드 3.

## 3. 입력 만들기

- 사례 하나 = (경로 `DI`/`GD`, 체 수 N, 질량, 정밀도 p, 차수 K, 상자 시드, 시스템 재사용 여부, 시스템을 만드는 정밀도). 전부 `fc1_cases.py`가 시드 `20261012`에서 **결정적으로** 만든다(정수 난수만 써서 노트북·서버에서 같은 사례). 사례 번호만 알면 `--replay 번호`로 그대로 재현한다.
- `GD` 경로는 같은 상자를 `lohner_step._seed`로 감싼 입력이다(블록 A는 DI와 GD가 **같은 상자**를 쓴다 — 계획 FC-2의 "FC-1과 같은 입력").
- 질량은 `m345`((3,4,5), 2체는 (3,4)/(4,5)/(3,5))와 `rand`(무작위 유리수: 정수, 작은 분수, 2의 거듭제곱 분모, 1/3·22/7 같은 비-이진 분수, 10²⁵, 10⁻²⁵, 정밀도 경계 정수 2^p±1, …). `rand`의 절반은 **되풀이되는 고정 풀**에서 뽑는다 — 같은 비-이진 질량이 여러 정밀도에서 다시 나와야 "상수를 첫 호출 정밀도로 캐시하는" 오류가 드러난다(§6 발견 4).

| 블록 | 내용 | `plan` 프로파일 사례 수 |
|---|---|---:|
| A | 계획 4절 FC-1/FC-2 격자: 경로 2 × 정밀도 {96,112,128,144,160,192,256} × K {1,5,20,30,40} × 질량 {m345, rand} × 체 수 {2,3}, **조합마다 50상자**(아래 A 범주 10개를 돌려 가며) | 14,000 (= 280조합 × 50) |
| B | 경계 위주: B 범주 16개 × 정밀도 {64,…,256} 8개 × K {0,1,2,5,20,30,40} × 체 수 {2,3} × 경로 2 × 반복(K별 3·3·3·3·2·1·1회) | 8,192 |
| O | **정밀도보다 긴 가수**(p+1, p+2, p+3, p+16, 2p 비트)를 가진 입력. 실제로는 `taylor_step`이 막는 입력이지만 오라클 동작은 정의돼 있다(별도 `overprec`로 센다) | 448 |
| X | 계획 밖 추가(`extra`): 체 수 1·4, 질량 0·음수, 정밀도 {8,24,53,63,65,127,129,255,257,512}, K {1,3,10,21,31,41}, 일반 기울기 입력 | 812 |
| 합계 | | **23,452** |

| 범주 | 만드는 입력 | 노리는 오류 |
|---|---|---|
| A: `pt` | 폭 0 점 구간(가수가 p비트 꽉 찬 무작위) | 점 구간 경로, 반올림 |
| A: `w30` `w10` `w3` `w01` | 반폭 10⁻³⁰, 10⁻¹⁰, 10⁻³, 0.1 (×1~2) | 폭에 따른 올림·내림 |
| A: `mixw` `fullrand` | 성분마다 폭·가수 모양(꽉 찬/성긴/짧은/전부 1/2^(p−1)+1/절반)을 무작위 | 일반 경로 |
| A: `zs` | 한 좌표의 차이가 0을 가로지르되 거리 제곱의 하한은 양수 + 0을 가로지르는 속도 | `sqr()`의 음이 아닌 하한, 부호 있는 곱 |
| A: `sgn` | 속도를 `[0,w]`, `[-w,0]`, `[-w,w]`, `(0,0)`, 음수 점으로 | 부호 경계, 0 끝점 |
| A: `hugev` | 속도 상자 반폭 최대 10⁶ | 아주 넓은 상자 |
| B: `allones` `pow2p1` `sparse` | 가수 2^p−1 / 2^(p−1)+1 / 켜진 비트 몇 개, 폭 1 ulp 안팎 | 자리올림, 가수 비트 경계(정확히 p비트), 정확한 결과 |
| B: `tie_sub` | 첫 연산(좌표 차)이 정확히 p+1비트가 되도록 짝 좌표를 만듦 | 1비트 버림 동률(§6 발견 3) |
| B: `zero_vel` `zero_end` | 속도 전부 0 / 끝점이 정확히 0인 상자 | 0 근처·0 끝점 |
| B: `hugep` | 좌표 ~10⁶, 위치 상자 반폭 최대 2¹⁷(충돌 아님) | 아주 넓은 위치 상자 |
| B: `exp_scaled` | 위치 ×2^E, 속도 ×2^(−E/2) (케플러 척도 불변), \|E\| ≤ 2²⁷. 차수 K에서 출력 지수 ≈ −1.5·K·E 가 **2³¹을 넘음** | 32비트 지수 오버플로 |
| B: `exp_mixed` | 성분마다 독립 지수 이동 ±300 | 지수 정렬(큰 자리 이동) |
| B: `near_coll` | 0번 천체 (0,0), 1번 천체 (2^−t₁, 2^−t₂), t ≤ 400 | 거리 제곱 하한이 아주 작은 양수 |
| B: `touch_x` `touch_coll` | 한 좌표만 접촉(비특이) / 두 좌표 모두 접촉 → 하한이 **정확히 0** | `_Singular` 판정 경계(`<=`) |
| B: `coll_full` `coll_point` `identical` | 충돌 배치(거리 구간이 0 포함) | FC-3: 같은 예외 |
| B: `pyth` | 피타고라스 3-4-5 시작 배치 (1,3),(−2,−1),(1,−1) — 거리 제곱 9·16·25가 **완전제곱**. 2의 거듭제곱 척도, 천체 순서 섞기, 점/1ulp/2^(−p/2)/10⁻³⁰ 폭 | sqrt의 '정확히 떨어짐' 분기(§6 발견 2) |
| O: `overprec` | 위 'fullrand'를 p보다 긴 정밀도로 만든 뒤 정밀도를 p로 낮춰 호출 | 입력을 반올림하지 않고 통과시키는 경로 |
| X: `gd_generic` | 기울기가 `None`·`DI(0)`(명시적 0)·일반 구간이 섞인 `_GD` 입력 | `None`/`DI(0)` 구분 |

- `K=0`은 오라클이 입력 객체를 그대로 돌려주는 경로다(계산 없음). `K=1`은 `taylor_step`의 사전 상자 검사(`F(B)`)가 쓰는 값, `K=차수`와 `K=차수+1`은 다항식·나머지 계산이 쓰는 값이다.
- 프로파일: `smoke` 223(노트북 확인용, 수 초) / `selftest` 1,516(자기시험) / `sens` 606(가짜별 민감도, 적발 후에도 계속) / **`plan` 23,452(FC-6 본 실행)** / `deep` 93,712(선택).


## 4. 판정과 종료 코드

- 대상별 판정: **PASS** = 모든 사례 일치(미실행·충돌 없음), **FAIL** = 불일치 또는 충돌 1건 이상, **UNRESOLVED** = 불일치는 없으나 미실행·`INVALID` 과다.
- 종료 코드: `0` 모두 PASS / `1` 불일치·충돌 있음 / `2` 불일치는 없으나 미실행 / `3` 검사기 오류 또는 오라클 봉인 해시 불일치.
- **PASS는 "시험한 범위 안에서 불일치 없음"이다. 형식 증명이 아니다.** 시험하지 않은 입력에서 다를 가능성은 남는다(§7).

## 5. 검사기 자기시험 (검수원이 믿을 만한가)

방법: 같은 사례에 **참조 대상 3개와 가짜 대상 44개, 동치 변형 3개**를 한꺼번에 돌린다(오라클은 사례마다 한 번만 계산). 가짜는 `fakes/_mutant.py`가 **원본 `dyad.py`·`taylor_step.py`·`lohner_step.py`의 소스 글자에서 정확히 한 군데**(anchor가 원본에 정확히 1번 나오는 문자열)만 바꾼 사본을 새 이름의 모듈로 실행해 만든다(anchor가 0번 또는 2번 이상이면 즉시 오류 — 원본이 바뀌었는데 가짜가 엉뚱한 곳을 망가뜨리는 일을 막는다). 사본의 결과는 원본 클래스(`dyad.DI`, `lohner_step._GD`)로 옮겨 담아 검사기 입장에서는 "원본 클래스를 돌려주는 `make_nbody_fast`"로 보인다. 가짜 정의 표는 `fakes/_gen_fakes.py`에 있다.

- 실행: 서버(GCP, Python 3.13.5), `--selftest --profile selftest --workers 2`, **1,516 사례**(plan 1,320 / extra 188 / overprec 8), 벽시계 493초, 작업자 최대 메모리 약 91MB. 적발된 가짜는 이후 사례를 건너뛰고(`--early-exit`), 참조 대상·동치 변형은 전 사례를 돌린다. 결과 파일: `FC1_CHECK_SELFTEST_RESULT.json`.
- 오라클 결과 분포: 정상 1,308 / `_Singular` 208 / 기타 0 (`INVALID` 0건).
- **`plan` 규모 오탐 시험:** `plan` 프로파일(23,452 사례)을 `ref_rebuilt`에 45분(2,703초) 동안 돌렸다 → **9,808 사례**(42%; 사례 순서가 섞여 있어 무작위 표본: plan 9,285 / extra 342 / overprec 181)에서 **불일치 0, 충돌 0, `INVALID` 0**, 오라클 해시 시작·끝 일치, 작업자 최대 메모리 약 179MB. 시간 제한 때문에 판정은 `UNRESOLVED`(종료 코드 2) — 설계된 동작이며 나머지 13,644 사례는 실행하지 않았다.
- 작업자 충돌 복구도 실제로 시험됐다: `process_exit` 가짜가 사례 37·152에서 프로세스를 죽였고, 검사기는 이를 `crash`로 기록하고 나머지를 이어 돌렸다.

### 5.1 참조 대상 (불일치 0 이어야 함)

| 대상 | 설명 | 실행 사례 | 불일치 | 충돌 | 판정 |
|---|---|---:|---:|---:|---|
| `ref_mutant_identity` | 패치 없는 가짜 공장 사본(공장이 결과를 바꾸지 않음을 보이는 참조 대상) | 1516 | 0 | 0 | PASS |
| `ref_original` | 원본 make_nbody 를 그대로 감싼 참조 대상(오탐 0 이어야 함) | 1516 | 0 | 0 | PASS |
| `ref_rebuilt` | 원본 결과를 새 DI/_GD/리스트 객체로 재구성해 돌려주는 참조 대상(오탐 0 이어야 함) | 1516 | 0 | 0 | PASS |

### 5.2 가짜 대상 (모두 적발해야 함)

| # | 가짜 | 무엇을 망가뜨렸나 | 적발 | 처음 적발한 사례 | 차이 종류 |
|---:|---|---|:-:|---|---|
| 1 | `assoc_order` | 곱셈 결합 순서 변경 (s*w)*c -> s*(w*c) | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 2 | `canon_off` | _canon 이 꼬리 0 비트를 안 지움 — 같은 값이지만 (가수, 지수) 표현이 다름 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_repr |
| 3 | `consts_at_make` | 질량 상수를 system 을 만들 때의 정밀도로 고정(결정 A 위반: 호출 시점 정밀도여야 함) | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 4 | `consts_cached` | 질량 상수를 첫 호출 때의 정밀도로 캐시하고 이후 호출에 재사용 | 예 | #322 A/pt/GD n=2 K=5 p=96 | endpoint_value |
| 5 | `container_tuple` | 결과를 리스트가 아닌 튜플로 돌려줌(모양 차이만; 값은 같음) | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | form |
| 6 | `conv_reversed` | 코시 곱 합산 순서를 j=k..0 으로 뒤집음 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 7 | `dim_wrong` | system.dim 이 틀림(4N 대신 4N+1) — API 모양 점검 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | api_dim |
| 8 | `div_int_recip` | a/(k+1) 를 a * DI(1/(k+1)) 로(역수를 먼저 반올림) | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 9 | `dn_nearest` | _dn 이 내림 대신 가장 가까운 값으로 반올림(반올림 동률 포함) | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 10 | `dn_trunc` | _dn(하한 내림)이 음수에서 0 쪽으로 자름(floor 대신 truncate) — GMP tdiv 를 쓴 이식 실수 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 11 | `exc_type` | 특이점에서 _Singular 대신 ZeroDivisionError | 예 | #30 B/touch_coll/GD n=3 K=5 p=64 | exception_type |
| 12 | `exc_wide` | 특이점에서 예외 대신 거대한 구간을 돌려줌 | 예 | #30 B/touch_coll/GD n=3 K=5 p=64 | result_vs_exception |
| 13 | `exp_wrap32` | 결과 끝점의 지수를 부호 있는 32비트로 감아 돌려줌(어댑터가 지수를 int32 로 저장하는 실수) — 지수가 ±2^31 을 넘는 입력에서만 틀림 | **아니오** |  |  |
| 14 | `gd_div_split` | _GD 나눗셈 미분 (da - q*db)/b 를 da/b - (q*db)/b 로 쪼갬 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 15 | `gd_mul_none` | _GD 곱셈에서 한쪽 기울기가 None 일 때 다른 쪽 항을 버림 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value,none_position |
| 16 | `gd_rsub_sign` | _GD 의 (상수 - _GD) 에서 기울기 부호를 뒤집지 않음 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 17 | `gd_rtruediv_sign` | _GD 의 (상수 / _GD) 미분에서 음수 부호를 빠뜨림 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 18 | `gd_sqr_factor` | _GD.sqr() 의 미분에서 2배를 빠뜨림 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 19 | `gd_sqrt_factor` | _GD.sqrt() 의 미분 ds/(2 sqrt s) 에서 2를 빠뜨림 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 20 | `k_short` | 재귀를 한 차수 덜 돌려 계수 목록이 짧음 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | length |
| 21 | `mass_hiprec` | 질량 상수 DI(m) 를 호출 정밀도보다 64비트 높게 만듦 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 22 | `mul_minmax` | 구간 곱의 최소·최대 선택이 지수 정렬 없이 가수만 비교 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 23 | `mutates_input` | 입력 리스트 X 를 제자리에서 바꿈(끝에 원소 추가) | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | input_mutated |
| 24 | `none_to_di0` | _GD 입력의 None(정확히 0)을 DI(0) 으로 바꿔 계산 — 값은 같지만 None 위치가 달라짐 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | none_position |
| 25 | `pairs_reversed` | 가속도 합산에서 쌍(pair) 순서를 뒤집음 (3체부터 반올림이 달라짐) | 예 | #16 X/w10/GD n=4 K=10 p=8 | endpoint_value |
| 26 | `pow_div_twice` | w_[k] 의 나눗셈을 한 번 대신 두 번으로 쪼갬(반올림 두 번) | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 27 | `prec_at_make` | 모든 계산을 system 을 만든 때의 정밀도로(호출 시점 정밀도 무시) | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 28 | `prec_clamp128` | 정밀도가 128비트를 넘으면 128로 잘라서 계산(고정 폭 한계 이식 실수) | 예 | #4 B/sparse/DI n=2 K=20 p=256 | endpoint_value |
| 29 | `prec_plus_one` | 내부 계산을 요구 정밀도보다 1비트 높게 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 30 | `prec_side_effect` | 계산 뒤에 전역 정밀도를 바꿔 놓음 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | prec_side_effect |
| 31 | `process_exit` | 특정 입력에서 프로세스가 죽음(세그폴트 모사: os._exit) | 예 | #152 B/zero_vel/DI n=2 K=20 p=64 | 충돌 |
| 32 | `qdn_ceil` | 나눗셈 하한 몫을 내림 대신 올림으로 계산(반올림 방향 효과가 드문 변형) | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 33 | `qdn_trunc` | 나눗셈 하한 몫을 내림 대신 0 쪽으로 자름(GMP tdiv 이식 실수) — 분자가 음수일 때만 틀림 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 34 | `qup_nobump` | 나눗셈 상한 몫에서 나머지가 있어도 +1 하지 않음 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 35 | `singular_lt` | 특이점 판정 경계 오류: lo <= 0 대신 lo < 0 (거리 하한이 정확히 0 인 접촉 배치) | 예 | #30 B/touch_coll/GD n=3 K=5 p=64 | exception_type |
| 36 | `sqconv_cauchy` | 제곱 코시 곱에서 대칭 쌍을 한 번 계산해 2배 하지 않고 두 번 곱해 모두 더함 | 예 | #2 B/exp_scaled/GD n=2 K=20 p=64 | endpoint_value |
| 37 | `sqconv_mid_mul` | _sq_conv 의 가운데 항을 sqr 대신 곱셈으로 | 예 | #10 B/pow2p1/DI n=3 K=20 p=256 | endpoint_value |
| 38 | `sqr_as_mul` | sqr() 를 x*x 로 대체 — 0 을 가로지르는 구간에서 하한이 음수로 나옴 | 예 | #10 B/pow2p1/DI n=3 K=20 p=256 | endpoint_value |
| 39 | `sqr_negative_branch` | x**2 에서 구간이 전부 음수일 때 끝점을 맞바꿈(하한이 상한보다 커짐) | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 40 | `sqrt_exact_ignored` | sqrt 의 '정확히 떨어짐' 판정을 빼먹음 — 완전제곱(예: 3-4-5 시작 배치)에서만 상한이 한 칸 넓어짐 | 예 | #20 B/allones/GD n=3 K=2 p=64 | endpoint_value |
| 41 | `sqrt_up_floor` | sqrt 상한이 완전제곱이 아닐 때 +1 을 빼먹음(올림 누락) | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 42 | `stateful` | 같은 system 의 2번째 호출부터 정밀도 +1 로 계산(내부 상태 오염 모사). 재사용 시스템의 재호출 점검으로만 적발됨 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 43 | `up_floor` | _up(상한 올림)이 올림 대신 내림 | 예 | #0 X/fullrand/DI n=3 K=1 p=64 | endpoint_value |
| 44 | `zero_to_none` | 결과 기울기의 정확한 0 인 DI(0) 을 None 으로 바꿈 | 예 | #8 B/zero_vel/GD n=3 K=20 p=128 | none_position |

**적발 43 / 44** (자기시험 프로파일 기준; 44번째 `exp_wrap32` 는 §5.4 에서 `plan` 프로파일로 적발)

'처음 적발한 사례'는 작업자 두 개가 각자 처음 적발한 사례 중 기록 순서상 첫 줄이다(사례 번호 최소값이 아닐 수 있다).

### 5.3 동치 변형 (비트 단위로 같은 결과 — 적발되지 않는 것이 정상)

| 변형 | 설명 | 불일치 | 판정 |
|---|---|---:|---|
| `eq_acc_plus_acc` | [동치] 2배를 acc*2 대신 acc+acc 로 (2의 거듭제곱 곱은 반올림과 교환되므로 비트 동일) | 0 | 정상(적발 없음) |
| `eq_acc_times_di2` | [동치] 2배를 정수 2 대신 DI(2) 로 | 0 | 정상(적발 없음) |
| `eq_add_commuted` | [동치] 덧셈의 피연산자 순서 교환(구간 덧셈은 교환법칙이 비트 단위로 성립) | 0 | 정상(적발 없음) |


### 5.4 `exp_wrap32`는 `plan` 프로파일에서만 적발된다

`exp_wrap32`(결과 끝점의 지수를 부호 있는 32비트로 감음)는 자기시험 프로파일에서 **적발되지 않았다.** 이 가짜는 출력 지수가 ±2³¹을 넘어야만 틀리는데, 그러려면 \|E\|=2²⁶이면 K≥30, \|E\|=2²⁷이면 K≥11 이어야 하고, `plan` 프로파일 23,452 사례 중 이런 `exp_scaled` 사례는 **14건**뿐이다. `selftest` 프로파일의 판정 코드는 이를 `기대와 다름`으로 표시하고 종료 코드 1을 낸다(그래서 `FC1_CHECK_SELFTEST_RESULT.json`의 `selftest_ok`는 `false`이며, 그 이유는 이 가짜 하나뿐이다). 검사기 코드는 이 때문에 바꾸지 않았다.

`plan` 프로파일을 이 가짜 하나에 돌리자 (`--selftest --only exp_wrap32 --profile plan`) 사례 **#252** (`B/exp_scaled/GD n=3 K=30 p=112`)에서 적발됐다: 오라클의 지수 `−2147483742`, 대상의 지수 `+2147483554` (32비트 감김), 그 사례에서 차이 2,820곳. 증거: `results/evidence/expwrap_plan2_shard0.jsonl`. 다른 작업자가 같은 가짜를 적발하기 전에 내가 이 실행을 중단했다(적발이 확인됐고 나머지는 같은 사실을 반복하므로).

### 5.5 민감도 — 어느 입력이 어느 가짜를 잡나

`--profile sens` (606 사례, **적발 후에도 계속**; `exp_scaled`는 제외) 를 가짜 44개 + 참조·동치에 돌려 가짜별 적발 비율을 쟀다(서버, 104초). 두 가짜(`exp_wrap32`, `process_exit`)는 이 프로파일에 필요한 입력(큰 지수 / K=20)이 없어 0건이다. 동치 변형 3개와 참조 3개는 0건(정상).

| 가짜 | 적발 사례 수 / 606 | 비율 | DI 경로 | GD 경로 |
|---|---:|---:|---:|---:|
| `eq_acc_plus_acc` | 0 | 0.0% | 0 | 0 |
| `eq_acc_times_di2` | 0 | 0.0% | 0 | 0 |
| `eq_add_commuted` | 0 | 0.0% | 0 | 0 |
| `exp_wrap32` | 0 | 0.0% | 0 | 0 |
| `process_exit` | 0 | 0.0% | 0 | 0 |
| `pairs_reversed` | 6 | 1.0% | 3 | 3 |
| `consts_cached` | 9 | 1.5% | 6 | 3 |
| `sqrt_exact_ignored` | 10 | 1.7% | 5 | 5 |
| `zero_to_none` | 21 | 3.5% | 0 | 21 |
| `sqconv_mid_mul` | 84 | 13.9% | 42 | 42 |
| `sqr_as_mul` | 84 | 13.9% | 42 | 42 |
| `sqrt_up_floor` | 95 | 15.7% | 44 | 51 |
| `exc_type` | 100 | 16.5% | 50 | 50 |
| `exc_wide` | 100 | 16.5% | 50 | 50 |
| `singular_lt` | 100 | 16.5% | 50 | 50 |
| `pow_div_twice` | 103 | 17.0% | 38 | 65 |
| `consts_at_make` | 108 | 17.8% | 53 | 55 |
| `gd_div_split` | 114 | 18.8% | 0 | 114 |
| `mass_hiprec` | 129 | 21.3% | 62 | 67 |
| `sqconv_cauchy` | 148 | 24.4% | 68 | 80 |
| `stateful` | 150 | 24.8% | 74 | 76 |
| `conv_reversed` | 155 | 25.6% | 76 | 79 |
| `div_int_recip` | 156 | 25.7% | 77 | 79 |
| `prec_clamp128` | 164 | 27.1% | 82 | 82 |
| `qdn_trunc` | 183 | 30.2% | 54 | 129 |
| `gd_rsub_sign` | 188 | 31.0% | 0 | 188 |
| `gd_mul_none` | 192 | 31.7% | 0 | 192 |
| `gd_rtruediv_sign` | 192 | 31.7% | 0 | 192 |
| `gd_sqr_factor` | 192 | 31.7% | 0 | 192 |
| `gd_sqrt_factor` | 192 | 31.7% | 0 | 192 |
| `assoc_order` | 215 | 35.5% | 94 | 121 |
| `none_to_di0` | 254 | 41.9% | 0 | 254 |
| `qdn_ceil` | 280 | 46.2% | 115 | 165 |
| `qup_nobump` | 285 | 47.0% | 115 | 170 |
| `sqr_negative_branch` | 311 | 51.3% | 154 | 157 |
| `prec_at_make` | 319 | 52.6% | 160 | 159 |
| `dn_trunc` | 371 | 61.2% | 179 | 192 |
| `canon_off` | 382 | 63.0% | 190 | 192 |
| `dn_nearest` | 382 | 63.0% | 190 | 192 |
| `mul_minmax` | 382 | 63.0% | 190 | 192 |
| `prec_plus_one` | 382 | 63.0% | 190 | 192 |
| `up_floor` | 382 | 63.0% | 190 | 192 |
| `k_short` | 418 | 69.0% | 208 | 210 |
| `container_tuple` | 506 | 83.5% | 252 | 254 |
| `prec_side_effect` | 506 | 83.5% | 252 | 254 |
| `dim_wrong` | 606 | 100.0% | 302 | 304 |
| `mutates_input` | 606 | 100.0% | 302 | 304 |


- **한 번에 걸리는 가짜**(적발 비율이 높음)는 반올림 방향·표준형·구간 곱·정밀도 취급 오류다. 거의 모든 입력에서 틀리므로 어떤 범주로도 잡힌다.
- **드문 효과 가짜**는 606 사례 중 1~3% 에서만 걸린다: `pairs_reversed` 6건(N=4만, §6-7), `consts_cached` 9건(고정 질량 풀 필요, §6-4), `sqrt_exact_ignored` 10건(`pyth` 8 + `allones` 2, §6-2), `zero_to_none` 21건(GD 경로만). `plan`(23,452 사례, 약 39배)에서는 이들이 수백 번 걸린다는 뜻이다(비례 추정).
- **어떤 범주도 혼자서만 잡는 가짜는 없었다** (이 프로파일에서). 범주마다 42개 중 31~37개를 잡았고, 충돌 범주(`coll_full`, `coll_point`, `identical`, `touch_coll`)는 계산이 예외로 끝나서 값 비교 가짜는 못 잡고, 예외 가짜 3개(`exc_type`, `exc_wide`, `singular_lt`)와 모양·API·상태 가짜 6개(`container_tuple`, `dim_wrong`, `k_short`, `mutates_input`, `none_to_di0`, `prec_side_effect` — K=0 사례 포함)만 잡는다. 즉 충돌 시험(FC-3)은 "같은 예외가 나는가"를 보는 별도 축이다. 범주를 많이 둔 것은 중복이지만, 한 범주가 빠져도 다른 범주가 잡는다는 뜻이다. 반면 범주별 사례 수가 적은(2~10건) 범주는 통계적으로 얇다.
- 가짜별·범주별 전체 표는 `analyze_run.py <폴더> --matrix` 로 다시 만들 수 있다.


## 6. 검사기를 만들면서 알게 된 것 (고친 것 포함)

1. **오라클이 지수를 어디까지 받는가 — 지수 극단 시험의 상한.** 처음에는 지수를 2⁷⁰까지 줬더니 오라클 자신이 `MemoryError`로 죽었다. `dyad.lo_q`가 `Fraction(m, 1 << -e)`를 만들어, 특이점 판정(`s0.lo_q <= 0`)에서 \|지수\|만큼의 비트 정수를 실제로 만들기 때문이다. 그래서 `exp_scaled` 범주는 \|E\| ≤ 2²⁷로 제한했다(오라클이 0.1~0.6초에 풀고 메모리도 작다). 그래도 차수 K에서 출력 지수는 ≈ −1.5·K·E 라서 E=2²⁷, K=20이면 약 4×10⁹ — **32비트 부호 있는 정수(±2.1×10⁹)를 넘는다.** 이 범위의 효과는 `exp_wrap32` 가짜로 확인했다(§5.4). 64비트 지수 오버플로는 오라클이 못 풀어 **시험 불가(UNKNOWN)** 다. 오라클이 `MemoryError`를 내는 사례가 생기면 비교하지 않고 `INVALID`로 따로 센다.
2. **sqrt의 "정확히 떨어짐" 분기는 무작위 입력으로는 거의 실행되지 않는다.** `_sqrt_up`은 `isqrt`가 완전제곱일 때와 아닐 때 다르게 동작한다. 일반 범주 4개(`fullrand`, `w10`, `sparse`, `pt`) × 24회 실행에서 이 분기 적중은 0회였다. 이 프로젝트의 실제 시작 배치 (1,3),(−2,−1),(1,−1) 는 거리 제곱이 9·16·25 라서 이 분기를 탄다(`pyth` 범주 24회 실행에서 27회). 그래서 `pyth` 범주를 추가했다. 이 분기만 망가뜨린 가짜(`sqrt_exact_ignored`)는 606 사례 민감도 실행에서 `pyth` 8건, `allones` 2건에서 적발됐다 — `pyth` 밖에서도 우연히 완전제곱이 나올 수는 있지만 드물다.
3. **"반올림 동률"은 표준형 때문에 1비트 버림뿐이다.** 끝점 가수는 홀수로 정규화되므로 덧셈·곱셈 결과에서 버려지는 비트열은 항상 1로 끝난다. 그래서 `100…0` 같은 긴 동률은 생길 수 없고 정확한 동률은 1비트 버림(= 홀수 결과)뿐이며, 이는 어느 범주에서나 `_dn` 호출의 약 15%다. 처음 설계한 `tie_sub`(긴 동률을 노린 구성)는 측정해 보니 무작위와 차이가 없었다. 그대로 두되 별도 이점을 주장하지 않는다. 반올림 방식 자체(내림/올림/가장 가까운 값/0 쪽 절단)는 `dn_trunc`, `up_floor`, `dn_nearest`, `qdn_trunc`, `qdn_ceil`, `qup_nobump` 가짜로 확인했다.
4. **질량이 매번 새 무작위 유리수이면 "상수 캐시" 오류를 못 잡는다.** `consts_cached` 가짜(질량 상수를 첫 호출 정밀도로 저장)는 처음 자기시험에서 **적발되지 않았다.** 정수·이진 분수 질량은 정밀도가 달라도 같은 `DI`가 되고, 무작위 질량은 두 번 다시 나오지 않기 때문이다. 같은 비-이진 질량 묶음(1/3, 2/7, 22/7, …)이 여러 정밀도에서 되풀이되는 **고정 풀**을 넣은 뒤 적발됐다. (606 사례 민감도 실행에서도 9건뿐 — 드문 효과다.)
5. **검사기 자신의 버그 하나를 자기시험이 드러냈다.** 출력 폴더를 재사용하면 이전 실행의 기록이 남아 있어(충돌 복구용으로 작업자는 기록된 사례를 건너뛴다) 대상 목록이 달라진 새 실행의 열과 어긋났다. 그래서 가짜 하나가 "적발 안 됨"으로 잘못 나왔다. 비어 있지 않은 출력 폴더는 거부하도록 고쳤다.
6. **작은 프로파일은 일부 가짜를 놓친다.** `smoke`(223 사례)는 `prec_clamp128`(128비트 초과 정밀도 필요), `process_exit`(K=20 필요)을 놓쳤다. `selftest`(1,516)부터 이 둘을 잡는다. `exp_wrap32`는 `selftest`로도 못 잡고 `plan`에서만 잡힌다(§5.4). 본 실행은 `plan`을 쓴다.
7. **N=3 이하에서는 쌍(pair) 합산 순서가 비트에 영향을 주지 않는다.** `pairs_reversed` 가짜는 606 사례 중 **N=4 사례 6건에서만** 적발됐고 N=2·3(592건)에서는 한 번도 안 걸렸다. 이유: N≤3이면 각 천체가 정확히 두 쌍에만 속하고, 가속도 합산은 `0 + 항1 + 항2` 의 두 항 합이라 순서를 바꿔도 비트가 같다(구간 덧셈은 교환법칙이 비트 단위로 성립). 즉 **계획의 N∈{2,3}만으로는 이 순서 오류를 구별할 수 없고(이 프로젝트에서는 무해), N≥4(`extra`)에서야 드러난다.** 그래서 `extra` 계층에 N=4 사례를 일부러 넣었다(`plan` 프로파일 200건).
8. **지수 오버플로 가짜를 처음엔 계산 중간에 넣었더니 쓸모없이 느려졌다.** 곱의 지수를 중간에 32비트로 감으면 이후 덧셈이 2³² 비트 정수 이동을 하느라 서버 메모리를 수 GB 쓸 뻔했다. 실제 int32 지수 버그는 이런 폭주 없이 틀린 값을 낼 뿐이므로, 가짜를 "결과 끝점의 지수를 돌려줄 때 32비트로 감는" 형태로 바꿨다(§5.4). 그 서버 실행은 내가 직접 중단했다.


## 7. 한계 — 확인 / 추정 / UNKNOWN

**확인(내가 직접 실행하거나 읽어서 본 것)**
- 계획서와 오라클 3파일의 sha256이 `GATE_FC1_SEAL.json`과 일치한다(노트북·서버 각각, 검사기가 실행 시작·끝에 다시 확인한다).
- 비교기·생성기·불러오기 단위 시험 43개 통과(노트북·서버).
- 참조 대상 3개는 §5.1의 사례에서 불일치 0, 가짜 44개는 전부 적발(43개는 자기시험 프로파일, 1개는 `plan` 프로파일), 동치 변형 3개는 불일치 0 (§5).
- 오라클은 지수가 약 2²⁷ 비트까지는 계산하지만 그보다 크면 `MemoryError`로 죽는다(`dyad._to_q`가 `1 << -e`를 만든다). 그래서 **지수 극단 시험은 \|E\| ≤ 2²⁷** 이다. 이 범위에서도 출력 지수는 최대 약 ±4×10⁹ 까지 가서 **32비트 지수 오버플로는 시험된다**(§6 발견 1, `exp_wrap32` 가짜로 확인).

**추정(근거는 있으나 직접 재지 않은 것)**
- `plan` 프로파일의 소요 시간: 오라클 시간 합계 ≈ 5,600초(1코어; 서버 자기시험 1,516 사례의 실측 오라클 시간에서 K=30·40을 K² 로 외삽한 값). 2코어 작업자면 오라클만 약 47분이고, 대상 호출·입력 생성·비교가 더해진다. 부분 실행 실측(서버 2코어, 서버는 다른 작업으로 거의 가득 참): `plan` 프로파일 9,808 사례(무작위 표본)에서 오라클 시간 합계 2,214초 → 전체 환산 약 5,300초(앞의 외삽 5,600초와 5% 이내), 대상 호출 2,238초, 입력 생성·비교 등 나머지 약 950초(작업자 합계 5,400초 중 약 18%). `ref_rebuilt`(호출마다 오라클을 한 번 더 불러 새 객체로 재구성하므로 오라클 약 2배 비용)는 분당 약 218 사례 → `plan` 전체 약 1.8시간. 빠른 C++ 대상이면 대상 호출 시간이 거의 0이므로 **약 1.1시간(추정; 오라클 약 0.77시간 + 오버헤드)** 으로 본다. 서버가 더 바쁘면 늘어난다.
- 실제 C++ 이식 오류가 내 가짜와 비슷한 부류(반올림 방향, 합산 순서, 정밀도 취급, 예외, None 처리, 상태 오염, 지수 폭)라는 가정. 전혀 다른 부류의 오류는 무작위·경계 입력이 우연히 건드려야 잡힌다.

**UNKNOWN(이 검사기로 알 수 없는 것)**
1. **실제 `fastcore/` 구현이 통과하는지** — 실행하지 않았다(NOT_TESTED).
2. \|지수\| > 2²⁷ 의 입력(특히 64비트 지수 오버플로): 오라클이 계산하지 못해 시험할 수 없다.
3. K > 41, 정밀도 > 512, 체 수 > 4, `DI`가 아닌 입력(정수·`float`), `X` 길이 불일치, 음수 K: 시험하지 않았다. 질량이 `float`·`bool`·문자열일 때의 오류 종류는 API 점검(§2)이 경고로만 알려준다.
4. 대상이 내부에서 오라클(`taylor_step`)을 불러 결과를 그대로 돌려주는 경우: 결과만 보면 구별할 수 없다(구현이 정말 C++를 쓰는지는 FC-5 속도와 코드 검토의 몫).
5. 수십만 번 호출 뒤의 메모리 누수·상태 오염: 작업자 최대 메모리를 기록할 뿐 판정하지 않는다. 스레드 안전성, 다른 CPU/컴파일러 환경은 시험하지 않았다.
6. **동률 반올림**: 표준형(홀수 가수) 때문에 정확한 동률은 1비트 버림에서만 생기고, 이는 모든 범주에서 반올림 호출의 약 15%를 차지한다(§6 발견 3). `tie_sub` 범주는 이를 첫 연산에 몰아주는 보조일 뿐이며 별도 이점은 측정되지 않았다.
7. 무작위·경계 입력으로는 "시험하지 않은 입력에서 다를 가능성"을 없앨 수 없다. **PASS는 형식 증명이 아니다.** 오라클(원본)의 버그는 못 잡는다 — 계획 5절이 정한 범위다.
8. `INVALID`(오라클이 자원 한계로 못 푼 사례)는 비교하지 않는다. 이 보고서의 모든 실행(자기시험, 민감도, `plan` 부분 실행)에서 `INVALID`는 0건이었다.


## 8. 실행 방법

저장소 루트에서(서버든 노트북이든 같다). 서버 CPU는 동시 2코어 이하(`--workers 2`)로 쓴다.

```bash
# (1) 검사기 자체의 단위 시험 (몇 초)
python3 -I research/pythagorean-3body/fastcore_check/check_fc1_unit_tests.py

# (2) 빠른 확인 — 223 사례, 수 초
python3 -I research/pythagorean-3body/fastcore_check/check_fc1.py \
    --target research/pythagorean-3body/fastcore/fast_system.py --profile smoke --workers 2 \
    --outdir research/pythagorean-3body/fastcore_check/results/fc1_smoke

# (3) FC-6 본 실행 — 23,452 사례 (소요 시간은 §7 참고)
python3 -I research/pythagorean-3body/fastcore_check/check_fc1.py \
    --target research/pythagorean-3body/fastcore/fast_system.py --profile plan --workers 2 \
    --outdir research/pythagorean-3body/fastcore_check/results/fc1_plan \
    --out research/pythagorean-3body/fastcore_check/results/FC1_CHECK_RESULT_plan.json

# (4) 불일치가 난 사례 하나를 자세히 (사례 번호는 위 출력의 '#번호')
python3 -I research/pythagorean-3body/fastcore_check/check_fc1.py \
    --target research/pythagorean-3body/fastcore/fast_system.py --profile plan --replay 12345

# (5) 검사기 자기시험 (참조·가짜 대상 전부, 2코어로 약 7~8분)
python3 -I research/pythagorean-3body/fastcore_check/check_fc1.py --selftest --profile selftest --workers 2 \
    --outdir research/pythagorean-3body/fastcore_check/results/selftest_x
```

- **출력 폴더는 비어 있어야 한다**(이전 실행과 섞이는 것을 막으려고 비어 있지 않으면 종료 코드 3). 매번 새 폴더 이름을 쓴다.
- `--target`은 여러 번 줄 수 있다(오라클은 사례마다 한 번만 계산해 대상 전부에 쓴다). 대상 폴더에 `__init__.py`가 있으면 패키지로(상대 import 가능), 없으면 파일 경로로 불러온다. 대상이 쓰는 확장 모듈(`.so`)은 대상 파일과 같은 폴더에 있으면 된다.
- 결과: 출력 폴더의 `summary.json`(대상별 판정·차이 종류·범주별 집계·첫 불일치 8건·API 경고·작업자 메모리·오라클 해시)과 `shard_*.jsonl`(사례별 기록). `--out`으로 같은 요약을 다른 경로에도 저장한다.
- 진행 상황은 60초마다 표준출력에 나온다. `--time-limit 초`를 주면 시간이 다 됐을 때 남은 사례는 `미실행`(종료 코드 2)이 된다 — 사례 순서가 섞여 있어 시간 제한 실행은 **무작위 표본**이 된다.
- 보조: `analyze_run.py <폴더> --coverage`(범주별 오라클 결과), `--matrix`(대상 × 범주 불일치 비율), `--estimate plan`(오라클 시간으로 소요 시간 추정).


## 9. 파일과 해시

| 파일 | 역할 | sha256 |
|---|---|---|
| `fastcore_check/check_fc1.py` | 검사기 본체(비교기, 실행기, 작업자·충돌 복구, 집계, 자기시험 판정) | `c0126d78e1af7abcf5cc971b6eb02956b52574dee1430c3a71ff3c4fb3a680f3` |
| `fastcore_check/fc1_cases.py` | 입력 생성기(결정적, 오라클 import 없음) | `bc1ffaad933263375492ebc405a51ffdeb5918c96db4dee81a48c9910f8d5a85` |
| `fastcore_check/check_fc1_unit_tests.py` | 비교기·생성기·불러오기 단위 시험 43개 | `a850e9ff7295088166ce2609f3b162da67ee4dce9144144638d6540bd6efef68` |
| `fastcore_check/analyze_run.py` | 실행 폴더 분석(민감도 표, 범주별 오라클 결과, 소요 시간 추정) | `9e2bf3857099178b1cc4973213e7c14dad22a81856a76873701fd92072426a2c` |
| `fastcore_check/fakes/_mutant.py` | 가짜 공장(원본 소스의 한 군데만 바꾼 사본을 만든다) | `af7cabdd1d96d35fe6deaa5f8c5ebd038ddb6890765a539ceab7462c3cb1718e` |
| `fastcore_check/fakes/_gen_fakes.py` | 가짜 44개(+동치 3개) 정의 표와 생성기 | `b6e20252ab540da8024543f0407976b65e97925a07ba5ca9b97ea6a762fb25e6` |
| `fastcore_check/ref/ref_original.py` | 참조 대상 1: 원본을 그대로 감쌈 | `2d866e985d6c1e37728b82c75091836c9b2a8bd7ac5c4aef8f79054bbbacce3b` |
| `fastcore_check/ref/ref_rebuilt.py` | 참조 대상 2: 결과를 새 객체로 재구성 | `a73abc751a0a1567114fc63e77328edf6b5090f49c93ff35fc2fb644bdf21122` |
| `fastcore_check/ref/ref_mutant_identity.py` | 참조 대상 3: 가짜 공장 무패치 사본 | `eefa8b4c0b507bdd8ba1fe873cabe4b657a658be1cc0197be3e6953ff8c81ef0` |
| `fastcore_check/FC1_CHECK_SELFTEST_RESULT.json` | 자기시험(selftest, 1,516 사례) 서버 실행 요약 | `a193e3cbafa6940eec6e54f7a1df87a494b5fafa9e49903cf93663eba7bf198c` |
| `fastcore_check/results/sens_a/summary.json` | 민감도 실행(sens, 606 사례) 요약; 사례별 기록 shard_*.jsonl 동봉 | `717ca2f8169270eb9de8372893f788ff1e7388c880849f203b4af26ae3ae1ddd` |
| `fastcore_check/results/evidence/expwrap_plan2_shard0.jsonl` | `exp_wrap32`의 `plan` 프로파일 적발 기록(사례 #252) | `af778c0faf2b64f0dec968007b8437dd81ceac87f5909c3599ed6db0f5b06104` |
| `fastcore_check/results/plan_ref_rebuilt/summary.json` | `plan` 프로파일 × `ref_rebuilt` 부분 실행 요약 | `c75bce9b18180d8797759a455e3fe399c6b21c746c5a83e55a64276774d6184c` |
| `fastcore_check/fakes/fake_*.py` | 가짜·동치 변형 47개 (생성기 출력 45개 + 수동 작성 2개: `dim_wrong`, `stateful`) | (생성기에서 재현) |

오라클 해시(봉인과 일치): `dyad.py` `4d49cf8f…abaab`, `taylor_step.py` `0183f0f8…c01`, `lohner_step.py` `936138ea…db6`, 계획서 `37276492…2ed`.

